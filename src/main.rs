extern crate clap;
extern crate csv;
extern crate flate2;
//extern crate serde;
mod aligner;
mod mating;
mod primerset;

use flate2::read::MultiGzDecoder;

use core::cmp::{max, min, Ordering};
use core::ops::Bound::{Excluded, Included};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use clap::Parser;

use aligner::Aligner;
use aligner::Alignment::{Forward, Reverse};
use mating::{mate, merge, rc};

use bio_streams::fasta::Fasta;
use bio_streams::fastq::Fastq;

use crate::primerset::{Primer, PrimerSet};
use crate::Amplicon::{Discarded, Merged, Paired};
use crate::Orientation::{F1R2, F2R1, R1F2, R2F1};

use store_interval_tree::{Interval, IntervalTree, IntervalTreeIterator};

#[derive(Parser)]
struct Cli {
    primers: PathBuf,
    reference: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
}

#[derive(Debug)]
struct Assembly {
    count: usize,
    //fwds: usize,
    //revs: usize,
    start: usize,
    end: usize,
}

enum Orientation {
    F1R2,
    F2R1,
    R1F2,
    R2F1,
}

enum Amplicon {
    Discarded,
    Merged(Orientation, usize, Vec<u8>),
    Paired(Orientation, usize, Vec<u8>, usize, Vec<u8>),
}

fn fix(template: &mut Vec<u8>, q: &Vec<u8>) {
    if template.len() != q.len() {
        return;
    }
    for (a, b) in template.iter_mut().zip(q) {
        if a != b {
            if *a == b'N' {
                *a = *b;
            }
        }
    }
}

fn merge_bin(bin: &HashMap<Vec<u8>, Assembly>) -> HashMap<Vec<u8>, Assembly> {
    let mut new = HashMap::new();
    let mut template: Vec<u8> = Vec::new();
    let mut assembly = Assembly {
        count: 0,
        start: 0,
        end: 0,
    };
    let mut max = 0;

    for (k, v) in bin {
        if v.count > max {
            max = v.count;
            template = k.to_vec();
            assembly.count = v.count;
            assembly.start = v.start;
            assembly.end = v.end;
        }
    }

    for (k, v) in bin {
        fix(&mut template, &k);
        assembly.count += v.count;
    }

    new.insert(template, assembly);
    new
}

#[inline]
fn merge_amplicon(p1: &Primer, r1: &[u8], p2: &Primer, r2: &[u8]) -> Amplicon {
    let start = min(p1.index, p2.index);
    let end = max(p1.index, p2.index);
    let max_indel = 80;
    let hint = ((end - start) / 2) - 1;
    match p1.index.cmp(&p2.index) {
        Ordering::Less => {
            if p1.forward && !p2.forward {
                // F1R2
                let r2rc = rc(r2);
                if hint < r1.len() - 30 {
                    match mate(r1, &r2rc, hint, max_indel) {
                        Some(seam) => Merged(F1R2, p1.index, merge(r1, &r2rc, seam)),
                        None => Paired(F1R2, p1.index, r1.to_vec(), p2.index, r2rc),
                    }
                } else {
                    Paired(F1R2, p1.index, r1.to_vec(), p2.index, r2rc)
                }
            } else if !p1.forward && p2.forward {
                // R1F2
                let r1rc = rc(r1);
                if hint < r1.len() - 30 {
                    match mate(&r1rc, r2, hint, max_indel) {
                        Some(seam) => Merged(R1F2, p1.index, merge(&r1rc, r2, seam)),
                        None => Paired(R1F2, p1.index, r1rc, p2.index, r2.to_vec()),
                    }
                } else {
                    Paired(R1F2, p1.index, r1rc, p2.index, r2.to_vec())
                }
            } else {
                Discarded
            }
        }
        _ => {
            if p1.forward && !p2.forward {
                // F2R1
                let r1rc = rc(r1);
                match mate(r2, &r1rc, hint, max_indel) {
                    Some(seam) => Merged(F2R1, p1.index, merge(r2, &r1rc, seam)),
                    None => Paired(F2R1, p2.index, r2.to_vec(), p1.index, r1rc),
                }
            } else if !p1.forward && p2.forward {
                // R2F1
                let r2rc = rc(r2);
                match mate(&r2rc, r1, hint, max_indel) {
                    Some(seam) => Merged(R2F1, p1.index, merge(&r2rc, r1, seam)),
                    None => Paired(R2F1, p2.index, r2rc, p1.index, r1.to_vec()),
                }
            } else {
                Discarded
            }
        }
    }
}

fn main() {
    let args = Cli::parse();
    let primers = PrimerSet::from_csv(&args.primers);

    //    let mut stats = Stats::new();

    let fq1: Fastq<BufReader<MultiGzDecoder<File>>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(File::open(&args.r1).unwrap()),
    ));
    let fq2: Fastq<BufReader<MultiGzDecoder<File>>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(File::open(&args.r2).unwrap()),
    ));

    let mut reference: Fasta<BufReader<File>> =
        Fasta::new(BufReader::new(File::open(&args.reference).unwrap()));

    let aligner = Aligner::new(&reference.next().unwrap().unwrap().seq);

    let mut f1r2 = 0;
    let mut f2r1 = 0;
    let mut r1f2 = 0;
    let mut r2f1 = 0;

    let mut total = 0;
    let mut merged = 0;
    let mut mapped = 0;

    let mut bins: HashMap<Vec<u8>, Assembly> = HashMap::new();
    let mut leftover: HashMap<Vec<u8>, usize> = HashMap::new();
    //    let mut tree: IntervalTree<usize, Vec<u8>> = IntervalTree::new();
    let mut tree: IntervalTree<usize, ()> = IntervalTree::new();
    let mut ibins: HashMap<Interval<usize>, HashMap<Vec<u8>, Assembly>> = HashMap::new();

    for (r1, r2) in fq1.zip(fq2) {
        total += 1;
        let amplicon = match (primers.get(&r1.unwrap().seq), primers.get(&r2.unwrap().seq)) {
            (Some(p1), Some(p2)) => {
                merged += 1;
                merge_amplicon(p1, &r1.unwrap().seq, p2, &r2.unwrap().seq)
            }
            _ => match (aligner.get(&r1.unwrap().seq), aligner.get(&r2.unwrap().seq)) {
                (Forward(r1pos), Reverse(r2pos)) => {
                    mapped += 1;
                    Paired(F1R2, r1pos, r1.seq, r2pos, rc(&r2.unwrap().seq))
                }
                (Reverse(r1pos), Forward(r2pos)) => {
                    mapped += 1;
                    Paired(R1F2, r1pos, rc(&r1.unwrap().seq), r2pos, r2.unwrap().seq)
                }
                _ => Discarded,
            },
        };
        match amplicon {
            Merged(orientation, pos, seq) => {
                match orientation {
                    F1R2 => f1r2 += 1,
                    F2R1 => f2r1 += 1,
                    R1F2 => r1f2 += 1,
                    R2F1 => r2f1 += 1,
                };
                let len = pos + &seq.len();
                let interval = Interval::new(Included(pos), Included(len));

                tree.insert(interval.clone(), ());
                ibins
                    .entry(interval)
                    .or_insert(HashMap::from([(
                        seq.clone(),
                        Assembly {
                            count: 0,
                            start: pos,
                            end: len,
                        },
                    )]))
                    .entry(seq)
                    .or_insert(Assembly {
                        count: 0,
                        start: pos,
                        end: len,
                    })
                    .count += 1;
            }
            Paired(orientation, start, r1, end, r2) => {
                match orientation {
                    F1R2 => f1r2 += 1,
                    F2R1 => f2r1 += 1,
                    R1F2 => r1f2 += 1,
                    R2F1 => r2f1 += 1,
                };

                /*
                bins.entry(r1)
                    .or_insert(Assembly {
                        count: 0,
                        start,
                        end,
                    })
                    .count += 1;
                bins.entry(r2)
                    .or_insert(Assembly {
                        count: 0,
                        start,
                        end,
                    })
                    .count += 1;
                    */
            }
            _ => (),
        }
    }

    /*
    println!(
        "{:?}",
        tree.intervals_between(&Interval::point(14), &Interval::point(29000))
    );
    */

    for interval in tree.intervals() {
        let bin = ibins.get(&interval).unwrap();

        for (k, v) in merge_bin(bin) {
            if v.count > 10 {
                println!(
                    ">{}-{},ref_length:{},count:{},lenth:{}\n{}",
                    v.start,
                    v.end,
                    v.end - v.start,
                    &v.count,
                    &k.len(),
                    String::from_utf8_lossy(&k[20..k.len() - 20]),
                );
            }
        }
    }
    eprintln!(
        "r1f2: {}\tf1r2: {}\tr2f1: {}\tf2r1: {}\tmerged: {}\tmapped: {}\ttotal: {}",
        r1f2, f1r2, r2f1, f2r1, merged, mapped, total
    );
}
