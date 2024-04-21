mod aligner;
mod mating;
mod primerset;

use flate2::read::MultiGzDecoder;

//use core::ops::Bound::{Excluded, Included};
use core::ops::Bound::Included;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use clap::Parser;

use aligner::edit_dist;

use bio_streams::fasta::Fasta;
use bio_streams::fastq::Fastq;

use bio_seq::prelude::*;

use crate::primerset::{
    Amplicon::{Discarded, Merged, Paired},
    Orientation::{F1R2, F2R1, R1F2, R2F1},
    Primer, PrimerSet,
};

use store_interval_tree::{Interval, IntervalTree}; //, IntervalTreeIterator};

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
fn merge_bin(bin: &HashMap<Seq<Dna>, Assembly>) -> HashMap<Seq<Dna>, Assembly> {
    let mut new = HashMap::new();
    let mut template: Seq<Dna> = Seq::new();
    let mut assembly = Assembly {
        count: 0,
        start: 0,
        end: 0,
    };
    let mut max = 0;

    for (k, v) in bin {
        if v.count > max {
            max = v.count;
            template = k.clone();
            assembly.count = v.count;
            assembly.start = v.start;
            assembly.end = v.end;
        }
    }

    for v in bin.values() {
        // Find majority allele?
        //fix(&mut template, &k);
        assembly.count += v.count;
    }

    new.insert(template, assembly);
    new
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

    //    let aligner = Aligner::new(&reference.next().unwrap().unwrap().seq);

    let mut f1r2 = 0;
    let mut f2r1 = 0;
    let mut r1f2 = 0;
    let mut r2f1 = 0;

    let mut total = 0;
    let mut merged = 0;
    let mut mapped = 0;

    let mut invalid_reads = 0;

    //let mut bins: HashMap<Seq<Dna>, Assembly> = HashMap::new();
    let mut bins: HashMap<(String, String), usize> = HashMap::new();
    //    let mut leftover: HashMap<Seq<Dna>, usize> = HashMap::new();
    //    let mut tree: IntervalTree<usize, Vec<u8>> = IntervalTree::new();
    let mut tree: IntervalTree<usize, ()> = IntervalTree::new();
    let mut ibins: HashMap<Interval<usize>, HashMap<Seq<Dna>, Assembly>> = HashMap::new();

    for (r1, r2) in fq1.zip(fq2) {
        match (r1, r2) {
            (Err(_e), _) | (_, Err(_e)) => {
                //println!("invalid sequence: {:?}", e);
                invalid_reads += 1;
                continue;
            }
            (Ok(r1), Ok(r2)) => {
                total += 1;
                if total > 1000 {
                    break;
                }
                match primers.get_amplicon(&r1.seq, &r2.seq) {
                    Merged(orientation, start, end, seq) => {
                        merged += 1;
                        match orientation {
                            F1R2 => f1r2 += 1,
                            F2R1 => f2r1 += 1,
                            R1F2 => r1f2 += 1,
                            R2F1 => r2f1 += 1,
                        };
                        let interval = Interval::new(Included(start), Included(end));

                        tree.insert(interval.clone(), ());
                        ibins
                            .entry(interval)
                            .or_insert(HashMap::from([(
                                seq.clone(),
                                Assembly {
                                    count: 0,
                                    start,
                                    end,
                                },
                            )]))
                            .entry(seq)
                            .or_insert(Assembly {
                                count: 0,
                                start,
                                end,
                            })
                            .count += 1;
                    }
                    Paired(orientation, _start, _end) => {
                        match orientation {
                            F1R2 => f1r2 += 1,
                            F2R1 => f2r1 += 1,
                            R1F2 => r1f2 += 1,
                            R2F1 => r2f1 += 1,
                        };
                    }
                    _ => (),
                }
            }
        }
    }

    println!(
        "{:?}",
        tree.intervals_between(&Interval::point(14), &Interval::point(29000))
    );

    for (k, v) in bins {
        if v > 100 {
            println!("{:?}: {}", k, v);
        }
    }

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
                    &k.to_string(),
                    //String::from_utf8_lossy(&k[20..k.len() - 20]),
                );
            }
        }
    }
    eprintln!(
        "r1f2: {}\tf1r2: {}\tr2f1: {}\tf2r1: {}\tmerged: {}\ttotal: {}\tinvalid: {}",
        r1f2, f1r2, r2f1, f2r1, merged, total, invalid_reads
    );
}
