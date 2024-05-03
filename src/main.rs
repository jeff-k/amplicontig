use flate2::read::MultiGzDecoder;

//use core::ops::Bound::{Excluded, Included};
use core::ops::Bound::Included;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use clap::Parser;
use store_interval_tree::{Interval, IntervalTree}; //, IntervalTreeIterator};

use ampliconlib::aligner::{edit_dist, merge_bin, pp, Assembly};

use bio_streams::fasta::Fasta;
use bio_streams::fastq::Fastq;

use bio_seq::prelude::*;

use ampliconlib::primerset::{
    Amplicon::{Merged, Paired},
    Orientation::{F1R2, F2R1, R1F2, R2F1},
    PrimerSet,
};

#[derive(Parser)]
struct Cli {
    primers: PathBuf,
    reference: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
}

fn main() {
    let args = Cli::parse();

    //    let mut stats = Stats::new();

    let fq1: Fastq<BufReader<MultiGzDecoder<File>>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(File::open(&args.r1).unwrap()),
    ));
    let fq2: Fastq<BufReader<MultiGzDecoder<File>>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(File::open(&args.r2).unwrap()),
    ));

    let mut reference: Fasta<BufReader<File>> =
        Fasta::new(BufReader::new(File::open(&args.reference).unwrap()));

    let ref_seq: Seq<Dna> = reference.next().unwrap().unwrap().seq; //    let aligner = Aligner::new(&reference.next().unwrap().unwrap().seq);

    let primers = PrimerSet::from_csv(&args.primers, Some(&ref_seq));
    let mut f1r2 = 0;
    let mut f2r1 = 0;
    let mut r1f2 = 0;
    let mut r2f1 = 0;

    let mut total = 0;
    let mut merged = 0;
    let mut _mapped = 0;

    let mut invalid_reads = 0;

    //let mut bins: HashMap<Seq<Dna>, Assembly> = HashMap::new();
    let bins: HashMap<(String, String), usize> = HashMap::new();
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
                    println!("ending early");
                    break;
                }
                match primers.get_amplicon(&r1.seq, &r2.seq) {
                    Merged(orientation, p1, p2, seq) => {
                        let start = p1.index;
                        let end = p2.index;
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
                    Paired(_orientation, _start, _end) => {
                        /*
                        match orientation {
                            F1R2 => f1r2 += 1,
                            F2R1 => f2r1 += 1,
                            R1F2 => r1f2 += 1,
                            R2F1 => r2f1 += 1,
                        };
                        */
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
            if v.end - v.start > 400 || k.len() > 400 {
                continue;
            }
            let ref_seg: &SeqSlice<Dna> = &ref_seq[v.start..v.end];
            if v.count > 5 {
                let (score, ops) = edit_dist(&ref_seg.to_string(), &k.to_string());
                if score > 0 {
                    println!("\n\nDISTANCE:\n{}\n{}\n{}\n\n", &ref_seg, pp(&ops), &k);
                }
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
