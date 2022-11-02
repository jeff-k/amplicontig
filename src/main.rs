extern crate clap;
extern crate csv;
extern crate flate2;
//extern crate serde;
mod primerset;

use flate2::read::MultiGzDecoder;

use core::cmp::{max, min};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::str;

use clap::Parser;

mod mating;
use mating::{mate, merge, rc};

use bio_streams::fastq::Fastq;
use bio_streams::Record;

use crate::primerset::{PrimerSet, Stats};

#[derive(Parser)]
struct Cli {
    primers: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
}

struct Assembly {
    count: usize,
    start: usize,
    end: usize,
}

fn main() {
    let args = Cli::parse();
    let primers = PrimerSet::from_csv(&args.primers);

    let mut stats = Stats::new();

    let fq1: Fastq<BufReader<MultiGzDecoder<File>>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(File::open(&args.r1).unwrap()),
    ));
    let fq2: Fastq<BufReader<MultiGzDecoder<File>>> = Fastq::new(BufReader::new(
        MultiGzDecoder::new(File::open(&args.r2).unwrap()),
    ));

    let mut mm = 0;
    let mut f1r2 = 0;
    let mut f2r1 = 0;
    let mut r1f2 = 0;
    let mut r2f1 = 0;

    let mut bins: HashMap<Vec<u8>, Assembly> = HashMap::new();
    for (r1, r2) in fq1.zip(fq2) {
        //let r1rc = rc(r1.seq);
        //let r2rc = rc(r2.seq);

        match (primers.get(&r1.seq), primers.get(&r2.seq)) {
            (Some(p1), Some(p2)) => {
                let start = min(p1.index, p2.index);
                let end = max(p1.index, p2.index);
                let max_indel = 80;
                let hint = ((end - start) / 2) - 1;
                if hint < r1.seq.len() - 30 {
                    let merged = if p1.index < p2.index {
                        if p1.forward && !p2.forward {
                            // F1R2
                            f1r2 += 1;
                            let r2rc = rc(&r2.seq);
                            match mate(&r1.seq, &r2rc, hint, max_indel) {
                                Some(seam) => Some(merge(&r1.seq, &r2rc, seam)),
                                None => None,
                            }
                        } else if (!p1.forward && p2.forward) {
                            // R1F2
                            r1f2 += 1;
                            let r1rc = rc(&r1.seq);
                            match mate(&r1rc, &r2.seq, hint, max_indel) {
                                Some(seam) => Some(merge(&r1rc, &r2.seq, seam)),
                                None => None,
                            }
                        } else {
                            None
                        }
                    } else if p2.index < p1.index {
                        if p1.forward && !p2.forward {
                            // F2R1
                            f2r1 += 1;
                            let r1rc = rc(&r1.seq);
                            match mate(&r2.seq, &r1rc, hint, max_indel) {
                                Some(seam) => Some(merge(&r2.seq, &r1rc, seam)),
                                None => None,
                            }
                        } else if (!p1.forward && p2.forward) {
                            // R2F1
                            r2f1 += 1;
                            let r2rc = rc(&r2.seq);
                            match mate(&r2rc, &r1.seq, hint, max_indel) {
                                Some(seam) => Some(merge(&r2rc, &r1.seq, seam)),
                                None => None,
                            }
                        } else {
                            None
                        }
                    } else {
                        // ???
                        None
                    };

                    if let Some(m) = merged {
                        bins.entry(m)
                            .or_insert(Assembly {
                                count: 0,
                                start,
                                end,
                            })
                            .count += 1;
                    }
                }
            }
            _ => {
                /*
                let r1rc = rc(&r1.seq);
                let r2rc = rc(&r2.seq);
                match (primers.get(&r1rc), primers.get(&r2rc)) {
                    (Some(p1), Some(p2)) => {
                        let start = min(p1.index, p2.index);
                        let end = max(p1.index, p2.index);
                        bins.entry(r1.seq)
                            .or_insert(Assembly {
                                count: 0,
                                start,
                                end,
                            })
                            .count += 1;
                    }
                    _ => mm += 1,
                }
                    */
                ()
            }
        }
        //primers.get(r1rc);
    }
    for (k, v) in bins {
        if v.count > 10 {
            println!(
                "{}-{}\t{}\tlen: {}\t{}:\t{}",
                v.start,
                v.end,
                v.end - v.start,
                &k.len(),
                String::from_utf8_lossy(&k),
                &v.count
            );
        }
    }
    println!(
        "r1f2: {}\tf1r2: {}\tr2f1: {}\tf2r1:\t{}\nvs. {}",
        r1f2, f1r2, r2f1, f2r1, mm
    );
}
