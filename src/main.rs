extern crate bio;
#[macro_use]
extern crate clap;
extern crate csv;
extern crate flate2;
extern crate fst;
extern crate serde;

use clap::{App, Arg};

use bio::io::fastq::{Reader, Record};

use flate2::bufread::MultiGzDecoder;

use std::fs::File;
use std::io::BufReader;
use std::process::exit;
use std::str;

mod mating;
mod primerset;

use mating::merge_records;
use primerset::{MatchedReads, PrimerSet, Stats};

fn printrec(r: &Record, pname: &str, start: usize, end: usize) {
    let desc = format!("{}:{}", pname, r.desc().unwrap());
    print!(
        "{}",
        Record::with_attrs(
            r.id(),
            Some(&desc),
            &r.seq()[start..end],
            &r.qual()[start..end]
        )
    );
}

fn main() {
    let args = clap_app!(myapp =>
        (version: "0.1.4")
        (author: "github.com/jeff-k")
        (about: "assemble reads from amplicon sequencing data")
        (@arg primers: +required "primer set")
        (@arg R1: +required "R1 reads")
        (@arg R2: +required "R2 reads")
        (@arg prefix: -p --prefix "set the output file(s) prefix")
        (@arg stats: -v "report primer match stats")
        (@subcommand test =>
            (about: "test reads against a set of primers")
        )
        (@subcommand match =>
            (about: "match reads against a primer set")
            (@arg trim: -t --trim "trim primer bases")
            (@arg merge: -m --merge "merge mated reads")
            (@arg filter: -f --filter "save unmatched reads")
        )
        (@subcommand assemble =>
            (about: "bin matched and merged read pairs into consensus")
        )
    )
    .get_matches();

    let primers = PrimerSet::from_csv(args.value_of("primers").unwrap().to_string());
    let mut stats = Stats {
        on_target: 0,
        off_target: 0,
        total_pairs: 0,
        matched: 0,
        mated: 0,
    };

    let trim: usize = match args.value_of("trim") {
        Some(s) => s.parse::<usize>().unwrap_or(0),
        None => 0,
    };

    let fq1 = Reader::new(MultiGzDecoder::new(BufReader::new(
        File::open(args.value_of("R1").unwrap()).unwrap(),
    )))
    .records();
    let fq2 = Reader::new(MultiGzDecoder::new(BufReader::new(
        File::open(args.value_of("R2").unwrap()).unwrap(),
    )))
    .records();

    let test_pset = args.is_present("test");
    let grep = !(args.is_present("stats"));
    let invert = args.is_present("invert");
    let excise = args.is_present("ex");

    for (read1, p1, read2, p2) in MatchedReads::new(Box::new(fq1.zip(fq2)), &primers) {
        match (p1, p2) {
            (Some(_), Some(_)) => stats.matched += 2,
            (Some(_), None) => stats.matched += 1,
            (None, Some(_)) => stats.matched += 1,
            _ => (),
        }
        stats.mated += 2;

        match (read1, read2) {
            (r1, r2) => {
                if let Some(r) = merge_records(&r1, &r2) {
                    let rlen = r.seq().len();
                    if rlen < 100 {
                        continue;
                    }

                    stats.mated += 1;
                }
            }
            _ => {
                eprintln!("Not proper fastq file pair");
                exit(1);
            }
        }
    }

    if test_pset {
        eprintln!("{:?},{:?}", stats.matched, stats.total_pairs);
        exit(0);
    }
    //for (primer, reads) in stats.readbins {
    //    print!(
    //        "{}\t{}\t{}",
    //        args.value_of("R1").unwrap(),
    //        primer,
    //        match stats.on_target.get(&primer) {
    //            Some(n) => *n,
    //            None => 0,
    //        }
    //    );
    //}
}
