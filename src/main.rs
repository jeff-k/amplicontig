extern crate bio;
extern crate clap;
extern crate csv;
extern crate flate2;
extern crate serde;

use clap::{App, Arg};

use bio::alignment::pairwise::*;
use bio::io::fastq::{Reader, Record};

use flate2::bufread::MultiGzDecoder;

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::BufReader;
use std::str;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct PrimerSet {
    name: String,
    primer: String,
    length: u32,
    index: u32,
    reference: String,
}

fn printrec(r: &Record, pname: &str) {
    let desc = format!("{}:{}", pname, r.desc().unwrap());
    print!(
        "{}",
        Record::with_attrs(r.id(), Some(&desc), r.seq(), r.qual())
    );
}

fn main() {
    let args = App::new("extract-amplicons")
        .version("0.1.0")
        .arg(
            Arg::with_name("primers")
                .index(1)
                .help("csv file for primer set"),
        )
        .arg(Arg::with_name("R1").index(2).help("first reads in pairs"))
        .arg(
            Arg::with_name("R2")
                .index(3)
                .help("second reads in pairs (reversed)"),
        )
        .arg(
            Arg::with_name("invert")
                .short("n")
                .required(false)
                .takes_value(false)
                .help("invert selection (show unmatching reads)"),
        )
        .arg(
            Arg::with_name("stats")
                .short("s")
                .required(false)
                .takes_value(false)
                .help("print stats about primer performance"),
        )
        .get_matches();

    let mut primers = HashMap::new();
    let mut off_target = HashMap::new();
    let mut on_target = HashMap::new();

    let mut csv = csv::Reader::from_reader(File::open(args.value_of("primers").unwrap()).unwrap());
    let mut readbins = HashMap::new();
    let mut seqs = HashMap::new();
    for result in csv.deserialize() {
        let record: PrimerSet = result.unwrap();
        let primer = String::from(&record.primer[..20]);
        let name = String::from(&record.name);
        primers.insert(primer, String::from(&record.name));
        readbins.insert(name, BTreeMap::new());
        seqs.insert(String::from(&record.name), String::from(&record.reference));
    }

    let fq1 = Reader::new(MultiGzDecoder::new(BufReader::new(
        File::open(args.value_of("R1").unwrap()).unwrap(),
    )))
    .records();
    let fq2 = Reader::new(MultiGzDecoder::new(BufReader::new(
        File::open(args.value_of("R2").unwrap()).unwrap(),
    )))
    .records();

    let mut total_pairs = 0;
    let mut matched = 0;
    let grep = !(args.is_present("stats"));
    let invert = args.is_present("invert");

    for (record1, record2) in fq1.zip(fq2) {
        match (record1, record2) {
            (Ok(r1), Ok(r2)) => {
                total_pairs += 1;
                let p1 = String::from_utf8(r1.seq()[..20].to_vec()).unwrap();
                let p2 = String::from_utf8(r2.seq()[..20].to_vec()).unwrap();
                match (primers.get(&p1), primers.get(&p2)) {
                    (Some(p1name), Some(p2name)) => {
                        if !(invert) & grep {
                            printrec(&r1, p1name);
                            printrec(&r2, p2name);
                        };
                        let seq = String::from_utf8(r1.seq()[20..].to_vec()).unwrap();
                        *(readbins.get_mut(&String::from(p1name)).unwrap())
                            .entry(seq)
                            .or_insert(0) += 1;
                        *on_target
                            .entry(format!("{}:{}", p1name, p2name))
                            .or_insert(0) += 1;
                        matched += 1
                    }
                    (Some(p1name), None) => {
                        *off_target.entry(p2).or_insert(0) += 1;
                        if invert & grep {
                            printrec(&r1, p1name);
                            print!("{}", r2);
                        }
                    }
                    (None, Some(p2name)) => {
                        *off_target.entry(p1).or_insert(0) += 1;
                        if invert & grep {
                            print!("{}", r1);
                            printrec(&r2, p2name);
                        }
                    }
                    (None, None) => {
                        *off_target.entry(p1).or_insert(0) += 1;
                        *off_target.entry(p2).or_insert(0) += 1;

                        if invert & grep {
                            print!("{}", r1);
                            print!("{}", r2);
                        }
                    }
                }
            }
            _ => eprintln!("Not proper fastq file pair"),
        }
    }
    if !grep {
        for (primer, count) in if invert { off_target } else { on_target } {
            if count > 1000 {
                println!("{}\t{}", primer, count);
            }
        }
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        for (primer, reads) in readbins {
            println!("{}", primer);
            for (read, count) in reads {
                if count > 500 {
                    //                    println!("\t{}", seqs.get(&primer).unwrap());
                    let r = seqs.get(&primer).unwrap().as_bytes();
                    let x = read.as_bytes();
                    let mut aligner = Aligner::with_capacity(r.len(), x.len(), -5, -1, &score);
                    let alignment = aligner.semiglobal(r, x);
                    println!("depth: {}, score: {}", count, alignment.score);
                    if alignment.score < 147 {
                        println!("{}", alignment.pretty(r, x));
                    }
                }
            }
        }
        eprintln!("{}/{} read pairs matched.", matched, total_pairs);
    }
}
