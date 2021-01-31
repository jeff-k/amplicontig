extern crate bio;
extern crate clap;
extern crate csv;
extern crate flate2;
extern crate serde;
extern crate trie;

use clap::{App, Arg};

use bio::io::fastq::{Reader, Record};

use flate2::bufread::MultiGzDecoder;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::str;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct PrimerSet {
    name: String,
    primer: String,
    index: u32,
    reference: String,
}

fn main() {
    let args = App::new("extract-amplicons")
        .version("0.1.0")
        .arg(Arg::with_name("primers").index(1))
        .arg(Arg::with_name("R1").index(2))
        .arg(Arg::with_name("R2").index(3))
        .get_matches();

    let mut primers = HashMap::new();
    let mut off_target = HashMap::new();

    let mut csv = csv::Reader::from_reader(File::open(args.value_of("primers").unwrap()).unwrap());
    for result in csv.deserialize() {
        let record: PrimerSet = result.unwrap();
        let primer = String::from(&record.primer[..20]);
        primers.insert(primer, String::from(&record.name));
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
    let grep = true;
    for (record1, record2) in fq1.zip(fq2) {
        match (record1, record2) {
            (Ok(r1), Ok(r2)) => {
                total_pairs += 1;
                let p1 = String::from_utf8(r1.seq()[..20].to_vec()).unwrap();
                let p2 = String::from_utf8(r2.seq()[..20].to_vec()).unwrap();
                match (primers.get(&p1), primers.get(&p2)) {
                    (Some(p1name), Some(p2name)) => {
                        if grep {
                            let desc1 = format!("{}:{}", p1name, r1.desc().unwrap());
                            print!(
                                "{}",
                                Record::with_attrs(r1.id(), Some(&desc1), r1.seq(), r1.qual())
                            );

                            let desc2 = format!("{}:{}", p2name, r2.desc().unwrap());
                            print!(
                                "{}",
                                Record::with_attrs(r2.id(), Some(&desc2), r2.seq(), r2.qual())
                            );
                        };

                        matched += 1
                    }
                    (Some(_), None) => {
                        *off_target.entry(p2).or_insert(0) += 1;
                    }
                    (None, Some(_)) => {
                        *off_target.entry(p1).or_insert(0) += 1;
                    }
                    (None, None) => {
                        *off_target.entry(p1).or_insert(0) += 1;
                        *off_target.entry(p2).or_insert(0) += 1;
                    }
                }
            }
            _ => eprintln!("Not proper fastq file pair"),
        }
    }
    println!("{}/{}", matched, total_pairs);
    println!("{:?}", off_target);
}
