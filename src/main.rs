extern crate bio;
extern crate clap;
extern crate csv;
extern crate flate2;
extern crate serde;

use clap::{App, Arg};

use bio::io::fastq::Reader;

use flate2::bufread::MultiGzDecoder;

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io;
use std::io::BufReader;
use std::str;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct Record {
    primer: String,
    left: String,
    right: String,
    index: String,
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

    let mut csv = csv::Reader::from_reader(File::open(args.value_of("primers").unwrap()).unwrap());
    for result in csv.deserialize() {
        let record: Record = result.unwrap();
        println!("{:?}", record);
        //y        let key = record.c[..21];
        let primer = (
            String::from(&record.left[..20]),
            String::from(&record.right[..20]),
        );
        primers.insert(primer, 0);
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
    for (record1, record2) in fq1.zip(fq2) {
        match (record1, record2) {
            (Ok(r1), Ok(r2)) => {
                total_pairs += 1;
                let primer = (
                    String::from_utf8(r1.seq()[..20].to_vec()).unwrap(),
                    String::from_utf8(r2.seq()[..20].to_vec()).unwrap(),
                );
                match primers.get(&primer) {
                    Some(_) => matched += 1,
                    None => (),
                }
            }
            _ => eprintln!("Not proper fastq file pair"),
        }
    }
    println!("{}/{}", matched, total_pairs);
}
