extern crate bio;
extern crate clap;
extern crate flate2;

use clap::{App, Arg};

use bio::io::fastq::Reader;

use flate2::bufread::MultiGzDecoder;

use std::fs::File;

use std::io::BufReader;

fn main() {
    let args = App::new("extract-amplicons")
        .version("0.1.0")
        .arg(Arg::with_name("R1").index(1))
        .arg(Arg::with_name("R2").index(2))
        .get_matches();
    let fq1 = Reader::new(MultiGzDecoder::new(BufReader::new(
        File::open(args.value_of("R1").unwrap()).unwrap(),
    )))
    .records();
    let fq2 = Reader::new(MultiGzDecoder::new(BufReader::new(
        File::open(args.value_of("R2").unwrap()).unwrap(),
    )))
    .records();

    let mut total_pairs = 0;
    for (_r1, _r2) in fq1.zip(fq2) {
        total_pairs += 1;
    }
    println!("{}", total_pairs);
}
