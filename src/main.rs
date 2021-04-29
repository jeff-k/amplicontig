extern crate bio;
extern crate clap;
extern crate csv;
extern crate flate2;
extern crate serde;

use clap::{App, Arg};

use bio::alphabets::dna;
use bio::io::fastq::{Error, Reader, Record};

use flate2::bufread::MultiGzDecoder;

use serde::Deserialize;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process::exit;
use std::str;

mod mating;
use mating::{mate, mend_consensus, merge, truncate};

#[derive(Debug, Deserialize, Clone)]
struct PrimerSet {
    name: String,
    primer: String,
    forward: bool,
    index: u32,
}

struct MatchedReads {
    zipfq: Box<dyn Iterator<Item = (Result<Record, Error>, Result<Record, Error>)>>,
    primers: HashMap<String, PrimerSet>,
}

impl Iterator for MatchedReads {
    type Item = (Record, Option<PrimerSet>, Record, Option<PrimerSet>);

    fn next(&mut self) -> Option<Self::Item> {
        match self.zipfq.next() {
            Some((rec1, rec2)) => match (rec1, rec2) {
                (Ok(r1), Ok(r2)) => {
                    let plen = 22;
                    let p1 = String::from_utf8(r1.seq()[..plen].to_vec()).unwrap();
                    let p2 = String::from_utf8(r2.seq()[..plen].to_vec()).unwrap();

                    let m1 = self.primers.get(&p1);
                    let m2 = self.primers.get(&p2);
                    Some((
                        r1,
                        match m1 {
                            None => None,
                            Some(p) => Some(p.clone()),
                        },
                        r2,
                        match m2 {
                            None => None,
                            Some(p) => Some(p.clone()),
                        },
                    ))
                }
                _ => {
                    eprintln!("Mismatched fastqs files (different lengths)");
                    exit(1);
                }
            },
            None => None,
        }
    }
}

impl MatchedReads {
    fn new(
        zipfq: Box<dyn Iterator<Item = (Result<Record, Error>, Result<Record, Error>)>>,
        primers: HashMap<String, PrimerSet>,
    ) -> MatchedReads {
        MatchedReads {
            zipfq: zipfq,
            primers: primers,
        }
    }
}

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

fn merge_records(r1: &Record, r2: &Record) -> Option<Record> {
    let r2_rc = dna::revcomp(r2.seq());
    let r1_rc = dna::revcomp(r1.seq());

    match mate(&r1.seq(), &r2_rc, 25, 20) {
        Some(overlap) => {
            let seq = merge(&r1.seq(), &r2_rc, overlap, mend_consensus);
            let qual = merge(&r1.qual(), &r2.qual(), overlap, max);
            Some(Record::with_attrs(r1.id(), None, &seq, &qual))
        }
        None => match mate(&r1_rc, &r2.seq(), 25, 20) {
            Some(overlap) => {
                let seq = truncate(&r1.seq(), &r2_rc, overlap, mend_consensus);
                let qual = truncate(&r1.qual(), &r2.qual(), overlap, max);
                Some(Record::with_attrs(r1.id(), None, &seq, &qual))
            }
            None => None,
        },
    }
}

fn main() {
    let args = App::new("amplicontig")
        .version("0.1.2")
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
        .arg(
            Arg::with_name("trim")
                .short("t")
                .required(false)
                .takes_value(true)
                .help("trim bases from 3' end"),
        )
        .arg(
            Arg::with_name("test")
                .required(false)
                .short("z")
                .takes_value(false)
                .help("test whether fastq(s) match primer set"),
        )
        .arg(
            Arg::with_name("ex")
                .short("x")
                .required(false)
                .takes_value(false)
                .help("excise primer sequence from reads"),
        )
        .get_matches();

    let mut primers: HashMap<String, PrimerSet> = HashMap::new();
    let mut off_target: HashMap<String, u32> = HashMap::new();
    let mut on_target: HashMap<String, u32> = HashMap::new();

    //    let mut readbins = HashMap::new();
    let mut plen = 100;
    let mut pmax = 0;

    let trim: usize = match args.value_of("trim") {
        Some(s) => s.parse::<usize>().unwrap_or(0),
        None => 0,
    };

    let mut primer_recs = Vec::new();
    for result in csv::Reader::from_reader(File::open(args.value_of("primers").unwrap()).unwrap())
        .deserialize()
    {
        let record: PrimerSet = result.unwrap();
        let l = record.primer.len();
        plen = min(l, plen);
        pmax = max(l, pmax);
        primer_recs.push(record);
    }

    println!("{:?}", plen);

    for record in primer_recs {
        let primer = String::from(&record.primer[..plen]);
        let name = String::from(&record.name);
        primers.insert(primer, record);
        //        readbins.insert(name, BTreeMap::new());
        //        primers.insert(String::from(dna::revcomp(primer)), record);
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
    let mut total_mated = 0;
    let mut matched = 0;
    let test_pset = args.is_present("test");
    let grep = !(args.is_present("stats"));
    let invert = args.is_present("invert");
    let excise = args.is_present("ex");

    //    for a,b,c,d in MatchedReads::new(fq1, fq2, primers) {

    for (a, b, c, d) in MatchedReads::new(Box::new(fq1.zip(fq2)), primers) {
        match (b, d) {
            (Some(_), Some(_)) => matched += 2,
            (Some(_), _) => matched += 1,
            (_, Some(_)) => matched += 1,
            _ => (),
        }
        total_mated += 1;

        println!("{:?}/{:?}", matched, total_mated);
    }

    //        total_pairs += 1;
    //        match (record1, record2) {
    //            (Ok(r1), Ok(r2)) => {
    //                if let Some(r) = merge_records(&r1, &r2) {
    //                    let rlen = r.seq().len();
    //                    if rlen < 100 {
    //                        continue;
    //                    }
    //                    let p1 = String::from_utf8(r.seq()[..plen].to_vec()).unwrap();
    //                    let p2 = String::from_utf8(r.seq()[..plen].to_vec()).unwrap();
    //
    //                    total_mated += 1;
    //                    match (primers.get(&p1), primers.get(&p2)) {
    //                        (Some((p1name, is_forward, p1len)), Some((p2name, _, p2len))) => {
    //                            matched += 2;
    //                            //eprintln!("{} {} {} {}", p1name, *p1len, p2name, *p2len);
    //                            if !(invert) & grep {
    //                                printrec(
    //                                    &r,
    //                                    &(format!("{}:{}", p1name, p2name)),
    //                                    *p1len,
    //                                    rlen - *p2len,
    //                                );
    //                            };
    //
    //                            let seq = if *is_forward {
    //                                String::from_utf8(r.seq()[*p1len..(rlen - *p2len)].to_vec())
    //                                    .unwrap()
    //                            } else {
    //                                String::from_utf8(dna::revcomp(
    //                                    r.seq()[(*p1len)..(rlen - *p2len)].to_vec(),
    //                                ))
    //                                .unwrap()
    //                            };
    //
    //                            *(readbins.get_mut(&String::from(p1name)).unwrap())
    //                                .entry(seq)
    //                                .or_insert(0) += 1;
    //
    //                            *full_match
    //                                .entry(format!("{}:{}", p1name, p2name))
    //                                .or_insert(0) += 1;
    //                            *on_target.entry(p1name.to_string()).or_insert(0) += 1;
    //                            *on_target.entry(p2name.to_string()).or_insert(0) += 1;
    //                        }
    //                        (Some((p1name, _, p1len)), None) => {
    //                            *off_target.entry(p2).or_insert(0) += 1;
    //                            *on_target.entry(p1name.to_string()).or_insert(0) += 1;
    //                            matched += 1;
    //                            if invert & grep {
    //                                //printrec(&r1, p1name, *p1len);
    //                                print!("{}", r2);
    //                            }
    //                        }
    //                        (None, Some((p2name, _, p2len))) => {
    //                            *off_target.entry(p1).or_insert(0) += 1;
    //                            *on_target.entry(p2name.to_string()).or_insert(0) += 1;
    //                            matched += 1;
    //                            if invert & grep {
    //                                print!("{}", r1);
    //                                //printrec(&r2, p2name, *p2len);
    //                            }
    //                        }
    //                        (None, None) => {
    //                            *off_target.entry(p1).or_insert(0) += 1;
    //                            *off_target.entry(p2).or_insert(0) += 1;
    //                        }
    //                    }
    //                }
    //            }
    //            _ => {
    //                eprintln!("Not proper fastq file pair");
    //                exit(1);
    //            }
    //        }
    //    }
    //    if test_pset {
    //        println!("{:?},{:?}", matched, total_pairs);
    //        exit(0);
    //    }
    //    if !grep {
    //        //  for (primer, count) in on_target {
    //        //     if count > 1000 {
    //        //        println!("{}\t{}", primer, count);
    //        //     }
    //        //  }
    //        // for (primer, count) in off_target {
    //        //     if count > 1000 {
    //        //        println!("{}\t{}", primer, count);
    //        //     }
    //        // }
    //
    //        for (primer, reads) in readbins {
    //            print!(
    //                "{}\t{}\t{}",
    //                args.value_of("R1").unwrap(),
    //                primer,
    //                match on_target.get(&primer) {
    //                    Some(n) => *n,
    //                    None => 0,
    //                }
    //            );
    //        }
    //    }
    //    eprintln!(
    //        "{}\t{}/{}\t({}) read pairs matched.",
    //        args.value_of("R1").unwrap(),
    //        matched,
    //        total_mated,
    //        total_pairs
    //    );
}
