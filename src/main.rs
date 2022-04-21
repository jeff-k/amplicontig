//extern crate bio;
//#[macro_use]
extern crate clap;
extern crate csv;
//extern crate flate2;
//extern crate serde;

use clap::Command;

//use bio::io::fastq::{Reader, Record};

//use flate2::bufread::MultiGzDecoder;

use std::str;

mod primerset;
//mod mating;

use crate::primerset::{PrimerSet, Stats};

fn main() {
    let args = Command::new("amplicontig")
        .version("0.3.0")
        .author("github.com/jeff-k")
        .about("amplicon aware contig assembly")
        .arg(clap::arg!(<primers> "primers").required(true))
        .arg(clap::arg!(<R1> "R1").required(true))
        .arg(clap::arg!(<R2> "R2").required(true))
        .subcommand_required(true)
        .subcommand(
            clap::command!("match")
                .arg(clap::arg!(--"trim"))
                .arg(clap::arg!(--"merge"))
                .arg(clap::arg!(--"filter")),
        )
        .subcommand(clap::command!("assemble"))
        .get_matches();

    match args.subcommand() {
        Some(("match", _match_args)) => {
            println!("doing match");
        }
        Some(("assemble", _)) => {
            println!("assembly");
        }
        _ => {
            println!("i owe you a coke");
        }
    }

    let primers = PrimerSet::from_csv(args.value_of("primers").unwrap().to_string());
    let mut stats = Stats {
        on_target: 0,
        off_target: 0,
        total_pairs: 0,
        matched: 0,
        mated: 0,
    };
}
/*
    let trim: usize = match args.value_of("trim") {
        Some(s) => s.parse::<usize>().unwrap_or(0),
        None => 0,
    };
*/
/*
 match matches.subcommand() {
     ("test", Some(test_args)) => {
         for readpair in reads {
             match (readpair.p1, readpair.p2) {
                 (Some(_), Some(_)) => stats.matched += 2,
                 (Some(_), None) => stats.matched += 1,
                 (None, Some(_)) => stats.matched += 1,
             }
         }

         eprintln!("{:?},{:?}", stats.matched, stats.total_pairs);
     }
     ("match", Some(match_args)) => {
         for readpair in reads {
             match (p1, p2) {
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
     }

     ("assemble", Some(assemble_args)) => {
         eprintln!("assemble reads");
         exit(1);
     }
 }
*/
