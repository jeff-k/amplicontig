extern crate csv;
extern crate fst;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

use fst::automaton::Levenshtein;

fn read_primers() {
    let _l = Levenshtein::new("foo", 3);
    //    let mut readbins = HashMap::new();
    let mut plen = 100;
    let mut pmax = 0;

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

    for record in primer_recs {
        let primer = String::from(&record.primer[..plen]);
        let name = String::from(&record.name);
        primers.insert(primer, record);
        //        readbins.insert(name, BTreeMap::new());
        //        primers.insert(String::from(dna::revcomp(primer)), record);
    }
}

#[derive(Debug, Deserialize, Clone)]
pub struct PrimerSet {
    name: String,
    primer: String,
    forward: bool,
    index: u32,
}

impl PrimerSet {
    fn new() -> PrimerSet {
        PrimerSet {
            name: "asdf",
            primer: "asdf",
            forward: false,
            index: 0,
        }
    }
}

pub struct MatchedReads<'a> {
    zipfq: Box<dyn Iterator<Item = (Result<Record, Error>, Result<Record, Error>)>>,
    primers: &'a HashMap<String, PrimerSet>,
}

impl<'a> Iterator for MatchedReads<'a> {
    type Item = (Record, Option<&'a PrimerSet>, Record, Option<&'a PrimerSet>);

    fn next(&mut self) -> Option<Self::Item> {
        match self.zipfq.next() {
            Some((rec1, rec2)) => match (rec1, rec2) {
                (Ok(r1), Ok(r2)) => {
                    let plen = 22;
                    let p1 = String::from_utf8(r1.seq()[..plen].to_vec()).unwrap();
                    let p2 = String::from_utf8(r2.seq()[..plen].to_vec()).unwrap();

                    Some((r1, self.primers.get(&p1), r2, self.primers.get(&p2)))
                }
                _ => {
                    eprintln!("Mismatched fastq files (different lengths)");
                    exit(1);
                }
            },
            None => None,
        }
    }
}

impl<'a> MatchedReads<'a> {
    fn new(
        zipfq: Box<dyn Iterator<Item = (Result<Record, Error>, Result<Record, Error>)>>,
        primers: &HashMap<String, PrimerSet>,
    ) -> MatchedReads {
        MatchedReads { zipfq, primers }
    }
}

pub struct Stats {
    on_target: u32,
    off_target: u32,
    total_pairs: u32,
    matched: u32,
    mated: u32,
}
