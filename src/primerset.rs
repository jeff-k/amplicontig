extern crate csv;

use std::process::exit;

use serde::Deserialize;

use bio::alphabets::dna;
use bio::io::fastq::Error;

use mating::Record;

#[derive(Debug, Deserialize, Clone)]
pub struct PrimerRow {
    name: String,
    seq: String,
    left: bool,
    forward: bool,
    index: u32,
    length: u8,
}

pub struct PrimerRow {
    seq: String,
    left: bool,
    index: u32,
    length: u8,
}

pub struct Amplicon {
    name: String,
    left: Primer,
    right: Primer,
}


#[derive(Debug, Deserialize, Clone)]
pub struct PrimerSet {
    name: String,
    forward: BTree<String, Primer>,
    reverse: BTree<String, Primer>,
}

impl PrimerSet {
    pub fn from_csv(src: Reader) -> PrimerSet {
        let mut readbins = HashMap::new();
        let plen = 100;
        let pmax = 0;

        let mut primer_recs = Vec::new();

        for result in csv::Reader::from_reader(src).deserialize() {
            let record: PrimerSet = result.unwrap();
            //let l = record.primer.len();
            //plen = min(l, plen);
            //pmax = max(l, pmax);
            primer_recs.push(record);
        }

        //for record in primer_recs {
        //    let primer = String::from(&record.primer[..plen]);
        //    let name = String::from(&record.name);
        //primers.insert(primer, record);
        //        readbins.insert(name, BTreeMap::new());
        //        primers.insert(String::from(dna::revcomp(primer)), record);
        //}
        PrimerSet {
            name: "asdf".to_string(),
        }
    }
    fn get(self: Self, p: &str) -> Option<&Primer> {
        //        Some(Primer {
        //    primer: "asdf".to_string(),
        //    forward: true,
        //    index: 0,
        //name: "asdf".to_string() })
        None
    }
}

pub struct MatchedReads<'a> {
    zipfq: Box<dyn Iterator<Item = (Result<Record, Error>, Result<Record, Error>)>>,
    primers: &'a PrimerSet,
}

impl<'a> Iterator for MatchedReads<'a> {
    type Item = (Record, Option<&'a Primer>, Record, Option<&'a Primer>);

    fn next(&mut self) -> Option<Self::Item> {
        match self.zipfq.next() {
            Some((rec1, rec2)) => match (rec1, rec2) {
                (Ok(r1), Ok(r2)) => {
                    let plen = 22;
                    let p1 = String::from_utf8(r1.seq()[..plen].to_vec()).unwrap();
                    let p2 = String::from_utf8(r2.seq()[..plen].to_vec()).unwrap();

                    //Some((r1, self.primers.get(&p1), r2, self.primers.get(&p2)))
                    Some((r1, None, r2, None))
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
    pub fn new(
        zipfq: Box<dyn Iterator<Item = (Result<Record, Error>, Result<Record, Error>)>>,
        primers: &PrimerSet,
    ) -> MatchedReads {
        MatchedReads { zipfq, primers }
    }

    pub fn from_fqs(fq1fn: String, fq2fn: String, primers: &PrimerSet) -> MatchedReads {
        let fq1 = Reader::new(MultiGzDecoder::new(BufReader::new(
            File::open(fq1fn).unwrap(),
        )))
        .records();
        let fq2 = Reader::new(MultiGzDecoder::new(BufReader::new(
            File::open(fq2fn).unwrap(),
        )))
        .records();

        MatchedReads::new(Box::new(fq1.zip(fq2)), primers)
    }
}

pub struct Stats {
    pub on_target: u32,
    pub off_target: u32,
    pub total_pairs: u32,
    pub matched: u32,
    pub mated: u32,
}
