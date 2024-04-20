use core::cmp::{max, min};
use std::fs::File;
use std::path::PathBuf;

use serde::Deserialize;

use std::collections::HashMap;

use bio_seq::prelude::*;

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct Primer {
    pub name: String,
    pub forward: bool,
    pub seq: String, //Vec<u8>,
    //    left: bool,
    pub index: usize,
    //    length: u8,
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PrimerSet {
    name: String,
    plen: usize,
    forward: HashMap<Seq<Dna>, Primer>,
    reverse: HashMap<Seq<Dna>, Primer>,
}

impl PrimerSet {
    pub fn from_csv(src: &PathBuf) -> PrimerSet {
        let mut rdr = csv::Reader::from_reader(File::open(src).unwrap());
        let mut forward = HashMap::new();
        let mut reverse = HashMap::new();
        let mut pmax = 0;
        let mut plen = 100;

        let mut primers = Vec::new();

        for result in rdr.deserialize() {
            let rx: Primer = result.unwrap();
            let index = if rx.forward {
                rx.index
            } else {
                rx.index + rx.seq.len()
            };
            let record = Primer {
                name: rx.name,
                forward: rx.forward,
                seq: rx.seq,
                index,
            };
            let l = record.seq.len();
            plen = min(l, plen);
            pmax = max(l, pmax);
            primers.push(record);
        }

        for primer in primers {
            if primer.forward {
                forward.insert(primer.seq[..plen].try_into().unwrap(), primer.clone());
            } else {
                reverse.insert(primer.seq[..plen].try_into().unwrap(), primer.clone());
            }
        }
        PrimerSet {
            name: "asdf".to_string(),
            plen,
            forward,
            reverse,
        }
    }
    pub fn get(&self, p: &SeqSlice<Dna>) -> Option<&Primer> {
        match self.forward.get(&p[..self.plen]) {
            Some(p) => Some(p),
            None => match self.reverse.get(&p[..self.plen]) {
                Some(p) => Some(p),
                None => None,
            },
        }
    }
}

#[allow(dead_code)]
pub struct Stats {
    pub on_target: u32,
    pub off_target: u32,
    pub total_pairs: u32,
    pub matched: u32,
    pub mated: u32,
}

#[allow(dead_code)]
impl Stats {
    pub fn new() -> Self {
        Stats {
            on_target: 0,
            off_target: 0,
            total_pairs: 0,
            matched: 0,
            mated: 0,
        }
    }
}
