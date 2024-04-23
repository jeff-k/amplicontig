use core::cmp::{max, min, Eq, Ordering};
use std::fs::File;
use std::path::PathBuf;

use serde::Deserialize;

use std::collections::HashMap;

use crate::mating::{mate, merge};
use bio_seq::prelude::*;

pub enum Orientation {
    F1R2,
    F2R1,
    R1F2,
    R2F1,
}

use Orientation::{F1R2, F2R1, R1F2, R2F1};

pub enum Amplicon {
    Discarded,
    Merged(Orientation, usize, usize, Seq<Dna>),
    Paired(Orientation, usize, usize),
}

use Amplicon::{Discarded, Merged, Paired};

#[derive(Debug, Deserialize, Clone, PartialEq, Eq, Hash)]
pub struct Primer {
    pub name: String,
    pub forward: bool,
    pub seq: String,
    //    left: bool,
    pub index: usize,
    //    length: u8,
}

#[derive(Debug, Clone)]
pub struct PrimerSet {
    plen: usize,
    forward: HashMap<Seq<Dna>, Primer>,
    reverse: HashMap<Seq<Dna>, Primer>,
}

#[inline]
fn merge_amplicon(p1: &Primer, r1: &SeqSlice<Dna>, p2: &Primer, r2: &SeqSlice<Dna>) -> Amplicon {
    let start = min(p1.index, p2.index);
    let end = max(p1.index, p2.index);
    let max_indel = 84;
    let hint = ((end - start) / 2) - 1;
    match (p1.index.cmp(&p2.index), p1.forward, p2.forward) {
        (Ordering::Less, true, false) => {
            // F1R2
            let r2rc = r2.revcomp();
            if hint < r1.len() - 30 {
                match mate(r1, &r2rc, hint, max_indel) {
                    Some(seam) => Merged(F1R2, start, end, merge(r1, &r2rc, seam)),
                    //None => Paired(F1R2, p1.index, r1.into(), p2.index, r2rc),
                    None => Paired(F1R2, start, end),
                }
            } else {
                Paired(F1R2, start, end)
            }
        }
        (Ordering::Less, false, true) => {
            // R1F2
            let r1rc = r1.revcomp();
            if hint < r1.len() - 30 {
                match mate(&r1rc, r2, hint, max_indel) {
                    Some(seam) => Merged(R1F2, start, end, merge(&r1rc, r2, seam)),
                    None => Paired(R1F2, start, end),
                }
            } else {
                Paired(R1F2, start, end)
            }
        }
        (Ordering::Greater, true, false) => {
            // F2R1
            let r2rc = r2.revcomp();
            match mate(&r2rc, r1, hint, max_indel) {
                Some(seam) => Merged(F2R1, start, end, merge(&r2rc, r1, seam)),
                None => Paired(F2R1, start, end),
            }
        }
        (Ordering::Greater, false, true) => {
            // R2F1
            let r1rc = r1.revcomp();
            match mate(r2, &r1rc, hint, max_indel) {
                Some(seam) => {
                    //                    println!("\tmerged: seam: {}, {}, {}, {}", seam, start, end, hint);
                    Merged(R2F1, start, end, merge(r2, &r1rc, seam))
                }
                None => {
                    //                    println!("\tpaired:\t{}\t{}\t{}\t{}", hint, end - start, start, end);
                    Paired(R2F1, start, end)
                }
            }
        }
        _ => Discarded,
    }
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
    pub fn get_amplicon(&self, r1: &SeqSlice<Dna>, r2: &SeqSlice<Dna>) -> Amplicon {
        match (self.get(r1), self.get(r2)) {
            (Some(p1), Some(p2)) => merge_amplicon(p1, r1, p2, r2),
            //                            *bins.entry((p1.name.clone(), p2.name.clone())).or_insert(1) += 1;
            _ => Amplicon::Discarded,
        }
    }
}

#[derive(Default)]
pub struct Stats {
    pub on_target: u32,
    pub off_target: u32,
    pub total_pairs: u32,
    pub matched: u32,
    pub mated: u32,
}
