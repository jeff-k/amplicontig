use bio_seq::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::aligner::Alignment::{Forward, Reverse, Unmapped};

pub enum Alignment {
    Forward(usize),
    Reverse(usize),
    Unmapped,
}

pub struct Aligner {
    k: usize,
    forward: HashMap<Seq<Dna>, usize>,
    reverse: HashMap<Seq<Dna>, usize>,
}

impl Aligner {
    pub fn new(reference: &SeqSlice<Dna>) -> Self {
        let k = 21;
        let mut forward = HashMap::new();
        let mut reverse = HashMap::new();
        let mut seen: HashSet<Seq<Dna>> = HashSet::new();

        for (index, window) in reference.windows(k).enumerate() {
            if !seen.contains(window) {
                forward.insert(window.into(), index);
                seen.insert(window.into());
            } else {
                forward.remove(window);
            }
        }
        for (index, window) in reference.revcomp().windows(k).enumerate() {
            if !seen.contains(window) {
                reverse.insert(window.into(), reference.len() - index);
                seen.insert(window.into());
            } else {
                forward.remove(window);
                reverse.remove(window);
            }
        }

        println!("{} forward, {} reverse", forward.len(), reverse.len());
        Aligner {
            k,
            forward,
            reverse,
        }
    }

    pub fn get(&self, query: &SeqSlice<Dna>) -> Alignment {
        match self.forward.get(&query[0..self.k]) {
            Some(pos) => Forward(*pos),
            None => match self.reverse.get(&query[0..self.k]) {
                Some(pos) => Reverse(*pos - query.len()),
                None => Unmapped,
            },
        }
    }

    /*
    pub fn align(&self, query: &[u8] ) -> Alignment {
        for i in range(0..query.len()) {
            let g1 = self.get(query);
            let g2 = self.get(query);
        }
    }
    */
}
