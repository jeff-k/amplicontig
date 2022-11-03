use crate::mating::rc;
use std::collections::HashMap;

use crate::aligner::Alignment::{Forward, Reverse, Unmapped};

pub enum Alignment {
    Forward(usize),
    Reverse(usize),
    Unmapped,
}

pub struct Aligner {
    k: usize,
    forward: HashMap<Vec<u8>, usize>,
    reverse: HashMap<Vec<u8>, usize>,
}

impl Aligner {
    pub fn new(reference: &[u8]) -> Self {
        let k = 21;
        let mut forward = HashMap::new();
        let mut reverse = HashMap::new();
        let mut seen: HashMap<Vec<u8>, bool> = HashMap::new();

        for (index, window) in reference.windows(k).enumerate() {
            if !seen.contains_key(window) {
                forward.insert(window.to_vec(), index);
                seen.insert(window.to_vec(), true);
            }
        }
        let ref_rc = rc(reference);
        for (index, window) in ref_rc.windows(k).enumerate() {
            if !seen.contains_key(window) {
                reverse.insert(window.to_vec(), ref_rc.len() - index);
                seen.insert(window.to_vec(), true);
            }
        }
        Aligner {
            k,
            forward,
            reverse,
        }
    }

    pub fn get(&self, query: &[u8]) -> Alignment {
        match self.forward.get(&query[0..self.k]) {
            Some(pos) => Forward(*pos),
            None => match self.reverse.get(&query[0..self.k]) {
                Some(pos) => Reverse(*pos - query.len()),
                None => Unmapped,
            },
        }
    }
}
