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
            } else {
                forward.remove(window);
            }
        }
        let ref_rc = rc(reference);
        for (index, window) in ref_rc.windows(k).enumerate() {
            if !seen.contains_key(window) {
                reverse.insert(window.to_vec(), ref_rc.len() - index);
                seen.insert(window.to_vec(), true);
            } else {
                forward.remove(window);
                reverse.remove(window);
            }
        }

        //println!("{} forward, {} reverse", forward.len(), reverse.len());
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

    /*
    pub fn align(&self, query: &[u8] ) -> Alignment {
        for i in range(0..query.len()) {
            let g1 = self.get(query);
            let g2 = self.get(query);
        }
    }
    */
}
