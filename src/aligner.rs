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
                Some(pos) => {
                    //println!("{} {}", pos, query.len());
                    //Reverse(*pos - query.len())
                    Reverse(*pos)
                }
                None => Unmapped,
            },
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Cigar {
    Match,
    Subs,
    Ins,
    Del,
}

pub fn edit_dist(a: &str, b: &str) -> (usize, Vec<Cigar>) {
    let a_len = a.len();
    let b_len = b.len();
    let mut dp = vec![vec![0; b_len + 1]; a_len + 1];
    let mut ops = vec![vec![Vec::new(); b_len + 1]; a_len + 1];

    for i in 0..=a_len {
        dp[i][0] = i;
        ops[i][0] = vec![Cigar::Del; i];
    }

    for j in 0..=b_len {
        dp[0][j] = j;
        ops[0][j] = vec![Cigar::Ins; j];
    }

    for i in 1..=a_len {
        for j in 1..=b_len {
            let base_a = a.chars().nth(i - 1).unwrap();
            let base_b = b.chars().nth(j - 1).unwrap();

            let ins_cost = dp[i][j - 1] + 1;
            let del_cost = dp[i - 1][j] + 1;
            let subs_cost = dp[i - 1][j - 1] + if base_a == base_b { 0 } else { 1 };

            let (penalty, best_op) = if subs_cost <= ins_cost && subs_cost <= del_cost {
                (
                    subs_cost,
                    if base_a == base_b {
                        Cigar::Match
                    } else {
                        Cigar::Subs
                    },
                )
            } else if ins_cost <= del_cost {
                (ins_cost, Cigar::Ins)
            } else {
                (del_cost, Cigar::Del)
            };

            dp[i][j] = penalty;

            ops[i][j] = {
                let prev_i = i - (best_op != Cigar::Ins) as usize;
                let prev_j = j - (best_op != Cigar::Del) as usize;

                let mut new_ops = ops[prev_i][prev_j].clone();
                new_ops.push(best_op);
                new_ops
            };
        }
    }
    (dp[a_len][b_len], ops[a_len][b_len].clone())
}

#[cfg(test)]
mod tests {
    use super::{edit_dist, Cigar};

    #[test]
    fn test_matches() {
        let r1 = "kittens";
        let r2 = "kittens";
        let (score, ops) = edit_dist(r1, r2);
        println!("{}\t{:?}", score, ops);
        assert_eq!(ops, vec![Cigar::Match; r1.len()]);
    }

    #[test]
    fn test_mm() {
        let r1 = "kittens";
        let r2 = "mitteny";
        let (score, ops) = edit_dist(r1, r2);
        println!("{}\t{:?}", score, ops);
        assert_eq!(
            ops,
            vec![
                Cigar::Subs,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
                Cigar::Subs
            ]
        );
    }
    #[test]
    fn test_ins() {
        let r1 = "ttttt";
        let r2 = "tttttttt";
        let (score, ops) = edit_dist(r1, r2);
        println!("{}\t{:?}", score, ops);
        assert_eq!(
            ops,
            vec![
                Cigar::Ins,
                Cigar::Ins,
                Cigar::Ins,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
            ]
        );
    }

    #[test]
    fn test_dels() {
        let r1 = "tttttttt";
        let r2 = "tttttt";
        let (score, ops) = edit_dist(r1, r2);
        println!("{}\t{:?}", score, ops);
        assert_eq!(
            ops,
            vec![
                Cigar::Del,
                Cigar::Del,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
            ]
        );
    }
    #[test]
    fn test_interior_ins() {
        let r1 = "cccaaa";
        let r2 = "ccctttaaa";
        let (score, ops) = edit_dist(r1, r2);
        println!("{}\t{:?}", score, ops);
        assert_eq!(
            ops,
            vec![
                Cigar::Match,
                Cigar::Match,
                Cigar::Match,
                Cigar::Ins,
                Cigar::Ins,
                Cigar::Ins,
                Cigar::Match,
                Cigar::Match,
                Cigar::Match
            ]
        );
    }
}
