[![Coverage Status](https://coveralls.io/repos/github/jeff-k/amplicontig/badge.svg?branch=main)](https://coveralls.io/github/jeff-k/amplicontig?branch=main)

# amplicontig

Identify paired reads that match primer sequences and assemble tiled amplicons.

This approach is valid for non-fragmented sequencing protocols.

### Synopsis
```
amplicontig 0.1.2

USAGE:
    amplicontig [FLAGS] [OPTIONS] [ARGS]

FLAGS:
    -x               excise primer sequence from reads
    -h, --help       Prints help information
    -n               invert selection (show unmatching reads)
    -s               print stats about primer performance
    -V, --version    Prints version information

OPTIONS:
    -t <trim>        trim bases from 3' end

ARGS:
    <primers>    csv file for primer set
    <R1>         first reads in pairs
    <R2>         second reads in pairs (reversed)
```

#### Primer spec

example:
```
name,forward,primer,position
nCoV-2019_1_LEFT,true,ACCAACCAACTTTCGATCTCTTGT,30
nCoV-2019_1_RIGHT,false,CATCTTTAAGATGTTGACGTGCCTC,230
```

### Installation

```
cargo build --release
```

### Examples

`target/release/amplicontig -s artic-v3.csv ERR4659819_1.fastq.gz ERR4659819_2.fastq.gz`

## Pipeline Description

### Primer identification

The 5' ends of individual reads are used to identify primer sequences.

* Exact string matching in O(n)
* Or lowest hamming distance up to threshold (default: 3)

This stage reports the most prevalent inexact matches and preserves unidentified reads.

### Mate polishing

Read pairs that belong to fragments that are shorter than twice the read length (eg. 500bp) will overlap at the 3' ends. These can be merged into single reads and sequencing errors are resolved as `N`s.

```
ACGTGTGTC->
   <-TCTCACGTCG
      |
ACGTGTNTCACGTCG
```

### Amplicon binning

Amplicons are binned and counted, `N`s are removed.
```
    ACGTGTNTCACGTCG
    ACGTGTGTCNCGTCG
    CCCTGGCTCACANCGC


result:
    ACGTGTGTCACGTCG, 2
    CCCTGGCTCACANCGC, 1
```

By default, bins are discarded for

* Fragment sizes < 50% shortest amplicon length
* Bins with fewer than 100 members

### Amplicon merging

Disagreements between forward and reverse fragments are resolved.

Amplicons that overlap according to primer position are assembled into contigs.

## Output

Fasta of contigs.

GFA support is planned.
