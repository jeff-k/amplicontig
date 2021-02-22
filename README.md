# extract-amplicons
Extract Illumina read pairs that exactly match primer sequences

### Synopsis
```
extract-amplicons 0.1.1

USAGE:
    extract-amplicons [FLAGS] [OPTIONS] [ARGS]

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

### Building

```
cargo build --release
```

### Examples

`target/release/extract-amplicons -s artic-v3.csv ERR4659819_1.fastq.gz ERR4659819_2.fastq.gz`
