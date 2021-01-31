# extract-amplicons
Extract Illumina read pairs that exactly match primer sequences

```
$ extract-amplicons --help
extract-amplicons 0.1.0

USAGE:
    extract-amplicons [FLAGS] [ARGS]

FLAGS:
    -h, --help       Prints help information
    -n               invert selection (show unmatching reads)
    -V, --version    Prints version information

ARGS:
    <primers>    csv file for primer set
    <R1>         first reads in pairs
    <R2>         second reads in pairs (reversed)
```
