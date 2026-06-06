# Byte-identity / benchmark test

`byte_identity_test.sh` runs grfmite-rs in `--rle-cigar` mode against stock GRF on the bundled
512 kbp fragments and reports both **byte-identity** and **timing**.

## Requirements (install yourself — GRF is not bundled)

On `PATH`:
- `grf-main`, `grf-filter` — from GRF / bioconda `genericrepeatfinder` (tested against **1.0.2**)
- `grf_rs` — this package

## Run

```shell
bash test/byte_identity_test.sh            # single-thread
THREADS=8 bash test/byte_identity_test.sh  # both tools at -t 8
```

Exit status is non-zero if any fragment differs. `--rle-cigar` emits the stock run-length cigar,
so a correct build is byte-identical to `grf-main` detection + `grf-filter`. (The DEFAULT mode
emits the integer TIR arm length instead — equivalent for TIR-Learner, not byte-identical.)

## Fragments

`fragments/*.fa` — 512 kbp prefixes of three public NCBI contigs (NC_081189.1, NC_081190.1,
NC_081234.1), chosen to span typical and repeat-denser composition.
