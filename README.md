# grfmite-rs

Rust port and acceleration of the GRFMite subprogram of
[Generic Repeat Finder](https://github.com/bioinfolabmu/GenericRepeatFinder), developed for
[TIR-Learner v4](https://github.com/KGerhardt/TIR-Learner).

TIR-Learner uses a genome fragmentation approach
([genomeSplitter](https://github.com/KGerhardt/genomeSplitter)) to produce ~5 Mbp sequence
chunks. Some chunks are pathological to repeat-finder programs like GRF: low-complexity sequence
causes GRF to emit a very large number of possible repeat matches, few of which are biologically
interesting targets. In centromeric and other satellite regions, the worst case involves chunks
that run hundreds to thousands of times slower than typical 5 Mbp fragments. Especially in
high-thread-count deployments, this dramatically harms TIR-Learner's performance, as a few
pathological chunks consume the vast majority of total runtime while most threads lie idle, 
having finished processing all non-pathological genome chunks in a fraction of the time.

This Rust port of GRFMite includes multiple algorithmic improvements that dramatically improve
performance while maintaining identical results (see
[Changes from stock GRF](docs/CHANGES_FROM_C.md)), plus a Rust parallel work-stealing model that
reduces the impact of pathological chunks by supplying additional worker threads to those chunks
when they appear.

**Scope:** only the MITE path (`grf-main -c 1`) — the one GRF mode TIR-Learner invokes — is
implemented. This is **not** a full GRF port (no `-c 0` / `-c 2`, no other GRF binaries); a
complete GRF overhaul is intended for a future release. Byte-identical to stock GRF in
`--rle-cigar` mode and **40–55× faster single-threaded** (see Benchmarks below).

## Benchmarks

Single-threaded, grfmite-rs vs stock GRF (`genericrepeatfinder` 1.0.2), on the bundled 512 kbp
fragments. grfmite-rs runs in `--rle-cigar` mode so its output is **byte-identical** to stock
GRF detection + `grf-filter`. Reproduce with `test/byte_identity_test.sh` (you install GRF).

| Fragment (512 kbp) | candidates | stock GRF | grfmite-rs | speedup | byte-identical |
|---|---:|---:|---:|---:|:---:|
| NC_081189.1 | 6,125 | 30.7 s | 0.7 s | 43.9× | ✓ |
| NC_081190.1 | 6,583 | 40.6 s | 0.9 s | 45.1× | ✓ |
| NC_081234.1 | 12,684 | 92.5 s | 1.7 s | 54.4× | ✓ |

The speedup grows with candidate density: stock GRF's cost rises super-linearly on repeat-dense
sequence, exactly where grfmite-rs's linear candidate finding wins most. (Single-thread, WSL2 dev
workstation; absolute times are machine-dependent — the ratio and byte-identity are the claims.)

### At scale (whole genome, parallel)

Measured on a ~2 Gbp genome (411 × ~5 Mbp chunks):

- Full ~5 Mbp chunk: **~10–15 s** single-thread (vs hundreds of seconds for stock GRF on the
  densest chunks).
- All 411 chunks, 10-way per-file parallel: **~98% parallel efficiency** — vs ~25% for the
  stock per-file pipeline, whose few repeat-dense stragglers strand cores.
- Worst straggler chunk, with intra-chunk work-stealing (`-t 20`): **95 s → 14 s**, so a single
  dense chunk no longer floors the wall on many-core machines.

See `docs/CHANGES_FROM_C.md` for the full list of algorithmic changes vs stock GRF.

## Install

Conda (once published):

```shell
conda install -c conda-forge -c bioconda grfmite-rs
```

From source (requires the Rust toolchain):

```shell
cargo build --release          # -> target/release/grf_rs
```

The binary is named `grf_rs` (the conda package also installs a `grfmite-rs` alias).

## Usage

Single FASTA (writes `<outdir>/candidate.fasta`):

```shell
grf_rs -i genome.fa -o outdir -c 1 -p 20 --min_space 10 --max_space 5000 \
  --max_indel 0 --min_tr 10 --min_spacer_len 10 --max_spacer_len 5000 -t 8
```

Batch — one work-stealing pool over many chunks (writes `<outdir>/<file_stem>/candidate.fasta`):

```shell
grf_rs --batch list_of_fasta_paths.txt -o outdir -t 16 -c 1 -p 20 \
  --min_space 10 --max_space 5000 --max_indel 0 --min_tr 10 \
  --min_spacer_len 10 --max_spacer_len 5000
```

### Output modes

- **default** — the cigar field of each candidate header is the integer TIR arm length
  (what TIR-Learner consumes). TIR-Learner-equivalent to stock GRF.
- **`--rle-cigar`** — emits the stock GRF run-length `m`/`M` cigar string, byte-identical to
  upstream `grf-main` (for general / non-TIR-Learner use).

Header format (1-based, inclusive): `>seqid:start:stop:<cigar|tir_len>:tsd`.

## License

GPL-3.0-or-later (see `LICENSE`). This is a clean-room reimplementation of GRF's MITE
detection; it produces equivalent output but shares no source with upstream GRF.
