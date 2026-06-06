#!/bin/bash
# Byte-identity + timing test: grfmite-rs (--rle-cigar) vs stock GRF, on the bundled fragments.
#
# Requires on PATH (install yourself; we do NOT bundle GRF):
#   grf-main, grf-filter   from GRF / bioconda `genericrepeatfinder` (tested against 1.0.2)
#   grf_rs                 this package
#
# --rle-cigar mode emits the stock run-length m/M cigar, so output is byte-identical to
# stock `grf-main` detection + `grf-filter`. (The DEFAULT mode instead emits the integer
# TIR arm length in the cigar field — TIR-Learner-equivalent, not byte-identical.)
#
# Env: THREADS (default 1).
set -uo pipefail

ARGS="-c 1 -p 20 --min_space 10 --max_space 5000 --max_indel 0 --min_tr 10 --min_spacer_len 10 --max_spacer_len 5000"
FRAGDIR="$(cd "$(dirname "$0")" && pwd)/fragments"
THREADS="${THREADS:-1}"

for b in grf-main grf-filter grf_rs; do
  command -v "$b" >/dev/null 2>&1 || { echo "ERROR: '$b' not on PATH"; exit 1; }
done

now() { date +%s.%N; }
fail=0
printf "%-22s %12s %12s %9s %9s %s\n" "fragment" "grf-main" "grfmite-rs" "speedup" "cands" "byte-identical"
for f in "$FRAGDIR"/*.fa; do
  name=$(basename "$f")
  ref=$(mktemp -d); rs=$(mktemp -d)
  # stock GRF reference: detection + explicit grf-filter (robust to grf-filter colocation)
  t0=$(now); grf-main -i "$f" -o "$ref" -t "$THREADS" $ARGS >/dev/null 2>&1
  grf-filter 10 2147483647 10 5000 "$ref/candidate.fasta" "$ref/filtered.fasta" >/dev/null 2>&1
  t1=$(now); gt=$(awk "BEGIN{printf \"%.1f\", $t1-$t0}")
  # grfmite-rs, byte-identical mode
  t0=$(now); grf_rs -i "$f" -o "$rs" -t "$THREADS" $ARGS --rle-cigar >/dev/null 2>&1
  t1=$(now); rt=$(awk "BEGIN{printf \"%.1f\", $t1-$t0}")
  spd=$(awk "BEGIN{ if ($rt>0) printf \"%.1fx\", $gt/$rt; else printf \"-\" }")
  n=$(grep -c '^>' "$rs/candidate.fasta" 2>/dev/null || echo 0)
  if diff -q "$rs/candidate.fasta" "$ref/filtered.fasta" >/dev/null 2>&1; then id="IDENTICAL"; else id="DIFFERS"; fail=1; fi
  printf "%-22s %11ss %11ss %9s %9s %s\n" "$name" "$gt" "$rt" "$spd" "$n" "$id"
  rm -rf "$ref" "$rs"
done
echo
if [ "$fail" -eq 0 ]; then echo "All fragments byte-identical to stock GRF."; else echo "MISMATCH detected."; fi
exit $fail
