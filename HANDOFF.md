# GRF → Rust port — session handoff (2026-06-05)

## Goal
Reimplement GRF `grf-main -c 1` (MITE detection) in Rust, **byte-identical** to the real
TIR-Learner GRF output, then optimize. Order of work: **correctness → single-thread
algorithmic → single-thread low-level (SIMD) → parallel**. (User principle: take the
HIGH-LEVEL algorithmic wins before low-level/SIMD; low-level only pays after that.)
Broader context in memory: `grf-rust-liftover-plan.md`, `cnn-phase-profiling-plan.md`.

## Current status
- **Port is byte-identical** to stock-detection + real `grf-filter` on 3 contigs
  (512 kbp samples): 8854 / 6682 / 6289 candidates. Now under git (`main` = safe baseline).
- `seqComplexity` = hand-rolled **suffix automaton** LZ76 (O(n), exact same count).
- **Perf (512 kbp sample, single-thread): ~27 s → ~4.3 s (~6.4×)**, vs C++ patched `-t1`
  **20.7 s** (so now ~4.8× faster than C++). All steps diff-validated byte-identical.
  Optimization journey (see `docs/PROFILING.md` for the supporting profile):
  1. **2-bit exact seed code** (commit `0800278`): the seed match was a lossy `abs_sum`
     skew proxy (~20% false-pass) over all 2.5e9 pairs + a byte two-pointer re-verify
     (~22 s). Replaced with precomputed `fcode`/`gcode` and an XOR + popcount test.
     `abs_sum` dropped entirely. ~27 s → ~11.4 s.
  2. **Code-keyed join** (commit `45a6085`): the O(1) test still visited all 2.5e9 pairs
     (~6.7 s). Replaced with a CSR index (fcode → ascending start list); per `end` look up
     the 1+3·(seed−1) neighbor codes of `gcode[end]` within the spacer window. Visits
     ~matches, not all pairs. ~11.4 s → ~4.3 s.
- **Next target:** `filter_low_complex` (gc / regex / `seqComplexity`) — now the dominant
  remaining cost (~3 s) since the scan collapsed. Then SIMD/parallel/PyO3 as before.
- Earlier micro-opts (still in place, were ~2%): homopolymer filter == C++ regex, integer
  `c_thresh` early-exit in the SAM, `target-cpu=native`, `getStem` reused thread-local.

## Key locations
- Rust port:            `~/tirlearner_test/grf_rs/`  (all logic in `src/main.rs`)
- Stock C++ (patched, same-results, 1.73× faster than upstream): `~/tirlearner_test/grf_src/`
  built binary: `~/tirlearner_test/grf_src/src/grf-main/grf-main`
- Iteration sample:     `~/tirlearner_test/grf_test/sample_512k.fa` (512 kbp, ~8854 cands)
- **Canonical reference**: `~/tirlearner_test/grf_test/REF_candidate.fasta`
  (= stock detection + REAL conda `grf-filter`; see GOTCHA below)
- More chunks for broader validation: `~/tirlearner_test/before/out/split_genome/*.fasta`
- Real grf-filter (what TIR-Learner uses): `/mnt/c/.../conda_envs/tirlearner/bin/grf-filter`
- C++ source read-reference: `/tmp/grf_inspect/src/grf-main/` (EPHEMERAL — persistent copy is `grf_src`)

## Locked params (TIR-Learner's `-c 1` invocation)
`seed=10, seed_mismatch=1, min_stem(--min_tr)=10, max_stem=INT_MAX, percent(-p)=20,`
`min_space=10, max_space=5000, min_tsd=2, max_tsd=10, max_indel=0, spacer(len3)∈[10,5000]`

## TESTING COMMANDS

```bash
DST=~/tirlearner_test
ARGS="-c 1 -p 20 --min_space 10 --max_space 5000 --max_indel 0 --min_tr 10 --min_spacer_len 10 --max_spacer_len 5000"
GRF_RS=$DST/grf_rs/target/release/grf_rs
GRF_STOCK=$DST/grf_src/src/grf-main/grf-main
GRFFILTER=/mnt/c/Users/kenji/Desktop/conda_envs/tirlearner/bin/grf-filter

# 1) BUILD the rust port
cd $DST/grf_rs && cargo build --release

# 2) RUN + VALIDATE (must be byte-identical to canonical REF)
od=$DST/grf_test/run; rm -rf $od; mkdir -p $od
$GRF_RS -i $DST/grf_test/sample_512k.fa -o $od $ARGS
diff -q $od/candidate.fasta $DST/grf_test/REF_candidate.fasta && echo IDENTICAL || \
  diff <(grep '^>' $od/candidate.fasta) <(grep '^>' $DST/grf_test/REF_candidate.fasta) | head

# 3) TIME (single-thread)
( time $GRF_RS -i $DST/grf_test/sample_512k.fa -o $od $ARGS >/dev/null 2>&1 ) 2>&1 | grep real

# 4) MAKE A 512 kbp SAMPLE from any chunk (handles wrapped lines)
sample() { python3 - "$1" "$2" <<'PY'
import sys; src,out=sys.argv[1],sys.argv[2]; N=512000; seq=[]; hdr=">s"
for line in open(src):
    if line.startswith('>'): hdr=line.rstrip()
    else: seq.append(line.strip())
s=''.join(seq)[:N]; open(out,'w').write(hdr+"\n"+"\n".join(s[i:i+80] for i in range(0,len(s),80))+"\n")
PY
}

# 5) BUILD a canonical reference for a NEW sample (stock detection + REAL grf-filter)
#    (stock grf-main's own filterByLen SILENTLY FAILS w/o the grf-filter sidecar -> filter manually)
mkref() { local s=$1 r=$2; rm -rf $r; mkdir -p $r
  $GRF_STOCK -i $s -o $r -t 4 $ARGS >/dev/null 2>&1
  $GRFFILTER 10 2147483647 10 5000 $r/candidate.fasta $r/filtered.fasta >/dev/null 2>&1; }

# 6) VALIDATE rust against a fresh contig end-to-end
validate() { local chunk=$1 tag=$2
  sample "$chunk" $DST/grf_test/s_$tag.fa
  mkref $DST/grf_test/s_$tag.fa $DST/grf_test/ref_$tag
  local u=$DST/grf_test/rs_$tag; rm -rf $u; mkdir -p $u
  $GRF_RS -i $DST/grf_test/s_$tag.fa -o $u $ARGS >/dev/null 2>&1
  diff -q $u/candidate.fasta $DST/grf_test/ref_$tag/filtered.fasta >/dev/null \
    && echo "[$tag] IDENTICAL ✓ ($(grep -c '^>' $u/candidate.fasta))" || echo "[$tag] DIFFERS"; }
# e.g.: validate $DST/before/out/split_genome/long_chunk_NC_081192.1_offset_0.fasta NC192

# QUICK COST PROBE (stub a function, build, time, REVERT): e.g. put `return 1.0;` at top of
# seq_complexity() -> isolates SAM cost (~4.6 s). Always revert + re-diff after.
```

## NEXT STEPS (algorithmic first)
DONE: instrument (localized ~22 s to seed verification, not abs_sum), 2-bit exact code
test (commit `0800278`), code-keyed join (commit `45a6085`). See `docs/PROFILING.md`.
The old "inverted index on S_TR" plan was WRONG — it targeted the ~1–2 s abs_sum
arithmetic, but the wall was the byte seed-verification downstream of it; the correct key
is the **2-bit seed code**, not S_TR. That is now built (the join).

1. **`filter_low_complex` (now the dominant ~3 s).** Runs on candidate survivors: gc_content,
   homopolymer/dinuc regex, and `seqComplexity` (suffix-automaton LZ76). Re-profile to split
   these three, then attack — likely the SAM build (per-candidate `HashMap` allocs) is the
   bulk; array transitions (5-symbol ACGTN) + buffer reuse should reclaim most.
2. **`get_stem` 2-bit** (~0.8 s): complement = bitwise-NOT, isPair ⟺ `enc(a)^enc(b)==0b11`;
   reuse the already-built `fcode`-style encoding instead of byte `is_pair`.
3. Only THEN low-level: SoA/SIMD on whatever hot loop remains.
4. Only THEN parallel: rayon over the spacer loop (was reverted — single-thread first),
   then PyO3 + fragment-aware orchestrator (dedup falls out structurally).

OPEN QUESTION (scale): the join only touches *matches*, so `abs_sum` as a pre-filter looks
moot. The real scale risk is giant `fcode` buckets in repetitive genomes inflating the
match count / downstream cost — watch at full-genome scale (per-bucket-size distribution).

## GOTCHAS / learnings (don't re-discover)
- **grf-filter sidecar**: stock `grf-main`'s `filterByLen` shells out to a separate
  `grf-filter` binary; if it's not next to grf-main the `system()` call fails SILENTLY and
  candidate.fasta is left UNFILTERED. Canonical REF must be filtered with the conda grf-filter.
- **getStem `compress(left.substr(0,pos[i]))`**: C++ `substr` CLAMPS the count to string
  length (pos[i] can exceed trimmed `left` from trailing-M positions) → rust needs `.min(left.len())`.
- **LZ76 ≠ LZ77**: the stock `seqComplexity` includes the breaking char in the factor and
  doesn't count a trailing partial factor → it's NOT the LPF/LZ77 a suffix array gives
  (e.g. "ABAB": code=2, LZ77=3). Use a suffix automaton (firstpos), not libsais/LPF.
- **`general-sam` crate** doesn't expose per-state first-occurrence (firstpos) → hand-rolled SAM.
- **seqComplexity early-exit**: caller only uses the `< 0.675` boolean and `c` only grows →
  return once `c >= ceil(0.675*denom)` (integer compare; c is integral so ceil is exact).
- Rust bounds checks come from index `s[i]`; elide via iterators (`zip`, proven in-bounds),
  not `unsafe`, unless a hot loop needs `get_unchecked` with a safety argument.
- `.cargo/config.toml` sets `target-cpu=native` (helped ~nothing here — scan is AoS, not vectorizing).
