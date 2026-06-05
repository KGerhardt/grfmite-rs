# Profiling & single-thread optimization journey (2026-06-05)

Sample: `~/tirlearner_test/grf_test/sample_512k.fa` (512 kbp, 8854 candidates).
Method: env-gated `GRF_PROBE` early-`continue` buckets + coarse region timers,
`std::hint::black_box` to defeat dead-code elimination. Instrumentation is added on a
throwaway branch, measured, then discarded; `probe=0` (default) is re-verified
byte-identical to REF each time. WSL timing is noisy (±1–2 s); structure is stable.

## Headline

| Stage | Time | vs original | Commit |
|---|---|---|---|
| Original baseline (abs_sum scan) | ~27–28 s | 1× | `127e078` |
| 2-bit exact seed code test | ~11.4 s | ~2.4× | `0800278` |
| Code-keyed join (CSR index) | ~4.3 s | ~6.4× | `45a6085` |
| (C++ patched `-t1` reference) | 20.7 s | — | — |

All steps byte-identical to `REF_candidate.fasta`.

## Step 0 — original profile (what localized the work)

`detect_str` was 100% of the wall (~27 s); fasta read / `s_tr` build / output were noise.
Within it:

| Region | Marginal |
|---|---|
| `abs_sum` arithmetic (2.5e9 iters) | ~1–2.5 s |
| threshold branch + `is_pair(start,end)` (~520M abs-passers) | ~8.5 s |
| unpair loop (≤9× `is_pair`, ~154M) | ~13.5 s |
| `detect_tsd` | ~0.5 s |
| `get_stem` | ~0.9 s |
| `filter_low_complex` (gc/regex/`seq_complexity`) | ~2.6 s |

Counts: 2.54e9 pairs, **520M (20.5%) pass `abs_sum ≤ 2`**, 154M pass the boundary pair,
295k pass the full seed check.

**Conclusion:** ~22 s was seed verification (`is_pair` + unpair), *downstream* of a weak
`abs_sum` filter. `abs_sum` is only a necessary proxy for the exact seed RC-match. This
overturned the old "inverted index on `S_TR`" plan (it targeted the ~1–2 s arithmetic).

## Step 1 — 2-bit exact code test (drop the proxy)

Encode bases 2 bits each (A=00 C=01 G=10 T=11; complement = XOR `0b11`). Precompute once:
`fcode[start]` = forward seed-window code, `gcode[end]` = reverse-complement-window code.
Seed match ⟺ field 0 of `fcode^gcode` is 0 AND ≤`seed_mismatch` fields differ
(popcount over the low bit of each 2-bit field). Exact, so `abs_sum` is dropped entirely;
both operands stream sequentially (no random `seq[]` chasing). The ~22 s seed-verification
block collapses. New breakdown of ~11.4 s:

| Region | Cost |
|---|---|
| F/G precompute | 0.06 s |
| code scan (2.5e9 pairs, O(1) each) | 6.7 s |
| `detect_tsd` | 0.6 s |
| `get_stem` | 0.8 s |
| `filter_low_complex` | 3.3 s |

The hybrid (keep `abs_sum` as a pre-filter in front of the code test) was rejected by
reasoning: the code test is already as cheap as `abs_sum`, so pre-filtering only adds
`s_tr` memory traffic. (Low-complexity-chunk value of `abs_sum` is a separate scale test.)

## Step 2 — code-keyed join (stop visiting non-matches)

The O(1) test still visited all 2.5e9 pairs (6.7 s). Build a CSR index `fcode value →
ascending start list`. For each `end` (ascending), enumerate the `1 + 3·(seed−1)` neighbor
codes of `gcode[end]` (itself + each inner field flipped to its 3 alternatives), look each
up, and keep starts in the spacer window `[end−max_off2, end−min_off2]`. Iterating `end`
ascending makes each spacer bucket fill in start-ascending order → byte-identical output.
Visits ~matches instead of all pairs. ~11.4 s → ~4.3 s.

## Next target

With the scan collapsed, `filter_low_complex` (gc / regex / `seq_complexity` suffix
automaton) is the dominant remaining cost (~3 s) and runs on candidate survivors. Re-profile
to split its three parts before attacking. Then `get_stem` (~0.8 s), then SIMD/parallel.

Open scale question: the join only touches matches, so `abs_sum` pre-filtering looks moot;
the real risk is giant `fcode` buckets in repetitive genomes inflating the match count.
