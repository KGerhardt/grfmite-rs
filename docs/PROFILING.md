# Profiling: where the single-thread wall lives (2026-06-05)

Sample: `~/tirlearner_test/grf_test/sample_512k.fa` (512 kbp, 8854 candidates).
Method: env-gated `GRF_PROBE` early-`continue` buckets + coarse region timers,
`std::hint::black_box` to defeat dead-code elimination. All instrumentation has since
been stripped; `probe=0` (default) was re-verified byte-identical to REF throughout.

## Region breakdown (single-thread, ~27s total `detect`)

| Region | Marginal cost | Notes |
|---|---|---|
| fasta read / `s_tr` build / output | ~0.02s | negligible — 100% of wall is in `detect_str` |
| `abs_sum` arithmetic (2.5e9 iters) | ~1–2.5s | the part the inverted-index-on-`S_TR` plan targeted |
| threshold branch + `is_pair(start,end)` | ~8.5s | runs on the ~520M abs-passers |
| unpair loop (up to 9× `is_pair`) | ~13.5s | **single biggest chunk** |
| `detect_tsd` | ~0.5s | |
| `get_stem` | ~0.9s | |
| `filter_low_complex` (gc/regex/`seq_complexity`) | ~2.6s | |

(WSL timing is noisy run-to-run; treat ±1–2s as noise. Structure is stable.)

## Conclusion — the queued plan targeted the wrong region

- **~22s of ~27s is seed verification** (`is_pair` + unpair loop), and it sits
  *downstream* of the `abs_sum` filter.
- `abs_sum ≤ 2` is a **weak filter**: ~20% of all positions (~520M) pass it, so the
  real selectivity (520M → ~295k true candidates) comes from the seed byte-comparisons,
  not the arithmetic.
- Therefore the **inverted-index-on-`S_TR`** step (old HANDOFF step 2) would emit exactly
  the 520M abs-passers and we'd still run `is_pair`+unpair on every one — net win ~1–2s,
  **not** the 22s. It is NOT the next move as written.

## Corrected next step — attack seed verification

Order still algorithmic-before-low-level:

- **A (algorithmic, lead):** index seed windows by their exact **2-bit seed code**
  (20 bits for seed=10) and look up reverse-complement partners directly (±1 mismatch for
  `seed_mismatch=1`). Emits ~only true seed-matches, attacking the 520M→295k *volume*
  instead of re-deriving the abs-pass set. The inverted index is not dead — it changes
  key from `S_TR` to the 2-bit seed code and does real selectivity work.
- **B (low-level, after A):** 2-bit pack `seq` once; seed check = packed-XOR +
  masked-popcount `≤ 1`, collapsing the 9-iteration unpair loop. Attacks per-comparison
  cost but still touches every abs-passer.

Both must stay diff-validated byte-identical against `REF_candidate.fasta`.
