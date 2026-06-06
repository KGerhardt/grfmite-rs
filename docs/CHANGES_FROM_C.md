# grfmite-rs vs stock GRF — algorithmic change log

This documents how grfmite-rs differs from stock GRF's MITE detection (`grf-main -c 1`,
tested against `genericrepeatfinder` 1.0.2). It is a **clean-room reimplementation**: it
produces equivalent results (byte-identical via `--rle-cigar`, see `test/byte_identity_test.sh`)
but shares no source with upstream GRF. Scope is the `-c 1` MITE path only.

Each entry: **Description** (what GRF does) · **Rationale** (why change it) ·
**Alteration** (what grfmite-rs does) · **Result** (effect + equivalence).

---

## Algorithmic

### 1. Seed-pair detection — skew `abs_sum` proxy + byte compare → exact 2-bit code test
- **Description.** GRF builds per-position (A−T, C−G) skew cumulative sums (`S_TR`), uses
  `|Δa| + |Δc| ≤ 2·seed_mismatch` as a necessary filter, then verifies each survivor by
  comparing bases pairwise.
- **Rationale.** The skew threshold passes ~20% of all positions — it is a weak filter; the
  real cost (~80% of runtime) was the byte verification, which chases two seq regions up to
  ~5 kb apart (cache-hostile).
- **Alteration.** Precompute a 2-bit code per seed window — forward `fcode`, reverse-complement
  `gcode` (complement = XOR `0b11`). A seed match is exact and branch-light: field 0 equal AND
  popcount of differing 2-bit fields ≤ `seed_mismatch`. This replaces *both* the skew proxy and
  the byte loop, and both operands stream sequentially.
- **Result.** Identical seed set; the dominant cost removed. Byte-identical output.

### 2. Seed scan — brute-force O(n × spacer_range) → code-keyed join (CSR index)
- **Description.** For each spacer length GRF scans every position, i.e. ~`n × (max_space−min_space)`
  ≈ 2.5×10⁹ pair evaluations on a 512 kbp window.
- **Rationale.** The overwhelming majority of those pairs are non-matches; enumerating them is wasted.
- **Alteration.** Build a CSR index mapping each 2-bit seed code → its ascending start positions.
  For each `end`, look up the `1 + 3·(seed−1)` neighbor codes of `gcode[end]` (exact + each inner
  field flipped) within the spacer-distance window. Iterate `end` ascending so per-spacer output
  stays start-ascending.
- **Result.** Same candidates, same order; cost drops from O(n·range) to ~O(n·matches).

### 3. `seqComplexity` (LZ76) — O(n²) `find()` → O(n) suffix automaton
- **Description.** GRF's complexity filter computes an LZ76 factorization via repeated substring
  search, O(n²) per candidate.
- **Rationale.** The O(n²) scan dominated the low-complexity (repeat-dense) tail.
- **Alteration.** Online suffix automaton tracking each state's first-occurrence end position;
  a factor extends while the current substring occurred earlier and commits otherwise.
- **Result.** Identical complexity value, O(n). (LZ76 ≠ LZ77 — see Correctness parity.)

### 4. `seqComplexity` — early-exit on an integer threshold
- **Description.** The caller only needs the boolean `complexity < 0.675`.
- **Rationale.** The factor count `c` only increases as the scan proceeds, so once it crosses the
  threshold the boolean is decided.
- **Alteration.** Return as soon as `c ≥ ceil(0.675 · denom)` (integer comparison; `c` is integral
  so the ceiling is exact).
- **Result.** Same boolean; the rest of the scan is skipped in the common (passing) case.

## Representation / data structures

### 5. SAM transitions — per-state hash maps → flat array + buffer reuse
- **Description.** A suffix automaton's transitions are naturally per-state maps.
- **Rationale.** After (3)/(4), per-candidate heap allocation + hashing + map-cloning of states
  became the dominant `seqComplexity` cost (the SAM is rebuilt per candidate).
- **Alteration.** Flat `trans[state·alpha + code]` table keyed by a dense byte→code map built once
  over the genome's actual alphabet (injective ⇒ identical SAM); the SAM buffers are reused across
  candidates via thread-local scratch.
- **Result.** Identical automaton, zero per-candidate allocation/hashing.

### 6. `getStem` extension — byte `is_pair` → 2-bit XOR pairing
- **Description.** Stem/arm extension compares bases by ASCII value.
- **Rationale.** Per-comparison cost over every candidate's full arms.
- **Alteration.** Per-base 2-bit encoding (ACGT→0..3, non-ACGT→4 sentinel); a pair is
  `enc[i] ^ enc[j] == 0b11`. The sentinel can never XOR to `0b11`, so non-ACGT never pairs.
- **Result.** Identical cigar, fewer instructions per base.

## Output

### 7. Cigar field — run-length `m`/`M` string → integer TIR arm length (default)
- **Description.** GRF writes an RLE cigar (e.g. `8m1M1m`) into each candidate header.
- **Rationale.** Downstream (TIR-Learner) uses only the **sum** of the cigar's run lengths
  (= the TIR arm length) and re-aligns the TIRs itself; the `m`/`M` structure is discarded.
  Building and parsing the string is pure overhead.
- **Alteration.** Default mode emits the integer arm length directly in the cigar field;
  `--rle-cigar` restores the exact stock RLE string. The emitted integer is `min(pos[idx], len2)`
  — the clamped arm length, matching what summing the stock cigar yields.
- **Result.** Default = downstream-equivalent (same candidates; the field carries the same value).
  `--rle-cigar` = byte-identical to stock GRF.

### 8. Length filter (`grf-filter`) — external sidecar → inlined
- **Description.** Stock `grf-main`'s `filterByLen` shells out to a separate `grf-filter` binary
  (TR-length and spacer-length bounds). If `grf-filter` is not colocated, the `system()` call
  fails **silently** and `candidate.fasta` is left unfiltered.
- **Rationale.** A fragile cross-process dependency for a trivial length check, with a silent
  failure mode.
- **Alteration.** Apply the TR-length / spacer-length filter inline in the output step.
- **Result.** Identical filtering; no sidecar process, no silent-fail mode.

### 9. Homopolymer / dinucleotide filter — regex → hand-rolled scan
- **Description.** The low-complexity pre-filter is the regex `(\w)\1{7,} | (\w\w)\2{3,}`.
- **Rationale.** Regex-engine overhead per candidate for a fixed pattern.
- **Alteration.** Direct scan: 8+ identical bases in a row, or 4+ consecutive copies of a 2-mer.
- **Result.** Identical decision, no regex engine.

## Parallelism

### 10. Scheduling — OpenMP static parallel-for → rayon work-stealing + global batch
- **Description.** GRF parallelizes the spacer loop with default (static) OpenMP scheduling, one
  process per chunk.
- **Rationale.** Static scheduling strands threads when heavy spacer-lengths cluster (satellites),
  and per-file parallelism floors the wall at the single longest chunk on many-core machines.
- **Alteration.** rayon work-stealing over end-range chunks within a contig (intra-file `-t N`),
  plus a `--batch` mode running many chunks through **one** global work-stealing pool — bulk runs
  inter-chunk, idle workers steal a dense chunk's sub-ranges at the tail. Output order preserved
  via in-order `collect`.
- **Result.** ~98% per-file efficiency on a 2 Gbp genome; dense stragglers absorbed (95 s → 14 s
  at 20 threads on the worst chunk); byte-identical at any thread count.

## Correctness parity (deliberately matched, not changed)

- **`getStem` substr clamp.** `compress(left.substr(0, pos[i]))` in C++ relies on `substr`
  clamping the count to the string length (`pos[i]` can exceed the trimmed arm). grfmite-rs uses
  `.min(left.len())` to reproduce this exactly.
- **LZ76 ≠ LZ77.** GRF's factorization includes the breaking character in the factor and does not
  count a trailing partial factor. The suffix automaton replicates *that* definition — it is not
  the LPF/LZ77 a suffix array would give (e.g. `ABAB` → 2, not 3).
- **Coordinates.** 1-based, inclusive `[start, stop]` in the header, as GRF emits.
