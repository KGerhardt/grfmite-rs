// Direct Rust port of GRF grf-main `-c 1` (MITE) detection, aiming for byte-identical
// candidate.fasta vs stock grf-main. Correctness first (byte-wise); 2-bit/SIMD later.
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Clone)]
pub struct Param {
    pub seed: i32,
    pub seed_mismatch: i32,
    pub min_stem: i32, // = --min_tr
    pub max_stem: i32,
    pub percent: i32, // = -p
    pub min_space: i32,
    pub max_space: i32,
    pub min_tsd: i32,
    pub max_tsd: i32,
    pub max_mismatch: i32,
    pub min_spacer_len: i32,
    pub max_spacer_len: i32,
    pub emit_cigar: bool, // --rle-cigar: emit stock GRF RLE cigar (byte-identical, general-tool);
                          // default false = emit the integer TIR arm length (TIR-Learner use-case)
    pub threads: usize,   // -t: intra-file parallelism over the end-loop (1 = serial)
    pub legacy_fasta: bool, // --legacy-fasta: emit the old candidate.fasta (header + full
                            // subsequence) for drop-in GRF compat. default false = emit
                            // candidate.json (coords only; TIR-Learner reslices seq/TSD by coord)
}
impl Default for Param {
    fn default() -> Self {
        Param { seed: 10, seed_mismatch: 1, min_stem: 0, max_stem: i32::MAX, percent: 10,
                min_space: -1, max_space: -1, min_tsd: 2, max_tsd: 10,
                max_mismatch: i32::MAX, min_spacer_len: 0, max_spacer_len: i32::MAX,
                emit_cigar: false, threads: 1, legacy_fasta: false }
    }
}

// arm = clamped TIR arm length (= what TIR-Learner's parse_cig recovers); cigar = the stock
// RLE m/M string, present only in --rle-cigar mode (None in the default integer mode).
struct Mite { start: usize, end: usize, arm: u32, tsd: String, cigar: Option<String> }

// 2-bit base code with complement == XOR 0b11 (A=00 C=01 G=10 T=11); None = non-ACGT.
#[inline]
fn code(b: u8) -> Option<u32> {
    match b { b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None }
}
const INVALID: u32 = u32::MAX; // sentinel: seed window contains a non-ACGT base

// run-length encode an m/M byte string -> e.g. "3m2M"
fn compress(s: &[u8]) -> String {
    if s.is_empty() { return String::new(); }
    let mut prev = s[0];
    let mut num: u32 = 0;
    let mut r = String::new();
    for &ch in s {
        if ch == prev { num += 1; }
        else { r.push_str(&num.to_string()); r.push(prev as char); prev = ch; num = 1; }
    }
    r.push_str(&num.to_string());
    r.push(prev as char);
    r
}

fn gc_content(s: &[u8]) -> f32 {
    let mut r = 0f32;
    for &i in s { if i == b'G' || i == b'C' { r += 1.0; } }
    r / s.len() as f32
}

// homopolymer / dinuc filter == regex (\w)\1{7,} | (\w\w)\2{3,}  (DNA: all word chars)
fn low_complex_regex(s: &[u8]) -> bool {
    let n = s.len();
    let mut run = 1usize; // (\w)\1{7,} : 8+ identical in a row
    for k in 1..n {
        if s[k] == s[k - 1] { run += 1; if run >= 8 { return true; } } else { run = 1; }
    }
    if n >= 8 { // (\w\w)\2{3,} : 4+ consecutive copies of a 2-mer
        for p in 0..=n - 8 {
            if s[p] == s[p + 2] && s[p] == s[p + 4] && s[p] == s[p + 6]
                && s[p + 1] == s[p + 3] && s[p + 1] == s[p + 5] && s[p + 1] == s[p + 7] {
                return true;
            }
        }
    }
    false
}

const SNONE: u32 = u32::MAX; // suffix-automaton "no transition" sentinel

// LZ76 complexity via a suffix automaton (O(n)), same result as the stock O(n^2) find().
// Build the SAM of seq tracking firstpos (first-occurrence END pos) per state; a factor
// extends while seq[send..=i] occurs earlier (firstpos < i) and commits otherwise.
// Transitions are a flat [state*alpha + code] table (a2c maps a byte to a dense code over
// the genome's actual alphabet) instead of per-state HashMaps, and the SAM buffers are
// reused across candidates via thread-local scratch — same SAM, no per-candidate alloc/hash.
fn seq_complexity(seq: &[u8], a2c: &[u8; 256], alpha: usize) -> f64 {
    let n = seq.len();
    let denom = n as f64 / ((n as f64).ln() / 4f64.ln());
    let c_thresh = (0.675 * denom).ceil() as u32;
    thread_local! {
        // (len, link, firstpos, trans)
        static SAM: std::cell::RefCell<(Vec<i32>, Vec<i32>, Vec<i32>, Vec<u32>)> =
            std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new()));
    }
    SAM.with(|cell| {
        let mut b = cell.borrow_mut();
        let (len, link, firstpos, trans) = &mut *b;
        len.clear(); link.clear(); firstpos.clear(); trans.clear();
        // state 0 = init. Invariant: after creating state s, trans.len() == (s+1)*alpha.
        len.push(0); link.push(-1); firstpos.push(-1);
        trans.resize(alpha, SNONE);
        let mut last = 0usize;
        for (pos, &ch) in seq.iter().enumerate() {
            let c = a2c[ch as usize] as usize;
            let cur = len.len();
            len.push(len[last] + 1);
            link.push(-1);
            firstpos.push(pos as i32);
            trans.resize((cur + 1) * alpha, SNONE);
            let mut p = last as i32;
            while p != -1 && trans[p as usize * alpha + c] == SNONE {
                trans[p as usize * alpha + c] = cur as u32;
                p = link[p as usize];
            }
            if p == -1 {
                link[cur] = 0;
            } else {
                let q = trans[p as usize * alpha + c] as usize;
                if len[p as usize] + 1 == len[q] {
                    link[cur] = q as i32;
                } else {
                    let clone = len.len();
                    len.push(len[p as usize] + 1);
                    link.push(link[q]);
                    firstpos.push(firstpos[q]); // clone keeps q's first occurrence
                    trans.extend_from_within(q * alpha..q * alpha + alpha); // copy q's edges
                    while p != -1 && trans[p as usize * alpha + c] == q as u32 {
                        trans[p as usize * alpha + c] = clone as u32;
                        p = link[p as usize];
                    }
                    link[q] = clone as i32;
                    link[cur] = clone as i32;
                }
            }
            last = cur;
        }
        // ---- LZ76 factorization ----
        let mut c: u32 = 1;
        let mut state = 0usize; // init state
        for i in 1..n {
            let t = trans[state * alpha + a2c[seq[i] as usize] as usize];
            let ext = if t != SNONE && (firstpos[t as usize] as usize) < i { Some(t as usize) } else { None };
            match ext {
                Some(nxt) => state = nxt,
                None => {
                    c += 1;
                    if c >= c_thresh { return c as f64 / denom; }
                    state = 0;
                }
            }
        }
        c as f64 / denom
    })
}

fn detect_tsd(seq: &[u8], start: usize, end: usize, min: i32, max: i32) -> String {
    let mut i = max;
    while i >= min {
        let iu = i as usize;
        if start < iu || end + iu > seq.len() - 1 { i -= 1; continue; }
        if seq[start - iu..start] == seq[end + 1..end + 1 + iu] {
            return String::from_utf8_lossy(&seq[start - iu..start]).into_owned();
        }
        i -= 1;
    }
    String::new()
}

// Default fast path: no-indel stem extension returning ONLY the clamped TIR arm length
// (= what stock GRF's cigar would sum to via TIR-Learner's parse_cig). No m/M buffer, no
// cigar string. `e` is the candidate's 2-bit encoding (ACGT->0..3, non-ACGT->4); pairing is
// e[i]^e[j]==0b11 (complement = XOR 0b11; the 4-sentinel can never XOR to 0b11). None = no arm.
fn stem_arm(e: &[u8], p: &Param) -> Option<u32> {
    thread_local! { static POS: std::cell::RefCell<Vec<i32>> = std::cell::RefCell::new(Vec::new()); }
    POS.with(|pc| {
        let pos = &mut *pc.borrow_mut();
        pos.clear();
        let len = e.len();
        let mut error = 0i32;
        let mut last_m = -1i32; // largest i that paired (== rposition of 'm')
        let (mut i, mut j) = (0usize, len - 1);
        while i < j {
            if e[i] ^ e[j] == 0b11 { last_m = i as i32; }
            else {
                error += 1;
                if error > p.max_mismatch { break; }
                pos.push(i as i32);
            }
            i += 1;
            j -= 1;
        }
        let len2 = last_m + 1; // index after last 'm', or 0
        pos.push(len2);
        for idx in (0..pos.len()).rev() {
            if pos[idx] >= p.min_stem && (idx as i32) * 100 <= p.percent * pos[idx] {
                return Some(pos[idx].min(len2) as u32); // clamp == C++ substr clamp == parse_cig value
            }
        }
        None
    })
}

// --rle-cigar path: same arm length, plus the stock GRF RLE m/M cigar (byte-identical output).
fn stem_cigar(e: &[u8], p: &Param) -> Option<(u32, String)> {
    thread_local! {
        static SCRATCH: std::cell::RefCell<(Vec<u8>, Vec<i32>)> =
            std::cell::RefCell::new((Vec::new(), Vec::new()));
    }
    SCRATCH.with(|sc| {
        let mut sc = sc.borrow_mut();
        let (left, pos) = &mut *sc;
        left.clear();
        pos.clear();
        let len = e.len();
        let mut error = 0i32;
        let (mut i, mut j) = (0usize, len - 1);
        while i < j {
            if e[i] ^ e[j] == 0b11 { left.push(b'm'); }
            else {
                error += 1;
                if error > p.max_mismatch { break; }
                left.push(b'M');
                pos.push(i as i32);
            }
            i += 1;
            j -= 1;
        }
        let len2 = left.iter().rposition(|&c| c == b'm').map(|x| x + 1).unwrap_or(0);
        left.truncate(len2);
        pos.push(len2 as i32);
        for idx in (0..pos.len()).rev() {
            if pos[idx] >= p.min_stem && (idx as i32) * 100 <= p.percent * pos[idx] {
                let arm = (pos[idx] as usize).min(left.len()); // C++ substr clamps
                return Some((arm as u32, compress(&left[0..arm])));
            }
        }
        None
    })
}

// arm == both TR lengths (no-indel => len1 == len2 == arm == old get_stem_len(cigar,_)).
fn filter_low_complex(seq: &[u8], arm: usize, a2c: &[u8; 256], alpha: usize) -> bool {
    let l = &seq[0..arm];
    let r = &seq[seq.len() - arm..];
    let gc1 = gc_content(l) as f64;
    let gc2 = gc_content(r) as f64;
    if gc1 < 0.2 || gc1 > 0.8 || gc2 < 0.2 || gc2 > 0.8 { return false; }
    if low_complex_regex(l) || low_complex_regex(r) { return false; }
    if seq_complexity(seq, a2c, alpha) < 0.675 { return false; }
    true
}

// Seed match via a code-keyed join instead of scanning all O(n*spacer_range) pairs.
// fcode[start] = forward window code; gcode[end] = RC-window code; a seed match (byte-
// identical to is_pair k=0 + <=seed_mismatch inner) means fcode[start] equals gcode[end]
// EXCEPT for <= seed_mismatch of the inner fields (field 0 must match exactly).
//
// For seed_mismatch==1 the set of fcode values matching a given gcode[end] is exactly:
//   { gcode[end] } U { gcode[end] with one inner field (1..seed) set to a different base }
// = 1 + 3*(seed-1) neighbor codes. We look each up in the CSR index (code -> ascending
// start list) and keep starts in the spacer window [end-max_off2, end-min_off2].
//
// Iterating `end` ascending makes each spacer bucket fill in ascending `start` order
// (start = end - off2 for that spacer), so output stays byte-identical to the old scan.
//
// Process one contiguous end-range [end_lo, end_hi) -> flat Vec<Mite> (end-ascending).
// Pure read of the immutable index/codes, so ranges run independently (rayon).
#[allow(clippy::too_many_arguments)]
fn detect_range(end_lo: usize, end_hi: usize, gcode: &[u32], offsets: &[u32],
                starts_by_code: &[u32], seq: &[u8], enc: &[u8], p: &Param,
                a2c: &[u8; 256], alpha: usize) -> Vec<Mite> {
    let seed = p.seed as usize;
    let min_off2 = (p.min_space + 2 * p.seed - 1) as usize; // end - start at min spacer
    let max_off2 = (p.max_space + 2 * p.seed - 1) as usize; // end - start at max spacer
    let mut out: Vec<Mite> = Vec::new();
    let mut codes: Vec<u32> = Vec::with_capacity(1 + 3 * (seed - 1));
    for end in end_lo.max(min_off2)..end_hi {
        let g = gcode[end];
        if g == INVALID { continue; }
        let hi = end - min_off2; // start <= hi  (i_sp >= min_space)
        let lo = end.saturating_sub(max_off2); // start >= lo  (i_sp <= max_space)
        // neighbor fcode values: g itself + each inner field flipped to its 3 alternatives
        codes.clear();
        codes.push(g);
        for k in 1..seed {
            let sh = 2 * k;
            let cur = (g >> sh) & 3;
            let base = g & !(3u32 << sh);
            for v in 0..4u32 { if v != cur { codes.push(base | (v << sh)); } }
        }
        for &c in &codes {
            let bucket = &starts_by_code[offsets[c as usize] as usize..offsets[c as usize + 1] as usize];
            let si = bucket.partition_point(|&x| (x as usize) < lo);
            for &st in &bucket[si..] {
                let start = st as usize;
                if start > hi { break; }
                let l = end - start + 1; // = i_sp + 2*seed
                let cand = &seq[start..start + l];
                let tsd = detect_tsd(seq, start, end, p.min_tsd, p.max_tsd);
                if tsd.is_empty() { continue; }
                let ecand = &enc[start..start + l]; // max_indel == 0
                let stem = if p.emit_cigar {
                    stem_cigar(ecand, p).map(|(a, c)| (a, Some(c)))
                } else {
                    stem_arm(ecand, p).map(|a| (a, None))
                };
                if let Some((arm, cigar)) = stem {
                    if filter_low_complex(cand, arm as usize, a2c, alpha) {
                        out.push(Mite { start, end, arm, tsd, cigar });
                    }
                }
            }
        }
    }
    out
}

// Driver: run detect_range over the contig, serial or rayon over fixed end-chunks.
// In-order `collect` keeps the flattened Mites end-ascending, so distributing them into
// per-spacer bins reproduces the serial start-ascending order -> byte-identical either way.
#[allow(clippy::too_many_arguments)]
fn join_detect(gcode: &[u32], offsets: &[u32], starts_by_code: &[u32],
               seq: &[u8], enc: &[u8], p: &Param, spacer_vecs: &mut [Vec<Mite>],
               a2c: &[u8; 256], alpha: usize) {
    let n = seq.len();
    let min_off2 = (p.min_space + 2 * p.seed - 1) as usize;
    let min_space = p.min_space as usize;
    let twoseed_m1 = (2 * p.seed - 1) as usize;
    let total = n.saturating_sub(min_off2);
    const CHUNK: usize = 16384; // ends per task; fine enough to spread dense regions
    let mites: Vec<Mite> = if p.threads <= 1 || total <= CHUNK {
        detect_range(min_off2, n, gcode, offsets, starts_by_code, seq, enc, p, a2c, alpha)
    } else {
        let nchunks = total.div_ceil(CHUNK);
        let parts: Vec<Vec<Mite>> = (0..nchunks).into_par_iter().map(|ci| {
            let lo = min_off2 + ci * CHUNK;
            let hi = (lo + CHUNK).min(n);
            detect_range(lo, hi, gcode, offsets, starts_by_code, seq, enc, p, a2c, alpha)
        }).collect();
        parts.into_iter().flatten().collect()
    };
    for m in mites {
        let i_sp = m.end - m.start - twoseed_m1;
        spacer_vecs[i_sp - min_space].push(m);
    }
}

// Process one input fasta end-to-end and write its candidate.fasta to out_path.
// In --batch mode this is the unit fed to the outer rayon pool; join_detect's inner
// chunk-parallelism nests in the same pool, so idle workers steal a dense file's chunks
// (global work-stealing: bulk = inter-file, tail = intra-file).
pub fn run_file(in_path: &str, out_path: &str, p: &Param) {
    // read fasta (uppercased); chrom name = first whitespace token of the header
    let mut chroms: Vec<String> = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut reader = parse_fastx_file(in_path).expect("open fasta");
    while let Some(rec) = reader.next() {
        let rec = rec.expect("record");
        let id = rec.id();
        let name = id.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(b"");
        chroms.push(String::from_utf8_lossy(name).into_owned());
        seqs.push(rec.seq().iter().map(|b| b.to_ascii_uppercase()).collect());
    }

    // Dense byte->code map over this file's actual alphabet (injective on every byte that
    // appears -> seq_complexity's array transitions are byte-identical to the HashMap version).
    let mut a2c = [0u8; 256];
    let mut alpha = 0usize;
    {
        let mut present = [false; 256];
        for s in &seqs { for &bb in s { present[bb as usize] = true; } }
        for bb in 0..256 { if present[bb] { a2c[bb] = alpha as u8; alpha += 1; } }
    }
    let alpha = alpha.max(1);

    let nspace = (p.max_space - p.min_space + 1) as usize;
    let mut candidates: Vec<Vec<Vec<Mite>>> = Vec::with_capacity(chroms.len());
    for seq in &seqs {
        let n = seq.len();
        let seed = p.seed as usize;
        let mut spacer_vecs: Vec<Vec<Mite>> = (0..nspace).map(|_| Vec::new()).collect();
        if n >= (2 * p.seed + p.max_space + 2 * p.max_tsd) as usize {
            // fcode[p] = forward 2-bit code of window [p, p+seed); INVALID if it has a non-ACGT.
            let mut fcode: Vec<u32> = vec![INVALID; n - seed + 1];
            for p_ in 0..=n - seed {
                let mut val = 0u32;
                let mut ok = true;
                for k in 0..seed {
                    match code(seq[p_ + k]) { Some(c) => val |= c << (2 * k), None => { ok = false; break; } }
                }
                if ok { fcode[p_] = val; }
            }
            // enc[pos] = per-base 2-bit code (ACGT->0..3) or 4 for non-ACGT; used by get_stem.
            let enc: Vec<u8> = seq.iter().map(|&b| code(b).map(|c| c as u8).unwrap_or(4)).collect();
            // gcode[end] = RC code of window [end-seed+1, end] = complement of each base, reversed,
            // so field k = comp(code(seq[end-k])). INVALID for end < seed-1 or any non-ACGT.
            let mut gcode: Vec<u32> = vec![INVALID; n];
            for end in (seed - 1)..n {
                let mut val = 0u32;
                let mut ok = true;
                for k in 0..seed {
                    match code(seq[end - k]) { Some(c) => val |= (c ^ 0b11) << (2 * k), None => { ok = false; break; } }
                }
                if ok { gcode[end] = val; }
            }
            // CSR index: fcode value -> ascending list of start positions (seed=10 => 2^20 codes).
            let ncode = 1usize << (2 * seed);
            let mut offsets = vec![0u32; ncode + 1];
            for &c in &fcode { if c != INVALID { offsets[c as usize + 1] += 1; } }
            for i in 0..ncode { offsets[i + 1] += offsets[i]; }
            let mut starts_by_code = vec![0u32; offsets[ncode] as usize];
            let mut cursor = offsets.clone();
            for (start, &c) in fcode.iter().enumerate() {
                if c != INVALID {
                    let cc = c as usize;
                    starts_by_code[cursor[cc] as usize] = start as u32;
                    cursor[cc] += 1;
                }
            }
            join_detect(&gcode, &offsets, &starts_by_code, seq, &enc, &p, &mut spacer_vecs, &a2c, alpha);
        }
        candidates.push(spacer_vecs);
    }

    // inline grf-filter (TR len + spacer len): a candidate is emitted iff this holds.
    let passes = |m: &Mite| {
        let (len1, len2) = (m.arm as i32, m.arm as i32); // no-indel: both TR lengths == arm
        let len3 = (m.end as i32) - (m.start as i32) + 1 - len1 - len2;
        len1 >= p.min_stem && len1 <= p.max_stem
            && len2 >= p.min_stem && len2 <= p.max_stem
            && len3 >= p.min_spacer_len && len3 <= p.max_spacer_len
    };
    let f = File::create(out_path).expect("create out");
    let mut w = BufWriter::new(f);
    if p.legacy_fasta {
        // Legacy candidate.fasta: per-candidate header + full subsequence (drop-in GRF compat).
        for j in 0..nspace {
            for ci in 0..chroms.len() {
                for m in &candidates[ci][j] {
                    if passes(m) {
                        // tir field: stock RLE cigar (--rle-cigar) or the integer arm length (default)
                        match &m.cigar {
                            Some(c) => writeln!(w, ">{}:{}:{}:{}:{}", chroms[ci], m.start + 1, m.end + 1, c, m.tsd).unwrap(),
                            None => writeln!(w, ">{}:{}:{}:{}:{}", chroms[ci], m.start + 1, m.end + 1, m.arm, m.tsd).unwrap(),
                        }
                        writeln!(w, "{}", std::str::from_utf8(&seqs[ci][m.start..=m.end]).unwrap()).unwrap();
                    }
                }
            }
        }
    } else {
        // Default candidate.json: chrom-keyed columnar coordinates, no sequence. All integers:
        // start/end (1-based inclusive), arm (TIR-arm length), tsd (TSD size). TIR-Learner reslices
        // the element sequence and TSD bases from the chunk fasta by these coordinates.
        write!(w, "{{").unwrap();
        let mut first_chrom = true;
        for ci in 0..chroms.len() {
            let mut sel: Vec<&Mite> = Vec::new();
            for j in 0..nspace {
                for m in &candidates[ci][j] {
                    if passes(m) { sel.push(m); }
                }
            }
            if sel.is_empty() { continue; }
            if !first_chrom { w.write_all(b",").unwrap(); }
            first_chrom = false;
            write!(w, "\"{}\":{{\"start\":[", chroms[ci]).unwrap();
            for (i, m) in sel.iter().enumerate() { if i > 0 { w.write_all(b",").unwrap(); } write!(w, "{}", m.start + 1).unwrap(); }
            write!(w, "],\"end\":[").unwrap();
            for (i, m) in sel.iter().enumerate() { if i > 0 { w.write_all(b",").unwrap(); } write!(w, "{}", m.end + 1).unwrap(); }
            write!(w, "],\"arm\":[").unwrap();
            for (i, m) in sel.iter().enumerate() { if i > 0 { w.write_all(b",").unwrap(); } write!(w, "{}", m.arm).unwrap(); }
            write!(w, "],\"tsd\":[").unwrap();
            for (i, m) in sel.iter().enumerate() { if i > 0 { w.write_all(b",").unwrap(); } write!(w, "{}", m.tsd.len()).unwrap(); }
            w.write_all(b"]}").unwrap();
        }
        w.write_all(b"}").unwrap();
    }
}

// Size the global rayon pool to `threads` (call once at startup; subsequent calls are no-ops).
pub fn set_threads(threads: usize) {
    rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build_global().ok();
}

// CLI --batch: outer par_iter over files + join_detect's inner chunk par_iter share one pool
// -> single global work-stealing pool (bulk = inter-file, tail = intra-file). Writes
// One flat file per input: <outdir>/<basename>.json (or <basename>.candidate.fasta for
// --legacy-fasta). No per-chunk directory, so batch output costs N inodes, not 2N — meaningful
// on HPC shared filesystems with file-count quotas. Returns elapsed seconds.
pub fn run_batch(paths: &[String], outdir: &str, p: &Param) -> f64 {
    let t0 = std::time::Instant::now();
    paths.par_iter().for_each(|ip| {
        let base = std::path::Path::new(ip).file_stem()
            .map(|s| s.to_string_lossy().into_owned()).unwrap_or_else(|| "out".into());
        let out_path = if p.legacy_fasta {
            format!("{}/{}.candidate.fasta", outdir, base)
        } else {
            format!("{}/{}.json", outdir, base)
        };
        run_file(ip, &out_path, p);
    });
    t0.elapsed().as_secs_f64()
}

// ---- in-memory fragment path (PyO3 orchestrator core) ----
// One candidate in GLOBAL 1-based coords: (seqid, start, stop, tir_arm_len, tsd_len).
// (used only by the `python` feature; kept compiled in both for parity/testing)
#[allow(dead_code)]
type Record = (String, u64, u64, u32, u32);

// Detect on one in-memory sequence (one contig/fragment), returning final length-filtered
// candidates in GLOBAL coords, with stride-ownership dedup: emit iff the candidate's anchor
// (local start) lands in this fragment's owned prefix [0, owned_len) -> overlap-tail dups are
// never emitted. `offset` = fragment's global 0-based start within the parent contig.
#[allow(dead_code)] // used by the `python` feature (PyO3 orchestrator entry)
fn process_fragment(seqid: &str, seq: &[u8], offset: usize, owned_len: usize, p: &Param) -> Vec<Record> {
    // per-fragment dense alphabet (injective over this seq's bytes -> seq_complexity identical)
    let mut a2c = [0u8; 256];
    let mut alpha = 0usize;
    { let mut present = [false; 256];
      for &bb in seq { present[bb as usize] = true; }
      for bb in 0..256 { if present[bb] { a2c[bb] = alpha as u8; alpha += 1; } } }
    let alpha = alpha.max(1);

    let n = seq.len();
    let seed = p.seed as usize;
    let nspace = (p.max_space - p.min_space + 1) as usize;
    let mut spacer_vecs: Vec<Vec<Mite>> = (0..nspace).map(|_| Vec::new()).collect();
    if n >= (2 * p.seed + p.max_space + 2 * p.max_tsd) as usize {
        let mut fcode: Vec<u32> = vec![INVALID; n - seed + 1];
        for p_ in 0..=n - seed {
            let mut val = 0u32; let mut ok = true;
            for k in 0..seed { match code(seq[p_ + k]) { Some(c) => val |= c << (2 * k), None => { ok = false; break; } } }
            if ok { fcode[p_] = val; }
        }
        let enc: Vec<u8> = seq.iter().map(|&b| code(b).map(|c| c as u8).unwrap_or(4)).collect();
        let mut gcode: Vec<u32> = vec![INVALID; n];
        for end in (seed - 1)..n {
            let mut val = 0u32; let mut ok = true;
            for k in 0..seed { match code(seq[end - k]) { Some(c) => val |= (c ^ 0b11) << (2 * k), None => { ok = false; break; } } }
            if ok { gcode[end] = val; }
        }
        let ncode = 1usize << (2 * seed);
        let mut offsets = vec![0u32; ncode + 1];
        for &c in &fcode { if c != INVALID { offsets[c as usize + 1] += 1; } }
        for i in 0..ncode { offsets[i + 1] += offsets[i]; }
        let mut starts_by_code = vec![0u32; offsets[ncode] as usize];
        let mut cursor = offsets.clone();
        for (start, &c) in fcode.iter().enumerate() {
            if c != INVALID { let cc = c as usize; starts_by_code[cursor[cc] as usize] = start as u32; cursor[cc] += 1; }
        }
        join_detect(&gcode, &offsets, &starts_by_code, seq, &enc, p, &mut spacer_vecs, &a2c, alpha);
    }

    let mut out: Vec<Record> = Vec::new();
    for j in 0..nspace {
        for m in &spacer_vecs[j] {
            let (len1, len2) = (m.arm as i32, m.arm as i32);
            let len3 = (m.end as i32) - (m.start as i32) + 1 - len1 - len2;
            if len1 >= p.min_stem && len1 <= p.max_stem
                && len2 >= p.min_stem && len2 <= p.max_stem
                && len3 >= p.min_spacer_len && len3 <= p.max_spacer_len
                && m.start < owned_len // stride-ownership dedup
            {
                out.push((seqid.to_string(),
                          (offset + m.start + 1) as u64, (offset + m.end + 1) as u64,
                          m.arm, m.tsd.len() as u32));
            }
        }
    }
    out
}

#[cfg(feature = "python")]
mod python {
    #![allow(unsafe_op_in_unsafe_fn)] // pyo3 0.22 macro-generated #[no_mangle] vs edition 2024
    use super::{process_fragment, Param, Record};
    use pyo3::prelude::*;

    /// Run GRF -c 1 (MITE) on a batch of in-memory fragments through one rayon work-stealing
    /// pool, with stride-ownership dedup. fragments = [(seqid, sequence_bytes, global_offset,
    /// owned_len), ...]. Returns [(seqid, global_start, global_stop, tir_len, tsd_len), ...].
    #[pyfunction]
    #[pyo3(signature = (fragments, threads=1, max_tir_len=5000))]
    fn do_rusty_grf(py: Python<'_>, fragments: Vec<(String, Vec<u8>, u64, u64)>,
                    threads: usize, max_tir_len: i32) -> PyResult<Vec<Record>> {
        // TIR-Learner's locked `-c 1` params; default integer-arm mode.
        let p = Param { seed: 10, seed_mismatch: 1, min_stem: 10, max_stem: i32::MAX, percent: 20,
                        min_space: 10, max_space: max_tir_len, min_tsd: 2, max_tsd: 10,
                        max_mismatch: i32::MAX, min_spacer_len: 10, max_spacer_len: max_tir_len,
                        emit_cigar: false, threads, legacy_fasta: false };
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build()
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        let recs = py.allow_threads(|| {
            pool.install(|| {
                use rayon::prelude::*;
                fragments.par_iter().flat_map(|(seqid, seq, off, owned)| {
                    let up: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
                    process_fragment(seqid, &up, *off as usize, *owned as usize, &p)
                }).collect::<Vec<_>>()
            })
        });
        Ok(recs)
    }

    #[pymodule]
    fn grf_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(do_rusty_grf, m)?)?;
        Ok(())
    }
}
