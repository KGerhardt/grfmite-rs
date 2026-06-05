// Direct Rust port of GRF grf-main `-c 1` (MITE) detection, aiming for byte-identical
// candidate.fasta vs stock grf-main. Correctness first (byte-wise); 2-bit/SIMD later.
use needletail::parse_fastx_file;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Clone)]
struct Param {
    seed: i32,
    seed_mismatch: i32,
    min_stem: i32, // = --min_tr
    max_stem: i32,
    percent: i32, // = -p
    min_space: i32,
    max_space: i32,
    min_tsd: i32,
    max_tsd: i32,
    max_mismatch: i32,
    min_spacer_len: i32,
    max_spacer_len: i32,
}
impl Default for Param {
    fn default() -> Self {
        Param { seed: 10, seed_mismatch: 1, min_stem: 0, max_stem: i32::MAX, percent: 10,
                min_space: -1, max_space: -1, min_tsd: 2, max_tsd: 10,
                max_mismatch: i32::MAX, min_spacer_len: 0, max_spacer_len: i32::MAX }
    }
}

struct Mite { start: usize, end: usize, tr1: u32, tr2: u32, tsd: String, tir: String }

#[inline]
fn is_pair(a: u8, b: u8) -> bool {
    matches!((a, b), (b'A', b'T') | (b'C', b'G') | (b'G', b'C') | (b'T', b'A'))
}

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

// sum of run-length counts whose symbol != c  (cigar parse)
fn get_stem_len(s: &str, c: char) -> i32 {
    let (mut num, mut count) = (0i32, 0i32);
    for ch in s.chars() {
        let d = ch as i32 - '0' as i32;
        if (0..=9).contains(&d) { count = count * 10 + d; }
        else { if ch != c { num += count; } count = 0; }
    }
    num
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

// no-indel stem extension -> (cigar, tr1, tr2); empty cigar = no valid arm.
// `e` is the candidate's 2-bit encoding (ACGT->0..3, non-ACGT->4): pairing is
// e[i]^e[j]==0b11 (complement = XOR 0b11; the 4-sentinel can never XOR to 0b11),
// exactly equivalent to byte is_pair, so the cigar stays byte-identical.
fn get_stem(e: &[u8], p: &Param) -> (String, u32, u32) {
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
        // trim trailing 'M' (len2 = index after last 'm', or 0 if none)
        let len2 = left.iter().rposition(|&c| c == b'm').map(|x| x + 1).unwrap_or(0);
        left.truncate(len2);
        // percent != 100 path
        pos.push(len2 as i32);
        for idx in (0..pos.len()).rev() {
            if pos[idx] >= p.min_stem && (idx as i32) * 100 <= p.percent * pos[idx] {
                let cigar = compress(&left[0..(pos[idx] as usize).min(left.len())]); // C++ substr clamps
                return (cigar, pos[idx] as u32, pos[idx] as u32);
            }
        }
        (String::new(), 0, 0)
    })
}

fn filter_low_complex(seq: &[u8], tir: &str, a2c: &[u8; 256], alpha: usize) -> bool {
    let len1 = get_stem_len(tir, 'D') as usize;
    let len2 = get_stem_len(tir, 'I') as usize;
    let l = &seq[0..len1];
    let r = &seq[seq.len() - len2..];
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
fn join_detect(gcode: &[u32], offsets: &[u32], starts_by_code: &[u32],
               seq: &[u8], enc: &[u8], p: &Param, spacer_vecs: &mut [Vec<Mite>],
               a2c: &[u8; 256], alpha: usize) {
    let seed = p.seed as usize;
    let n = seq.len();
    let min_off2 = (p.min_space + 2 * p.seed - 1) as usize; // end - start at min spacer
    let max_off2 = (p.max_space + 2 * p.seed - 1) as usize; // end - start at max spacer
    let twoseed_m1 = (2 * p.seed - 1) as usize;
    let min_space = p.min_space as usize;
    let mut codes: Vec<u32> = Vec::with_capacity(1 + 3 * (seed - 1));
    for end in min_off2..n {
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
                let i_sp = end - start - twoseed_m1; // in [min_space, max_space]
                let cand = &seq[start..start + l];
                let tsd = detect_tsd(seq, start, end, p.min_tsd, p.max_tsd);
                if tsd.is_empty() { continue; }
                let (cigar, tr1, tr2) = get_stem(&enc[start..start + l], p); // max_indel == 0
                if !cigar.is_empty() && filter_low_complex(cand, &cigar, a2c, alpha) {
                    spacer_vecs[i_sp - min_space].push(Mite { start, end, tr1, tr2, tsd, tir: cigar });
                }
            }
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut p = Param::default();
    let (mut input, mut outdir) = (String::new(), String::from("."));
    let mut k = 1;
    while k < args.len() {
        let a = args[k].clone();
        let mut val = |k: &mut usize| { *k += 1; args[*k].clone() };
        match a.as_str() {
            "-i" => input = val(&mut k),
            "-o" => outdir = val(&mut k),
            "-c" => { val(&mut k); }
            "-t" => { val(&mut k); }
            "-p" => p.percent = val(&mut k).parse().unwrap(),
            "-s" => p.seed = val(&mut k).parse().unwrap(),
            "--seed_mismatch" => p.seed_mismatch = val(&mut k).parse().unwrap(),
            "--min_tr" => p.min_stem = val(&mut k).parse().unwrap(),
            "--max_indel" => { val(&mut k); }
            "--min_space" => p.min_space = val(&mut k).parse().unwrap(),
            "--max_space" => p.max_space = val(&mut k).parse().unwrap(),
            "--min_tsd" => p.min_tsd = val(&mut k).parse().unwrap(),
            "--max_tsd" => p.max_tsd = val(&mut k).parse().unwrap(),
            "--min_spacer_len" => p.min_spacer_len = val(&mut k).parse().unwrap(),
            "--max_spacer_len" => p.max_spacer_len = val(&mut k).parse().unwrap(),
            _ => {}
        }
        k += 1;
    }

    // read fasta (uppercased); chrom name = first whitespace token of the header
    let mut chroms: Vec<String> = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut reader = parse_fastx_file(&input).expect("open fasta");
    while let Some(rec) = reader.next() {
        let rec = rec.expect("record");
        let id = rec.id();
        let name = id.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(b"");
        chroms.push(String::from_utf8_lossy(name).into_owned());
        seqs.push(rec.seq().iter().map(|b| b.to_ascii_uppercase()).collect());
    }

    // Dense byte->code map over the genome's actual alphabet (injective on every byte that
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

    // output (spacer, chrom, candidate order) + inline grf-filter (TR len + spacer len)
    let f = File::create(format!("{}/candidate.fasta", outdir)).expect("create out");
    let mut w = BufWriter::new(f);
    for j in 0..nspace {
        for ci in 0..chroms.len() {
            for m in &candidates[ci][j] {
                let len1 = get_stem_len(&m.tir, 'D');
                let len2 = get_stem_len(&m.tir, 'I');
                let len3 = (m.end as i32) - (m.start as i32) + 1 - len1 - len2;
                if len1 >= p.min_stem && len1 <= p.max_stem
                    && len2 >= p.min_stem && len2 <= p.max_stem
                    && len3 >= p.min_spacer_len && len3 <= p.max_spacer_len {
                    writeln!(w, ">{}:{}:{}:{}:{}", chroms[ci], m.start + 1, m.end + 1, m.tir, m.tsd).unwrap();
                    writeln!(w, "{}", std::str::from_utf8(&seqs[ci][m.start..=m.end]).unwrap()).unwrap();
                }
            }
        }
    }
}
