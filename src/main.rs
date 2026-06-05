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

// LZ76 complexity via a suffix automaton (O(n)), same result as the stock O(n^2) find().
// Build the SAM of seq tracking firstpos (first-occurrence END pos) per state; a factor
// extends while seq[send..=i] occurs earlier (firstpos < i) and commits otherwise.
fn seq_complexity(seq: &[u8]) -> f64 {
    use std::collections::HashMap;
    let n = seq.len();
    let denom = n as f64 / ((n as f64).ln() / 4f64.ln());
    let c_thresh = (0.675 * denom).ceil() as u32;
    // ---- online suffix automaton with firstpos ----
    let mut len: Vec<i32> = vec![0];
    let mut link: Vec<i32> = vec![-1];
    let mut next: Vec<HashMap<u8, u32>> = vec![HashMap::new()];
    let mut firstpos: Vec<i32> = vec![-1];
    let mut last = 0usize;
    for (pos, &ch) in seq.iter().enumerate() {
        let cur = len.len();
        len.push(len[last] + 1);
        link.push(-1);
        next.push(HashMap::new());
        firstpos.push(pos as i32);
        let mut p = last as i32;
        while p != -1 && !next[p as usize].contains_key(&ch) {
            next[p as usize].insert(ch, cur as u32);
            p = link[p as usize];
        }
        if p == -1 {
            link[cur] = 0;
        } else {
            let q = next[p as usize][&ch] as usize;
            if len[p as usize] + 1 == len[q] {
                link[cur] = q as i32;
            } else {
                let clone = len.len();
                len.push(len[p as usize] + 1);
                link.push(link[q]);
                next.push(next[q].clone());
                firstpos.push(firstpos[q]); // clone keeps q's first occurrence
                while p != -1 && next[p as usize].get(&ch) == Some(&(q as u32)) {
                    next[p as usize].insert(ch, clone as u32);
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
    let mut cur = 0usize; // init state
    for i in 1..n {
        let ext = match next[cur].get(&seq[i]) {
            Some(&nxt) if (firstpos[nxt as usize] as usize) < i => Some(nxt as usize),
            _ => None,
        };
        match ext {
            Some(nxt) => cur = nxt,
            None => {
                c += 1;
                if c >= c_thresh { return c as f64 / denom; }
                cur = 0;
            }
        }
    }
    c as f64 / denom
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

// no-indel stem extension -> (cigar, tr1, tr2); empty cigar = no valid arm
fn get_stem(s: &[u8], p: &Param) -> (String, u32, u32) {
    thread_local! {
        static SCRATCH: std::cell::RefCell<(Vec<u8>, Vec<i32>)> =
            std::cell::RefCell::new((Vec::new(), Vec::new()));
    }
    SCRATCH.with(|sc| {
        let mut sc = sc.borrow_mut();
        let (left, pos) = &mut *sc;
        left.clear();
        pos.clear();
        let len = s.len();
        let mut error = 0i32;
        let (mut i, mut j) = (0usize, len - 1);
        while i < j {
            if is_pair(s[i], s[j]) { left.push(b'm'); }
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

fn filter_low_complex(seq: &[u8], tir: &str) -> bool {
    let len1 = get_stem_len(tir, 'D') as usize;
    let len2 = get_stem_len(tir, 'I') as usize;
    let l = &seq[0..len1];
    let r = &seq[seq.len() - len2..];
    let gc1 = gc_content(l) as f64;
    let gc2 = gc_content(r) as f64;
    if gc1 < 0.2 || gc1 > 0.8 || gc2 < 0.2 || gc2 > 0.8 { return false; }
    if low_complex_regex(l) || low_complex_regex(r) { return false; }
    if seq_complexity(seq) < 0.675 { return false; }
    true
}

fn detect_str(i_sp: i32, s_tr: &[(i32, i32)], seq: &[u8], p: &Param, v: &mut Vec<Mite>) {
    let seed = p.seed;
    let num_i = s_tr.len() as i64 - i_sp as i64 - seed as i64;
    if num_i <= 0 { return; }
    let num = num_i as usize;
    let l = (i_sp + 2 * seed) as usize;
    let max_abs = p.seed_mismatch * 2;
    let off = (seed + i_sp) as usize;
    // zip elides bounds checks on the 2.5e9 hot path (off+num == s_tr.len())
    for (j, (a, b)) in s_tr[..num].iter().zip(&s_tr[off..]).enumerate() {
        let abs_sum = (a.0 + b.0).abs() + (a.1 + b.1).abs();
        if abs_sum > max_abs { continue; }
        let (start, end) = (j, j + l - 1);
        if !is_pair(seq[start], seq[end]) { continue; }
        let mut unpair = 0i32;
        for k in 1..seed as usize {
            if !is_pair(seq[start + k], seq[end - k]) { unpair += 1; }
        }
        if unpair > p.seed_mismatch { continue; }
        let cand = &seq[start..start + l];
        let tsd = detect_tsd(seq, start, end, p.min_tsd, p.max_tsd);
        if tsd.is_empty() { continue; }
        let (cigar, tr1, tr2) = get_stem(cand, p); // max_indel == 0
        if !cigar.is_empty() && filter_low_complex(cand, &cigar) {
            v.push(Mite { start, end, tr1, tr2, tsd, tir: cigar });
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

    let nspace = (p.max_space - p.min_space + 1) as usize;
    let mut candidates: Vec<Vec<Vec<Mite>>> = Vec::with_capacity(chroms.len());
    for seq in &seqs {
        let n = seq.len();
        let seed = p.seed as usize;
        let mut spacer_vecs: Vec<Vec<Mite>> = (0..nspace).map(|_| Vec::new()).collect();
        if n >= (2 * p.seed + p.max_space + 2 * p.max_tsd) as usize {
            let mut cum: Vec<(i32, i32)> = Vec::with_capacity(n);
            let (mut a, mut b) = (0i32, 0i32);
            for &ch in seq {
                let (x, y) = match ch { b'A' => (1, 0), b'C' => (0, 1), b'G' => (0, -1), b'T' => (-1, 0), _ => (100, 0) };
                a += x; b += y;
                cum.push((a, b));
            }
            let mut s_tr: Vec<(i32, i32)> = vec![(0, 0); n - seed + 1];
            s_tr[0] = cum[seed - 1];
            for i in seed..n {
                s_tr[i - seed + 1] = (cum[i].0 - cum[i - seed].0, cum[i].1 - cum[i - seed].1);
            }
            for i_sp in p.min_space..=p.max_space {
                detect_str(i_sp, &s_tr, seq, &p, &mut spacer_vecs[(i_sp - p.min_space) as usize]);
            }
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
