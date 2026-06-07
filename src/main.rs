// CLI front-end for the grf_rs library (GRF -c 1 / MITE detection).
// All logic lives in lib.rs; this just parses args and dispatches.
use grf_rs::{run_batch, run_file, set_threads, Param};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut p = Param::default();
    let (mut input, mut outdir) = (String::new(), String::from("."));
    let mut batch: Option<String> = None; // --batch <listfile>: one input fasta path per line
    let mut k = 1;
    while k < args.len() {
        let a = args[k].clone();
        let mut val = |k: &mut usize| { *k += 1; args[*k].clone() };
        match a.as_str() {
            "-i" => input = val(&mut k),
            "-o" => outdir = val(&mut k),
            "-c" => { val(&mut k); }
            "-t" => p.threads = val(&mut k).parse().unwrap_or(1).max(1),
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
            "--rle-cigar" => p.emit_cigar = true,
            "--legacy-fasta" => p.legacy_fasta = true,
            "--batch" => batch = Some(val(&mut k)),
            _ => {}
        }
        k += 1;
    }

    set_threads(p.threads);

    match batch {
        Some(listfile) => {
            let list = std::fs::read_to_string(&listfile).expect("read batch list");
            let paths: Vec<String> = list.lines().map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty()).collect();
            let secs = run_batch(&paths, &outdir, &p);
            eprintln!("[batch] files={} threads={} wall={:.1}s", paths.len(), p.threads, secs);
        }
        None => {
            let fname = if p.legacy_fasta { "candidate.fasta" } else { "candidate.json" };
            run_file(&input, &format!("{}/{}", outdir, fname), &p);
        }
    }
}
