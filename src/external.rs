/// External tool integration for repeat filtering and indexing.
///
/// Calls ntcard, nthits, physlr-makebf, and indexlr as external processes,
/// matching the pipeline used in profile_pipeline.sh.
use anyhow::{bail, Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Check that an external tool is available in PATH.
fn require_tool(name: &str) -> Result<()> {
    let output = Command::new("which")
        .arg(name)
        .output()
        .with_context(|| format!("Failed to check for {}", name))?;
    if !output.status.success() {
        bail!(
            "'{}' not found in PATH. Install btllib (provides ntcard, nthits, indexlr) \
             and ensure physlr-makebf is built and in PATH.",
            name
        );
    }
    Ok(())
}

/// Run ntcard to produce a k-mer frequency histogram.
///
/// Equivalent to: `ntcard -t THREADS -k K -o OUTPREFIX.histogram INPUT...`
pub fn run_ntcard(
    input_files: &[String],
    k: usize,
    threads: usize,
    output_histogram: &Path,
) -> Result<()> {
    require_tool("ntcard")?;

    log::info!(
        "Running ntcard (k={}, threads={}) on {} file(s)...",
        k,
        threads,
        input_files.len()
    );

    let mut cmd = Command::new("ntcard");
    cmd.arg("-t")
        .arg(threads.to_string())
        .arg("-k")
        .arg(k.to_string())
        .arg("-o")
        .arg(output_histogram.to_str().unwrap());

    for f in input_files {
        cmd.arg(f);
    }

    let status = cmd.status().with_context(|| "Failed to run ntcard")?;

    if !status.success() {
        bail!("ntcard failed with exit code {:?}", status.code());
    }

    // ntcard appends _kK to the output filename
    let actual_path = format!("{}_k{}", output_histogram.to_str().unwrap(), k);
    if Path::new(&actual_path).exists() && !output_histogram.exists() {
        std::fs::rename(&actual_path, output_histogram).with_context(|| {
            format!("Failed to rename {} to {:?}", actual_path, output_histogram)
        })?;
    }

    log::info!("ntcard histogram written to {:?}", output_histogram);
    Ok(())
}

/// Parse ntcard histogram to find the mode (first mode after minimum).
///
/// Reimplements find-ntcard-mode.py in Rust.
pub fn find_ntcard_mode(histogram_path: &Path) -> Result<u64> {
    let content = std::fs::read_to_string(histogram_path)
        .with_context(|| format!("Failed to read histogram {:?}", histogram_path))?;

    let freq_counts: Vec<u64> = content
        .lines()
        .filter(|line| !line.starts_with('k') && !line.starts_with('F') && !line.is_empty())
        .filter_map(|line| {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                parts[2].parse::<u64>().ok()
            } else {
                // Try space-separated
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 3 {
                    parts[2].parse::<u64>().ok()
                } else {
                    None
                }
            }
        })
        .collect();

    if freq_counts.is_empty() {
        bail!("No frequency data found in histogram {:?}", histogram_path);
    }

    // Find the first minimum
    let mut min_idx = 0;
    let mut min_val = freq_counts[0];
    for (idx, &freq) in freq_counts.iter().enumerate() {
        if freq > min_val {
            min_idx = if idx > 0 { idx - 1 } else { 0 };
            break;
        }
        min_val = freq;
    }

    // Find the mode after the minimum
    let slice = &freq_counts[min_idx..];
    let max_val = *slice.iter().max().unwrap();
    let mode_offset = slice.iter().position(|&v| v == max_val).unwrap();
    let mode = (mode_offset + 1 + min_idx) as u64;

    log::info!("K-mer histogram mode: {}", mode);
    Ok(mode)
}

/// Run nthits to find repetitive k-mers above a threshold.
///
/// Equivalent to: `nthits -t THREADS -k K -c THRESHOLD -p PREFIX INPUT...`
/// Produces PREFIX_kK.rep file containing repetitive k-mers.
pub fn run_nthits(
    input_files: &[String],
    k: usize,
    threshold: u64,
    threads: usize,
    output_prefix: &str,
) -> Result<PathBuf> {
    require_tool("nthits")?;

    log::info!(
        "Running nthits (k={}, threshold={}, threads={})...",
        k,
        threshold,
        threads
    );

    let mut cmd = Command::new("nthits");
    cmd.arg("-t")
        .arg(threads.to_string())
        .arg("-k")
        .arg(k.to_string())
        .arg("-c")
        .arg(threshold.to_string())
        .arg("-p")
        .arg(output_prefix);

    for f in input_files {
        cmd.arg(f);
    }

    let status = cmd.status().with_context(|| "Failed to run nthits")?;

    if !status.success() {
        bail!("nthits failed with exit code {:?}", status.code());
    }

    // nthits outputs PREFIX_kK.rep
    let rep_path = PathBuf::from(format!("{}_k{}.rep", output_prefix, k));
    if !rep_path.exists() {
        // Try without _kK suffix
        let alt_path = PathBuf::from(format!("{}.rep", output_prefix));
        if alt_path.exists() {
            return Ok(alt_path);
        }
        bail!(
            "nthits output not found at {:?} or {:?}",
            rep_path,
            alt_path
        );
    }

    log::info!("nthits repetitive k-mers written to {:?}", rep_path);
    Ok(rep_path)
}

/// Run physlr-makebf to build a btllib Bloom filter from nthits output.
///
/// Equivalent to: `physlr-makebf -k K -b BF_SIZE -t THREADS -o OUTPUT INPUT.rep`
pub fn run_makebf(
    rep_file: &Path,
    k: usize,
    bf_size: u64,
    threads: usize,
    output_bf: &Path,
    makebf_path: Option<&str>,
) -> Result<()> {
    let makebf = makebf_path.unwrap_or("physlr-makebf");
    require_tool(makebf)?;

    log::info!(
        "Running physlr-makebf (k={}, bf_size={}, threads={})...",
        k,
        bf_size,
        threads
    );

    let status = Command::new(makebf)
        .arg("-k")
        .arg(k.to_string())
        .arg("-b")
        .arg(bf_size.to_string())
        .arg("-t")
        .arg(threads.to_string())
        .arg("-o")
        .arg(output_bf.to_str().unwrap())
        .arg(rep_file.to_str().unwrap())
        .status()
        .with_context(|| format!("Failed to run {}", makebf))?;

    if !status.success() {
        bail!("physlr-makebf failed with exit code {:?}", status.code());
    }

    log::info!("Bloom filter written to {:?}", output_bf);
    Ok(())
}

/// Run indexlr with repeat filtering via a btllib Bloom filter.
///
/// Equivalent to: `indexlr --bx -t5 -k K -w W -r BF_PATH -o OUTPUT INPUT...`
/// Returns the parsed barcode → minimizer set map.
pub fn run_indexlr(
    input_files: &[String],
    k: usize,
    w: usize,
    threads: usize,
    repeat_bf: Option<&Path>,
    long_reads: bool,
) -> Result<FxHashMap<String, FxHashSet<u64>>> {
    require_tool("indexlr")?;

    log::info!(
        "Running indexlr (k={}, w={}, threads={}, repeat_filter={}, long={})...",
        k,
        w,
        threads,
        repeat_bf.is_some(),
        long_reads
    );

    let mut cmd = Command::new("indexlr");

    if long_reads {
        cmd.arg("--long");
    } else {
        cmd.arg("--bx");
    }

    cmd.arg("-k")
        .arg(k.to_string())
        .arg("-w")
        .arg(w.to_string())
        .arg("-t")
        .arg(threads.to_string());

    if let Some(bf_path) = repeat_bf {
        cmd.arg("-r").arg(bf_path.to_str().unwrap());
    }

    for f in input_files {
        cmd.arg(f);
    }

    cmd.stdout(Stdio::piped()).stderr(Stdio::piped());

    let child = cmd.spawn().with_context(|| "Failed to spawn indexlr")?;

    let stdout = child
        .stdout
        .ok_or_else(|| anyhow::anyhow!("Failed to capture indexlr stdout"))?;

    let reader = BufReader::new(stdout);
    let mut bx_to_mxs: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();
    let mut n_lines = 0u64;

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }

        let mut parts = line.splitn(2, '\t');
        let barcode = match parts.next() {
            Some(b) if !b.is_empty() && b != "NA" => b.to_string(),
            _ => continue,
        };

        let mxs_str = match parts.next() {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };

        // For long reads without --bx, the first column is the read name
        let entry = bx_to_mxs.entry(barcode).or_default();
        for token in mxs_str.split_whitespace() {
            // indexlr may output hash:pos:strand:seq — take only the hash
            let hash_str = token.split(':').next().unwrap_or(token);
            if let Ok(hash) = hash_str.parse::<u64>() {
                entry.insert(hash);
            }
        }

        n_lines += 1;
        if n_lines % 1_000_000 == 0 {
            log::info!(
                "Parsed {} lines from indexlr, {} barcodes",
                n_lines,
                bx_to_mxs.len()
            );
        }
    }

    // Remove empty barcodes
    bx_to_mxs.retain(|_, mxs| !mxs.is_empty());

    log::info!(
        "indexlr: parsed {} lines into {} barcodes",
        n_lines,
        bx_to_mxs.len()
    );

    Ok(bx_to_mxs)
}

/// Write a homopolymer-compressed version of FASTQ/FASTA files.
///
/// For long reads, collapses homopolymer runs before passing to external tools.
/// Returns the path to the compressed file.
fn write_homopolymer_compressed_fastq(input_files: &[String], output_path: &Path) -> Result<()> {
    use std::io::Write;

    let mut writer = std::io::BufWriter::new(
        std::fs::File::create(output_path)
            .with_context(|| format!("Failed to create {:?}", output_path))?,
    );

    let mut total_reads = 0u64;
    let mut total_bases = 0u64;
    let mut compressed_bases = 0u64;

    for input_file in input_files {
        let mut reader = needletail::parse_fastx_file(input_file)
            .map_err(|e| anyhow::anyhow!("Cannot open {}: {}", input_file, e))?;

        while let Some(record) = reader.next() {
            let record = record?;
            let header = std::str::from_utf8(record.id()).unwrap_or("");
            let seq = record.seq();
            total_bases += seq.len() as u64;

            let (compressed_seq, coord_map) = crate::minimizer::homopolymer_compress(&seq);
            compressed_bases += compressed_seq.len() as u64;

            // Compress quality string to match: keep quality at each run's start position
            let compressed_qual: Vec<u8> = match record.qual() {
                Some(qual) if !qual.is_empty() => coord_map.iter().map(|&pos| qual[pos]).collect(),
                _ => vec![b'I'; compressed_seq.len()], // No quality available, use placeholder
            };

            writeln!(writer, "@{}", header)?;
            writer.write_all(&compressed_seq)?;
            writeln!(writer)?;
            writeln!(writer, "+")?;
            writer.write_all(&compressed_qual)?;
            writeln!(writer)?;

            total_reads += 1;
            if total_reads % 100_000 == 0 {
                log::info!(
                    "Compressed {} reads ({:.1}% compression ratio)",
                    total_reads,
                    compressed_bases as f64 / total_bases as f64 * 100.0
                );
            }
        }
    }

    log::info!(
        "Homopolymer compression: {} reads, {}/{} bases ({:.1}%)",
        total_reads,
        compressed_bases,
        total_bases,
        compressed_bases as f64 / total_bases as f64 * 100.0
    );

    Ok(())
}

/// Full repeat-filter + index pipeline using external tools.
///
/// Runs: ntcard → find-ntcard-mode → nthits → physlr-makebf → indexlr -r
///
/// This matches the pipeline in profile_pipeline.sh and produces results
/// consistent with the published benchmarks.
///
/// For long reads, a homopolymer-compressed FASTQ is generated first and
/// used as input to the external tools.
pub fn repeat_filter_and_index(
    input_files: &[String],
    k: usize,
    w: usize,
    threads: usize,
    work_dir: &Path,
    prefix: &str,
    bf_size: u64,
    multiplier: u64,
    long_reads: bool,
    makebf_path: Option<&str>,
) -> Result<FxHashMap<String, FxHashSet<u64>>> {
    std::fs::create_dir_all(work_dir)?;

    // For long reads, pre-compress FASTQ with homopolymer compression
    let effective_inputs: Vec<String>;
    let _compressed_file_guard: Option<PathBuf>; // keep path alive for cleanup

    if long_reads {
        let compressed_path = work_dir.join(format!("{}.hpc.fq", prefix));
        log::info!("Homopolymer-compressing input for long-read mode...");
        write_homopolymer_compressed_fastq(input_files, &compressed_path)?;
        effective_inputs = vec![compressed_path.to_str().unwrap().to_string()];
        _compressed_file_guard = Some(compressed_path);
    } else {
        effective_inputs = input_files.to_vec();
        _compressed_file_guard = None;
    }

    // Step 1: ntcard
    let histogram_path = work_dir.join(format!("{}_k{}.histogram", prefix, k));
    run_ntcard(&effective_inputs, k, threads, &histogram_path)?;

    // Step 2: find mode
    let mode = find_ntcard_mode(&histogram_path)?;
    let threshold = mode * multiplier;
    log::info!(
        "Repeat threshold: {} (mode={} x {})",
        threshold,
        mode,
        multiplier
    );

    // Step 3: nthits
    let nthits_prefix = work_dir.join(prefix).to_str().unwrap().to_string();
    let rep_path = run_nthits(&effective_inputs, k, threshold, threads, &nthits_prefix)?;

    // Step 4: physlr-makebf
    let bf_path = work_dir.join(format!("{}.k{}.bf", prefix, k));
    run_makebf(&rep_path, k, bf_size, threads, &bf_path, makebf_path)?;

    // Step 5: indexlr with repeat filter
    let bx_to_mxs = run_indexlr(&effective_inputs, k, w, threads, Some(&bf_path), long_reads)?;

    Ok(bx_to_mxs)
}

/// Run indexlr on a reference FASTA to get ordered minimizer lists per sequence.
///
/// Uses `--long --id` mode (no `--bx`). Returns a map from sequence name to
/// an ordered Vec of minimizer hashes, compatible with `map_paf`.
///
/// This ensures the reference uses the same ntHash function as the reads
/// (which were also indexed by indexlr).
pub fn run_indexlr_reference_ordered(
    fasta_path: &str,
    k: usize,
    w: usize,
    threads: usize,
) -> Result<FxHashMap<String, Vec<u64>>> {
    require_tool("indexlr")?;

    log::info!(
        "Running indexlr on reference (k={}, w={}, threads={})...",
        k, w, threads
    );

    let mut cmd = Command::new("indexlr");
    cmd.arg("--long")
        .arg("--id")
        .arg("-k").arg(k.to_string())
        .arg("-w").arg(w.to_string())
        .arg("-t").arg(threads.to_string())
        .arg(fasta_path)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    let child = cmd.spawn().with_context(|| "Failed to spawn indexlr for reference")?;
    let stdout = child
        .stdout
        .ok_or_else(|| anyhow::anyhow!("Failed to capture indexlr stdout"))?;

    let reader = BufReader::new(stdout);
    let mut result: FxHashMap<String, Vec<u64>> = FxHashMap::default();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }

        let mut parts = line.splitn(2, '\t');
        let name = match parts.next() {
            Some(n) if !n.is_empty() => n.to_string(),
            _ => continue,
        };

        let mxs_str = match parts.next() {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };

        let mut mxs = Vec::new();
        let mut seen = FxHashSet::default();
        for token in mxs_str.split_whitespace() {
            let hash_str = token.split(':').next().unwrap_or(token);
            if let Ok(hash) = hash_str.parse::<u64>() {
                // Deduplicate while preserving order (matches extract_minimizers behavior)
                if seen.insert(hash) {
                    mxs.push(hash);
                }
            }
        }

        log::info!("{}: {} minimizers", name, mxs.len());
        result.insert(name, mxs);
    }

    log::info!("indexlr reference: indexed {} sequences", result.len());
    Ok(result)
}

/// Run indexlr on a reference FASTA with `--pos` to build a position map.
///
/// Writes a TSV with columns: chr, minimizer_index, bp_position, total_minimizers, seq_length.
/// The `step` parameter controls subsampling (every Nth minimizer, plus first and last).
///
/// For HPC mode, pass the HPC-compressed FASTA and an `HpcIndex` to translate
/// positions back to original coordinates.
pub fn run_indexlr_reference_positions(
    fasta_path: &str,
    k: usize,
    w: usize,
    threads: usize,
    step: usize,
    hpc_index: Option<&crate::minimizer::HpcIndex>,
    writer: &mut dyn std::io::Write,
) -> Result<()> {
    require_tool("indexlr")?;

    // First pass: get sequence lengths from the FASTA
    let mut seq_lengths: FxHashMap<String, usize> = FxHashMap::default();
    let mut reader = needletail::parse_fastx_file(fasta_path)
        .map_err(|e| anyhow::anyhow!("Cannot open {}: {}", fasta_path, e))?;
    while let Some(record) = reader.next() {
        let record = record?;
        let name = std::str::from_utf8(record.id())
            .unwrap_or("")
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        seq_lengths.insert(name, record.seq().len());
    }

    log::info!(
        "Running indexlr --pos on reference (k={}, w={}, threads={})...",
        k, w, threads
    );

    let mut cmd = Command::new("indexlr");
    cmd.arg("--long")
        .arg("--id")
        .arg("--pos")
        .arg("-k").arg(k.to_string())
        .arg("-w").arg(w.to_string())
        .arg("-t").arg(threads.to_string())
        .arg(fasta_path)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    let child = cmd.spawn().with_context(|| "Failed to spawn indexlr for positions")?;
    let stdout = child
        .stdout
        .ok_or_else(|| anyhow::anyhow!("Failed to capture indexlr stdout"))?;

    let buf_reader = BufReader::new(stdout);

    for line in buf_reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }

        let mut parts = line.splitn(2, '\t');
        let name = match parts.next() {
            Some(n) if !n.is_empty() => n.to_string(),
            _ => continue,
        };

        let mxs_str = match parts.next() {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };

        // Parse hash:pos tokens
        let mut positions: Vec<usize> = Vec::new();
        let mut seen_hashes = FxHashSet::default();
        for token in mxs_str.split_whitespace() {
            let fields: Vec<&str> = token.split(':').collect();
            if fields.len() >= 2 {
                if let (Ok(hash), Ok(pos)) = (
                    fields[0].parse::<u64>(),
                    fields[1].parse::<usize>(),
                ) {
                    // Deduplicate by hash, preserving order
                    if seen_hashes.insert(hash) {
                        positions.push(pos);
                    }
                }
            }
        }

        let total = positions.len();
        if total == 0 {
            continue;
        }

        // Get the appropriate sequence length
        let (seq_len, hpc_seq_map) = if let Some(idx) = hpc_index {
            // HPC mode: use original sequence length
            let sm = idx.get(&name);
            let orig_len = sm.map_or(
                *seq_lengths.get(&name).unwrap_or(&0),
                |s| s.original_len,
            );
            (orig_len, sm)
        } else {
            (*seq_lengths.get(&name).unwrap_or(&0), None)
        };

        for (i, &bp) in positions.iter().enumerate() {
            if i == 0 || i == total - 1 || i % step == 0 {
                // Translate HPC position to original if needed
                let out_bp = if let Some(sm) = hpc_seq_map {
                    if bp < sm.coord_map.len() {
                        sm.coord_map[bp]
                    } else {
                        sm.original_len
                    }
                } else {
                    bp
                };
                writeln!(writer, "{}\t{}\t{}\t{}\t{}", name, i, out_bp, total, seq_len)?;
            }
        }

        log::info!("{}: {} bp, {} minimizers", name, seq_len, total);
    }

    Ok(())
}

/// Write a minimizer list TSV from ordered minimizer data (from indexlr).
///
/// Each line: sequence_name\thash1 hash2 hash3 ...
pub fn write_minimizer_list_tsv_ordered(
    ref_mxs: &FxHashMap<String, Vec<u64>>,
    writer: &mut dyn std::io::Write,
) -> Result<()> {
    for (name, mxs) in ref_mxs {
        let mx_strs: Vec<String> = mxs.iter().map(|h| h.to_string()).collect();
        writeln!(writer, "{}\t{}", name, mx_strs.join(" "))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_ntcard_mode_basic() {
        // Create a temporary histogram file
        let dir = std::env::temp_dir().join("physlr_test_ntcard");
        std::fs::create_dir_all(&dir).unwrap();
        let hist_path = dir.join("test.histogram");

        // Simulate ntcard histogram: columns are multiplicity, f, count
        // Typical shape: high at f=1 (errors), dip around f=5-6, then peak at coverage mode ~11
        // The "k" and "F0" lines are header/summary lines that get filtered out.
        let content = "\
k\t32\t0
F0\t1000000\t0
1\t500000\t500000
2\t200000\t200000
3\t100000\t100000
4\t80000\t80000
5\t70000\t70000
6\t60000\t60000
7\t65000\t65000
8\t80000\t80000
9\t100000\t100000
10\t120000\t120000
11\t130000\t130000
12\t125000\t125000
13\t110000\t110000
14\t90000\t90000
15\t70000\t70000
";

        std::fs::write(&hist_path, content).unwrap();

        let mode = find_ntcard_mode(&hist_path).unwrap();
        // The F0 line (count=0) is filtered by the "k" prefix check but F0 isn't.
        // F0 has count 0, then f=1 has 500000, so the minimum is at F0 (idx=0, val=0).
        // Since 500000 > 0, min_idx = 0-1 clamped to 0. Then mode is max of full slice = f=11.
        assert_eq!(mode, 11, "mode={} expected 11", mode);

        std::fs::remove_dir_all(&dir).ok();
    }
}
