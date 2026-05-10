/// Minimizer extraction and filtering.
///
/// Minimizers are small k-mer hashes used as fingerprints for sequences.
/// Each barcode's linked reads produce a set of minimizers that characterize
/// the genomic region covered by that barcode's molecules.
///
/// This module handles:
/// - Extracting (k,w)-minimizers from FASTA/FASTQ sequences
/// - Filtering barcodes by minimizer count
/// - Removing singleton and high-frequency minimizers
/// - Building the inverted index (minimizer → barcodes)
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::{BufRead, Write};

/// Filter barcodes by minimizer count.
/// Removes barcodes with fewer than `min_n` or >= `max_n` minimizers.
/// Returns (removed_too_few, removed_too_many).
pub fn filter_barcodes(
    bx_to_mxs: &mut FxHashMap<String, FxHashSet<u64>>,
    min_n: usize,
    max_n: usize,
) -> (usize, usize) {
    let initial = bx_to_mxs.len();
    let mut too_few = 0;
    let mut too_many = 0;

    bx_to_mxs.retain(|_, mxs| {
        if mxs.len() < min_n {
            too_few += 1;
            false
        } else if mxs.len() >= max_n {
            too_many += 1;
            false
        } else {
            true
        }
    });

    log::info!(
        "Filtered barcodes: {} total, {} too few (<{}), {} too many (>={}), {} remaining",
        initial,
        too_few,
        min_n,
        too_many,
        max_n,
        bx_to_mxs.len()
    );
    (too_few, too_many)
}

/// Count the frequency of each minimizer across all barcodes.
pub fn count_minimizers(bx_to_mxs: &FxHashMap<String, FxHashSet<u64>>) -> FxHashMap<u64, u32> {
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    for mxs in bx_to_mxs.values() {
        for &mx in mxs {
            *counts.entry(mx).or_insert(0) += 1;
        }
    }
    counts
}

/// Remove singleton minimizers (those occurring in only one barcode).
/// Returns the number of singletons removed.
///
/// Uses a two-layer cascading Bloom filter for memory efficiency:
/// - Layer 1 (`first_occ`): tracks minimizers seen at least once
/// - Layer 2 (`second_occ`): tracks minimizers seen at least twice
///
/// A minimizer in `first_occ` but NOT in `second_occ` is a singleton.
/// This uses ~600 MB for human-scale data instead of ~6 GB for a HashMap.
///
/// False positives in `second_occ` mean some singletons are kept (harmless —
/// they just won't contribute to overlaps). No non-singletons are lost.
pub fn remove_singletons(bx_to_mxs: &mut FxHashMap<String, FxHashSet<u64>>) -> usize {
    // Count total minimizer occurrences to size the BFs
    let mut total_mxs = 0u64;
    for mxs in bx_to_mxs.values() {
        total_mxs += mxs.len() as u64;
    }

    if total_mxs == 0 {
        return 0;
    }

    // Size BFs for ~1% FPR: ~10 bits per element, 3 hash functions
    // Each layer needs to hold up to `total_mxs` distinct elements
    let bits_per_element = 10u64;
    let n_hashes = 3u32;
    let bf_size_bytes = (total_mxs * bits_per_element / 8).max(1024);

    log::info!(
        "Cascading BF singleton removal: {} total minimizer occurrences, BF size={:.1} MB each",
        total_mxs,
        bf_size_bytes as f64 / 1e6
    );

    // Pass 1: populate the two BF layers
    let mut first_occ = crate::repeat::BloomFilter::new(bf_size_bytes, n_hashes);
    let mut second_occ = crate::repeat::BloomFilter::new(bf_size_bytes, n_hashes);

    for mxs in bx_to_mxs.values() {
        for &mx in mxs {
            if first_occ.contains(mx) {
                second_occ.insert(mx);
            } else {
                first_occ.insert(mx);
            }
        }
    }

    log::info!(
        "Cascading BF: first_occ FPR={:.4}%, second_occ FPR={:.4}%",
        first_occ.fpr() * 100.0,
        second_occ.fpr() * 100.0
    );

    // Pass 2: remove minimizers NOT in second_occ (singletons)
    let mut n_removed = 0usize;
    for mxs in bx_to_mxs.values_mut() {
        let before = mxs.len();
        mxs.retain(|&mx| second_occ.contains(mx));
        n_removed += before - mxs.len();
    }

    log::info!("Removed {} singleton minimizers", n_removed);
    n_removed
}

/// Remove singleton minimizers using exact HashMap counting.
/// This is the original implementation — accurate but uses more memory.
/// Kept as a fallback for small datasets or when exact counts are needed.
#[allow(dead_code)]
pub fn remove_singletons_exact(bx_to_mxs: &mut FxHashMap<String, FxHashSet<u64>>) -> usize {
    let counts = count_minimizers(bx_to_mxs);
    let singletons: FxHashSet<u64> = counts
        .iter()
        .filter(|(_, &count)| count < 2)
        .map(|(&mx, _)| mx)
        .collect();

    let n_singletons = singletons.len();
    for mxs in bx_to_mxs.values_mut() {
        mxs.retain(|mx| !singletons.contains(mx));
    }

    log::info!(
        "Removed {} singleton minimizers of {} total (exact)",
        n_singletons,
        counts.len()
    );
    n_singletons
}

/// Remove high-frequency (repetitive) minimizers.
/// If `max_count` is None, uses Q3 + 1.5*IQR as the threshold.
/// Returns the number of minimizers removed.
pub fn remove_repetitive(
    bx_to_mxs: &mut FxHashMap<String, FxHashSet<u64>>,
    max_count: Option<u32>,
) -> usize {
    let counts = count_minimizers(bx_to_mxs);
    let mut values: Vec<u32> = counts.values().copied().collect();
    values.sort_unstable();

    let threshold = match max_count {
        Some(c) => c,
        None => {
            if values.is_empty() {
                return 0;
            }
            let q1 = values[values.len() / 4];
            let q3 = values[3 * values.len() / 4];
            let iqr = q3 - q1;
            q3 + (1.5 * iqr as f64) as u32
        }
    };

    let repetitive: FxHashSet<u64> = counts
        .iter()
        .filter(|(_, &count)| count >= threshold)
        .map(|(&mx, _)| mx)
        .collect();

    let n_removed = repetitive.len();
    for mxs in bx_to_mxs.values_mut() {
        mxs.retain(|mx| !repetitive.contains(mx));
    }

    // Remove empty barcodes
    bx_to_mxs.retain(|_, mxs| !mxs.is_empty());

    log::info!(
        "Removed {} repetitive minimizers (threshold={}), {} barcodes remaining",
        n_removed,
        threshold,
        bx_to_mxs.len()
    );
    n_removed
}

/// Build an inverted index: minimizer → set of barcodes containing it.
pub fn build_inverted_index(
    bx_to_mxs: &FxHashMap<String, FxHashSet<u64>>,
) -> FxHashMap<u64, Vec<String>> {
    let mut mx_to_bxs: FxHashMap<u64, Vec<String>> = FxHashMap::default();
    for (bx, mxs) in bx_to_mxs {
        for &mx in mxs {
            mx_to_bxs.entry(mx).or_default().push(bx.clone());
        }
    }
    log::info!("Indexed {} minimizers", mx_to_bxs.len());
    mx_to_bxs
}

/// Build an inverted index using integer barcode IDs for memory efficiency.
pub fn build_inverted_index_ids(
    bx_to_mxs: &FxHashMap<String, FxHashSet<u64>>,
    bx_to_id: &FxHashMap<String, u32>,
) -> FxHashMap<u64, Vec<u32>> {
    let mut mx_to_ids: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
    for (bx, mxs) in bx_to_mxs {
        if let Some(&id) = bx_to_id.get(bx) {
            for &mx in mxs {
                mx_to_ids.entry(mx).or_default().push(id);
            }
        }
    }
    mx_to_ids
}

/// Assign integer IDs to barcodes. Returns (bx→id, id→bx) maps.
pub fn assign_barcode_ids(
    bx_to_mxs: &FxHashMap<String, FxHashSet<u64>>,
) -> (FxHashMap<String, u32>, Vec<String>) {
    let mut bx_to_id = FxHashMap::default();
    let mut id_to_bx = Vec::new();
    let mut sorted_bxs: Vec<&String> = bx_to_mxs.keys().collect();
    sorted_bxs.sort();
    for bx in sorted_bxs {
        let id = id_to_bx.len() as u32;
        bx_to_id.insert(bx.clone(), id);
        id_to_bx.push(bx.clone());
    }
    (bx_to_id, id_to_bx)
}

// ─── Minimizer extraction from sequences ─────────────────────────────────────

/// Hash a k-mer using the invertible hash function from the original Physlr.
/// This is a bijective hash for 64-bit integers.
fn hash_kmer(mut key: u64) -> u64 {
    key = (!key).wrapping_add(key << 21);
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8);
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4);
    key ^= key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

/// Encode a DNA base as a 2-bit value. Returns None for non-ACGT bases.
pub fn encode_base(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Complement of a 2-bit encoded base.
pub fn complement_2bit(b: u64) -> u64 {
    b ^ 3
}

/// Hash a k-mer byte slice to u64 using FNV-1a. Used for k > 32.
#[inline]
fn hash_kmer_bytes(data: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Compute canonical hash of a k-mer byte slice (for k > 32).
/// Canonical = min(hash(fwd), hash(revcomp)).
#[inline]
fn canonical_kmer_hash_bytes(seq: &[u8]) -> u64 {
    let fwd = hash_kmer_bytes(seq);
    let rc: Vec<u8> = seq
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect();
    let rc_hash = hash_kmer_bytes(&rc);
    fwd.min(rc_hash)
}

/// Extract the canonical k-mer hash at each position of a sequence.
/// Returns a vector of (position, hash) pairs.
/// For k <= 32, uses 2-bit encoding. For k > 32, uses byte-level hashing.
fn kmer_hashes(seq: &[u8], k: usize) -> Vec<(usize, u64)> {
    if seq.len() < k || k == 0 {
        return Vec::new();
    }

    if k <= 32 {
        // Fast path: 2-bit encoding
        let mask = if k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        let mut hashes = Vec::with_capacity(seq.len() - k + 1);
        let mut fwd: u64 = 0;
        let mut rev: u64 = 0;
        let mut valid = 0;

        for (i, &base) in seq.iter().enumerate() {
            if let Some(b) = encode_base(base) {
                fwd = ((fwd << 2) | b) & mask;
                rev = (rev >> 2) | (complement_2bit(b) << (2 * (k - 1)));
                valid += 1;
            } else {
                valid = 0;
                fwd = 0;
                rev = 0;
            }
            if valid >= k {
                let canonical = fwd.min(rev);
                hashes.push((i + 1 - k, hash_kmer(canonical)));
            }
        }
        hashes
    } else {
        // k > 32: slide window, hash each k-mer from bytes
        let mut hashes = Vec::with_capacity(seq.len() - k + 1);
        let mut run_start = 0usize;
        let mut in_run = false;

        for (i, &base) in seq.iter().enumerate() {
            let valid = matches!(base, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't');
            if valid {
                if !in_run {
                    run_start = i;
                    in_run = true;
                }
            } else {
                if in_run {
                    let run = &seq[run_start..i];
                    if run.len() >= k {
                        for j in 0..=(run.len() - k) {
                            let h = canonical_kmer_hash_bytes(&run[j..j + k]);
                            hashes.push((run_start + j, hash_kmer(h)));
                        }
                    }
                }
                in_run = false;
            }
        }
        if in_run {
            let run = &seq[run_start..];
            if run.len() >= k {
                for j in 0..=(run.len() - k) {
                    let h = canonical_kmer_hash_bytes(&run[j..j + k]);
                    hashes.push((run_start + j, hash_kmer(h)));
                }
            }
        }
        hashes
    }
}

/// Extract (k,w)-minimizers from a sequence.
/// A minimizer is the minimum hash in each window of w consecutive k-mers.
/// Returns a deduplicated, ordered list of minimizer hashes.
pub fn extract_minimizers(seq: &[u8], k: usize, w: usize) -> Vec<u64> {
    let hashes = kmer_hashes(seq, k);
    if hashes.is_empty() || w == 0 {
        return Vec::new();
    }

    let mut minimizers = Vec::new();
    let mut last_min: Option<(usize, u64)> = None;

    for window_start in 0..hashes.len().saturating_sub(w - 1) {
        let window_end = (window_start + w).min(hashes.len());
        let window = &hashes[window_start..window_end];

        let min_entry = window.iter().min_by_key(|(_, h)| *h).unwrap();

        if last_min.is_none_or(|(pos, hash)| pos != min_entry.0 || hash != min_entry.1) {
            minimizers.push(min_entry.1);
            last_min = Some(*min_entry);
        }
    }

    // Deduplicate while preserving order
    let mut seen = FxHashSet::default();
    minimizers.retain(|h| seen.insert(*h));
    minimizers
}

/// Extract minimizers from a FASTA/FASTQ file, grouping by barcode.
///
/// For linked reads, the barcode is extracted from the read header (BX:Z: tag).
/// For plain FASTA (e.g., contigs), each sequence name is its own "barcode".
///
/// Returns a map from barcode/name → set of minimizer hashes.
pub fn index_file(
    path: &str,
    k: usize,
    w: usize,
) -> anyhow::Result<FxHashMap<String, FxHashSet<u64>>> {
    let mut bx_to_mxs: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();

    let mut reader = needletail::parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot open {}: {}", path, e))?;

    let mut n_reads = 0u64;
    let mut n_minimizers = 0u64;

    while let Some(record) = reader.next() {
        let record = record?;
        let header = std::str::from_utf8(record.id()).unwrap_or("");

        // Extract barcode from BX:Z: tag, or use sequence name
        let barcode = extract_barcode(header).unwrap_or_else(|| {
            header
                .split_whitespace()
                .next()
                .unwrap_or(header)
                .to_string()
        });

        let seq = record.seq();
        let mxs = extract_minimizers(&seq, k, w);

        n_reads += 1;
        n_minimizers += mxs.len() as u64;

        let entry = bx_to_mxs.entry(barcode).or_default();
        for mx in mxs {
            entry.insert(mx);
        }

        if n_reads.is_multiple_of(1_000_000) {
            log::info!("Processed {} reads, {} barcodes", n_reads, bx_to_mxs.len());
        }
    }

    log::info!(
        "Indexed {} reads into {} barcodes with {} total minimizers",
        n_reads,
        bx_to_mxs.len(),
        n_minimizers
    );

    Ok(bx_to_mxs)
}

/// Extract minimizers from a FASTA file, preserving order per sequence.
/// Used for mapping contigs where position matters.
pub fn index_file_ordered(
    path: &str,
    k: usize,
    w: usize,
) -> anyhow::Result<FxHashMap<String, Vec<u64>>> {
    let mut result: FxHashMap<String, Vec<u64>> = FxHashMap::default();

    let mut reader = needletail::parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot open {}: {}", path, e))?;

    while let Some(record) = reader.next() {
        let record = record?;
        let name = std::str::from_utf8(record.id())
            .unwrap_or("")
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();

        let seq = record.seq();
        let mxs = extract_minimizers(&seq, k, w);
        result.insert(name, mxs);
    }

    Ok(result)
}

// ─── Homopolymer compression ──────────────────────────────────────────────────

/// Compress homopolymer runs in a sequence.
///
/// Collapses consecutive identical bases (e.g., `AAACCCGG` → `ACG`).
/// Returns the compressed sequence and a coordinate map where
/// `coord_map[compressed_pos] = original_pos` for reversibility.
pub fn homopolymer_compress(seq: &[u8]) -> (Vec<u8>, Vec<usize>) {
    if seq.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let mut compressed = Vec::with_capacity(seq.len() / 2);
    let mut coord_map = Vec::with_capacity(seq.len() / 2);

    let mut prev = seq[0];
    compressed.push(prev);
    coord_map.push(0);

    for (i, &base) in seq.iter().enumerate().skip(1) {
        if base != prev {
            compressed.push(base);
            coord_map.push(i);
            prev = base;
        }
    }

    compressed.shrink_to_fit();
    coord_map.shrink_to_fit();
    (compressed, coord_map)
}

/// Index a FASTQ/FASTA file for long reads.
///
/// Differences from linked-read indexing:
/// - Read name is always used as the barcode (no BX:Z: extraction)
/// - Homopolymer compression is applied before minimizer extraction
/// - Each read is its own "barcode"
pub fn index_file_long(
    path: &str,
    k: usize,
    w: usize,
) -> anyhow::Result<FxHashMap<String, FxHashSet<u64>>> {
    index_file_long_inner(path, k, w, None)
}

/// Index a FASTQ/FASTA file for long reads with optional repeat filtering.
pub fn index_file_long_with_repeat_filter(
    path: &str,
    k: usize,
    w: usize,
    repeat_bf: &crate::repeat::BloomFilter,
) -> anyhow::Result<FxHashMap<String, FxHashSet<u64>>> {
    index_file_long_inner(path, k, w, Some(repeat_bf))
}

/// Inner implementation for long-read indexing with optional repeat filter.
fn index_file_long_inner(
    path: &str,
    k: usize,
    w: usize,
    repeat_bf: Option<&crate::repeat::BloomFilter>,
) -> anyhow::Result<FxHashMap<String, FxHashSet<u64>>> {
    let mut bx_to_mxs: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();

    let mut reader = needletail::parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot open {}: {}", path, e))?;

    let mut n_reads = 0u64;
    let mut n_minimizers = 0u64;
    let mut total_bases = 0u64;
    let mut compressed_bases = 0u64;

    while let Some(record) = reader.next() {
        let record = record?;
        let header = std::str::from_utf8(record.id()).unwrap_or("");

        // Always use read name as barcode for long reads
        let read_name = header
            .split_whitespace()
            .next()
            .unwrap_or(header)
            .to_string();

        let seq = record.seq();
        total_bases += seq.len() as u64;

        // Apply homopolymer compression
        let (compressed_seq, _coord_map) = homopolymer_compress(&seq);
        compressed_bases += compressed_seq.len() as u64;

        let mxs = if let Some(bf) = repeat_bf {
            extract_minimizers_filtered(&compressed_seq, k, w, bf)
        } else {
            extract_minimizers(&compressed_seq, k, w)
        };

        n_reads += 1;
        n_minimizers += mxs.len() as u64;

        let entry = bx_to_mxs.entry(read_name).or_default();
        for mx in mxs {
            entry.insert(mx);
        }

        if n_reads.is_multiple_of(100_000) {
            log::info!(
                "Processed {} reads, {} unique read IDs, compression ratio: {:.2}%",
                n_reads,
                bx_to_mxs.len(),
                if total_bases > 0 {
                    compressed_bases as f64 / total_bases as f64 * 100.0
                } else {
                    100.0
                }
            );
        }
    }

    log::info!(
        "Long-read indexing: {} reads, {} minimizers, compression {}/{} bases ({:.1}%), repeat_filter={}",
        n_reads,
        n_minimizers,
        compressed_bases,
        total_bases,
        if total_bases > 0 {
            compressed_bases as f64 / total_bases as f64 * 100.0
        } else {
            100.0
        },
        repeat_bf.is_some()
    );

    Ok(bx_to_mxs)
}

/// Extract the barcode from a FASTQ header.
/// Looks for BX:Z: tag (10x Genomics format) or #barcode (stLFR format).
fn extract_barcode(header: &str) -> Option<String> {
    // 10x Genomics format: BX:Z:ACGTACGT-1
    for part in header.split_whitespace() {
        if let Some(bx) = part.strip_prefix("BX:Z:") {
            return Some(bx.to_string());
        }
    }

    // stLFR format: readname#barcode (barcode is e.g. 123_456_789)
    // #0 and #0_0_0 are unassigned barcodes
    if let Some(pos) = header.find('#') {
        let bx = &header[pos + 1..];
        let bx = bx.split_whitespace().next().unwrap_or(bx);
        if !bx.is_empty() && bx != "0" && bx != "0_0_0" {
            return Some(bx.to_string());
        }
    }

    None
}

/// Write minimizer index to TSV format.
pub fn write_minimizer_tsv(
    bx_to_mxs: &FxHashMap<String, FxHashSet<u64>>,
    writer: &mut dyn Write,
) -> std::io::Result<()> {
    let mut sorted_bxs: Vec<&String> = bx_to_mxs.keys().collect();
    sorted_bxs.sort();
    for bx in sorted_bxs {
        let mxs = &bx_to_mxs[bx];
        if mxs.is_empty() {
            continue;
        }
        let mut sorted_mxs: Vec<u64> = mxs.iter().copied().collect();
        sorted_mxs.sort_unstable();
        write!(writer, "{}\t", bx)?;
        for (i, mx) in sorted_mxs.iter().enumerate() {
            if i > 0 {
                write!(writer, " ")?;
            }
            write!(writer, "{}", mx)?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Write ordered minimizer index to TSV format.
pub fn write_minimizer_list_tsv(
    bx_to_mxs: &FxHashMap<String, Vec<u64>>,
    writer: &mut dyn Write,
) -> std::io::Result<()> {
    let mut sorted_bxs: Vec<&String> = bx_to_mxs.keys().collect();
    sorted_bxs.sort();
    for bx in sorted_bxs {
        let mxs = &bx_to_mxs[bx];
        if mxs.is_empty() {
            continue;
        }
        write!(writer, "{}\t", bx)?;
        for (i, mx) in mxs.iter().enumerate() {
            if i > 0 {
                write!(writer, " ")?;
            }
            write!(writer, "{}", mx)?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Extract minimizers using btllib's indexlr external tool.
///
/// Calls `indexlr --bx -k K -w W -t T <input>` and parses the output.
/// indexlr output format with --bx: `barcode\thash1 hash2 hash3...`
/// Each hash may include colon-separated sub-fields (pos, strand, seq);
/// we take only the first colon-delimited component as the hash value.
pub fn index_file_btllib(
    path: &str,
    k: usize,
    w: usize,
    threads: usize,
) -> anyhow::Result<FxHashMap<String, FxHashSet<u64>>> {
    // Verify indexlr is available
    let which = std::process::Command::new("which").arg("indexlr").output();
    if which.is_err() || !which.unwrap().status.success() {
        anyhow::bail!("btllib indexlr not found in PATH. Install btllib or use --indexer builtin.");
    }

    log::info!(
        "Running btllib indexlr (k={}, w={}, threads={}) on {}",
        k,
        w,
        threads,
        path
    );

    let child = std::process::Command::new("indexlr")
        .arg("--bx")
        .arg("-k")
        .arg(k.to_string())
        .arg("-w")
        .arg(w.to_string())
        .arg("-t")
        .arg(threads.to_string())
        .arg(path)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| anyhow::anyhow!("Failed to spawn indexlr: {}", e))?;

    let stdout = child
        .stdout
        .ok_or_else(|| anyhow::anyhow!("Failed to capture indexlr stdout"))?;

    let reader = std::io::BufReader::new(stdout);
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

        let entry = bx_to_mxs.entry(barcode).or_default();
        for token in mxs_str.split_whitespace() {
            // indexlr may output hash:pos:strand:seq — take only the hash
            let hash_str = token.split(':').next().unwrap_or(token);
            if let Ok(hash) = hash_str.parse::<u64>() {
                entry.insert(hash);
            }
        }

        n_lines += 1;
        if n_lines.is_multiple_of(1_000_000) {
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
        "btllib indexlr: parsed {} lines into {} barcodes",
        n_lines,
        bx_to_mxs.len()
    );

    Ok(bx_to_mxs)
}

// ---------------------------------------------------------------------------
// Bloom-filter-aware minimizer extraction (repeat filtering)
// ---------------------------------------------------------------------------

/// Compute k-mer hashes, skipping k-mers present in the repeat Bloom filter.
/// Returns (position, hash) pairs for non-repetitive k-mers only.
fn kmer_hashes_filtered(
    seq: &[u8],
    k: usize,
    repeat_bf: &crate::repeat::BloomFilter,
) -> Vec<(usize, u64)> {
    if k == 0 || seq.len() < k {
        return Vec::new();
    }

    if k <= 32 {
        let mask: u64 = if k < 32 {
            (1u64 << (2 * k)) - 1
        } else {
            u64::MAX
        };
        let mut fwd: u64 = 0;
        let mut rev: u64 = 0;
        let mut valid = 0usize;
        let mut result = Vec::new();

        for (i, &base) in seq.iter().enumerate() {
            if let Some(b) = encode_base(base) {
                fwd = ((fwd << 2) | b) & mask;
                rev = (rev >> 2) | (complement_2bit(b) << (2 * (k - 1)));
                valid += 1;
                if valid >= k {
                    let canonical = fwd.min(rev);
                    if repeat_bf.contains(canonical) {
                        continue;
                    }
                    result.push((i + 1 - k, hash_kmer(canonical)));
                }
            } else {
                valid = 0;
                fwd = 0;
                rev = 0;
            }
        }
        result
    } else {
        // k > 32: slide window, hash from bytes, filter via BF
        let mut result = Vec::new();
        let mut run_start = 0usize;
        let mut in_run = false;

        for (i, &base) in seq.iter().enumerate() {
            let valid = matches!(base, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't');
            if valid {
                if !in_run {
                    run_start = i;
                    in_run = true;
                }
            } else {
                if in_run {
                    let run = &seq[run_start..i];
                    if run.len() >= k {
                        for j in 0..=(run.len() - k) {
                            let canonical = canonical_kmer_hash_bytes(&run[j..j + k]);
                            if repeat_bf.contains(canonical) {
                                continue;
                            }
                            result.push((run_start + j, hash_kmer(canonical)));
                        }
                    }
                }
                in_run = false;
            }
        }
        if in_run {
            let run = &seq[run_start..];
            if run.len() >= k {
                for j in 0..=(run.len() - k) {
                    let canonical = canonical_kmer_hash_bytes(&run[j..j + k]);
                    if repeat_bf.contains(canonical) {
                        continue;
                    }
                    result.push((run_start + j, hash_kmer(canonical)));
                }
            }
        }
        result
    }
}

/// Extract minimizers from a sequence, skipping repetitive k-mers.
pub fn extract_minimizers_filtered(
    seq: &[u8],
    k: usize,
    w: usize,
    repeat_bf: &crate::repeat::BloomFilter,
) -> Vec<u64> {
    let hashes = kmer_hashes_filtered(seq, k, repeat_bf);
    if hashes.is_empty() || w == 0 {
        return Vec::new();
    }

    let mut minimizers = Vec::new();
    let mut prev_min: Option<u64> = None;

    for window in hashes.windows(w) {
        let min_hash = window.iter().map(|(_, h)| *h).min().unwrap();
        if prev_min != Some(min_hash) {
            minimizers.push(min_hash);
            prev_min = Some(min_hash);
        }
    }

    // If fewer hashes than window size, take the minimum of all
    if hashes.len() < w && !hashes.is_empty() {
        let min_hash = hashes.iter().map(|(_, h)| *h).min().unwrap();
        if prev_min != Some(min_hash) {
            minimizers.push(min_hash);
        }
    }

    minimizers
}

/// Index a FASTQ/FASTA file with repeat filtering via Bloom filter.
///
/// This is the equivalent of `indexlr -r repeat.bf` in the original pipeline.
pub fn index_file_with_repeat_filter(
    path: &str,
    k: usize,
    w: usize,
    repeat_bf: &crate::repeat::BloomFilter,
) -> anyhow::Result<FxHashMap<String, FxHashSet<u64>>> {
    log::info!(
        "Indexing minimizers from {} (k={}, w={}) with repeat filter...",
        path,
        k,
        w
    );

    let mut bx_to_mxs: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();
    let mut n_reads = 0u64;
    let mut n_minimizers = 0u64;
    let _n_skipped = 0u64;

    let mut fq_reader = needletail::parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot parse {}: {}", path, e))?;

    while let Some(record) = fq_reader.next() {
        let record = record?;
        let header = std::str::from_utf8(record.id()).unwrap_or("");

        let barcode = extract_barcode(header).unwrap_or_else(|| {
            header
                .split_whitespace()
                .next()
                .unwrap_or(header)
                .to_string()
        });

        let seq = record.seq();
        let mxs = extract_minimizers_filtered(&seq, k, w, repeat_bf);

        n_reads += 1;
        n_minimizers += mxs.len() as u64;

        let entry = bx_to_mxs.entry(barcode).or_default();
        for mx in mxs {
            entry.insert(mx);
        }

        if n_reads.is_multiple_of(1_000_000) {
            log::info!(
                "Processed {} reads, {} barcodes ({} minimizers)",
                n_reads,
                bx_to_mxs.len(),
                n_minimizers
            );
        }
    }

    log::info!(
        "Indexed {} reads into {} barcodes: {} minimizers (repetitive k-mers filtered)",
        n_reads,
        bx_to_mxs.len(),
        n_minimizers
    );

    Ok(bx_to_mxs)
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // extract_minimizers
    // -----------------------------------------------------------------------

    #[test]
    fn test_extract_minimizers_basic() {
        let seq = b"ACGTACGTACGTACGT";
        let mxs = extract_minimizers(seq, 8, 4);
        assert!(!mxs.is_empty());
    }

    #[test]
    fn test_extract_minimizers_deterministic() {
        let seq = b"ACGTACGTACGTACGT";
        let m1 = extract_minimizers(seq, 8, 4);
        let m2 = extract_minimizers(seq, 8, 4);
        assert_eq!(m1, m2);
    }

    #[test]
    fn test_extract_minimizers_short_seq() {
        // Sequence shorter than k → no minimizers
        let seq = b"ACGT";
        let mxs = extract_minimizers(seq, 8, 4);
        assert!(mxs.is_empty());
    }

    #[test]
    fn test_extract_minimizers_empty() {
        let mxs = extract_minimizers(b"", 8, 4);
        assert!(mxs.is_empty());
    }

    #[test]
    fn test_extract_minimizers_w1() {
        // w=1: every k-mer is a minimizer (before dedup)
        let seq = b"ACGTACGTACGT";
        let mxs = extract_minimizers(seq, 4, 1);
        assert!(!mxs.is_empty());
    }

    #[test]
    fn test_extract_minimizers_canonical() {
        // Reverse complement should give same minimizers
        let fwd = b"ACGTACGT";
        let rev = b"ACGTACGT"; // ACGTACGT is its own reverse complement
        let m1 = extract_minimizers(fwd, 4, 2);
        let m2 = extract_minimizers(rev, 4, 2);
        assert_eq!(m1, m2);
    }

    #[test]
    fn test_extract_minimizers_with_n() {
        // N bases should break k-mer chains
        let seq = b"ACGTNNNNACGT";
        let mxs = extract_minimizers(seq, 4, 1);
        // Should still get minimizers from the valid parts
        assert!(!mxs.is_empty());
    }

    // -----------------------------------------------------------------------
    // kmer_hashes
    // -----------------------------------------------------------------------

    #[test]
    fn test_kmer_hashes_basic() {
        let hashes = kmer_hashes(b"ACGTACGT", 4);
        assert_eq!(hashes.len(), 5); // 8 - 4 + 1 = 5
    }

    #[test]
    fn test_kmer_hashes_too_short() {
        let hashes = kmer_hashes(b"ACG", 4);
        assert!(hashes.is_empty());
    }

    #[test]
    fn test_kmer_hashes_k_zero() {
        let hashes = kmer_hashes(b"ACGT", 0);
        assert!(hashes.is_empty());
    }

    #[test]
    fn test_kmer_hashes_k_too_large() {
        let hashes = kmer_hashes(b"ACGT", 33);
        assert!(hashes.is_empty());
    }

    #[test]
    fn test_kmer_hashes_positions() {
        let hashes = kmer_hashes(b"ACGTAC", 4);
        // Positions should be 0, 1, 2
        let positions: Vec<usize> = hashes.iter().map(|(p, _)| *p).collect();
        assert_eq!(positions, vec![0, 1, 2]);
    }

    // -----------------------------------------------------------------------
    // extract_barcode
    // -----------------------------------------------------------------------

    #[test]
    fn test_barcode_10x() {
        let bc = extract_barcode("read1 BX:Z:ACGTACGT-1 other:stuff");
        assert_eq!(bc, Some("ACGTACGT-1".to_string()));
    }

    #[test]
    fn test_barcode_stlfr() {
        let bc = extract_barcode("read1#123_456_789");
        assert_eq!(bc, Some("123_456_789".to_string()));
    }

    #[test]
    fn test_barcode_stlfr_unassigned_0() {
        let bc = extract_barcode("read1#0");
        assert_eq!(bc, None);
    }

    #[test]
    fn test_barcode_stlfr_unassigned_000() {
        let bc = extract_barcode("read1#0_0_0");
        assert_eq!(bc, None);
    }

    #[test]
    fn test_barcode_none() {
        let bc = extract_barcode("plain_read_name");
        assert_eq!(bc, None);
    }

    #[test]
    fn test_barcode_bx_priority() {
        // BX:Z: should take priority even if # is present
        let bc = extract_barcode("read#123 BX:Z:ACGT-1");
        assert_eq!(bc, Some("ACGT-1".to_string()));
    }

    // -----------------------------------------------------------------------
    // filter_barcodes
    // -----------------------------------------------------------------------

    #[test]
    fn test_filter_barcodes() {
        let mut bx: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();
        bx.insert("a".into(), [1, 2, 3].into_iter().collect());
        bx.insert("b".into(), [1].into_iter().collect());
        bx.insert("c".into(), (1..=100).collect());

        let (too_few, too_many) = filter_barcodes(&mut bx, 2, 50);
        assert_eq!(too_few, 1); // "b" removed
        assert_eq!(too_many, 1); // "c" removed
        assert_eq!(bx.len(), 1); // only "a" remains
    }

    // -----------------------------------------------------------------------
    // count_minimizers
    // -----------------------------------------------------------------------

    #[test]
    fn test_count_minimizers() {
        let mut bx: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();
        bx.insert("a".into(), [10, 20, 30].into_iter().collect());
        bx.insert("b".into(), [20, 30, 40].into_iter().collect());

        let counts = count_minimizers(&bx);
        assert_eq!(counts[&10], 1);
        assert_eq!(counts[&20], 2);
        assert_eq!(counts[&30], 2);
        assert_eq!(counts[&40], 1);
    }

    // -----------------------------------------------------------------------
    // remove_singletons
    // -----------------------------------------------------------------------

    #[test]
    fn test_remove_singletons() {
        let mut bx: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();
        bx.insert("a".into(), [10, 20].into_iter().collect());
        bx.insert("b".into(), [10, 30].into_iter().collect());
        // minimizer 20 appears only in "a", 30 only in "b", 10 in both

        let removed = remove_singletons(&mut bx);
        assert!(removed > 0);
        // After removal, only minimizer 10 should remain (appears in 2 barcodes)
        for mxs in bx.values() {
            assert!(mxs.contains(&10));
            assert!(!mxs.contains(&20));
            assert!(!mxs.contains(&30));
        }
    }

    #[test]
    fn test_remove_singletons_bloom_vs_exact() {
        // Verify cascading BF produces same results as exact method on a larger dataset
        let mut bx_bloom: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();
        let mut bx_exact: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();

        // Create 100 barcodes with overlapping minimizers
        // Minimizers 0-49 appear in multiple barcodes, 1000+ are singletons
        for i in 0..100u64 {
            let mut mxs = FxHashSet::default();
            // Shared minimizers (appear in ~10 barcodes each)
            for j in 0..5 {
                mxs.insert((i / 10) * 5 + j); // 0-49
            }
            // Singleton minimizer unique to this barcode
            mxs.insert(1000 + i);
            bx_bloom.insert(format!("bx_{}", i), mxs.clone());
            bx_exact.insert(format!("bx_{}", i), mxs);
        }

        let removed_bloom = remove_singletons(&mut bx_bloom);
        let removed_exact = remove_singletons_exact(&mut bx_exact);

        // Both should remove the same singletons
        assert_eq!(removed_bloom, removed_exact,
            "bloom removed {} vs exact removed {}", removed_bloom, removed_exact);

        // Verify same minimizers remain
        for key in bx_bloom.keys() {
            assert_eq!(bx_bloom[key], bx_exact[key],
                "mismatch for barcode {}", key);
        }
    }

    #[test]
    fn test_remove_singletons_empty() {
        let mut bx: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();
        let removed = remove_singletons(&mut bx);
        assert_eq!(removed, 0);
    }

    // -----------------------------------------------------------------------
    // hash_kmer
    // -----------------------------------------------------------------------

    #[test]
    fn test_hash_kmer_deterministic() {
        assert_eq!(hash_kmer(42), hash_kmer(42));
    }

    #[test]
    fn test_hash_kmer_different() {
        assert_ne!(hash_kmer(0), hash_kmer(1));
    }

    // -----------------------------------------------------------------------
    // encode_base / complement_2bit
    // -----------------------------------------------------------------------

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'A'), Some(0));
        assert_eq!(encode_base(b'C'), Some(1));
        assert_eq!(encode_base(b'G'), Some(2));
        assert_eq!(encode_base(b'T'), Some(3));
        assert_eq!(encode_base(b'a'), Some(0));
        assert_eq!(encode_base(b'N'), None);
    }

    #[test]
    fn test_complement_2bit() {
        assert_eq!(complement_2bit(0), 3); // A -> T
        assert_eq!(complement_2bit(1), 2); // C -> G
        assert_eq!(complement_2bit(2), 1); // G -> C
        assert_eq!(complement_2bit(3), 0); // T -> A
    }

    // -----------------------------------------------------------------------
    // homopolymer_compress
    // -----------------------------------------------------------------------

    #[test]
    fn test_homopolymer_compress_basic() {
        let (compressed, coord_map) = homopolymer_compress(b"AAACCCGGG");
        assert_eq!(compressed, b"ACG");
        assert_eq!(coord_map, vec![0, 3, 6]);
    }

    #[test]
    fn test_homopolymer_compress_no_runs() {
        let (compressed, coord_map) = homopolymer_compress(b"ACGT");
        assert_eq!(compressed, b"ACGT");
        assert_eq!(coord_map, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_homopolymer_compress_single_base() {
        let (compressed, coord_map) = homopolymer_compress(b"AAAA");
        assert_eq!(compressed, b"A");
        assert_eq!(coord_map, vec![0]);
    }

    #[test]
    fn test_homopolymer_compress_empty() {
        let (compressed, coord_map) = homopolymer_compress(b"");
        assert!(compressed.is_empty());
        assert!(coord_map.is_empty());
    }

    #[test]
    fn test_homopolymer_compress_mixed() {
        // AATTCCGG -> ATCG, coords [0, 2, 4, 6]
        let (compressed, coord_map) = homopolymer_compress(b"AATTCCGG");
        assert_eq!(compressed, b"ATCG");
        assert_eq!(coord_map, vec![0, 2, 4, 6]);
    }

    #[test]
    fn test_homopolymer_compress_reversible() {
        // Verify we can map compressed positions back to original
        let original = b"AAACGTTTAAC";
        let (compressed, coord_map) = homopolymer_compress(original);
        assert_eq!(compressed, b"ACGTAC");
        // Each compressed position maps to the first base of its run
        for (comp_pos, &orig_pos) in coord_map.iter().enumerate() {
            assert_eq!(compressed[comp_pos], original[orig_pos]);
        }
    }

    #[test]
    fn test_homopolymer_compress_with_n() {
        // N bases are treated like any other character
        let (compressed, _) = homopolymer_compress(b"AAANNNACGT");
        assert_eq!(compressed, b"ANACGT");
    }
}
