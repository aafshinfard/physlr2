//! Repetitive k-mer detection and filtering.
//!
//! Uses a Count-Min Sketch for streaming k-mer frequency estimation,
//! avoiding the need to store all k-mers in memory. This allows processing
//! human-scale datasets (250GB+ FASTQ) in bounded memory (~4GB for the CMS).

use std::io::Write;

// --- Bloom filter ---

/// A Bloom filter for k-mer membership queries.
pub struct BloomFilter {
    bits: Vec<u64>,
    n_bits: u64,
    n_hashes: u32,
}

impl BloomFilter {
    pub fn new(size_bytes: u64, n_hashes: u32) -> Self {
        let n_bits = size_bytes * 8;
        let n_words = n_bits.div_ceil(64) as usize;
        Self {
            bits: vec![0u64; n_words],
            n_bits,
            n_hashes,
        }
    }

    pub fn insert(&mut self, kmer_hash: u64) {
        for i in 0..self.n_hashes {
            let bit = self.hash_to_bit(kmer_hash, i);
            self.bits[(bit / 64) as usize] |= 1u64 << (bit % 64);
        }
    }

    pub fn contains(&self, kmer_hash: u64) -> bool {
        for i in 0..self.n_hashes {
            let bit = self.hash_to_bit(kmer_hash, i);
            if self.bits[(bit / 64) as usize] & (1u64 << (bit % 64)) == 0 {
                return false;
            }
        }
        true
    }

    #[inline]
    fn hash_to_bit(&self, kmer_hash: u64, i: u32) -> u64 {
        let h1 = kmer_hash;
        let h2 = kmer_hash.wrapping_mul(0x9E3779B97F4A7C15).wrapping_shr(17);
        (h1.wrapping_add((i as u64).wrapping_mul(h2))) % self.n_bits
    }

    pub fn popcount(&self) -> u64 {
        self.bits.iter().map(|w| w.count_ones() as u64).sum()
    }

    pub fn fpr(&self) -> f64 {
        let p = self.popcount() as f64 / self.n_bits as f64;
        p.powi(self.n_hashes as i32)
    }

    pub fn size_bytes(&self) -> u64 {
        self.n_bits / 8
    }

    pub fn save(&self, path: &str) -> std::io::Result<()> {
        let mut f = std::fs::File::create(path)?;
        f.write_all(&self.n_bits.to_le_bytes())?;
        f.write_all(&self.n_hashes.to_le_bytes())?;
        f.write_all(&0u32.to_le_bytes())?;
        for word in &self.bits {
            f.write_all(&word.to_le_bytes())?;
        }
        Ok(())
    }

    pub fn load(path: &str) -> std::io::Result<Self> {
        use std::io::Read;
        let mut f = std::fs::File::open(path)?;
        let mut buf8 = [0u8; 8];
        let mut buf4 = [0u8; 4];
        f.read_exact(&mut buf8)?;
        let n_bits = u64::from_le_bytes(buf8);
        f.read_exact(&mut buf4)?;
        let n_hashes = u32::from_le_bytes(buf4);
        f.read_exact(&mut buf4)?;
        let n_words = n_bits.div_ceil(64) as usize;
        let mut bits = vec![0u64; n_words];
        for word in &mut bits {
            f.read_exact(&mut buf8)?;
            *word = u64::from_le_bytes(buf8);
        }
        Ok(Self {
            bits,
            n_bits,
            n_hashes,
        })
    }
}

// --- Count-Min Sketch ---

/// Count-Min Sketch for approximate k-mer frequency estimation.
/// Uses fixed memory regardless of distinct k-mer count.
/// Estimates are always >= true count (one-sided error).
pub struct CountMinSketch {
    counters: Vec<u16>,
    width: u64,
    depth: u32,
}

impl CountMinSketch {
    /// Create a CMS with the given memory budget and number of hash rows.
    pub fn new(size_bytes: u64, depth: u32) -> Self {
        let width = size_bytes / 2 / depth as u64; // u16 = 2 bytes
        log::info!(
            "Count-Min Sketch: {:.1}GB, depth={}, width={}",
            size_bytes as f64 / 1_073_741_824.0,
            depth,
            width
        );
        Self {
            counters: vec![0u16; (width * depth as u64) as usize],
            width,
            depth,
        }
    }

    /// Create a CMS without logging (for use in parallel threads).
    pub fn new_quiet(size_bytes: u64, depth: u32) -> Self {
        let width = size_bytes / 2 / depth as u64;
        Self {
            counters: vec![0u16; (width * depth as u64) as usize],
            width,
            depth,
        }
    }

    #[inline]
    pub fn insert(&mut self, kmer: u64) -> u16 {
        let mut min_count = u16::MAX;
        for row in 0..self.depth {
            let idx = row as u64 * self.width + self.hash(kmer, row);
            let c = &mut self.counters[idx as usize];
            *c = c.saturating_add(1);
            min_count = min_count.min(*c);
        }
        min_count
    }

    #[inline]
    pub fn query(&self, kmer: u64) -> u16 {
        let mut min_count = u16::MAX;
        for row in 0..self.depth {
            let idx = row as u64 * self.width + self.hash(kmer, row);
            min_count = min_count.min(self.counters[idx as usize]);
        }
        min_count
    }

    /// Merge another CMS into this one (element-wise saturating add).
    pub fn merge(&mut self, other: &CountMinSketch) {
        assert_eq!(self.width, other.width);
        assert_eq!(self.depth, other.depth);
        for (a, b) in self.counters.iter_mut().zip(other.counters.iter()) {
            *a = a.saturating_add(*b);
        }
    }

    #[inline]
    fn hash(&self, kmer: u64, row: u32) -> u64 {
        let mixed = match row {
            0 => kmer.wrapping_mul(0x9E3779B97F4A7C15),
            1 => kmer.wrapping_mul(0x517CC1B727220A95),
            2 => kmer.wrapping_mul(0x6C62272E07BB0142),
            3 => kmer.wrapping_mul(0xBF58476D1CE4E5B9),
            _ => kmer
                .wrapping_mul(0x94D049BB133111EB)
                .wrapping_add(row as u64),
        };
        let mixed = mixed ^ (mixed >> 33);
        let mixed = mixed.wrapping_mul(0xFF51AFD7ED558CCD);
        let mixed = mixed ^ (mixed >> 33);
        mixed % self.width
    }
}

// --- Streaming k-mer helpers ---

/// Hash a byte slice (k-mer sequence) to u64 using a fast, well-distributed hash.
/// Used for k > 32 where 2-bit encoding doesn't fit in u64.
#[inline]
fn hash_bytes(data: &[u8]) -> u64 {
    // FNV-1a 64-bit
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Compute a canonical hash for a k-mer sequence (for k > 32).
/// Canonical = min(hash(fwd), hash(revcomp)).
#[inline]
fn canonical_kmer_hash(seq: &[u8]) -> u64 {
    let fwd_hash = hash_bytes(seq);
    // Build reverse complement
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
    let rc_hash = hash_bytes(&rc);
    fwd_hash.min(rc_hash)
}

/// Stream canonical k-mers from a FASTQ/FASTA file, calling `f` for each.
/// For k <= 32, uses efficient 2-bit encoding. For k > 32, uses byte-level hashing.
fn stream_kmers<F>(path: &str, k: usize, mut f: F) -> anyhow::Result<u64>
where
    F: FnMut(u64),
{
    let mut reader = needletail::parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Failed to parse {}: {}", path, e))?;

    let mut total_kmers = 0u64;

    if k <= 32 {
        // Fast path: 2-bit encoding fits in u64
        use crate::minimizer::{complement_2bit, encode_base};
        let mask: u64 = if k < 32 {
            (1u64 << (2 * k)) - 1
        } else {
            u64::MAX
        };

        while let Some(record) = reader.next() {
            let record = record.map_err(|e| anyhow::anyhow!("Read error: {}", e))?;
            let seq = record.seq();
            if seq.len() < k {
                continue;
            }

            let mut fwd: u64 = 0;
            let mut rev: u64 = 0;
            let mut valid = 0usize;

            for &base in seq.iter() {
                if let Some(b) = encode_base(base) {
                    fwd = ((fwd << 2) | b) & mask;
                    rev = (rev >> 2) | (complement_2bit(b) << (2 * (k - 1)));
                    valid += 1;
                    if valid >= k {
                        f(fwd.min(rev));
                        total_kmers += 1;
                    }
                } else {
                    valid = 0;
                    fwd = 0;
                    rev = 0;
                }
            }
        }
    } else {
        // k > 32: slide a window over the sequence, hash each k-mer
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| anyhow::anyhow!("Read error: {}", e))?;
            let seq = record.seq();
            if seq.len() < k {
                continue;
            }

            // Find runs of valid bases (no N) and extract k-mers from each run
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
                        // Process k-mers in the run [run_start, i)
                        let run = &seq[run_start..i];
                        if run.len() >= k {
                            for j in 0..=(run.len() - k) {
                                f(canonical_kmer_hash(&run[j..j + k]));
                                total_kmers += 1;
                            }
                        }
                    }
                    in_run = false;
                }
            }
            // Handle last run
            if in_run {
                let run = &seq[run_start..];
                if run.len() >= k {
                    for j in 0..=(run.len() - k) {
                        f(canonical_kmer_hash(&run[j..j + k]));
                        total_kmers += 1;
                    }
                }
            }
        }
    }

    Ok(total_kmers)
}

/// Write histogram in ntcard format.
pub fn write_histogram(hist: &[u64], path: &str) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "k\tcount\tfrequency")?;
    for (i, &freq) in hist.iter().enumerate() {
        if i == 0 {
            continue;
        }
        if freq > 0 {
            writeln!(f, "\t{}\t{}", i, freq)?;
        }
    }
    Ok(())
}

// --- Histogram mode ---

/// Find the first mode after the minimum in a k-mer frequency histogram.
///
/// Matches the original Python implementation:
///   1. Walk from left until frequency increases (find the error/noise valley)
///   2. From that valley, find the maximum (the coverage mode)
pub fn find_histogram_mode(hist: &[u64]) -> Option<u32> {
    if hist.len() <= 1 {
        return None;
    }
    let freq_count: Vec<u64> = hist[1..].to_vec();
    if freq_count.is_empty() {
        return None;
    }

    // Find the first local minimum
    let mut min_idx = 0;
    let mut min_val = freq_count[0];
    for (idx, &freq) in freq_count.iter().enumerate() {
        if freq > min_val {
            min_idx = if idx > 0 { idx - 1 } else { 0 };
            break;
        }
        min_val = freq;
    }

    // From the minimum, find the maximum (the mode)
    let after_min = &freq_count[min_idx..];
    if after_min.is_empty() {
        return None;
    }

    let max_pos = after_min
        .iter()
        .enumerate()
        .max_by_key(|(_, &v)| v)
        .map(|(i, _)| i)?;

    Some((max_pos + min_idx + 1) as u32)
}

/// Compute the repeat threshold from the histogram mode.
pub fn repeat_threshold(mode: u32, multiplier: u32) -> u32 {
    mode * multiplier
}

// --- Main detect_repeats ---

/// Run the full repeat-detection pipeline on one or more FASTQ/FASTA files.
///
/// Optimized two-pass streaming approach:
///   Pass 1: Stream k-mers into CMS + build histogram inline (parallel across files)
///   Pass 2: Re-stream to collect repetitive k-mers into Bloom filter (parallel across files)
///
/// Files are processed in parallel using std::thread (one CMS per file, merged after).
/// `cms_size_bytes`: CMS memory budget per file (default 4GB)
pub fn detect_repeats(
    paths: &[&str],
    k: usize,
    multiplier: u32,
    bloom_size_bytes: u64,
    cms_size_bytes: u64,
) -> anyhow::Result<(BloomFilter, Vec<u64>)> {
    log::info!(
        "Detecting repetitive k-mers from {} file(s) (k={}, multiplier={}, {:.1}GB CMS/file)",
        paths.len(),
        k,
        multiplier,
        cms_size_bytes as f64 / 1_073_741_824.0
    );

    let max_hist: u16 = 10000;
    let hist_cap = max_hist as usize + 1;
    let cms_depth = 4u32;

    // --- Pass 1: Count k-mers + build histogram (parallel across files) ---
    log::info!(
        "Pass 1: counting k-mers across {} file(s) in parallel...",
        paths.len()
    );

    let owned_paths: Vec<String> = paths.iter().map(|s| s.to_string()).collect();

    // When processing multiple files, we use one shared CMS so that counts
    // from R1 and R2 accumulate correctly. Files are processed sequentially
    // in pass 1 (CMS requires &mut), but pass 2 is fully parallel.
    // For a single file, this is the same as before.
    let mut cms = CountMinSketch::new(cms_size_bytes, cms_depth);
    let mut hist = vec![0u64; hist_cap];
    let mut total_kmers = 0u64;

    for path in paths {
        log::info!("Pass 1: {}", path);
        let n = stream_kmers(path, k, |kmer| {
            let prev = cms.query(kmer);
            let new_est = cms.insert(kmer);
            // Track histogram transitions: move this k-mer from hist[prev] to hist[new_est].
            // At the end, each distinct k-mer sits in hist[its_final_count].
            let new_idx = (new_est as usize).min(hist_cap - 1);
            hist[new_idx] += 1;
            if prev > 0 {
                let old_idx = (prev as usize).min(hist_cap - 1);
                hist[old_idx] = hist[old_idx].saturating_sub(1);
            }
        })?;
        total_kmers += n;
        log::info!("Pass 1 done: {} — {} k-mers", path, n);
    }
    log::info!(
        "Pass 1 complete: {} total k-mers across all files",
        total_kmers
    );

    // Find mode and threshold
    let mode = find_histogram_mode(&hist)
        .ok_or_else(|| anyhow::anyhow!("Could not find mode in k-mer histogram"))?;
    log::info!("K-mer histogram mode: {}", mode);

    let threshold = repeat_threshold(mode, multiplier);
    log::info!(
        "Repeat threshold: {} (mode={} x {})",
        threshold,
        mode,
        multiplier
    );

    // --- Pass 2: Collect repetitive k-mers into Bloom filter (parallel across files) ---
    log::info!("Pass 2: collecting repetitive k-mers into Bloom filter...");
    let threshold_u16 = threshold.min(u16::MAX as u32) as u16;

    // Each thread builds its own BF, then merge via bitwise OR
    let bf_results: Vec<anyhow::Result<BloomFilter>> = std::thread::scope(|s| {
        let handles: Vec<_> = owned_paths
            .iter()
            .map(|path| {
                let cms_ref = &cms;
                let bf_size = bloom_size_bytes;
                s.spawn(move || {
                    let mut bf = BloomFilter::new(bf_size, 1);
                    log::info!("Pass 2: {}", path);
                    stream_kmers(path, k, |kmer| {
                        if cms_ref.query(kmer) >= threshold_u16 {
                            bf.insert(kmer);
                        }
                    })?;
                    log::info!("Pass 2 done: {}", path);
                    Ok(bf)
                })
            })
            .collect();

        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    // Merge BFs via bitwise OR
    let mut final_bf = BloomFilter::new(bloom_size_bytes, 1);
    for result in bf_results {
        let bf = result?;
        for (a, b) in final_bf.bits.iter_mut().zip(bf.bits.iter()) {
            *a |= *b;
        }
    }

    log::info!(
        "Repeat BF: popcount={}, size={} bytes, FPR={:.4}%",
        final_bf.popcount(),
        final_bf.size_bytes(),
        final_bf.fpr() * 100.0
    );

    drop(cms);
    Ok((final_bf, hist))
}

// --- Tests ---

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bloom_filter_basic() {
        let mut bf = BloomFilter::new(1024, 3);
        bf.insert(42);
        bf.insert(100);
        assert!(bf.contains(42));
        assert!(bf.contains(100));
        assert!(!bf.contains(999));
    }

    #[test]
    fn test_bloom_filter_empty() {
        let bf = BloomFilter::new(1024, 3);
        assert!(!bf.contains(42));
        assert_eq!(bf.popcount(), 0);
    }

    #[test]
    fn test_bloom_filter_save_load() {
        let mut bf = BloomFilter::new(256, 2);
        bf.insert(1);
        bf.insert(2);
        bf.insert(3);
        let path = "/tmp/test_bf_cms.bin";
        bf.save(path).unwrap();
        let bf2 = BloomFilter::load(path).unwrap();
        assert!(bf2.contains(1));
        assert!(bf2.contains(2));
        assert!(bf2.contains(3));
        assert_eq!(bf.popcount(), bf2.popcount());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_bloom_filter_fpr() {
        let bf = BloomFilter::new(1024, 3);
        assert_eq!(bf.fpr(), 0.0);
    }

    #[test]
    fn test_cms_basic() {
        let mut cms = CountMinSketch::new(1_000_000, 4);
        for _ in 0..50 {
            cms.insert(42);
        }
        let est = cms.query(42);
        assert!(est >= 50, "Expected >= 50, got {}", est);
        assert!(est < 60, "Estimate too high: {}", est);
        let est_absent = cms.query(999999);
        assert!(est_absent < 5, "Absent k-mer too high: {}", est_absent);
    }

    #[test]
    fn test_cms_multiple_kmers() {
        let mut cms = CountMinSketch::new(1_000_000, 4);
        for _ in 0..100 {
            cms.insert(1);
        }
        for _ in 0..10 {
            cms.insert(2);
        }
        cms.insert(3);
        assert!(cms.query(1) >= 100);
        assert!(cms.query(2) >= 10);
        assert!(cms.query(3) >= 1);
    }

    #[test]
    fn test_find_histogram_mode_simple() {
        let hist = vec![0, 1000, 500, 200, 100, 150, 300, 500, 400, 200, 100];
        assert_eq!(find_histogram_mode(&hist), Some(7));
    }

    #[test]
    fn test_find_histogram_mode_no_valley() {
        let hist = vec![0, 1000, 500, 200, 100, 50];
        assert_eq!(find_histogram_mode(&hist), Some(1));
    }

    #[test]
    fn test_find_histogram_mode_empty() {
        assert_eq!(find_histogram_mode(&[0]), None);
    }

    #[test]
    fn test_repeat_threshold() {
        assert_eq!(repeat_threshold(10, 3), 30);
        assert_eq!(repeat_threshold(25, 3), 75);
    }
}
