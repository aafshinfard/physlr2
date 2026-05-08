use rustc_hash::{FxHashMap, FxHashSet};

/// Configuration for scaffolding.
#[derive(Debug, Clone)]
pub struct ScaffoldConfig {
    /// Gap size (number of Ns) between contigs in scaffolds.
    pub gap_size: usize,
    /// Minimum mapping score to use a contig.
    pub min_score: u32,
    /// Minimum scaffold length to output.
    pub min_length: usize,
}

impl Default for ScaffoldConfig {
    fn default() -> Self {
        Self {
            gap_size: 100,
            min_score: 10,
            min_length: 0,
        }
    }
}

/// Complement table for DNA sequences.
const COMPLEMENT: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table[b'N' as usize] = b'N';
    table[b'n' as usize] = b'n';
    table
};

/// Reverse complement a DNA sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| COMPLEMENT[b as usize]).collect()
}

/// Get a sequence in the specified orientation.
fn get_oriented_sequence(
    seqs: &FxHashMap<String, Vec<u8>>,
    name: &str,
    orientation: char,
) -> Option<Vec<u8>> {
    let seq = seqs.get(name)?;
    match orientation {
        '+' => Some(seq.clone()),
        '-' => Some(reverse_complement(seq)),
        _ => Some(seq.clone()), // Unknown orientation — use forward
    }
}

/// Scaffold a draft assembly using the physical map.
///
/// Returns scaffolded sequences as (name, sequence) pairs.
pub fn scaffold_assembly(
    draft_seqs: &FxHashMap<String, Vec<u8>>,
    paths: &[Vec<(String, char)>],
    config: &ScaffoldConfig,
) -> Vec<(String, Vec<u8>)> {
    let gap: Vec<u8> = vec![b'N'; config.gap_size];
    let mut scaffolds = Vec::new();
    let mut used_seqs: FxHashSet<String> = FxHashSet::default();
    let mut scaffold_num = 0;

    for path in paths {
        if path.is_empty() {
            continue;
        }

        // Skip paths where ALL contigs are unoriented (matching original)
        if path.iter().all(|(_, ori)| *ori == '.') {
            continue;
        }

        let mut seq = Vec::new();
        let mut n_contigs = 0;

        for (name, orientation) in path {
            if let Some(contig_seq) = get_oriented_sequence(draft_seqs, name, *orientation) {
                if !seq.is_empty() {
                    seq.extend_from_slice(&gap);
                }
                seq.extend_from_slice(&contig_seq);
                used_seqs.insert(name.clone());
                n_contigs += 1;
            } else if *orientation == '.' {
                // Unoriented — try to include as Ns
                if let Some(contig_seq) = draft_seqs.get(name) {
                    if !seq.is_empty() {
                        seq.extend_from_slice(&gap);
                    }
                    seq.extend(std::iter::repeat_n(b'N', contig_seq.len()));
                    used_seqs.insert(name.clone());
                    n_contigs += 1;
                }
            }
        }

        if seq.len() >= config.min_length && n_contigs > 0 {
            scaffold_num += 1;
            let name = format!("{:07} LN:i:{} xn:i:{}", scaffold_num, seq.len(), n_contigs);
            scaffolds.push((name, seq));
        }
    }

    // Add unused sequences as singletons
    let mut unused: Vec<(&String, &Vec<u8>)> = draft_seqs
        .iter()
        .filter(|(name, _)| !used_seqs.contains(*name))
        .collect();
    unused.sort_by_key(|(name, _)| (*name).clone());

    for (_name, seq) in unused {
        if seq.len() >= config.min_length {
            scaffold_num += 1;
            let scaffold_name = format!("{:07} LN:i:{} xn:i:1", scaffold_num, seq.len());
            scaffolds.push((scaffold_name, seq.clone()));
        }
    }

    log::info!(
        "Produced {} scaffolds ({} from paths, {} singletons)",
        scaffolds.len(),
        paths.len(),
        scaffolds.len() - paths.len()
    );
    scaffolds
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ACGTACGT"), b"ACGTACGT");
    }

    #[test]
    fn test_reverse_complement_single() {
        assert_eq!(reverse_complement(b"A"), b"T");
        assert_eq!(reverse_complement(b"C"), b"G");
    }

    #[test]
    fn test_reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn test_reverse_complement_lowercase() {
        assert_eq!(reverse_complement(b"acgt"), b"acgt");
    }

    #[test]
    fn test_reverse_complement_with_n() {
        assert_eq!(reverse_complement(b"ACNGT"), b"ACNGT");
    }

    #[test]
    fn test_reverse_complement_asymmetric() {
        assert_eq!(reverse_complement(b"AACG"), b"CGTT");
    }

    #[test]
    fn test_scaffold_config_default() {
        let config = ScaffoldConfig::default();
        assert_eq!(config.gap_size, 100);
        assert_eq!(config.min_score, 10);
        assert_eq!(config.min_length, 0);
    }

    #[test]
    fn test_scaffold_assembly_empty() {
        let contigs = FxHashMap::default();
        let paths: Vec<Vec<(String, char)>> = Vec::new();
        let config = ScaffoldConfig::default();
        let result = scaffold_assembly(&contigs, &paths, &config);
        assert!(result.is_empty());
    }

    #[test]
    fn test_scaffold_assembly_single_contig() {
        let mut contigs = FxHashMap::default();
        contigs.insert("ctg1".to_string(), b"ACGTACGT".to_vec());
        let paths = vec![vec![("ctg1".to_string(), '+')]];
        let config = ScaffoldConfig {
            gap_size: 5,
            min_score: 0,
            min_length: 0,
        };
        let result = scaffold_assembly(&contigs, &paths, &config);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].1, b"ACGTACGT");
    }

    #[test]
    fn test_scaffold_assembly_two_contigs_fwd() {
        let mut contigs = FxHashMap::default();
        contigs.insert("a".to_string(), b"AAAA".to_vec());
        contigs.insert("b".to_string(), b"CCCC".to_vec());
        let paths = vec![vec![("a".to_string(), '+'), ("b".to_string(), '+')]];
        let config = ScaffoldConfig {
            gap_size: 3,
            min_score: 0,
            min_length: 0,
        };
        let result = scaffold_assembly(&contigs, &paths, &config);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].1, b"AAAANNNCCCC");
    }

    #[test]
    fn test_scaffold_assembly_reverse() {
        let mut contigs = FxHashMap::default();
        contigs.insert("a".to_string(), b"AACG".to_vec());
        let paths = vec![vec![("a".to_string(), '-')]];
        let config = ScaffoldConfig {
            gap_size: 3,
            min_score: 0,
            min_length: 0,
        };
        let result = scaffold_assembly(&contigs, &paths, &config);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].1, reverse_complement(b"AACG"));
    }
}
