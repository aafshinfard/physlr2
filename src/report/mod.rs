use serde::Serialize;
use std::io::Write;

/// Assembly metrics.
#[derive(Debug, Clone, Serialize)]
pub struct AssemblyMetrics {
    pub num_sequences: usize,
    pub total_length: u64,
    pub max_length: u64,
    pub min_length: u64,
    pub mean_length: f64,
    pub n50: u64,
    pub n75: u64,
    pub ng50: Option<u64>,
    pub ng75: Option<u64>,
    pub l50: usize,
    pub l75: usize,
}

/// Physical map metrics.
#[derive(Debug, Clone, Serialize)]
pub struct PhysicalMapMetrics {
    pub num_paths: usize,
    pub total_molecules: usize,
    pub max_path_length: usize,
    pub min_path_length: usize,
    pub mean_path_length: f64,
    pub ng50: Option<usize>,
    pub ng75: Option<usize>,
}

/// Compute NGxx metric. `lengths` must be sorted descending.
pub fn compute_ngxx(lengths: &[u64], genome_size: u64, proportion: f64) -> u64 {
    let target = (genome_size as f64 * proportion) as u64;
    let mut running_sum = 0u64;
    for &len in lengths {
        running_sum += len;
        if running_sum >= target {
            return len;
        }
    }
    0
}

/// Compute Nxx metric (using total assembly size as genome size).
pub fn compute_nxx(lengths: &[u64], proportion: f64) -> u64 {
    let total: u64 = lengths.iter().sum();
    compute_ngxx(lengths, total, proportion)
}

/// Compute Lxx metric (number of sequences needed to reach xx% of genome).
pub fn compute_lxx(lengths: &[u64], genome_size: u64, proportion: f64) -> usize {
    let target = (genome_size as f64 * proportion) as u64;
    let mut running_sum = 0u64;
    for (i, &len) in lengths.iter().enumerate() {
        running_sum += len;
        if running_sum >= target {
            return i + 1;
        }
    }
    lengths.len()
}

/// Compute assembly metrics from a set of sequences.
pub fn compute_assembly_metrics(
    seqs: &[(String, Vec<u8>)],
    genome_size: Option<u64>,
) -> AssemblyMetrics {
    let mut lengths: Vec<u64> = seqs.iter().map(|(_, s)| s.len() as u64).collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a));

    let total: u64 = lengths.iter().sum();
    let n50 = compute_nxx(&lengths, 0.5);
    let n75 = compute_nxx(&lengths, 0.75);

    let ng50 = genome_size.map(|g| compute_ngxx(&lengths, g, 0.5));
    let ng75 = genome_size.map(|g| compute_ngxx(&lengths, g, 0.75));

    let l50 = compute_lxx(&lengths, total, 0.5);
    let l75 = compute_lxx(&lengths, total, 0.75);

    AssemblyMetrics {
        num_sequences: lengths.len(),
        total_length: total,
        max_length: *lengths.first().unwrap_or(&0),
        min_length: *lengths.last().unwrap_or(&0),
        mean_length: if lengths.is_empty() {
            0.0
        } else {
            total as f64 / lengths.len() as f64
        },
        n50,
        n75,
        ng50,
        ng75,
        l50,
        l75,
    }
}

/// Compute physical map metrics from backbone paths.
pub fn compute_physical_map_metrics(
    paths: &[Vec<String>],
    expected_molecules: Option<usize>,
) -> PhysicalMapMetrics {
    let mut lengths: Vec<usize> = paths.iter().map(|p| p.len()).collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a));

    let total: usize = lengths.iter().sum();

    let ng50 = expected_molecules.map(|g| {
        let target = (g as f64 * 0.5) as usize;
        let mut running = 0;
        for &len in &lengths {
            running += len;
            if running >= target {
                return len;
            }
        }
        0
    });

    let ng75 = expected_molecules.map(|g| {
        let target = (g as f64 * 0.75) as usize;
        let mut running = 0;
        for &len in &lengths {
            running += len;
            if running >= target {
                return len;
            }
        }
        0
    });

    PhysicalMapMetrics {
        num_paths: lengths.len(),
        total_molecules: total,
        max_path_length: *lengths.first().unwrap_or(&0),
        min_path_length: *lengths.last().unwrap_or(&0),
        mean_path_length: if lengths.is_empty() {
            0.0
        } else {
            total as f64 / lengths.len() as f64
        },
        ng50,
        ng75,
    }
}

/// Write metrics as a TSV table.
pub fn write_metrics_tsv(
    metrics: &AssemblyMetrics,
    label: &str,
    writer: &mut dyn Write,
) -> std::io::Result<()> {
    writeln!(writer, "Assembly\t{}", label)?;
    writeln!(writer, "Sequences\t{}", metrics.num_sequences)?;
    writeln!(writer, "Total length\t{}", metrics.total_length)?;
    writeln!(writer, "Max length\t{}", metrics.max_length)?;
    writeln!(writer, "Min length\t{}", metrics.min_length)?;
    writeln!(writer, "Mean length\t{:.0}", metrics.mean_length)?;
    writeln!(writer, "N50\t{}", metrics.n50)?;
    writeln!(writer, "N75\t{}", metrics.n75)?;
    if let Some(ng50) = metrics.ng50 {
        writeln!(writer, "NG50\t{}", ng50)?;
    }
    if let Some(ng75) = metrics.ng75 {
        writeln!(writer, "NG75\t{}", ng75)?;
    }
    writeln!(writer, "L50\t{}", metrics.l50)?;
    writeln!(writer, "L75\t{}", metrics.l75)?;
    Ok(())
}

/// Write physical map metrics as TSV.
pub fn write_physical_map_metrics_tsv(
    metrics: &PhysicalMapMetrics,
    writer: &mut dyn Write,
) -> std::io::Result<()> {
    writeln!(writer, "Paths\t{}", metrics.num_paths)?;
    writeln!(writer, "Total molecules\t{}", metrics.total_molecules)?;
    writeln!(writer, "Max path length\t{}", metrics.max_path_length)?;
    writeln!(writer, "Min path length\t{}", metrics.min_path_length)?;
    writeln!(writer, "Mean path length\t{:.1}", metrics.mean_path_length)?;
    if let Some(ng50) = metrics.ng50 {
        writeln!(writer, "NG50\t{}", ng50)?;
    }
    if let Some(ng75) = metrics.ng75 {
        writeln!(writer, "NG75\t{}", ng75)?;
    }
    Ok(())
}

/// Generate a DOT (GraphViz) representation of backbone paths for visualization.
pub fn backbone_to_dot(paths: &[Vec<String>], writer: &mut dyn Write) -> std::io::Result<()> {
    writeln!(writer, "graph backbone {{")?;
    writeln!(writer, "  graph [rankdir=LR]")?;
    writeln!(writer, "  node [shape=point width=0.1]")?;
    writeln!(writer, "  edge [color=steelblue]")?;

    for (tid, path) in paths.iter().enumerate() {
        if path.len() < 2 {
            continue;
        }
        writeln!(writer, "  subgraph cluster_{} {{", tid)?;
        writeln!(
            writer,
            "    label=\"Path {} ({} molecules)\"",
            tid,
            path.len()
        )?;
        for i in 0..path.len() - 1 {
            writeln!(writer, "    \"{}\" -- \"{}\"", path[i], path[i + 1])?;
        }
        writeln!(writer, "  }}")?;
    }

    writeln!(writer, "}}")?;
    Ok(())
}

/// Generate a JSON report combining all metrics.
pub fn write_json_report(
    physical_map_metrics: &PhysicalMapMetrics,
    assembly_metrics_before: Option<&AssemblyMetrics>,
    assembly_metrics_after: Option<&AssemblyMetrics>,
    writer: &mut dyn Write,
) -> std::io::Result<()> {
    #[derive(Serialize)]
    struct Report<'a> {
        physical_map: &'a PhysicalMapMetrics,
        #[serde(skip_serializing_if = "Option::is_none")]
        assembly_before: Option<&'a AssemblyMetrics>,
        #[serde(skip_serializing_if = "Option::is_none")]
        assembly_after: Option<&'a AssemblyMetrics>,
    }

    let report = Report {
        physical_map: physical_map_metrics,
        assembly_before: assembly_metrics_before,
        assembly_after: assembly_metrics_after,
    };

    let json = serde_json::to_string_pretty(&report).map_err(std::io::Error::other)?;
    writeln!(writer, "{}", json)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n50_simple() {
        // Input must be sorted descending. [30, 20, 10], total=60, need >= 30
        // Cumulative: 30 → N50 = 30
        let n50 = compute_nxx(&[30, 20, 10], 0.5);
        assert_eq!(n50, 30);
    }

    #[test]
    fn test_n50_equal() {
        let n50 = compute_nxx(&[100, 100, 100], 0.5);
        assert_eq!(n50, 100);
    }

    #[test]
    fn test_n50_single() {
        let n50 = compute_nxx(&[42], 0.5);
        assert_eq!(n50, 42);
    }

    #[test]
    fn test_n50_empty() {
        let n50 = compute_nxx(&[], 0.5);
        assert_eq!(n50, 0);
    }

    #[test]
    fn test_ng50() {
        // Sorted desc: [300, 200, 100], genome_size=1000
        // Need cumulative >= 500. 300, 500 → NG50 = 200
        let ng50 = compute_ngxx(&[300, 200, 100], 1000, 0.5);
        assert_eq!(ng50, 200);
    }

    #[test]
    fn test_ng50_genome_smaller() {
        // Sorted desc: [300, 200, 100], genome_size=100
        // Need cumulative >= 50. 300 >= 50 → NG50 = 300
        let ng50 = compute_ngxx(&[300, 200, 100], 100, 0.5);
        assert_eq!(ng50, 300);
    }

    #[test]
    fn test_ng50_empty() {
        let ng50 = compute_ngxx(&[], 1000, 0.5);
        assert_eq!(ng50, 0);
    }

    #[test]
    fn test_lxx() {
        // Sorted desc: [300, 200, 100], genome_size=1000
        // Need cumulative >= 500. 300(1), 500(2) → L50 = 2
        let l50 = compute_lxx(&[300, 200, 100], 1000, 0.5);
        assert_eq!(l50, 2);
    }

    #[test]
    fn test_lxx_empty() {
        let l50 = compute_lxx(&[], 1000, 0.5);
        assert_eq!(l50, 0);
    }

    #[test]
    fn test_n90() {
        // Sorted desc: [50, 40, 30, 20, 10], total=150
        // Need cumulative >= 135. 50, 90, 120, 140 → N90 = 20
        let n90 = compute_nxx(&[50, 40, 30, 20, 10], 0.9);
        assert_eq!(n90, 20);
    }

    #[test]
    fn test_backbone_to_dot() {
        let paths = vec![
            vec!["a".to_string(), "b".to_string(), "c".to_string()],
            vec!["d".to_string(), "e".to_string()],
        ];
        let mut buf = Vec::new();
        backbone_to_dot(&paths, &mut buf).unwrap();
        let dot = String::from_utf8(buf).unwrap();
        assert!(
            dot.contains("graph backbone"),
            "Expected 'graph backbone' in DOT output"
        );
        assert!(dot.contains("a"), "Expected node 'a' in DOT output");
        assert!(dot.contains("--"), "Expected '--' edges in DOT output");
    }
}
