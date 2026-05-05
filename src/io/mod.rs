/// I/O routines for reading/writing Physlr data formats.
///
/// Supported formats:
/// - TSV: The Physlr overlap graph format (vertices + edges)
/// - Minimizer TSV: barcode → minimizer hash list
/// - FASTA/FASTQ: sequence files (via needletail)
/// - Path files: one path per line, space-separated vertex names
use crate::graph::NamedGraph;
use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

/// Open a file for reading, transparently handling .gz compression.
pub fn open_reader(path: &str) -> Result<Box<dyn BufRead>> {
    if path == "-" || path == "/dev/stdin" {
        return Ok(Box::new(BufReader::new(io::stdin())));
    }
    let file = File::open(path).with_context(|| format!("Cannot open {}", path))?;
    if path.ends_with(".gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Open a file for writing.
pub fn open_writer(path: &str) -> Result<Box<dyn Write>> {
    if path == "-" || path == "/dev/stdout" {
        return Ok(Box::new(BufWriter::new(io::stdout())));
    }
    let file = File::create(path).with_context(|| format!("Cannot create {}", path))?;
    Ok(Box::new(BufWriter::new(file)))
}

// ─── Minimizer TSV ───────────────────────────────────────────────────────────

/// Read a minimizer TSV file: `barcode\thash1 hash2 hash3 ...`
/// Returns a map from barcode name to set of minimizer hashes.
pub fn read_minimizers(path: &str) -> Result<FxHashMap<String, FxHashSet<u64>>> {
    let reader = open_reader(path)?;
    let mut bx_to_mxs: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();

    for line in reader.lines() {
        let line = line?;
        let mut parts = line.splitn(2, '\t');
        let bx = match parts.next() {
            Some(b) if !b.is_empty() => b.to_string(),
            _ => continue,
        };
        let mxs_str = match parts.next() {
            Some(s) => s,
            None => continue,
        };
        let mxs: FxHashSet<u64> = mxs_str
            .split_whitespace()
            .filter_map(|s| {
                // Handle "hash:pos" format — take only the hash part
                let hash_str = s.split(':').next().unwrap_or(s);
                hash_str.parse::<u64>().ok()
            })
            .collect();
        bx_to_mxs
            .entry(bx)
            .and_modify(|existing| existing.extend(&mxs))
            .or_insert(mxs);
    }
    Ok(bx_to_mxs)
}

/// Read a minimizer TSV file preserving order (list, not set).
pub fn read_minimizers_list(path: &str) -> Result<FxHashMap<String, Vec<u64>>> {
    let reader = open_reader(path)?;
    let mut bx_to_mxs: FxHashMap<String, Vec<u64>> = FxHashMap::default();

    for line in reader.lines() {
        let line = line?;
        let mut parts = line.splitn(2, '\t');
        let bx = match parts.next() {
            Some(b) if !b.is_empty() => b.to_string(),
            _ => continue,
        };
        let mxs_str = match parts.next() {
            Some(s) => s,
            None => continue,
        };
        let mxs: Vec<u64> = mxs_str
            .split_whitespace()
            .filter_map(|s| {
                let hash_str = s.split(':').next().unwrap_or(s);
                hash_str.parse::<u64>().ok()
            })
            .collect();
        bx_to_mxs.insert(bx, mxs);
    }
    Ok(bx_to_mxs)
}

/// Write minimizers in TSV format.
pub fn write_minimizers(
    bx_to_mxs: &FxHashMap<String, FxHashSet<u64>>,
    writer: &mut dyn Write,
) -> Result<()> {
    let mut sorted_bxs: Vec<&String> = bx_to_mxs.keys().collect();
    sorted_bxs.sort();
    for bx in sorted_bxs {
        let mxs = &bx_to_mxs[bx];
        write!(writer, "{}", bx)?;
        let mut first = true;
        for mx in mxs {
            if first {
                write!(writer, "\t{}", mx)?;
                first = false;
            } else {
                write!(writer, " {}", mx)?;
            }
        }
        writeln!(writer)?;
    }
    Ok(())
}

// ─── Overlap Graph TSV ───────────────────────────────────────────────────────

/// Read a Physlr overlap graph in TSV format.
///
/// Format:
/// ```text
/// U\tm
/// barcode1\t150
/// barcode2\t200
///
/// U\tV\tm
/// barcode1\tbarcode2\t45
/// ```
pub fn read_graph_tsv(path: &str) -> Result<NamedGraph> {
    let reader = open_reader(path)?;
    let mut g = NamedGraph::new();
    let mut reading_vertices = true;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Skip headers
        if trimmed == "U\tm" || trimmed == "U\tm\tmol" {
            reading_vertices = true;
            continue;
        }
        if trimmed == "U\tV\tm" {
            reading_vertices = false;
            continue;
        }
        if trimmed.is_empty() {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();

        if reading_vertices {
            if fields.len() >= 2 {
                let name = fields[0];
                let m: u32 = fields[1].parse().unwrap_or(0);
                let idx = g.add_vertex(name, m);
                if fields.len() >= 3 {
                    if let Ok(mol) = fields[2].parse::<u32>() {
                        g.graph[idx].mol = Some(mol);
                    }
                }
            }
        } else if fields.len() >= 3 {
            let u_name = fields[0];
            let v_name = fields[1];
            let m: u32 = fields[2].parse().unwrap_or(0);
            let u = g.add_vertex(u_name, 0);
            let v = g.add_vertex(v_name, 0);
            g.add_edge(u, v, m);
        }
    }
    Ok(g)
}

/// Write a Physlr overlap graph in TSV format.
pub fn write_graph_tsv(g: &NamedGraph, writer: &mut dyn Write) -> Result<()> {
    let has_mol = g.graph.node_indices().any(|n| g.graph[n].mol.is_some());

    if has_mol {
        writeln!(writer, "U\tm\tmol")?;
    } else {
        writeln!(writer, "U\tm")?;
    }

    // Sort vertices by name for deterministic output
    let mut vertices: Vec<_> = g
        .graph
        .node_indices()
        .filter_map(|idx| g.names.get_name(idx).map(|name| (name.to_string(), idx)))
        .collect();
    vertices.sort_by(|a, b| a.0.cmp(&b.0));

    for (name, idx) in &vertices {
        if has_mol {
            writeln!(
                writer,
                "{}\t{}\t{}",
                name,
                g.graph[*idx].m,
                g.graph[*idx].mol.unwrap_or(0)
            )?;
        } else {
            writeln!(writer, "{}\t{}", name, g.graph[*idx].m)?;
        }
    }

    writeln!(writer)?;
    writeln!(writer, "U\tV\tm")?;

    // Sort edges for deterministic output
    let mut edges: Vec<_> = g
        .graph
        .edge_indices()
        .filter_map(|e| {
            let (u, v) = g.graph.edge_endpoints(e).unwrap();
            let u_name = g.names.get_name(u)?;
            let v_name = g.names.get_name(v)?;
            let (u_name, v_name) = if u_name <= v_name {
                (u_name.to_string(), v_name.to_string())
            } else {
                (v_name.to_string(), u_name.to_string())
            };
            Some((u_name, v_name, g.graph[e].m))
        })
        .collect();
    edges.sort();

    for (u, v, m) in &edges {
        writeln!(writer, "{}\t{}\t{}", u, v, m)?;
    }
    Ok(())
}

// ─── Path files ──────────────────────────────────────────────────────────────

/// Read path files. Each line is a space-separated list of vertex names.
pub fn read_paths(path: &str) -> Result<Vec<Vec<String>>> {
    let reader = open_reader(path)?;
    let mut paths = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let path: Vec<String> = trimmed.split_whitespace().map(String::from).collect();
        if !path.is_empty() {
            paths.push(path);
        }
    }
    Ok(paths)
}

/// Write paths to a file.
pub fn write_paths(paths: &[Vec<String>], writer: &mut dyn Write) -> Result<()> {
    for path in paths {
        writeln!(writer, "{}", path.join(" "))?;
    }
    Ok(())
}

// ─── FASTA ───────────────────────────────────────────────────────────────────

/// Read a FASTA file into a map of name → sequence.
pub fn read_fasta(path: &str) -> Result<FxHashMap<String, Vec<u8>>> {
    let mut seqs = FxHashMap::default();
    let mut reader = needletail::parse_fastx_file(path)
        .with_context(|| format!("Cannot open FASTA file: {}", path))?;

    while let Some(record) = reader.next() {
        let record = record?;
        let name = std::str::from_utf8(record.id())
            .unwrap_or("unknown")
            .to_string();
        seqs.insert(name, record.seq().to_vec());
    }
    Ok(seqs)
}

/// Read a FASTA file, returning sequences in order with their names.
pub fn read_fasta_ordered(path: &str) -> Result<Vec<(String, Vec<u8>)>> {
    let mut seqs = Vec::new();
    let mut reader = needletail::parse_fastx_file(path)
        .with_context(|| format!("Cannot open FASTA file: {}", path))?;

    while let Some(record) = reader.next() {
        let record = record?;
        let name = std::str::from_utf8(record.id())
            .unwrap_or("unknown")
            .to_string();
        seqs.push((name, record.seq().to_vec()));
    }
    Ok(seqs)
}

/// Write sequences in FASTA format.
pub fn write_fasta(seqs: &[(String, Vec<u8>)], writer: &mut dyn Write) -> Result<()> {
    for (name, seq) in seqs {
        writeln!(writer, ">{}", name)?;
        writer.write_all(seq)?;
        writeln!(writer)?;
    }
    Ok(())
}

// ─── BED ─────────────────────────────────────────────────────────────────────

/// A BED record for mapping results.
#[derive(Debug, Clone)]
pub struct BedRecord {
    pub tname: String,
    pub tstart: u64,
    pub tend: u64,
    pub qname: String,
    pub score: u32,
    pub orientation: char,
}

/// Read a BED file.
pub fn read_bed(path: &str) -> Result<Vec<BedRecord>> {
    let reader = open_reader(path)?;
    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() >= 5 {
            let orientation = if fields.len() >= 6 {
                fields[5].chars().next().unwrap_or('.')
            } else {
                '.'
            };
            records.push(BedRecord {
                tname: fields[0].to_string(),
                tstart: fields[1].parse().unwrap_or(0),
                tend: fields[2].parse().unwrap_or(0),
                qname: fields[3].to_string(),
                score: fields[4].parse().unwrap_or(0),
                orientation,
            });
        }
    }
    Ok(records)
}
