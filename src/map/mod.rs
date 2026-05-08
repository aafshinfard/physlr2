/// Mapping sequences to the physical map.
///
/// Matches the original Physlr `physlr_map`, `physlr_map_paf`, and `physlr_bed_to_path` logic.
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::Write;

/// A mapping of a query sequence to a backbone position.
#[derive(Debug, Clone)]
pub struct Mapping {
    pub tid: usize,
    pub tpos: usize,
    pub qname: String,
    pub score: u32,
    pub orientation: char,
}

/// Index minimizer positions along backbone paths.
///
/// Matches original `index_minimizers_in_backbones`:
///   for tid, path in enumerate(backbones):
///     for tpos, bx in enumerate(path):
///       for mx in bxtomxs[bx]:
///         mxtopos[mx].append((tid, tpos))
pub fn index_backbone_minimizers(
    backbones: &[Vec<String>],
    mol_to_mxs: &FxHashMap<String, FxHashSet<u64>>,
) -> FxHashMap<u64, Vec<(usize, usize)>> {
    let mut mx_to_pos: FxHashMap<u64, Vec<(usize, usize)>> = FxHashMap::default();

    for (tid, path) in backbones.iter().enumerate() {
        for (pos, mol_name) in path.iter().enumerate() {
            // Try the molecule name directly, then strip _N suffixes
            let mxs = mol_to_mxs
                .get(mol_name)
                .or_else(|| {
                    mol_name
                        .rsplit_once('_')
                        .and_then(|(base, _)| mol_to_mxs.get(base))
                })
                .or_else(|| {
                    mol_name
                        .rsplit_once('_')
                        .and_then(|(base, _)| base.rsplit_once('_'))
                        .and_then(|(base, _)| mol_to_mxs.get(base))
                });

            if let Some(mxs) = mxs {
                for &mx in mxs {
                    mx_to_pos.entry(mx).or_default().push((tid, pos));
                }
            }
        }
    }

    log::info!("Indexed {} minimizers from backbone paths", mx_to_pos.len());
    mx_to_pos
}

/// Determine orientation from three median query positions.
///
/// Matches original `determine_orientation(x, y, z)`:
///   if x is not None and z is not None:
///     return "." if x == y == z else "+" if x <= y <= z else "-" if x >= y >= z else "."
///   if x is not None:
///     return "." if x == y else "+" if x < y else "-"
///   if z is not None:
///     return "." if y == z else "+" if y < z else "-"
///   return "."
fn determine_orientation(x: Option<i64>, y: Option<i64>, z: Option<i64>) -> char {
    let y = match y {
        Some(v) => v,
        None => return '.',
    };

    match (x, z) {
        (Some(xv), Some(zv)) => {
            if xv == y && y == zv {
                '.'
            } else if xv <= y && y <= zv {
                '+'
            } else if xv >= y && y >= zv {
                '-'
            } else {
                '.'
            }
        }
        (Some(xv), None) => {
            if xv == y {
                '.'
            } else if xv < y {
                '+'
            } else {
                '-'
            }
        }
        (None, Some(zv)) => {
            if y == zv {
                '.'
            } else if y < zv {
                '+'
            } else {
                '-'
            }
        }
        (None, None) => '.',
    }
}

/// Compute median_low (lower median) of a slice.
/// Matches Python's `statistics.median_low`.
fn median_low(values: &[i64]) -> i64 {
    if values.is_empty() {
        return 0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    sorted[(sorted.len() - 1) / 2]
}

/// Map query sequences to backbone paths.
///
/// Matches original `physlr_map` with `--map-pos 1` (default):
///   For each query:
///     1. Build tidpos_to_qpos: for each (qpos, mx), for each tidpos in mxtopos[mx],
///        collect qpos values, then take median_low
///     2. Count tidpos_to_n: number of minimizers at each tidpos
///     3. For each tidpos with score >= min_score:
///        orientation = determine_orientation(qpos[tid,tpos-1], qpos[tid,tpos], qpos[tid,tpos+1])
pub fn map_to_backbone(
    query_mxs: &FxHashMap<String, Vec<u64>>,
    mx_to_pos: &FxHashMap<u64, Vec<(usize, usize)>>,
    min_score: u32,
) -> Vec<Mapping> {
    let mut mappings = Vec::new();
    let mut num_mapped = 0;

    for (qname, mxs) in query_mxs {
        // Map each target position to a list of query positions
        let mut tidpos_to_qpos_list: FxHashMap<(usize, usize), Vec<i64>> = FxHashMap::default();
        for (qpos, &mx) in mxs.iter().enumerate() {
            if let Some(positions) = mx_to_pos.get(&mx) {
                for &tidpos in positions {
                    tidpos_to_qpos_list
                        .entry(tidpos)
                        .or_default()
                        .push(qpos as i64);
                }
            }
        }

        // Compute median_low for each tidpos
        let tidpos_to_qpos: FxHashMap<(usize, usize), i64> = tidpos_to_qpos_list
            .iter()
            .map(|(&tidpos, qpos_list)| (tidpos, median_low(qpos_list)))
            .collect();

        // Count minimizers at each tidpos
        let mut tidpos_to_n: FxHashMap<(usize, usize), u32> = FxHashMap::default();
        for &mx in mxs {
            if let Some(positions) = mx_to_pos.get(&mx) {
                for &tidpos in positions {
                    *tidpos_to_n.entry(tidpos).or_insert(0) += 1;
                }
            }
        }

        let mut mapped = false;
        for (&(tid, tpos), &score) in &tidpos_to_n {
            if score >= min_score {
                mapped = true;

                // Determine orientation using map_pos=1 (default)
                let x = tidpos_to_qpos.get(&(tid, tpos.wrapping_sub(1))).copied();
                let y = tidpos_to_qpos.get(&(tid, tpos)).copied();
                let z = tidpos_to_qpos.get(&(tid, tpos + 1)).copied();
                let orientation = determine_orientation(x, y, z);

                mappings.push(Mapping {
                    tid,
                    tpos,
                    qname: qname.clone(),
                    score,
                    orientation,
                });
            }
        }
        if mapped {
            num_mapped += 1;
        }
    }

    log::info!(
        "Mapped {} of {} query sequences",
        num_mapped,
        query_mxs.len()
    );
    mappings
}

/// Convert BED mappings to scaffold paths.
///
/// Matches original `physlr_bed_to_path`:
///   1. Group by qname -> {tname -> [(tstart, orientation), ...]}
///   2. For each qname, pick tname with most positions
///   3. tstart = median_low of all tstart values for that tname
///   4. orientation = most_common orientation; if tied, "."
///   5. Sort by (tname, tstart)
///   6. Group consecutive same-tname entries into scaffold paths
pub fn bed_to_scaffold_paths(mappings: &[Mapping], min_score: u32) -> Vec<Vec<(String, char)>> {
    // Step 1: Group by qname -> {tname -> [(tstart, orientation)]}
    let mut qnames: FxHashMap<String, FxHashMap<usize, Vec<(usize, char)>>> = FxHashMap::default();
    for m in mappings {
        if m.score < min_score {
            continue;
        }
        qnames
            .entry(m.qname.clone())
            .or_default()
            .entry(m.tid)
            .or_default()
            .push((m.tpos, m.orientation));
    }

    // Step 2-4: For each qname, pick best target, compute median position and voted orientation
    let mut scaffolds: Vec<(usize, i64, char, String)> = Vec::new();
    for (qname, targets) in &qnames {
        // Pick target with most positions (matching: max by len(positions))
        let (best_tid, positions) = targets
            .iter()
            .max_by_key(|(_, positions)| positions.len())
            .unwrap();

        // median_low of tstart values
        let tstarts: Vec<i64> = positions.iter().map(|(ts, _)| *ts as i64).collect();
        let tstart = median_low(&tstarts);

        // Orientation voting: most_common; if tied, "."
        let mut ori_counts: FxHashMap<char, usize> = FxHashMap::default();
        for (_, ori) in positions {
            *ori_counts.entry(*ori).or_insert(0) += 1;
        }
        let mut ori_vec: Vec<(char, usize)> = ori_counts.into_iter().collect();
        ori_vec.sort_by_key(|x| std::cmp::Reverse(x.1));

        let orientation = if ori_vec.len() == 1 || ori_vec[0].1 > ori_vec[1].1 {
            ori_vec[0].0
        } else {
            '.'
        };

        scaffolds.push((*best_tid, tstart, orientation, qname.clone()));
    }

    // Step 5: Sort by (tname, tstart)
    scaffolds.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    // Step 6: Group consecutive same-tname entries into paths
    let mut paths: Vec<Vec<(String, char)>> = Vec::new();
    let mut current_path: Vec<(String, char)> = Vec::new();
    let mut prev_tid: Option<usize> = None;

    for (tid, _, orientation, qname) in &scaffolds {
        if prev_tid.is_some() && prev_tid != Some(*tid) && !current_path.is_empty() {
            paths.push(std::mem::take(&mut current_path));
        }
        current_path.push((qname.clone(), *orientation));
        prev_tid = Some(*tid);
    }
    if !current_path.is_empty() {
        paths.push(current_path);
    }

    log::info!(
        "Produced {} scaffold paths from {} contigs",
        paths.len(),
        scaffolds.len()
    );
    paths
}

// ---------------------------------------------------------------------------
// PAF mapping (backbone vs reference visualization)
// ---------------------------------------------------------------------------

/// A PAF record for backbone-to-reference mapping.
#[derive(Debug)]
pub struct PafRecord {
    pub qname: String,
    pub qlength: usize,
    pub qstart: usize,
    pub qend: usize,
    pub orientation: char,
    pub tname: String,
    pub tlength: usize,
    pub tstart: usize,
    pub tend: usize,
    pub score: u32,
    pub length: usize,
    pub mapq: u32,
}

impl PafRecord {
    pub fn write_to<W: Write + ?Sized>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.qname,
            self.qlength,
            self.qstart,
            self.qend,
            self.orientation,
            self.tname,
            self.tlength,
            self.tstart,
            self.tend,
            self.score,
            self.length,
            self.mapq
        )
    }
}

/// Compute quantile values. Uses nearest-rank method.
fn quantile(probs: &[f64], data: &[usize]) -> Vec<usize> {
    if data.is_empty() {
        return probs.iter().map(|_| 0).collect();
    }
    let mut sorted = data.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    probs
        .iter()
        .map(|&p| {
            let idx = ((p * (n as f64 - 1.0)).round() as usize).min(n - 1);
            sorted[idx]
        })
        .collect()
}

/// Map query sequences to backbone paths, output PAF records.
/// Faithful reimplementation of original `physlr_map_paf`.
pub fn map_paf(
    query_mxs: &FxHashMap<String, Vec<u64>>,
    mx_to_pos: &FxHashMap<u64, Vec<(usize, usize)>>,
    backbones: &[Vec<String>],
    min_score: u32,
    coef: f64,
) -> Vec<PafRecord> {
    let mut records = Vec::new();
    let mut num_mapped = 0u64;

    for (qname, mxs) in query_mxs {
        let mut tidpos_to_qpos: FxHashMap<(usize, usize), Vec<usize>> = FxHashMap::default();
        for (qpos, &mx) in mxs.iter().enumerate() {
            if let Some(positions) = mx_to_pos.get(&mx) {
                for &tidpos in positions {
                    tidpos_to_qpos.entry(tidpos).or_default().push(qpos);
                }
            }
        }

        let mut tidpos_bounds: FxHashMap<(usize, usize), (usize, usize, usize)> =
            FxHashMap::default();
        for (&tidpos, qpos_list) in &tidpos_to_qpos {
            let q = quantile(&[0.0, 0.25, 0.5, 0.75, 1.0], qpos_list);
            let iqr = q[3].saturating_sub(q[1]);
            let wr = (coef * iqr as f64) as usize;
            let low = q[0].max(q[1].saturating_sub(wr));
            let high = q[4].min(q[3] + wr);
            tidpos_bounds.insert(tidpos, (low, q[2], high));
        }

        let mut tidpos_to_n: FxHashMap<(usize, usize), u32> = FxHashMap::default();
        for &mx in mxs {
            if let Some(positions) = mx_to_pos.get(&mx) {
                for &tidpos in positions {
                    *tidpos_to_n.entry(tidpos).or_insert(0) += 1;
                }
            }
        }

        let qlength = mxs.len();
        let mut mapped = false;

        for (&(tid, tpos), &score) in &tidpos_to_n {
            if score >= min_score {
                mapped = true;
                let (qstart, _qmed, qend) = tidpos_bounds[&(tid, tpos)];
                let qmed_before = tidpos_bounds
                    .get(&(tid, tpos.wrapping_sub(1)))
                    .map(|b| b.1 as i64);
                let qmed = Some(tidpos_bounds[&(tid, tpos)].1 as i64);
                let qmed_after = tidpos_bounds.get(&(tid, tpos + 1)).map(|b| b.1 as i64);
                let orientation = determine_orientation(qmed_before, qmed, qmed_after);
                let tlength = backbones.get(tid).map_or(0, |p| p.len());
                let length = if qend > qstart { qend - qstart } else { 1 };
                let mapq = (100 * score as usize / length).min(255) as u32;

                records.push(PafRecord {
                    qname: qname.clone(),
                    qlength,
                    qstart,
                    qend,
                    orientation,
                    tname: format!("{}", tid),
                    tlength,
                    tstart: tpos,
                    tend: tpos + 1,
                    score,
                    length,
                    mapq,
                });
            }
        }
        if mapped {
            num_mapped += 1;
        }
    }

    log::info!(
        "Mapped {} of {} query sequences to backbone",
        num_mapped,
        query_mxs.len()
    );
    records
}

/// Read backbone paths from a .backbone.tsv file.
pub fn read_backbone_paths(path: &str) -> anyhow::Result<Vec<Vec<String>>> {
    use std::io::BufRead;
    let reader = crate::io::open_reader(path)?;
    let mut paths = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let nodes: Vec<String> = line.split_whitespace().map(|s| s.to_string()).collect();
        if !nodes.is_empty() {
            paths.push(nodes);
        }
    }
    log::info!("Read {} backbone paths from {}", paths.len(), path);
    Ok(paths)
}

/// Filter PAF records to keep only the top N backbone paths by total score.
pub fn filter_paf_top_n(records: &[PafRecord], n: usize) -> Vec<usize> {
    let mut path_scores: FxHashMap<&str, u32> = FxHashMap::default();
    for r in records {
        *path_scores.entry(&r.tname).or_insert(0) += r.score;
    }
    let mut scored: Vec<(&str, u32)> = path_scores.into_iter().collect();
    scored.sort_by_key(|&(_, s)| std::cmp::Reverse(s));
    let top: FxHashSet<&str> = scored.iter().take(n).map(|&(name, _)| name).collect();
    records
        .iter()
        .enumerate()
        .filter(|(_, r)| top.contains(r.tname.as_str()))
        .map(|(i, _)| i)
        .collect()
}
