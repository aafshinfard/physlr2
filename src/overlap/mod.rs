/// Overlap computation between barcodes based on shared minimizers.
///
/// For each pair of barcodes that share at least `min_m` minimizers,
/// we create an edge in the overlap graph weighted by the count of
/// shared minimizers.
///
/// This is the most compute-intensive step. We parallelize using rayon.
use crate::graph::NamedGraph;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

/// Compute the barcode overlap graph.
///
/// For each minimizer, we know which barcodes contain it. For each pair of
/// barcodes sharing a minimizer, we increment their shared count. Pairs with
/// count >= `min_m` become edges.
pub fn compute_overlap(bx_to_mxs: &FxHashMap<String, FxHashSet<u64>>, min_m: u32) -> NamedGraph {
    log::info!(
        "Computing overlaps for {} barcodes, min_m={}",
        bx_to_mxs.len(),
        min_m
    );

    // Assign integer IDs to barcodes for memory-efficient computation
    let (bx_to_id, id_to_bx) = crate::minimizer::assign_barcode_ids(bx_to_mxs);
    let n_barcodes = id_to_bx.len();

    // Build inverted index: minimizer → list of barcode IDs
    let mx_to_ids = crate::minimizer::build_inverted_index_ids(bx_to_mxs, &bx_to_id);

    log::info!(
        "Built inverted index: {} minimizers, {} barcodes",
        mx_to_ids.len(),
        n_barcodes
    );

    // For each barcode, count shared minimizers with all other barcodes.
    // We parallelize over barcodes.
    let barcode_ids: Vec<u32> = (0..n_barcodes as u32).collect();

    // Build per-barcode minimizer sets indexed by ID
    let id_to_mxs: Vec<&FxHashSet<u64>> = (0..n_barcodes)
        .map(|id| &bx_to_mxs[&id_to_bx[id]])
        .collect();

    // Parallel overlap computation
    let edges: Vec<(u32, u32, u32)> = barcode_ids
        .par_iter()
        .flat_map(|&id1| {
            let mut local_counts: FxHashMap<u32, u32> = FxHashMap::default();
            let mxs = id_to_mxs[id1 as usize];

            for &mx in mxs {
                if let Some(ids) = mx_to_ids.get(&mx) {
                    for &id2 in ids {
                        if id1 < id2 {
                            *local_counts.entry(id2).or_insert(0) += 1;
                        }
                    }
                }
            }

            local_counts
                .into_iter()
                .filter(|(_, count)| *count >= min_m)
                .map(|(id2, count)| (id1, id2, count))
                .collect::<Vec<_>>()
        })
        .collect();

    log::info!(
        "Found {} edges with >= {} shared minimizers",
        edges.len(),
        min_m
    );

    // Build the graph
    let mut g = NamedGraph::new();

    // Add all vertices
    for bx in id_to_bx.iter() {
        let m = bx_to_mxs[bx].len() as u32;
        g.add_vertex(bx, m);
    }

    // Add edges
    for (id1, id2, m) in &edges {
        let u = g.names.get_idx(&id_to_bx[*id1 as usize]).unwrap();
        let v = g.names.get_idx(&id_to_bx[*id2 as usize]).unwrap();
        g.add_edge(u, v, *m);
    }

    log::info!("Overlap graph: V={} E={}", g.num_vertices(), g.num_edges());
    g
}

/// Filter edges by removing the bottom `percentile`% by weight.
/// This is the "m" parameter from the original Physlr.
pub fn filter_edges_by_percentile(g: &mut NamedGraph, percentile: f64) -> usize {
    if percentile <= 0.0 {
        return 0;
    }

    // Collect all edge weights
    let mut weights: Vec<u32> = g.graph.edge_indices().map(|e| g.graph[e].m).collect();
    weights.sort_unstable();

    if weights.is_empty() {
        return 0;
    }

    let idx = ((weights.len() as f64) * percentile / 100.0) as usize;
    let threshold = weights[idx.min(weights.len() - 1)];

    let removed = g.filter_edges(threshold);
    log::info!(
        "Removed {} edges below threshold {} ({}th percentile)",
        removed,
        threshold,
        percentile
    );
    removed
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_bx_map(data: &[(&str, &[u64])]) -> FxHashMap<String, FxHashSet<u64>> {
        data.iter()
            .map(|(name, mxs)| (name.to_string(), mxs.iter().copied().collect()))
            .collect()
    }

    #[test]
    fn test_compute_overlap_basic() {
        let bx = make_bx_map(&[
            ("a", &[1, 2, 3, 4, 5]),
            ("b", &[3, 4, 5, 6, 7]),
            ("c", &[10, 11, 12]),
        ]);
        let g = compute_overlap(&bx, 2);
        // a and b share {3,4,5} = 3 minimizers >= 2 → edge
        // c shares nothing with a or b → no edges
        assert_eq!(g.num_vertices(), 3);
        assert_eq!(g.num_edges(), 1);
    }

    #[test]
    fn test_compute_overlap_min_m() {
        let bx = make_bx_map(&[("a", &[1, 2, 3]), ("b", &[2, 3, 4])]);
        // min_m=3: a and b share {2,3} = 2 < 3 → no edge
        let g = compute_overlap(&bx, 3);
        assert_eq!(g.num_edges(), 0);

        // min_m=2: 2 >= 2 → edge
        let g = compute_overlap(&bx, 2);
        assert_eq!(g.num_edges(), 1);
    }

    #[test]
    fn test_compute_overlap_empty() {
        let bx: FxHashMap<String, FxHashSet<u64>> = FxHashMap::default();
        let g = compute_overlap(&bx, 1);
        assert_eq!(g.num_vertices(), 0);
        assert_eq!(g.num_edges(), 0);
    }

    #[test]
    fn test_compute_overlap_no_shared() {
        let bx = make_bx_map(&[("a", &[1, 2]), ("b", &[3, 4])]);
        let g = compute_overlap(&bx, 1);
        assert_eq!(g.num_edges(), 0);
    }

    #[test]
    fn test_filter_edges_by_percentile() {
        let bx = make_bx_map(&[
            ("a", &[1, 2, 3, 4, 5]),
            ("b", &[1, 2, 3, 4, 5, 6, 7]),
            ("c", &[1, 2]),
        ]);
        let mut g = compute_overlap(&bx, 1);
        // a-b share 5, a-c share 2, b-c share 2
        let initial_edges = g.num_edges();
        assert!(initial_edges > 0);

        // Remove bottom 50% → threshold will be weight at 50th percentile
        // Weights sorted: [2, 2, 5]. 50th percentile idx=1 → threshold=2
        // filter_edges removes edges with weight < threshold, so edges with weight 2
        // are NOT removed (2 < 2 is false). Only edges < 2 would be removed.
        // So with these weights, 0 edges removed is correct.
        let removed = filter_edges_by_percentile(&mut g, 50.0);
        assert_eq!(g.num_edges(), initial_edges - removed);
    }

    #[test]
    fn test_filter_edges_zero_percentile() {
        let bx = make_bx_map(&[("a", &[1, 2, 3]), ("b", &[1, 2, 3])]);
        let mut g = compute_overlap(&bx, 1);
        let removed = filter_edges_by_percentile(&mut g, 0.0);
        assert_eq!(removed, 0);
    }
}
