/// Backbone path extraction from the overlap graph.
///
/// Matches the original Physlr `determine_backbones` flow:
///   1. Copy graph
///   2. If prune_bridges > 0: remove_bridges iteratively
///      - Each iteration: compute MST+prune, find bridges, remove from original graph
///   3. Compute MST+prune
///   4. While MST not empty:
///      a. Extract backbones from trees (with junction splitting)
///      b. Keep paths >= prune_branches
///      c. Remove those nodes from MST
///      d. Repeat
///   5. Sort by length descending
use crate::graph::{
    connected_components, diameter_path, maximum_spanning_tree, measure_branch_lengths,
    prune_branches, NamedGraph, OverlapGraph,
};
use petgraph::stable_graph::NodeIndex;
use petgraph::visit::EdgeRef;
use rustc_hash::{FxHashMap, FxHashSet};

/// Configuration for backbone extraction.
#[derive(Debug, Clone)]
pub struct BackboneConfig {
    /// Minimum branch size to keep during pruning.
    pub prune_branch_size: usize,
    /// Minimum bridge size to keep.
    pub prune_bridge_size: usize,
    /// Minimum branch size at junctions to trigger splitting.
    pub prune_junction_size: usize,
    /// Minimum path length to report.
    pub min_path_size: usize,
}

impl Default for BackboneConfig {
    fn default() -> Self {
        Self {
            prune_branch_size: 10,
            prune_bridge_size: 10,
            prune_junction_size: 200,
            min_path_size: 50,
        }
    }
}

/// Compute a pruned MST of the given graph.
/// Matches original `determine_pruned_mst`: MST + iterative branch pruning.
fn determine_pruned_mst(
    graph: &OverlapGraph,
    prune_branch_size: usize,
) -> (OverlapGraph, FxHashMap<NodeIndex, NodeIndex>) {
    let (mut mst, node_map) = maximum_spanning_tree(graph);
    prune_branches(&mut mst, prune_branch_size);
    (mst, node_map)
}

/// Identify bridges in the MST.
///
/// Matches original `identify_bridges`:
///   1. Find junctions (degree >= 3)
///   2. Remove junctions, find connected components
///   3. Bridges = short paths where all nodes have degree 2 in the MST
///   4. Also include edges between junctions
fn identify_bridges(mst: &OverlapGraph, bridge_length: usize) -> Vec<Vec<NodeIndex>> {
    // Find junctions
    let junctions: FxHashSet<NodeIndex> = mst
        .node_indices()
        .filter(|&n| mst.neighbors(n).count() >= 3)
        .collect();

    // Remove junctions and find connected components (contiguous paths)
    let remaining: FxHashSet<NodeIndex> = mst
        .node_indices()
        .filter(|n| !junctions.contains(n))
        .collect();

    let mut visited = FxHashSet::default();
    let mut bridges = Vec::new();

    for &node in &remaining {
        if visited.contains(&node) {
            continue;
        }
        let mut path = Vec::new();
        let mut stack = vec![node];
        while let Some(n) = stack.pop() {
            if !visited.insert(n) {
                continue;
            }
            path.push(n);
            for neighbor in mst.neighbors(n) {
                if remaining.contains(&neighbor) && !visited.contains(&neighbor) {
                    stack.push(neighbor);
                }
            }
        }

        // A bridge is a short path where all nodes have degree 2 in the original MST
        if path.len() < bridge_length && path.iter().all(|&n| mst.neighbors(n).count() == 2) {
            bridges.push(path);
        }
    }

    // Also include edges between junctions (matching original: `bridges += g.subgraph(junctions).edges`)
    for &j1 in &junctions {
        for neighbor in mst.neighbors(j1) {
            if junctions.contains(&neighbor) && j1.index() < neighbor.index() {
                bridges.push(vec![j1, neighbor]);
            }
        }
    }

    bridges
}

/// Remove bridges iteratively from the graph.
///
/// Matches original `remove_bridges`:
///   while True:
///     gmst = determine_pruned_mst(g)
///     bridges = identify_bridges(gmst, bridge_length)
///     if no bridges: break
///     g.remove_nodes_from(bridge_nodes)
fn remove_bridges(graph: &mut OverlapGraph, bridge_length: usize, prune_branch_size: usize) {
    let mut iterations = 0;
    let mut total_removed = 0;

    loop {
        let (mst, node_map) = determine_pruned_mst(graph, prune_branch_size);
        let bridges = identify_bridges(&mst, bridge_length);

        if bridges.is_empty() {
            break;
        }

        // Build reverse map: MST index -> original graph index
        let reverse_map: FxHashMap<NodeIndex, NodeIndex> = node_map
            .iter()
            .map(|(&orig, &mst_idx)| (mst_idx, orig))
            .collect();

        // Remove bridge nodes from the ORIGINAL graph (not the MST)
        let mut removed = 0;
        for bridge in &bridges {
            for &mst_node in bridge {
                if let Some(&orig_node) = reverse_map.get(&mst_node) {
                    if graph.node_weight(orig_node).is_some() {
                        graph.remove_node(orig_node);
                        removed += 1;
                    }
                }
            }
        }

        total_removed += removed;
        iterations += 1;

        if removed == 0 {
            break;
        }
    }

    if iterations > 0 {
        log::info!(
            "Removed {} vertices in bridges over {} iterations",
            total_removed,
            iterations
        );
    }
}

/// Extract backbone paths as named paths (lists of vertex names).
///
/// Matches original `determine_backbones`:
///   1. Copy graph
///   2. remove_bridges if configured
///   3. gmst = determine_pruned_mst(g)
///   4. while gmst not empty:
///      paths = determine_backbones_of_trees(gmst, prune_junctions)
///      backbones += [p for p in paths if len(p) >= prune_branches]
///      gmst.remove_nodes_from(all path nodes)
///   5. sort by length descending
pub fn extract_named_backbones(g: &NamedGraph, config: &BackboneConfig) -> Vec<Vec<String>> {
    log::info!("Extracting backbone paths...");

    // Work on a copy of the graph
    let mut work_graph = g.graph.clone();

    // Step 1: Remove bridges iteratively
    if config.prune_bridge_size > 0 {
        remove_bridges(
            &mut work_graph,
            config.prune_bridge_size,
            config.prune_branch_size,
        );
    }

    // Step 2: Compute pruned MST
    let (mut mst, node_map) = determine_pruned_mst(&work_graph, config.prune_branch_size);
    log::info!("MST: V={} E={}", mst.node_count(), mst.edge_count());

    // Build reverse map: MST index -> original graph index
    let reverse_map: FxHashMap<NodeIndex, NodeIndex> = node_map
        .iter()
        .map(|(&orig, &mst_idx)| (mst_idx, orig))
        .collect();

    // Step 3: Iteratively extract backbones
    let mut all_paths: Vec<Vec<String>> = Vec::new();

    loop {
        if mst.node_count() == 0 {
            break;
        }

        // Extract backbones from all trees in the MST
        let paths = determine_backbones_of_trees(&mst, config.prune_junction_size);

        if paths.is_empty() {
            break;
        }

        let mut any_added = false;
        for path in &paths {
            // Convert MST indices to original names
            let named_path: Vec<String> = path
                .iter()
                .filter_map(|&mst_idx| {
                    let orig_idx = reverse_map.get(&mst_idx)?;
                    g.names.get_name(*orig_idx).map(String::from)
                })
                .collect();

            if named_path.len() >= config.prune_branch_size {
                all_paths.push(named_path);
            }

            // Remove path nodes from MST
            for &node in path {
                mst.remove_node(node);
            }
            any_added = true;
        }

        if !any_added {
            break;
        }
    }

    all_paths.sort_by_key(|p| std::cmp::Reverse(p.len()));

    let total_nodes: usize = all_paths.iter().map(|p| p.len()).sum();
    log::info!(
        "Extracted {} backbone paths containing {} molecules",
        all_paths.len(),
        total_nodes
    );

    // Merge paths that are connected through the original graph
    let merged = merge_adjacent_named_paths(&all_paths, g, config.min_path_size);
    log::info!(
        "After merging: {} paths ({} molecules)",
        merged.len(),
        merged.iter().map(|p| p.len()).sum::<usize>()
    );

    merged
}

/// Merge backbone paths that are connected through the original molecule graph.
///
/// Two paths can be merged if an endpoint of one path shares a neighbor
/// in the original graph with an endpoint of another path, and that
/// neighbor is not already in any backbone path.
///
/// This bridges gaps caused by pruning or junction splitting.
fn merge_adjacent_named_paths(
    paths: &[Vec<String>],
    g: &NamedGraph,
    min_path_size: usize,
) -> Vec<Vec<String>> {
    if paths.len() <= 1 {
        return paths.to_vec();
    }

    // Build set of all backbone node names
    let mut backbone_names: FxHashSet<&str> = FxHashSet::default();
    for path in paths {
        for name in path {
            backbone_names.insert(name.as_str());
        }
    }

    // Map endpoint names to (path_index, is_last_endpoint)
    let mut endpoint_to_path: FxHashMap<&str, (usize, bool)> = FxHashMap::default();
    for (i, path) in paths.iter().enumerate() {
        if path.len() >= 2 {
            endpoint_to_path.insert(path[0].as_str(), (i, false));
            endpoint_to_path.insert(path.last().unwrap().as_str(), (i, true));
        }
    }

    // Find merge candidates: pairs of paths whose endpoints are connected
    // through the original graph (directly or via one non-backbone bridge node).
    let mut merge_edges: Vec<(usize, usize, bool, bool)> = Vec::new();

    for (&ep_name, &(path_idx, is_last)) in &endpoint_to_path {
        let ep_node = match g.names.get_idx(ep_name) {
            Some(idx) => idx,
            None => continue,
        };
        if g.graph.node_weight(ep_node).is_none() {
            continue;
        }

        for neighbor in g.graph.neighbors(ep_node) {
            let neighbor_name = match g.names.get_name(neighbor) {
                Some(n) => n,
                None => continue,
            };

            // Direct: endpoint -> other_endpoint
            if let Some(&(other_path, other_is_last)) = endpoint_to_path.get(neighbor_name) {
                if other_path != path_idx {
                    merge_edges.push((path_idx, other_path, is_last, other_is_last));
                }
            }

            // One-hop bridge: endpoint -> bridge -> other_endpoint
            if !backbone_names.contains(neighbor_name) {
                for neighbor2 in g.graph.neighbors(neighbor) {
                    let n2_name = match g.names.get_name(neighbor2) {
                        Some(n) => n,
                        None => continue,
                    };
                    if let Some(&(other_path, other_is_last)) = endpoint_to_path.get(n2_name) {
                        if other_path != path_idx {
                            merge_edges.push((path_idx, other_path, is_last, other_is_last));
                        }
                    }
                }
            }
        }
    }

    if merge_edges.is_empty() {
        let mut result = paths.to_vec();
        result.retain(|p| p.len() >= min_path_size);
        result.sort_by_key(|p| std::cmp::Reverse(p.len()));
        return result;
    }

    log::info!(
        "Found {} merge candidates between path endpoints",
        merge_edges.len()
    );

    // Union-find to group paths
    let n = paths.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut merge_info: FxHashMap<(usize, usize), (bool, bool)> = FxHashMap::default();

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    // Prefer merging larger paths first
    merge_edges.sort_by_key(|&(a, b, _, _)| std::cmp::Reverse(paths[a].len().min(paths[b].len())));

    for (a, b, a_is_last, b_is_last) in merge_edges {
        let ra = find(&mut parent, a);
        let rb = find(&mut parent, b);
        if ra != rb {
            parent[ra] = rb;
            merge_info.insert((a, b), (a_is_last, b_is_last));
        }
    }

    // Group paths by root
    let mut groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for i in 0..n {
        let root = find(&mut parent, i);
        groups.entry(root).or_default().push(i);
    }

    let mut result: Vec<Vec<String>> = Vec::new();
    for group in groups.values() {
        if group.len() == 1 {
            result.push(paths[group[0]].clone());
            continue;
        }

        // Build adjacency among paths in this group
        let group_set: FxHashSet<usize> = group.iter().copied().collect();
        let mut path_adj: FxHashMap<usize, Vec<(usize, bool, bool)>> = FxHashMap::default();

        for (&(a, b), &(a_last, b_last)) in &merge_info {
            if group_set.contains(&a) && group_set.contains(&b) {
                path_adj.entry(a).or_default().push((b, a_last, b_last));
                path_adj.entry(b).or_default().push((a, b_last, a_last));
            }
        }

        // Start from a path with degree ≤ 1 in the path adjacency graph
        let start = group
            .iter()
            .copied()
            .find(|&p| path_adj.get(&p).map_or(0, |v| v.len()) <= 1)
            .unwrap_or(group[0]);

        let mut chain: Vec<String> = Vec::new();
        let mut visited_paths: FxHashSet<usize> = FxHashSet::default();
        let mut current = start;
        let mut need_reverse = false;

        loop {
            visited_paths.insert(current);
            let mut p = paths[current].clone();
            if need_reverse {
                p.reverse();
            }
            chain.extend(p);

            let next = path_adj
                .get(&current)
                .and_then(|adjs| {
                    adjs.iter()
                        .find(|(other, _, _)| !visited_paths.contains(other))
                })
                .copied();

            match next {
                Some((next_path, _my_is_last, other_is_last)) => {
                    need_reverse = other_is_last;
                    current = next_path;
                }
                None => break,
            }
        }

        // Add unchained paths from this group
        for &p in group {
            if !visited_paths.contains(&p) {
                result.push(paths[p].clone());
            }
        }

        result.push(chain);
    }

    result.retain(|p| p.len() >= min_path_size);
    result.sort_by_key(|p| std::cmp::Reverse(p.len()));
    result
}

/// Extract backbone paths from all trees in the MST.
///
/// Matches original `determine_backbones_of_trees`:
///   For each connected component:
///     If has junctions with 3+ large branches: split at junctions
///     Else: find diameter path
fn determine_backbones_of_trees(
    mst: &OverlapGraph,
    prune_junction_size: usize,
) -> Vec<Vec<NodeIndex>> {
    let components = connected_components(mst);
    let mut all_paths = Vec::new();

    for comp in &components {
        if comp.len() < 2 {
            all_paths.push(comp.clone());
            continue;
        }

        // Detect junctions
        if prune_junction_size > 0 {
            let junctions = detect_junctions(mst, comp, prune_junction_size);
            if !junctions.is_empty() {
                let paths = split_at_junctions(mst, comp, &junctions);
                all_paths.extend(paths);
                continue;
            }
        }

        // No junctions — find diameter path
        all_paths.push(diameter_path(mst, comp));
    }

    all_paths
}

/// Detect junction nodes in a tree component.
///
/// Matches original `detect_junctions_of_tree`:
///   candidate_junctions = [u for u, deg in g.degree() if deg >= 3]
///   For each candidate: count branches >= minor_branch_size
///   If 3+ large branches: it's a junction
fn detect_junctions(
    tree: &OverlapGraph,
    component: &[NodeIndex],
    min_branch_size: usize,
) -> Vec<NodeIndex> {
    let branch_lengths = measure_branch_lengths(tree, component);
    let mut junctions = Vec::new();

    for &node in component {
        let degree = tree.neighbors(node).count();
        if degree < 3 {
            continue;
        }

        let mut large_branches = 0;
        for neighbor in tree.neighbors(node) {
            if let Some(&len) = branch_lengths.get(&(node, neighbor)) {
                if len >= min_branch_size {
                    large_branches += 1;
                }
            }
        }

        if large_branches >= 3 {
            junctions.push(node);
        }
    }
    junctions
}

/// Split a tree at junction nodes.
///
/// Matches original `split_junctions_of_tree` with `keep_largest=0` (default):
///   For each junction: remove ALL edges (disconnecting the junction entirely)
///   Then find connected components and extract diameter paths
fn split_at_junctions(
    tree: &OverlapGraph,
    component: &[NodeIndex],
    junctions: &[NodeIndex],
) -> Vec<Vec<NodeIndex>> {
    // Remove ALL edges at junctions (matching original keep_largest=0 default)
    let mut edges_to_skip: FxHashSet<(NodeIndex, NodeIndex)> = FxHashSet::default();

    for &junction in junctions {
        for edge in tree.edges(junction) {
            let other = if edge.source() == junction {
                edge.target()
            } else {
                edge.source()
            };
            edges_to_skip.insert((junction, other));
            edges_to_skip.insert((other, junction));
        }
    }

    // Find connected components after removing junction edges
    let comp_set: FxHashSet<NodeIndex> = component.iter().copied().collect();
    let mut visited = FxHashSet::default();
    let mut paths = Vec::new();

    for &start in component {
        if visited.contains(&start) {
            continue;
        }

        let mut comp = Vec::new();
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            if !visited.insert(node) {
                continue;
            }
            comp.push(node);
            for neighbor in tree.neighbors(node) {
                if comp_set.contains(&neighbor)
                    && !visited.contains(&neighbor)
                    && !edges_to_skip.contains(&(node, neighbor))
                {
                    stack.push(neighbor);
                }
            }
        }

        if !comp.is_empty() {
            let path = diameter_path(tree, &comp);
            paths.push(path);
        }
    }

    paths
}
/// Merge adjacent backbone paths using split-minimizer bridge evidence.
///
/// Finds non-backbone molecules whose split minimizers overlap with
/// endpoint molecules from two different backbone paths. These bridge
/// molecules indicate that the two paths are adjacent on the same
/// chromosome and can be merged.
///
/// This is an optional post-processing step for the physical map.
/// Configuration for path merging.
#[derive(Debug, Clone)]
pub struct MergePathsConfig {
    /// Number of molecules from each path end to use as endpoints.
    pub endpoint_depth: usize,
    /// Minimum shared minimizers between a bridge molecule and an endpoint.
    pub min_shared_mx: usize,
    /// Minimum bridge molecules to consider a link.
    pub min_bridges: usize,
    /// Maximum number of different paths a bridge molecule can connect to
    /// (specificity filter — higher values allow more noise).
    pub max_path_connections: usize,
    /// Minimum path length (molecules) to include in merging.
    pub min_path_size: usize,
    /// Maximum candidate links per endpoint. Endpoints appearing in more
    /// links than this are likely non-specific and are excluded.
    pub max_links_per_endpoint: usize,
    /// Minimum bridge density: bridge_count / min(len_A, len_B).
    /// Filters out links with high absolute bridge count but low relative
    /// evidence. Set to 0.0 to disable.
    pub min_bridge_density: f64,
    /// Minimum number of endpoint molecules a bridge must connect to on
    /// EACH side of a link. A bridge molecule must share >= min_shared_mx
    /// minimizers with >= min_endpoint_hits distinct endpoint molecules on
    /// both sides. Higher values require stronger neighborhood evidence.
    /// Set to 1 for original behavior.
    pub min_endpoint_hits: usize,
}

impl Default for MergePathsConfig {
    fn default() -> Self {
        Self {
            endpoint_depth: 25,
            min_shared_mx: 3,
            min_bridges: 2,
            max_path_connections: 2,
            min_path_size: 50,
            max_links_per_endpoint: 1,
            min_bridge_density: 0.01,
            min_endpoint_hits: 4,
        }
    }
}

/// An endpoint of a backbone path.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Endpoint {
    path_idx: usize,
    is_end: bool, // false = start, true = end
}

/// A merge link between two path endpoints.
#[derive(Debug)]
struct MergeLink {
    ep1: Endpoint,
    ep2: Endpoint,
    bridge_count: usize,
}

/// Merge adjacent backbone paths using split-minimizer bridge molecules.
///
/// Returns a new set of paths where some have been merged based on
/// bridge evidence. This is an optional post-processing step.
pub fn merge_paths(
    paths: &[Vec<String>],
    split_mxs: &FxHashMap<String, FxHashSet<u64>>,
    config: &MergePathsConfig,
) -> Vec<Vec<String>> {
    if paths.is_empty() {
        return Vec::new();
    }

    // Build backbone molecule set and path membership
    let mut mol_to_path: FxHashMap<&str, usize> = FxHashMap::default();
    for (i, path) in paths.iter().enumerate() {
        for mol in path {
            mol_to_path.insert(mol.as_str(), i);
        }
    }
    let backbone_mols: FxHashSet<&str> = mol_to_path.keys().copied().collect();

    // Identify endpoint molecules. Each gets a unique ID for per-molecule
    // hit tracking. endpoint_info maps molecule name -> Endpoint (path-level).
    // endpoint_mol_ids maps molecule name -> unique ID.
    let mut endpoint_info: FxHashMap<&str, Endpoint> = FxHashMap::default();
    let mut endpoint_mol_ids: FxHashMap<&str, usize> = FxHashMap::default();
    let mut endpoint_mol_to_ep: Vec<Endpoint> = Vec::new(); // mol_id -> Endpoint
    let mut next_mol_id = 0usize;

    for (i, path) in paths.iter().enumerate() {
        let depth = config.endpoint_depth.min(path.len());
        for mol in &path[..depth] {
            let ep = Endpoint {
                path_idx: i,
                is_end: false,
            };
            endpoint_info.insert(mol.as_str(), ep);
            endpoint_mol_ids.insert(mol.as_str(), next_mol_id);
            endpoint_mol_to_ep.push(ep);
            next_mol_id += 1;
        }
        let start = if path.len() > depth {
            path.len() - depth
        } else {
            0
        };
        for mol in &path[start..] {
            let ep = Endpoint {
                path_idx: i,
                is_end: true,
            };
            endpoint_info.insert(mol.as_str(), ep);
            endpoint_mol_ids.insert(mol.as_str(), next_mol_id);
            endpoint_mol_to_ep.push(ep);
            next_mol_id += 1;
        }
    }

    log::info!(
        "Scaffold: {} paths, {} endpoint molecules (depth={})",
        paths.len(),
        endpoint_info.len(),
        config.endpoint_depth
    );

    // Build minimizer index: mx -> list of (endpoint_mol_id)
    let mut mx_to_mol_ids: FxHashMap<u64, Vec<usize>> = FxHashMap::default();
    for (&mol_name, &mol_id) in &endpoint_mol_ids {
        if let Some(mxs) = split_mxs.get(mol_name) {
            for &mx in mxs {
                mx_to_mol_ids.entry(mx).or_default().push(mol_id);
            }
        }
    }

    log::info!(
        "Indexed {} minimizers from endpoint molecules",
        mx_to_mol_ids.len()
    );

    // Scan all non-backbone molecules for bridge evidence.
    // For each bridge molecule, count shared minimizers with each individual
    // endpoint molecule, then aggregate to path-endpoints.
    //
    // Bridge evidence is tracked at the path-endpoint pair level:
    // (Endpoint, Endpoint) -> set of bridge molecule names
    let mut bridge_evidence: FxHashMap<(Endpoint, Endpoint), FxHashSet<String>> =
        FxHashMap::default();

    let mut scanned = 0usize;
    for (mol_name, mxs) in split_mxs {
        // Skip backbone molecules
        if backbone_mols.contains(mol_name.as_str()) {
            continue;
        }

        // Count shared minimizers with each individual endpoint molecule
        let mut mol_hits: FxHashMap<usize, usize> = FxHashMap::default();
        for &mx in mxs {
            if let Some(mol_ids) = mx_to_mol_ids.get(&mx) {
                for &mid in mol_ids {
                    *mol_hits.entry(mid).or_insert(0) += 1;
                }
            }
        }

        // Group by path-endpoint: count how many endpoint molecules this bridge
        // shares >= min_shared_mx minimizers with, per path-end
        let mut ep_mol_counts: FxHashMap<Endpoint, usize> = FxHashMap::default();
        for (&mid, &shared_count) in &mol_hits {
            if shared_count >= config.min_shared_mx {
                let ep = endpoint_mol_to_ep[mid];
                *ep_mol_counts.entry(ep).or_insert(0) += 1;
            }
        }

        // Filter: only keep path-endpoints where the bridge connects to
        // enough endpoint molecules (min_endpoint_hits)
        let strong_eps: Vec<Endpoint> = ep_mol_counts
            .iter()
            .filter(|(_, &count)| count >= config.min_endpoint_hits)
            .map(|(&ep, _)| ep)
            .collect();

        // Check specificity: how many different paths does this molecule connect to?
        let connected_paths: FxHashSet<usize> = strong_eps.iter().map(|ep| ep.path_idx).collect();

        if connected_paths.len() < 2 || connected_paths.len() > config.max_path_connections {
            scanned += 1;
            continue;
        }

        // Record bridge evidence for each pair of path-endpoints from different paths
        for i in 0..strong_eps.len() {
            for j in (i + 1)..strong_eps.len() {
                let ep1 = strong_eps[i];
                let ep2 = strong_eps[j];
                if ep1.path_idx == ep2.path_idx {
                    continue;
                }
                let key = if ep1 < ep2 { (ep1, ep2) } else { (ep2, ep1) };
                bridge_evidence
                    .entry(key)
                    .or_default()
                    .insert(mol_name.clone());
            }
        }

        scanned += 1;
    }

    log::info!(
        "Scanned {} non-backbone molecules, found {} endpoint pairs with bridges",
        scanned,
        bridge_evidence.len()
    );

    // Collect merge links above threshold, filtered by path size and bridge density
    let mut links: Vec<MergeLink> = bridge_evidence
        .into_iter()
        .filter_map(|((ep1, ep2), bridges)| {
            let count = bridges.len();
            if count < config.min_bridges {
                return None;
            }
            let len_a = paths[ep1.path_idx].len();
            let len_b = paths[ep2.path_idx].len();
            // Both paths must be large enough
            if len_a < config.min_path_size || len_b < config.min_path_size {
                return None;
            }
            // Bridge density filter: bridges / min(len_A, len_B)
            if config.min_bridge_density > 0.0 {
                let min_len = len_a.min(len_b) as f64;
                let density = count as f64 / min_len;
                if density < config.min_bridge_density {
                    return None;
                }
            }
            Some(MergeLink {
                ep1,
                ep2,
                bridge_count: count,
            })
        })
        .collect();

    // Sort by bridge count descending (strongest evidence first)
    links.sort_by_key(|l| std::cmp::Reverse(l.bridge_count));

    log::info!(
        "Found {} merge links (min_bridges={}, min_path_size={}, min_density={:.3}, min_ep_hits={})",
        links.len(),
        config.min_bridges,
        config.min_path_size,
        config.min_bridge_density,
        config.min_endpoint_hits
    );

    // Filter out promiscuous endpoints that appear in too many links.
    // These are likely non-specific (repetitive minimizers).
    let mut endpoint_link_count: FxHashMap<Endpoint, usize> = FxHashMap::default();
    for link in &links {
        *endpoint_link_count.entry(link.ep1).or_insert(0) += 1;
        *endpoint_link_count.entry(link.ep2).or_insert(0) += 1;
    }

    let before_filter = links.len();
    links.retain(|link| {
        let c1 = endpoint_link_count.get(&link.ep1).copied().unwrap_or(0);
        let c2 = endpoint_link_count.get(&link.ep2).copied().unwrap_or(0);
        c1 <= config.max_links_per_endpoint && c2 <= config.max_links_per_endpoint
    });

    if links.len() < before_filter {
        log::info!(
            "Filtered {} links with promiscuous endpoints (max_links_per_endpoint={}), {} remaining",
            before_filter - links.len(),
            config.max_links_per_endpoint,
            links.len()
        );
    }

    for link in &links {
        log::info!(
            "  Link: path {} ({}) <-> path {} ({}) — {} bridge molecules",
            link.ep1.path_idx,
            if link.ep1.is_end { "end" } else { "start" },
            link.ep2.path_idx,
            if link.ep2.is_end { "end" } else { "start" },
            link.bridge_count
        );
    }

    if links.is_empty() {
        return paths.to_vec();
    }

    // Greedy path merging using union-find
    // Each endpoint can participate in at most one link
    let n = paths.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut used_endpoints: FxHashSet<Endpoint> = FxHashSet::default();

    // For each path, track which other path is joined at each end
    // join_at[path_idx] = (start_join, end_join)
    // where each join is Option<(other_path_idx, other_is_end)>
    let mut join_at_start: FxHashMap<usize, (usize, bool)> = FxHashMap::default();
    let mut join_at_end: FxHashMap<usize, (usize, bool)> = FxHashMap::default();

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut [usize], x: usize, y: usize) {
        let rx = find(parent, x);
        let ry = find(parent, y);
        if rx != ry {
            parent[rx] = ry;
        }
    }

    let mut accepted = 0;
    for link in &links {
        // Skip if either endpoint is already used
        if used_endpoints.contains(&link.ep1) || used_endpoints.contains(&link.ep2) {
            continue;
        }

        // Skip if paths are already in the same scaffold (would create a cycle)
        if find(&mut parent, link.ep1.path_idx) == find(&mut parent, link.ep2.path_idx) {
            continue;
        }

        // Accept this link
        used_endpoints.insert(link.ep1);
        used_endpoints.insert(link.ep2);
        union(&mut parent, link.ep1.path_idx, link.ep2.path_idx);

        if link.ep1.is_end {
            join_at_end.insert(link.ep1.path_idx, (link.ep2.path_idx, link.ep2.is_end));
        } else {
            join_at_start.insert(link.ep1.path_idx, (link.ep2.path_idx, link.ep2.is_end));
        }

        if link.ep2.is_end {
            join_at_end.insert(link.ep2.path_idx, (link.ep1.path_idx, link.ep1.is_end));
        } else {
            join_at_start.insert(link.ep2.path_idx, (link.ep1.path_idx, link.ep1.is_end));
        }

        accepted += 1;
    }

    log::info!("Accepted {} merge links", accepted);

    if accepted == 0 {
        return paths.to_vec();
    }

    // Build scaffold chains by following the join links
    let mut used_in_scaffold: FxHashSet<usize> = FxHashSet::default();
    let mut merged_paths: Vec<Vec<String>> = Vec::new();

    for start_path in 0..n {
        if used_in_scaffold.contains(&start_path) {
            continue;
        }
        if paths[start_path].len() < config.min_path_size {
            continue;
        }

        // Check if this path has any joins
        let has_start_join = join_at_start.contains_key(&start_path);
        let has_end_join = join_at_end.contains_key(&start_path);

        if !has_start_join && !has_end_join {
            continue; // No joins, will be added as-is later
        }

        // Find the leftmost path in this chain by walking backwards
        // through start-joins, with cycle detection
        let mut chain_start = start_path;
        let mut visited_backward: FxHashSet<usize> = FxHashSet::default();
        visited_backward.insert(chain_start);

        loop {
            // If chain_start has a start-join, the previous path connects here
            if let Some(&(prev_path, _prev_is_end)) = join_at_start.get(&chain_start) {
                if !visited_backward.contains(&prev_path) && !used_in_scaffold.contains(&prev_path)
                {
                    visited_backward.insert(prev_path);
                    chain_start = prev_path;
                    continue;
                }
            }
            break;
        }

        // Walk forward from chain_start, building the merged path
        let mut merged = Vec::new();
        let mut current = chain_start;
        let mut entering_from_end = false;

        // Determine initial orientation: if we reached chain_start via its end-join,
        // we need to check how we got here
        // Actually, chain_start is the leftmost — we always start from its natural beginning
        // unless its start has no join (which is why it's the chain start)

        let mut chain_visited: FxHashSet<usize> = FxHashSet::default();

        loop {
            if chain_visited.contains(&current) || used_in_scaffold.contains(&current) {
                break;
            }
            chain_visited.insert(current);
            used_in_scaffold.insert(current);

            let path_mols = &paths[current];
            if entering_from_end {
                // We entered from the end, so reverse the path
                for mol in path_mols.iter().rev() {
                    merged.push(mol.clone());
                }
                // Exit from the start — check start join
                if let Some(&(next_path, next_is_end)) = join_at_start.get(&current) {
                    if !chain_visited.contains(&next_path) && !used_in_scaffold.contains(&next_path)
                    {
                        entering_from_end = next_is_end;
                        current = next_path;
                        continue;
                    }
                }
            } else {
                // Normal direction
                merged.extend(path_mols.iter().cloned());
                // Exit from the end — check end join
                if let Some(&(next_path, next_is_end)) = join_at_end.get(&current) {
                    if !chain_visited.contains(&next_path) && !used_in_scaffold.contains(&next_path)
                    {
                        entering_from_end = next_is_end;
                        current = next_path;
                        continue;
                    }
                }
            }
            break;
        }

        if !merged.is_empty() {
            merged_paths.push(merged);
        }
    }
    // Add unmerged paths
    for (i, path) in paths.iter().enumerate() {
        if !used_in_scaffold.contains(&i) {
            merged_paths.push(path.clone());
        }
    }

    // Sort by length descending
    merged_paths.sort_by_key(|p| std::cmp::Reverse(p.len()));

    let total_mols: usize = merged_paths.iter().map(|p| p.len()).sum();
    log::info!(
        "After merging: {} paths ({} molecules)",
        merged_paths.len(),
        total_mols
    );

    merged_paths
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_scaffold_empty() {
        let paths: Vec<Vec<String>> = Vec::new();
        let mxs = FxHashMap::default();
        let config = MergePathsConfig::default();
        let result = merge_paths(&paths, &mxs, &config);
        assert!(result.is_empty());
    }

    #[test]
    fn test_scaffold_no_bridges() {
        let paths = vec![
            (0..100).map(|i| format!("a_{}", i)).collect::<Vec<_>>(),
            (0..100).map(|i| format!("b_{}", i)).collect::<Vec<_>>(),
        ];
        let mxs = FxHashMap::default();
        let config = MergePathsConfig::default();
        let result = merge_paths(&paths, &mxs, &config);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_scaffold_with_bridge() {
        // Two paths of 60 molecules each
        let path_a: Vec<String> = (0..60).map(|i| format!("a_{}", i)).collect();
        let path_b: Vec<String> = (0..60).map(|i| format!("b_{}", i)).collect();
        let paths = vec![path_a, path_b];

        // Endpoint molecules: a_50..a_59 (end of path A), b_0..b_9 (start of path B)
        // Bridge molecules: bridge_0, bridge_1 share minimizers with both endpoints
        let mut mxs = FxHashMap::default();

        // Give endpoint molecules some minimizers
        for i in 50..60 {
            let mut s = FxHashSet::default();
            s.insert(1000 + i as u64);
            s.insert(2000 + i as u64);
            s.insert(9000); // shared with bridge
            s.insert(9001);
            s.insert(9002);
            mxs.insert(format!("a_{}", i), s);
        }
        for i in 0..10 {
            let mut s = FxHashSet::default();
            s.insert(3000 + i as u64);
            s.insert(4000 + i as u64);
            s.insert(8000); // shared with bridge
            s.insert(8001);
            s.insert(8002);
            mxs.insert(format!("b_{}", i), s);
        }

        // Bridge molecules share minimizers with both endpoints
        for b in 0..3 {
            let mut s = FxHashSet::default();
            s.insert(9000); // shared with path A endpoint
            s.insert(9001);
            s.insert(9002);
            s.insert(8000); // shared with path B endpoint
            s.insert(8001);
            s.insert(8002);
            mxs.insert(format!("bridge_{}", b), s);
        }

        let config = MergePathsConfig {
            endpoint_depth: 10,
            min_shared_mx: 3,
            min_bridges: 2,
            max_path_connections: 4,
            min_path_size: 50,
            max_links_per_endpoint: 3,
            min_bridge_density: 0.0,
            min_endpoint_hits: 1,
        };

        let result = merge_paths(&paths, &mxs, &config);
        // Should merge into 1 path
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].len(), 120);
        // First 60 should be path A, next 60 should be path B
        assert_eq!(result[0][0], "a_0");
        assert_eq!(result[0][59], "a_59");
        assert_eq!(result[0][60], "b_0");
        assert_eq!(result[0][119], "b_59");
    }
}
