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
fn identify_bridges(
    mst: &OverlapGraph,
    bridge_length: usize,
) -> Vec<Vec<NodeIndex>> {
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
fn remove_bridges(
    graph: &mut OverlapGraph,
    bridge_length: usize,
    prune_branch_size: usize,
) {
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
        remove_bridges(&mut work_graph, config.prune_bridge_size, config.prune_branch_size);
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

    log::info!("Found {} merge candidates between path endpoints", merge_edges.len());

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
    merge_edges.sort_by_key(|&(a, b, _, _)| {
        std::cmp::Reverse(paths[a].len().min(paths[b].len()))
    });

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
    for (_root, group) in &groups {
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
        let start = group.iter().copied()
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

            let next = path_adj.get(&current).and_then(|adjs| {
                adjs.iter().find(|(other, _, _)| !visited_paths.contains(other))
            }).copied();

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
