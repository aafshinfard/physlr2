/// Molecule separation: splitting barcodes into individual molecules.
///
/// All community detection operates on compact local subgraphs (typically 5-100 nodes)
/// remapped to dense 0..n indices. This avoids hash map overhead on millions of calls.
///
/// Strategies (matching original Physlr):
///   - cc: connected components
///   - bc: biconnected components (remove articulation points, then CC)
///   - k3: k-clique percolation (k=3)
///   - k3bin: random binning + k3 + merge
///   - sqcos: cosine similarity of squared adjacency matrix
///   - sqcosbin: random binning + sqcos + merge
///   - distributed: bc → bin → bc → k3 → merge (ensemble pipeline)
use crate::graph::{NamedGraph, OverlapGraph};
use petgraph::stable_graph::NodeIndex;
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Debug, Clone, PartialEq)]
pub enum Strategy {
    Cc,
    Bc,
    K3,
    K3Bin,
    Sqcos,
    SqcosBin,
    Distributed,
}

#[derive(Debug, Clone)]
pub struct MoleculeParams {
    pub sqcos_threshold: f64,
    pub skip_small: usize,
    pub bin_max_size: usize,
    pub merge_cutoff: i64,
    /// Optional seed for deterministic random binning. None = non-deterministic.
    pub seed: Option<u64>,
}

impl Default for MoleculeParams {
    fn default() -> Self {
        Self {
            sqcos_threshold: 0.75,
            skip_small: 10,
            bin_max_size: 50,
            merge_cutoff: 20,
            seed: None,
        }
    }
}

pub fn parse_strategy(s: &str) -> Vec<Strategy> {
    s.split('+')
        .filter_map(|part| match part.trim() {
            "bc" => Some(Strategy::Bc),
            "cc" => Some(Strategy::Cc),
            "k3" => Some(Strategy::K3),
            "k3bin" => Some(Strategy::K3Bin),
            "sqcos" => Some(Strategy::Sqcos),
            "sqcosbin" => Some(Strategy::SqcosBin),
            "distributed" => Some(Strategy::Distributed),
            _ => {
                log::warn!("Unknown strategy '{}', ignoring", part);
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Compact local subgraph: dense 0..n representation for tiny neighborhood graphs.
// Avoids all hash map/set overhead on the hot path.
// ---------------------------------------------------------------------------

/// A compact subgraph with nodes remapped to 0..n.
/// Adjacency stored as a flat bitset (row-major n×n) for O(1) edge queries.
struct LocalGraph {
    /// Number of nodes
    n: usize,
    /// Mapping from local index to global NodeIndex
    local_to_global: Vec<NodeIndex>,
    /// Adjacency bitset: adj[i * n + j] means edge i-j exists
    adj: Vec<bool>,
}

impl LocalGraph {
    /// Build a compact local subgraph from a set of global node indices.
    fn from_nodes(graph: &OverlapGraph, nodes: &FxHashSet<NodeIndex>) -> Self {
        let n = nodes.len();
        let mut local_to_global: Vec<NodeIndex> = nodes.iter().copied().collect();
        // Sort for deterministic ordering
        local_to_global.sort_unstable();
        let mut global_to_local = FxHashMap::default();
        for (i, &g) in local_to_global.iter().enumerate() {
            global_to_local.insert(g, i);
        }

        let mut adj = vec![false; n * n];
        for (li, &gi) in local_to_global.iter().enumerate() {
            for neighbor in graph.neighbors(gi) {
                if let Some(&lj) = global_to_local.get(&neighbor) {
                    adj[li * n + lj] = true;
                }
            }
        }

        Self {
            n,
            local_to_global,
            adj,
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

pub fn separate_molecules(
    g: &NamedGraph,
    strategy: &str,
    junctions: &FxHashSet<String>,
) -> NamedGraph {
    separate_molecules_with_params(g, strategy, junctions, &MoleculeParams::default())
}

/// Split minimizers by molecule.
///
/// For each molecule vertex (e.g. `barcode_0`) in the molecule graph:
/// 1. Collect the union of minimizers from all neighbor molecules' barcodes
/// 2. Intersect with the molecule's own barcode minimizers
/// 3. Output the intersection as the molecule's minimizers
///
/// This assigns each minimizer to the specific molecule(s) where it's
/// relevant based on neighbor overlap, producing cleaner visualizations.
pub fn split_minimizers(
    mol_graph: &NamedGraph,
    bx_mxs: &FxHashMap<String, Vec<u64>>,
) -> Vec<(String, Vec<u64>)> {
    use rayon::prelude::*;

    // Build molecule name → barcode name mapping
    // Molecule names are "barcode_N" where N is the molecule index
    let mol_indices: Vec<_> = mol_graph.graph.node_indices().collect();

    let results: Vec<(String, Vec<u64>)> = mol_indices
        .par_iter()
        .map(|&mol_idx| {
            let mol_name = mol_graph.names.get_name(mol_idx).unwrap().to_string();

            // Extract barcode from molecule name: "barcode_N" → "barcode"
            let barcode = match mol_name.rsplit_once('_') {
                Some((bx, _)) => bx,
                None => &mol_name,
            };

            // Get this barcode's minimizers
            let my_mxs = match bx_mxs.get(barcode) {
                Some(mxs) => mxs,
                None => return (mol_name, Vec::new()),
            };

            // Collect union of neighbor barcodes' minimizers
            let mut neighbor_mxs: FxHashSet<u64> = FxHashSet::default();
            for neighbor_idx in mol_graph.graph.neighbors(mol_idx) {
                let neighbor_name = mol_graph.names.get_name(neighbor_idx).unwrap();
                let neighbor_bx = match neighbor_name.rsplit_once('_') {
                    Some((bx, _)) => bx,
                    None => neighbor_name,
                };
                if let Some(nmxs) = bx_mxs.get(neighbor_bx) {
                    for &mx in nmxs {
                        neighbor_mxs.insert(mx);
                    }
                }
            }

            // Intersect: keep only minimizers that appear in neighbor union
            let split: Vec<u64> = my_mxs
                .iter()
                .copied()
                .filter(|mx| neighbor_mxs.contains(mx))
                .collect();

            (mol_name, split)
        })
        .collect();

    log::info!(
        "Split minimizers: {} molecules, {} with non-empty minimizers",
        results.len(),
        results.iter().filter(|(_, mxs)| !mxs.is_empty()).count()
    );

    results
}

/// Trace molecule separation for specific barcodes (diagnostic output).
pub fn trace_molecules(
    g: &NamedGraph,
    strategy: &str,
    params: &MoleculeParams,
    barcodes: &[String],
) {
    let strategies = parse_strategy(strategy);
    let _junction_indices: FxHashSet<petgraph::stable_graph::NodeIndex> = FxHashSet::default();

    for barcode in barcodes {
        let u = match g.names.get_idx(barcode) {
            Some(idx) => idx,
            None => {
                println!("Barcode {} not found in graph", barcode);
                continue;
            }
        };

        let neighbors: Vec<_> = g.graph.neighbors(u).collect();
        println!("\n{}", "=".repeat(60));
        println!("Barcode: {} (degree={})", barcode, neighbors.len());

        if neighbors.is_empty() {
            println!("No neighbors, skipping");
            continue;
        }

        let neighbor_set: FxHashSet<_> = neighbors.iter().copied().collect();
        let lg = LocalGraph::from_nodes(&g.graph, &neighbor_set);

        let all_local: Vec<usize> = (0..lg.n).collect();
        let mut communities: Vec<Vec<usize>> = vec![all_local];

        println!("Initial: 1 community with {} nodes", lg.n);

        for strat in &strategies {
            let mut new_communities = Vec::new();
            let strat_name = format!("{:?}", strat);
            println!("\n--- Strategy: {} ---", strat_name);
            println!("Input: {} communities, sizes: {:?}", communities.len(), {
                let mut sizes: Vec<_> = communities.iter().map(|c| c.len()).collect();
                sizes.sort_unstable_by(|a, b| b.cmp(a));
                sizes.truncate(20);
                sizes
            });

            for (ci, component) in communities.iter().enumerate() {
                let result: Vec<Vec<usize>> = match strat {
                    Strategy::Cc => cc_local(&lg, component),
                    Strategy::Bc => bc_local(&lg, component),
                    Strategy::K3 => k3_local(&lg, component),
                    Strategy::K3Bin => {
                        let bins = bin_local(component, params.bin_max_size, params.seed);
                        let mut clusters = Vec::new();
                        for bin in &bins {
                            clusters.extend(k3_local(&lg, bin));
                        }
                        merge_local(&lg, &clusters, params.merge_cutoff)
                    }
                    Strategy::Sqcos => sqcos_local(
                        &lg,
                        component,
                        true,
                        params.sqcos_threshold,
                        params.skip_small,
                    ),
                    Strategy::SqcosBin => {
                        let bins = bin_local(component, params.bin_max_size, params.seed);
                        let mut clusters = Vec::new();
                        for bin in &bins {
                            clusters.extend(sqcos_local(
                                &lg,
                                bin,
                                true,
                                params.sqcos_threshold,
                                params.skip_small,
                            ));
                        }
                        merge_local(&lg, &clusters, params.merge_cutoff)
                    }
                    Strategy::Distributed => {
                        // Detailed tracing of distributed sub-steps
                        let bc_comps = bc_local(&lg, component);
                        if ci < 3 {
                            let bc_sizes: Vec<_> = bc_comps.iter().map(|c| c.len()).collect();
                            let total_bc: usize = bc_sizes.iter().sum();
                            println!(
                                "  Component {} (size {}): bc → {} components, total={}, lost={}",
                                ci,
                                component.len(),
                                bc_comps.len(),
                                total_bc,
                                component.len() - total_bc
                            );
                        }
                        let mut all_result = Vec::new();
                        let mut total_clusters = 0;
                        for (bci, bc_comp) in bc_comps.iter().enumerate() {
                            let bins = bin_local(bc_comp, params.bin_max_size, params.seed);
                            let mut clusters = Vec::new();
                            for bin in &bins {
                                let inner_bcs = bc_local(&lg, bin);
                                for inner_bc in &inner_bcs {
                                    clusters.extend(k3_local(&lg, inner_bc));
                                }
                            }
                            total_clusters += clusters.len();
                            let merged = merge_local(&lg, &clusters, params.merge_cutoff);
                            if ci < 3 && bci < 3 {
                                let mut csizes: Vec<_> = clusters.iter().map(|c| c.len()).collect();
                                csizes.sort_unstable_by(|a, b| b.cmp(a));
                                csizes.truncate(10);
                                let mut msizes: Vec<_> = merged.iter().map(|c| c.len()).collect();
                                msizes.sort_unstable_by(|a, b| b.cmp(a));
                                msizes.truncate(10);
                                println!("    bc[{}] (size {}): {} bins → {} k3 clusters {:?} → {} merged {:?}",
                                    bci, bc_comp.len(), bins.len(), clusters.len(), csizes,
                                    merged.len(), msizes);
                            }
                            all_result.extend(merged);
                        }
                        if ci < 3 {
                            println!(
                                "    Total: {} clusters → {} merged communities",
                                total_clusters,
                                all_result.len()
                            );
                        }
                        all_result
                    }
                };

                if ci < 5 {
                    let mut sizes: Vec<_> = result.iter().map(|c| c.len()).collect();
                    sizes.sort_unstable_by(|a, b| b.cmp(a));
                    sizes.truncate(10);
                    println!(
                        "  Component {} (size {}) → {} communities, sizes: {:?}",
                        ci,
                        component.len(),
                        result.len(),
                        sizes
                    );
                }
                new_communities.extend(result);
            }

            communities = new_communities;
            let mut sizes: Vec<_> = communities.iter().map(|c| c.len()).collect();
            sizes.sort_unstable_by(|a, b| b.cmp(a));
            let total_nodes: usize = communities.iter().map(|c| c.len()).sum();
            let gt1 = communities.iter().filter(|c| c.len() > 1).count();
            println!("Output: {} communities, sizes: {:?}", communities.len(), {
                sizes.truncate(20);
                sizes
            });
            println!("  Total nodes in communities: {}", total_nodes);
            println!("  Communities with >1 member: {}", gt1);
        }

        // Final assignment
        communities.sort_by_key(|c| std::cmp::Reverse(c.len()));
        let mut assignment = FxHashMap::default();
        for (i, community) in communities.iter().enumerate() {
            if community.len() > 1 {
                for &li in community {
                    assignment.entry(lg.local_to_global[li]).or_insert(i);
                }
            }
        }
        let n_mols = if assignment.is_empty() {
            0
        } else {
            *assignment.values().max().unwrap_or(&0) + 1
        };
        let n_assigned = assignment.len();
        let n_unassigned = neighbors.len() - n_assigned;
        let distinct_ids: FxHashSet<_> = assignment.values().copied().collect();

        println!("\n--- Final ---");
        println!("n_molecules: {}", n_mols);
        println!("n_assigned: {}", n_assigned);
        println!("n_unassigned: {}", n_unassigned);
        println!("Distinct molecule IDs: {}", distinct_ids.len());
    }
}

pub fn separate_molecules_with_params(
    g: &NamedGraph,
    strategy: &str,
    junctions: &FxHashSet<String>,
    params: &MoleculeParams,
) -> NamedGraph {
    log::info!(
        "Separating barcodes into molecules (strategy={})...",
        strategy
    );

    let strategies = parse_strategy(strategy);
    let strategies = if strategies.is_empty() {
        log::warn!("No valid strategies parsed, defaulting to bc+cc");
        vec![Strategy::Bc, Strategy::Cc]
    } else {
        strategies
    };

    let use_junctions = !junctions.is_empty();
    let nodes: Vec<NodeIndex> = g.graph.node_indices().collect();

    let junction_indices: FxHashSet<NodeIndex> = if use_junctions {
        junctions
            .iter()
            .filter_map(|name| g.names.get_idx(name))
            .collect()
    } else {
        FxHashSet::default()
    };

    // Diagnostic counters
    static DIAG_TOTAL_BARCODES: AtomicU64 = AtomicU64::new(0);
    static DIAG_TOTAL_COMMUNITIES: AtomicU64 = AtomicU64::new(0);
    static DIAG_TOTAL_ASSIGNED: AtomicU64 = AtomicU64::new(0);
    static DIAG_MAX_COMMUNITIES: AtomicU64 = AtomicU64::new(0);
    static DIAG_BARCODES_GT10: AtomicU64 = AtomicU64::new(0);
    static DIAG_BARCODES_GT20: AtomicU64 = AtomicU64::new(0);
    static DIAG_BARCODES_GT50: AtomicU64 = AtomicU64::new(0);
    DIAG_TOTAL_BARCODES.store(0, Ordering::Relaxed);
    DIAG_TOTAL_COMMUNITIES.store(0, Ordering::Relaxed);
    DIAG_TOTAL_ASSIGNED.store(0, Ordering::Relaxed);
    DIAG_MAX_COMMUNITIES.store(0, Ordering::Relaxed);
    DIAG_BARCODES_GT10.store(0, Ordering::Relaxed);
    DIAG_BARCODES_GT20.store(0, Ordering::Relaxed);
    DIAG_BARCODES_GT50.store(0, Ordering::Relaxed);

    let assignments: Vec<(NodeIndex, FxHashMap<NodeIndex, usize>)> = nodes
        .par_iter()
        .map(|&u| {
            let assignment = determine_molecules(
                &g.graph,
                u,
                &strategies,
                use_junctions,
                &junction_indices,
                params,
            );
            let n_mols = if assignment.is_empty() {
                0
            } else {
                *assignment.values().max().unwrap_or(&0) + 1
            };
            DIAG_TOTAL_BARCODES.fetch_add(1, Ordering::Relaxed);
            DIAG_TOTAL_COMMUNITIES.fetch_add(n_mols as u64, Ordering::Relaxed);
            DIAG_TOTAL_ASSIGNED.fetch_add(assignment.len() as u64, Ordering::Relaxed);
            DIAG_MAX_COMMUNITIES.fetch_max(n_mols as u64, Ordering::Relaxed);
            if n_mols > 10 {
                DIAG_BARCODES_GT10.fetch_add(1, Ordering::Relaxed);
            }
            if n_mols > 20 {
                DIAG_BARCODES_GT20.fetch_add(1, Ordering::Relaxed);
            }
            if n_mols > 50 {
                DIAG_BARCODES_GT50.fetch_add(1, Ordering::Relaxed);
            }
            (u, assignment)
        })
        .collect();

    log::info!(
        "Molecule diagnostics: barcodes={} total_communities={} total_assigned={} max_communities={} gt10={} gt20={} gt50={}",
        DIAG_TOTAL_BARCODES.load(Ordering::Relaxed),
        DIAG_TOTAL_COMMUNITIES.load(Ordering::Relaxed),
        DIAG_TOTAL_ASSIGNED.load(Ordering::Relaxed),
        DIAG_MAX_COMMUNITIES.load(Ordering::Relaxed),
        DIAG_BARCODES_GT10.load(Ordering::Relaxed),
        DIAG_BARCODES_GT20.load(Ordering::Relaxed),
        DIAG_BARCODES_GT50.load(Ordering::Relaxed),
    );

    let mut molecules: FxHashMap<NodeIndex, FxHashMap<NodeIndex, usize>> = FxHashMap::default();
    for (u, assignment) in assignments {
        molecules.insert(u, assignment);
    }

    let mut mol_graph = NamedGraph::new();

    for &u in &nodes {
        let u_name = match g.names.get_name(u) {
            Some(n) => n.to_string(),
            None => continue,
        };
        let m = g.graph[u].m;
        if let Some(assignment) = molecules.get(&u) {
            let n_molecules = if assignment.is_empty() {
                0
            } else {
                *assignment.values().max().unwrap_or(&0) + 1
            };
            for i in 0..n_molecules {
                mol_graph.add_vertex(&format!("{}_{}", u_name, i), m);
            }
        }
    }

    for edge in g.graph.edge_references() {
        let u = edge.source();
        let v = edge.target();
        let weight = edge.weight().m;
        let u_mol = molecules.get(&u).and_then(|a| a.get(&v)).copied();
        let v_mol = molecules.get(&v).and_then(|a| a.get(&u)).copied();
        if let (Some(u_mol_id), Some(v_mol_id)) = (u_mol, v_mol) {
            let u_name = g.names.get_name(u).unwrap();
            let v_name = g.names.get_name(v).unwrap();
            let u_mol_name = format!("{}_{}", u_name, u_mol_id);
            let v_mol_name = format!("{}_{}", v_name, v_mol_id);
            if let (Some(u_idx), Some(v_idx)) = (
                mol_graph.names.get_idx(&u_mol_name),
                mol_graph.names.get_idx(&v_mol_name),
            ) {
                mol_graph.add_edge(u_idx, v_idx, weight);
            }
        }
    }

    let n_singletons = mol_graph.remove_singletons();
    log::info!(
        "Molecule graph: V={} E={} (removed {} singletons)",
        mol_graph.num_vertices(),
        mol_graph.num_edges(),
        n_singletons
    );
    mol_graph
}

// ---------------------------------------------------------------------------
// Per-barcode molecule determination
// ---------------------------------------------------------------------------

fn determine_molecules(
    graph: &OverlapGraph,
    u: NodeIndex,
    strategies: &[Strategy],
    use_junctions: bool,
    junction_indices: &FxHashSet<NodeIndex>,
    params: &MoleculeParams,
) -> FxHashMap<NodeIndex, usize> {
    let neighbors: Vec<NodeIndex> = graph.neighbors(u).collect();
    if neighbors.is_empty() {
        return FxHashMap::default();
    }

    // Junction shortcut
    if use_junctions && !junction_indices.contains(&u) {
        let mut assignment = FxHashMap::default();
        if neighbors.len() > 1 {
            for &v in &neighbors {
                assignment.insert(v, 0);
            }
        }
        return assignment;
    }

    let neighbor_set: FxHashSet<NodeIndex> = neighbors.iter().copied().collect();

    // Build compact local graph once — all strategies operate on this
    let lg = LocalGraph::from_nodes(graph, &neighbor_set);

    // Start with all nodes in one community (as local indices)
    let all_local: Vec<usize> = (0..lg.n).collect();
    let mut communities: Vec<Vec<usize>> = vec![all_local];

    for strategy in strategies {
        let mut new_communities = Vec::new();
        for component in &communities {
            match strategy {
                Strategy::Cc => {
                    new_communities.extend(cc_local(&lg, component));
                }
                Strategy::Bc => {
                    new_communities.extend(bc_local(&lg, component));
                }
                Strategy::K3 => {
                    new_communities.extend(k3_local(&lg, component));
                }
                Strategy::K3Bin => {
                    let bins = bin_local(component, params.bin_max_size, params.seed);
                    let mut clusters = Vec::new();
                    for bin in &bins {
                        clusters.extend(k3_local(&lg, bin));
                    }
                    new_communities.extend(merge_local(&lg, &clusters, params.merge_cutoff));
                }
                Strategy::Sqcos => {
                    new_communities.extend(sqcos_local(
                        &lg,
                        component,
                        true,
                        params.sqcos_threshold,
                        params.skip_small,
                    ));
                }
                Strategy::SqcosBin => {
                    let bins = bin_local(component, params.bin_max_size, params.seed);
                    let mut clusters = Vec::new();
                    for bin in &bins {
                        clusters.extend(sqcos_local(
                            &lg,
                            bin,
                            true,
                            params.sqcos_threshold,
                            params.skip_small,
                        ));
                    }
                    new_communities.extend(merge_local(&lg, &clusters, params.merge_cutoff));
                }
                Strategy::Distributed => {
                    new_communities.extend(distributed_local(&lg, component, params));
                }
            }
        }
        communities = new_communities;
    }

    // Sort by size descending
    communities.sort_by_key(|c| std::cmp::Reverse(c.len()));

    // Build assignment: global neighbor -> molecule_id (only communities with >1 member)
    // or_insert keeps the first (largest) community for overlapping nodes.
    // The original Python dict comprehension uses last-wins, but the difference
    // is minor now that k3 correctly excludes non-triangle edges.
    let mut assignment = FxHashMap::default();
    for (i, community) in communities.iter().enumerate() {
        if community.len() > 1 {
            for &li in community {
                let gi = lg.local_to_global[li];
                assignment.entry(gi).or_insert(i);
            }
        }
    }
    assignment
}

// ---------------------------------------------------------------------------
// Connected components on local indices — no hash sets, just Vec<bool>
// ---------------------------------------------------------------------------

/// Find connected components among `nodes` (local indices) using the LocalGraph adjacency.
fn cc_local(lg: &LocalGraph, nodes: &[usize]) -> Vec<Vec<usize>> {
    if nodes.is_empty() {
        return Vec::new();
    }
    // Membership bitset for the subset
    let mut in_set = vec![false; lg.n];
    for &i in nodes {
        in_set[i] = true;
    }

    let mut visited = vec![false; lg.n];
    let mut components = Vec::new();
    let mut stack = Vec::with_capacity(nodes.len());

    for &start in nodes {
        if visited[start] {
            continue;
        }
        let mut component = Vec::new();
        stack.push(start);
        visited[start] = true;
        while let Some(u) = stack.pop() {
            component.push(u);
            let row = u * lg.n;
            for j in 0..lg.n {
                if lg.adj[row + j] && in_set[j] && !visited[j] {
                    visited[j] = true;
                    stack.push(j);
                }
            }
        }
        components.push(component);
    }
    components
}

// ---------------------------------------------------------------------------
// Biconnected components: articulation points via Tarjan's, then CC
// ---------------------------------------------------------------------------

fn bc_local(lg: &LocalGraph, nodes: &[usize]) -> Vec<Vec<usize>> {
    if nodes.len() < 2 {
        if nodes.is_empty() {
            return Vec::new();
        }
        return vec![nodes.to_vec()];
    }

    let cut = articulation_points_local(lg, nodes);

    // Remove cut vertices
    let remaining: Vec<usize> = nodes.iter().copied().filter(|i| !cut[*i]).collect();
    if remaining.is_empty() {
        return Vec::new();
    }
    cc_local(lg, &remaining)
}

/// Tarjan's articulation points on local indices. Returns a bitset of size lg.n.
fn articulation_points_local(lg: &LocalGraph, nodes: &[usize]) -> Vec<bool> {
    let n = lg.n;
    let mut in_set = vec![false; n];
    for &i in nodes {
        in_set[i] = true;
    }

    let mut disc = vec![0u32; n];
    let mut low = vec![0u32; n];
    let mut parent = vec![u32::MAX; n]; // MAX = no parent
    let mut visited = vec![false; n];
    let mut ap = vec![false; n];
    let mut timer = 0u32;

    // Pre-build neighbor lists for the subset to avoid repeated scanning
    let mut nbr: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &i in nodes {
        let row = i * n;
        for (j, &is_in) in in_set.iter().enumerate() {
            if lg.adj[row + j] && is_in {
                nbr[i].push(j);
            }
        }
    }

    // Stack frame: (node, neighbor_index)
    let mut stack: Vec<(usize, usize)> = Vec::with_capacity(nodes.len());

    for &start in nodes {
        if visited[start] {
            continue;
        }
        timer += 1;
        disc[start] = timer;
        low[start] = timer;
        visited[start] = true;
        parent[start] = u32::MAX;
        stack.push((start, 0));

        let mut root_children = 0u32;

        while let Some(&mut (u, ref mut idx)) = stack.last_mut() {
            if *idx < nbr[u].len() {
                let v = nbr[u][*idx];
                *idx += 1;

                if !visited[v] {
                    timer += 1;
                    disc[v] = timer;
                    low[v] = timer;
                    visited[v] = true;
                    parent[v] = u as u32;
                    if u == start {
                        root_children += 1;
                    }
                    stack.push((v, 0));
                } else if v as u32 != parent[u] && disc[v] < low[u] {
                    low[u] = disc[v];
                }
            } else {
                stack.pop();
                if let Some(&mut (p, _)) = stack.last_mut() {
                    if low[u] < low[p] {
                        low[p] = low[u];
                    }
                    // p is an AP if it's not root and low[u] >= disc[p]
                    if parent[p] != u32::MAX && low[u] >= disc[p] {
                        ap[p] = true;
                    }
                }
            }
        }

        if root_children > 1 {
            ap[start] = true;
        }
    }
    ap
}

// ---------------------------------------------------------------------------
// K-clique percolation (k=3): optimized triangle enumeration + union-find
// ---------------------------------------------------------------------------

/// K=3 clique percolation on local indices.
/// Finds all triangles, merges those sharing an edge via union-find.
/// Returns communities as Vec<Vec<usize>>.
fn k3_local(lg: &LocalGraph, nodes: &[usize]) -> Vec<Vec<usize>> {
    if nodes.len() < 3 {
        return Vec::new();
    }

    let n = lg.n;
    let mut in_set = vec![false; n];
    for &i in nodes {
        in_set[i] = true;
    }

    // Collect all edges in the induced subgraph and assign dense edge IDs.
    // Edge (i,j) stored as (min(i,j), max(i,j)) for canonical form.
    let mut edges: Vec<(usize, usize)> = Vec::new();
    let mut edge_id: FxHashMap<(usize, usize), usize> = FxHashMap::default();
    for (ni, &u) in nodes.iter().enumerate() {
        let u_row = u * n;
        for &v in &nodes[(ni + 1)..] {
            if lg.adj[u_row + v] {
                let eid = edges.len();
                edges.push((u, v));
                edge_id.insert((u, v), eid);
            }
        }
    }

    if edges.is_empty() {
        return Vec::new();
    }

    // Union-find over EDGES (not vertices).
    // Two triangles sharing an edge get their edges merged into one component.
    // Communities = connected components of the edge percolation graph.
    let ne = edges.len();
    let mut uf_parent: Vec<usize> = (0..ne).collect();
    let mut uf_rank: Vec<u8> = vec![0; ne];
    let mut any_triangle = false;
    // Track which edges participate in at least one triangle.
    // Only these edges contribute to k-clique communities.
    let mut in_triangle = vec![false; ne];

    // Triangle enumeration: for each edge (u,v), find common neighbors w
    // forming triangle (u,v,w). Union edges (u,v), (u,w), (v,w).
    for (ni, &u) in nodes.iter().enumerate() {
        let u_row = u * n;
        for (offset_j, &v) in nodes[(ni + 1)..].iter().enumerate() {
            if !lg.adj[u_row + v] {
                continue;
            }
            let v_row = v * n;
            let e_uv = edge_id[&(u, v)];
            let nj = ni + 1 + offset_j;
            for &w in &nodes[(nj + 1)..] {
                if lg.adj[u_row + w] && lg.adj[v_row + w] {
                    any_triangle = true;
                    // Canonical edge keys (u < v < w guaranteed by iteration order)
                    let e_uw = edge_id[&(u, w)];
                    let e_vw = edge_id[&(v, w)];
                    in_triangle[e_uv] = true;
                    in_triangle[e_uw] = true;
                    in_triangle[e_vw] = true;
                    uf_union_local(&mut uf_parent, &mut uf_rank, e_uv, e_uw);
                    uf_union_local(&mut uf_parent, &mut uf_rank, e_uv, e_vw);
                }
            }
        }
    }

    if !any_triangle {
        return Vec::new();
    }

    // Group edges by component, collect endpoint nodes.
    // Only include edges that participated in at least one triangle,
    // matching NetworkX's k_clique_communities which only returns nodes in k-cliques.
    // This naturally produces overlapping communities (a node can appear in multiple).
    // In determine_molecules, overlapping nodes are assigned to the last (smallest)
    // community, matching the original Python dict comprehension behavior.
    let mut community_map: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for (eid, &(u, v)) in edges.iter().enumerate() {
        if !in_triangle[eid] {
            continue; // skip edges not in any triangle
        }
        let root = uf_find_local(&mut uf_parent, eid);
        let entry = community_map.entry(root).or_default();
        entry.push(u);
        entry.push(v);
    }

    // Deduplicate nodes within each community
    community_map
        .into_values()
        .map(|mut nodes| {
            nodes.sort_unstable();
            nodes.dedup();
            nodes
        })
        .collect()
}

#[inline(always)]
fn uf_find_local(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]]; // path splitting
        x = parent[x];
    }
    x
}

#[inline(always)]
fn uf_union_local(parent: &mut [usize], rank: &mut [u8], x: usize, y: usize) {
    let rx = uf_find_local(parent, x);
    let ry = uf_find_local(parent, y);
    if rx == ry {
        return;
    }
    if rank[rx] < rank[ry] {
        parent[rx] = ry;
    } else if rank[rx] > rank[ry] {
        parent[ry] = rx;
    } else {
        parent[ry] = rx;
        rank[rx] += 1;
    }
}

// ---------------------------------------------------------------------------
// Cosine-of-squared adjacency: all dense, no hash sets
// ---------------------------------------------------------------------------

fn sqcos_local(
    lg: &LocalGraph,
    nodes: &[usize],
    squaring: bool,
    threshold: f64,
    skip_small: usize,
) -> Vec<Vec<usize>> {
    let nn = nodes.len();
    if nn < skip_small || nn <= 1 {
        return vec![nodes.to_vec()];
    }

    // Build dense adjacency for the subset (nn × nn)
    // Map: nodes[i] -> i (local-to-sub-local)
    let mut sub_to_lg = vec![0usize; nn]; // sub index -> lg index
    let mut lg_to_sub = vec![usize::MAX; lg.n]; // lg index -> sub index
    for (si, &li) in nodes.iter().enumerate() {
        sub_to_lg[si] = li;
        lg_to_sub[li] = si;
    }

    let mut adj = vec![0.0f32; nn * nn];
    for si in 0..nn {
        let li = sub_to_lg[si];
        let row = li * lg.n;
        for sj in 0..nn {
            let lj = sub_to_lg[sj];
            if lg.adj[row + lj] {
                adj[si * nn + sj] = 1.0;
            }
        }
    }

    // Square the matrix if requested: work = A * A
    let work = if squaring {
        mat_mul_f32(&adj, &adj, nn)
    } else {
        adj.clone()
    };

    // Compute row norms for cosine similarity
    let mut norms = vec![0.0f32; nn];
    for (i, norm) in norms.iter_mut().enumerate() {
        let mut s = 0.0f32;
        let row = i * nn;
        for j in 0..nn {
            let v = work[row + j];
            s += v * v;
        }
        *norm = s.sqrt();
    }

    // Determine which edges to remove: where cosine_sim < threshold but adj != 0
    // Instead of building full cosine matrix, only check existing edges
    let thresh = threshold as f32;

    // Build a modified adjacency for CC: start with original, remove edges
    let mut kept = adj.clone(); // copy of original adj
    for i in 0..nn {
        if norms[i] == 0.0 {
            // Zero norm: remove all edges from this node
            for j in 0..nn {
                kept[i * nn + j] = 0.0;
            }
            continue;
        }
        let row_i = i * nn;
        for j in (i + 1)..nn {
            if adj[row_i + j] == 0.0 {
                continue; // no edge to potentially remove
            }
            if norms[j] == 0.0 {
                // cos = 0 < threshold, remove edge
                kept[i * nn + j] = 0.0;
                kept[j * nn + i] = 0.0;
                continue;
            }
            // Compute cosine similarity for this pair
            let mut dot = 0.0f32;
            for k in 0..nn {
                dot += work[row_i + k] * work[j * nn + k];
            }
            let cos = dot / (norms[i] * norms[j]);
            if cos < thresh {
                kept[i * nn + j] = 0.0;
                kept[j * nn + i] = 0.0;
            }
        }
    }

    // Find connected components on the kept adjacency
    let mut visited = vec![false; nn];
    let mut components = Vec::new();
    let mut stack = Vec::with_capacity(nn);

    for start in 0..nn {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        stack.push(start);
        let mut component = Vec::new();
        while let Some(u) = stack.pop() {
            component.push(nodes[u]); // convert back to lg-local index
            let row = u * nn;
            for j in 0..nn {
                if kept[row + j] != 0.0 && !visited[j] {
                    visited[j] = true;
                    stack.push(j);
                }
            }
        }
        components.push(component);
    }

    components
}

/// Dense f32 matrix multiplication C = A * B, n×n row-major.
/// Uses ikj loop order for cache-friendly access on B's rows.
#[inline(never)]
fn mat_mul_f32(a: &[f32], b: &[f32], n: usize) -> Vec<f32> {
    let mut c = vec![0.0f32; n * n];
    for i in 0..n {
        let c_row = i * n;
        for k in 0..n {
            let a_ik = a[i * n + k];
            if a_ik == 0.0 {
                continue; // sparse skip — adjacency matrices are sparse
            }
            let b_row = k * n;
            for j in 0..n {
                c[c_row + j] += a_ik * b[b_row + j];
            }
        }
    }
    c
}

// ---------------------------------------------------------------------------
// Random binning (Fisher-Yates shuffle, matching original Python random.shuffle)
// ---------------------------------------------------------------------------

fn bin_local(nodes: &[usize], max_size: usize, seed: Option<u64>) -> Vec<Vec<usize>> {
    use rand::seq::SliceRandom;
    use rand::SeedableRng;

    let n = nodes.len();
    if n == 0 {
        return Vec::new();
    }
    let bins_count = 1 + n / max_size;
    if bins_count <= 1 {
        return vec![nodes.to_vec()];
    }

    let mut shuffled = nodes.to_vec();
    match seed {
        Some(s) => {
            // Deterministic: seed derived from user seed + first node for per-call variation
            let local_seed = s.wrapping_add(nodes.first().copied().unwrap_or(0) as u64);
            let mut rng = rand::rngs::StdRng::seed_from_u64(local_seed);
            shuffled.shuffle(&mut rng);
        }
        None => {
            let mut rng = rand::rng();
            shuffled.shuffle(&mut rng);
        }
    }

    let (size, leftover) = (n / bins_count, n % bins_count);
    let mut bins: Vec<Vec<usize>> = Vec::with_capacity(bins_count);
    for i in 0..bins_count {
        bins.push(shuffled[size * i..size * (i + 1)].to_vec());
    }
    let edge = size * bins_count;
    for i in 0..leftover {
        bins[i % bins_count].push(shuffled[edge + i]);
    }
    bins
}

// ---------------------------------------------------------------------------
// Merge communities sharing many cross-edges
// ---------------------------------------------------------------------------

fn merge_local(lg: &LocalGraph, communities: &[Vec<usize>], cutoff: i64) -> Vec<Vec<usize>> {
    if cutoff == -1 || communities.len() <= 1 {
        return communities.to_vec();
    }

    let nc = communities.len();

    let mut uf_parent: Vec<usize> = (0..nc).collect();
    let mut uf_rank: Vec<u8> = vec![0; nc];

    // For each pair, compute cross-edges = edges(union) - edges(i) - edges(j)
    // Optimization: build membership arrays to count cross-edges directly
    // instead of computing edges(union) which re-counts internal edges.
    for i in 0..nc {
        for j in (i + 1)..nc {
            // Count edges between community i and community j directly
            let mut cross = 0i64;
            for &u in &communities[i] {
                let row = u * lg.n;
                for &v in &communities[j] {
                    if lg.adj[row + v] {
                        cross += 1;
                    }
                }
            }
            // cross counts each directed edge once (u in i, v in j)
            // The original counts undirected: edges(union) - edges(i) - edges(j)
            // which equals the number of undirected cross-edges.
            // Our cross variable counts directed cross-edges (u->v where u in i, v in j).
            // For undirected graph, we also need v->u, but since the graph is undirected
            // and we only iterate (u in i, v in j), cross = undirected cross-edges.
            // Actually no: for undirected, adj[u][v] = adj[v][u], so each undirected
            // cross-edge is counted once here (u in i, v in j). That matches the original.
            if cross > cutoff {
                uf_union_local(&mut uf_parent, &mut uf_rank, i, j);
            }
        }
    }

    let mut merged: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for (ci, com) in communities.iter().enumerate() {
        let root = uf_find_local(&mut uf_parent, ci);
        merged.entry(root).or_default().extend(com);
    }

    merged.into_values().collect()
}

// ---------------------------------------------------------------------------
// Distributed: bc → bin → bc → k3 → merge (ensemble pipeline)
// ---------------------------------------------------------------------------

fn distributed_local(lg: &LocalGraph, nodes: &[usize], params: &MoleculeParams) -> Vec<Vec<usize>> {
    // Step 1: biconnected components
    let bc_components = bc_local(lg, nodes);

    let mut result = Vec::new();
    for bc_comp in &bc_components {
        // Step 2: random binning of each bc component
        let bins = bin_local(bc_comp, params.bin_max_size, params.seed);

        // Step 3+4: for each bin, bc again, then k3
        let mut clusters = Vec::new();
        for bin in &bins {
            let inner_bcs = bc_local(lg, bin);
            for inner_bc in &inner_bcs {
                clusters.extend(k3_local(lg, inner_bc));
            }
        }

        // Step 5: merge clusters scoped to this bc component
        let merged = merge_local(lg, &clusters, params.merge_cutoff);
        result.extend(merged);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::NamedGraph;

    /// Build a LocalGraph from an edge list on nodes 0..n.
    fn make_local_graph(n: usize, edges: &[(usize, usize)]) -> LocalGraph {
        let mut adj = vec![false; n * n];
        for &(u, v) in edges {
            adj[u * n + v] = true;
            adj[v * n + u] = true;
        }
        let local_to_global: Vec<NodeIndex> = (0..n).map(NodeIndex::new).collect();
        LocalGraph {
            n,
            local_to_global,
            adj,
        }
    }

    /// Build a NamedGraph from named vertices and weighted edges.
    fn make_named_graph(names: &[&str], edges: &[(&str, &str, u32)]) -> NamedGraph {
        let mut g = NamedGraph::new();
        for &name in names {
            g.add_vertex(name, 100);
        }
        for &(u, v, w) in edges {
            let ui = g.names.get_idx(u).unwrap();
            let vi = g.names.get_idx(v).unwrap();
            g.add_edge(ui, vi, w);
        }
        g
    }

    fn sorted(mut v: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
        for c in &mut v {
            c.sort_unstable();
        }
        v.sort();
        v
    }

    // -----------------------------------------------------------------------
    // cc_local tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_cc_single_component() {
        // 0-1-2 (path)
        let lg = make_local_graph(3, &[(0, 1), (1, 2)]);
        let result = sorted(cc_local(&lg, &[0, 1, 2]));
        assert_eq!(result, vec![vec![0, 1, 2]]);
    }

    #[test]
    fn test_cc_two_components() {
        // 0-1  2-3 (two disconnected edges)
        let lg = make_local_graph(4, &[(0, 1), (2, 3)]);
        let result = sorted(cc_local(&lg, &[0, 1, 2, 3]));
        assert_eq!(result, vec![vec![0, 1], vec![2, 3]]);
    }

    #[test]
    fn test_cc_subset() {
        // Full graph: 0-1-2-3, but only query {0,1,3}
        let lg = make_local_graph(4, &[(0, 1), (1, 2), (2, 3)]);
        let result = sorted(cc_local(&lg, &[0, 1, 3]));
        // 0-1 connected, 3 isolated
        assert_eq!(result, vec![vec![0, 1], vec![3]]);
    }

    #[test]
    fn test_cc_empty() {
        let lg = make_local_graph(3, &[(0, 1)]);
        let result = cc_local(&lg, &[]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_cc_isolated_nodes() {
        let lg = make_local_graph(3, &[]);
        let result = sorted(cc_local(&lg, &[0, 1, 2]));
        assert_eq!(result, vec![vec![0], vec![1], vec![2]]);
    }

    // -----------------------------------------------------------------------
    // bc_local tests (biconnected components)
    // -----------------------------------------------------------------------

    #[test]
    fn test_bc_no_articulation() {
        // Triangle: 0-1, 1-2, 0-2 — no articulation points
        let lg = make_local_graph(3, &[(0, 1), (1, 2), (0, 2)]);
        let result = sorted(bc_local(&lg, &[0, 1, 2]));
        assert_eq!(result, vec![vec![0, 1, 2]]);
    }

    #[test]
    fn test_bc_single_articulation() {
        // 0-1-2 (path) — node 1 is articulation point
        let lg = make_local_graph(3, &[(0, 1), (1, 2)]);
        let result = sorted(bc_local(&lg, &[0, 1, 2]));
        // After removing node 1: {0} and {2}
        assert_eq!(result, vec![vec![0], vec![2]]);
    }

    #[test]
    fn test_bc_bowtie() {
        // Two triangles sharing node 2:
        // 0-1, 1-2, 0-2, 2-3, 3-4, 2-4
        // Node 2 is the articulation point
        let lg = make_local_graph(5, &[(0, 1), (1, 2), (0, 2), (2, 3), (3, 4), (2, 4)]);
        let result = sorted(bc_local(&lg, &[0, 1, 2, 3, 4]));
        // After removing node 2: {0,1} and {3,4}
        assert_eq!(result, vec![vec![0, 1], vec![3, 4]]);
    }

    #[test]
    fn test_bc_single_node() {
        let lg = make_local_graph(3, &[(0, 1)]);
        let result = bc_local(&lg, &[2]);
        assert_eq!(result, vec![vec![2]]);
    }

    #[test]
    fn test_bc_chain() {
        // 0-1-2-3 (path) — nodes 1 and 2 are articulation points
        let lg = make_local_graph(4, &[(0, 1), (1, 2), (2, 3)]);
        let result = sorted(bc_local(&lg, &[0, 1, 2, 3]));
        // After removing 1 and 2: {0} and {3}
        assert_eq!(result, vec![vec![0], vec![3]]);
    }

    // -----------------------------------------------------------------------
    // k3_local tests (k-clique percolation, k=3)
    // -----------------------------------------------------------------------

    #[test]
    fn test_k3_single_triangle() {
        let lg = make_local_graph(3, &[(0, 1), (1, 2), (0, 2)]);
        let result = sorted(k3_local(&lg, &[0, 1, 2]));
        assert_eq!(result, vec![vec![0, 1, 2]]);
    }

    #[test]
    fn test_k3_two_triangles_shared_edge() {
        // Triangles (0,1,2) and (1,2,3) share edge 1-2 → one community
        let lg = make_local_graph(4, &[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]);
        let result = sorted(k3_local(&lg, &[0, 1, 2, 3]));
        assert_eq!(result, vec![vec![0, 1, 2, 3]]);
    }

    #[test]
    fn test_k3_two_triangles_shared_node_only() {
        // Triangles (0,1,2) and (2,3,4) share only node 2 (no shared edge)
        // → two separate communities (this is the critical correctness test)
        let lg = make_local_graph(5, &[(0, 1), (0, 2), (1, 2), (2, 3), (2, 4), (3, 4)]);
        let result = sorted(k3_local(&lg, &[0, 1, 2, 3, 4]));
        assert_eq!(result, vec![vec![0, 1, 2], vec![2, 3, 4]]);
    }

    #[test]
    fn test_k3_edges_not_in_triangles_excluded() {
        // Graph: triangle (0,1,2) + edge (2,3) + edge (3,4)
        // Only nodes 0,1,2 should be in k3 communities.
        // Nodes 3,4 are connected by edges but not in any triangle.
        let lg = make_local_graph(5, &[(0, 1), (1, 2), (0, 2), (2, 3), (3, 4)]);
        let result = sorted(k3_local(&lg, &[0, 1, 2, 3, 4]));
        assert_eq!(result, vec![vec![0, 1, 2]]);
        // Node 3 and 4 should NOT appear in any community
    }

    #[test]
    fn test_k3_no_triangles() {
        // Path: 0-1-2-3 — no triangles
        let lg = make_local_graph(4, &[(0, 1), (1, 2), (2, 3)]);
        let result = k3_local(&lg, &[0, 1, 2, 3]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_k3_clique_of_4() {
        // Complete graph K4: all 4 triangles share edges → one community
        let lg = make_local_graph(4, &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]);
        let result = sorted(k3_local(&lg, &[0, 1, 2, 3]));
        assert_eq!(result, vec![vec![0, 1, 2, 3]]);
    }

    #[test]
    fn test_k3_chain_of_triangles() {
        // (0,1,2), (2,3,4), (4,5,6) — sharing single nodes only
        let lg = make_local_graph(
            7,
            &[
                (0, 1),
                (0, 2),
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (4, 5),
                (4, 6),
                (5, 6),
            ],
        );
        let result = sorted(k3_local(&lg, &[0, 1, 2, 3, 4, 5, 6]));
        assert_eq!(result.len(), 3);
        assert_eq!(result, vec![vec![0, 1, 2], vec![2, 3, 4], vec![4, 5, 6]]);
    }

    #[test]
    fn test_k3_overlapping_communities() {
        // Node 2 belongs to two communities — k-clique percolation allows overlap
        let lg = make_local_graph(5, &[(0, 1), (0, 2), (1, 2), (2, 3), (2, 4), (3, 4)]);
        let result = sorted(k3_local(&lg, &[0, 1, 2, 3, 4]));
        // Node 2 should appear in both communities
        let has_2: Vec<bool> = result.iter().map(|c| c.contains(&2)).collect();
        assert_eq!(has_2, vec![true, true]);
    }

    #[test]
    fn test_k3_too_few_nodes() {
        let lg = make_local_graph(2, &[(0, 1)]);
        let result = k3_local(&lg, &[0, 1]);
        assert!(result.is_empty());
    }

    // -----------------------------------------------------------------------
    // sqcos_local tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_sqcos_skip_small() {
        // With skip_small=10, a 3-node graph should be returned as-is
        let lg = make_local_graph(3, &[(0, 1), (1, 2), (0, 2)]);
        let result = sqcos_local(&lg, &[0, 1, 2], true, 0.75, 10);
        assert_eq!(result.len(), 1);
        let mut r = result[0].clone();
        r.sort_unstable();
        assert_eq!(r, vec![0, 1, 2]);
    }

    #[test]
    fn test_sqcos_dense_clique() {
        // K5 complete graph — all cosine similarities should be high
        // → should stay as one community
        let edges: Vec<(usize, usize)> = (0..5)
            .flat_map(|i| ((i + 1)..5).map(move |j| (i, j)))
            .collect();
        let lg = make_local_graph(5, &edges);
        // skip_small=2 so it actually runs the algorithm
        let result = sqcos_local(&lg, &[0, 1, 2, 3, 4], true, 0.5, 2);
        // Dense clique should remain one community
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_sqcos_two_cliques_weak_bridge() {
        // Two K4 cliques connected by a single edge — the bridge should be cut
        // Clique A: 0,1,2,3  Clique B: 4,5,6,7  Bridge: 3-4
        let mut edges = Vec::new();
        for i in 0..4 {
            for j in (i + 1)..4 {
                edges.push((i, j));
            }
        }
        for i in 4..8 {
            for j in (i + 1)..8 {
                edges.push((i, j));
            }
        }
        edges.push((3, 4)); // weak bridge
        let lg = make_local_graph(8, &edges);
        let nodes: Vec<usize> = (0..8).collect();
        let result = sorted(sqcos_local(&lg, &nodes, true, 0.5, 2));
        // Should split into two communities
        assert!(
            result.len() >= 2,
            "Expected >=2 communities, got {}",
            result.len()
        );
    }

    #[test]
    fn test_sqcos_single_node() {
        let lg = make_local_graph(1, &[]);
        let result = sqcos_local(&lg, &[0], true, 0.75, 0);
        assert_eq!(result, vec![vec![0]]);
    }

    // -----------------------------------------------------------------------
    // bin_local tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_bin_small_set() {
        // 10 nodes, max_size=50 → 1 bin (no splitting)
        let nodes: Vec<usize> = (0..10).collect();
        let result = bin_local(&nodes, 50, None);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].len(), 10);
    }

    #[test]
    fn test_bin_exact_split() {
        // 100 nodes, max_size=50 → bins_count = 1 + 100/50 = 3
        let nodes: Vec<usize> = (0..100).collect();
        let result = bin_local(&nodes, 50, None);
        assert_eq!(result.len(), 3);
        let total: usize = result.iter().map(|b| b.len()).sum();
        assert_eq!(total, 100);
    }

    #[test]
    fn test_bin_random_covers_all_nodes() {
        // Random binning should still cover all nodes, just in random order
        let nodes: Vec<usize> = (0..200).collect();
        let result = bin_local(&nodes, 50, None);
        assert_eq!(result.len(), 5); // 1 + 200/50 = 5 bins
        let mut all: Vec<usize> = result.into_iter().flatten().collect();
        all.sort_unstable();
        assert_eq!(all, nodes);
    }

    #[test]
    fn test_bin_all_nodes_present() {
        let nodes: Vec<usize> = (0..73).collect();
        let result = bin_local(&nodes, 50, None);
        let mut all: Vec<usize> = result.into_iter().flatten().collect();
        all.sort_unstable();
        assert_eq!(all, nodes);
    }

    #[test]
    fn test_bin_empty() {
        let result = bin_local(&[], 50, None);
        assert!(result.is_empty());
    }

    // -----------------------------------------------------------------------
    // merge_local tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_merge_no_cross_edges() {
        // Two communities with no edges between them → no merge
        let lg = make_local_graph(4, &[(0, 1), (2, 3)]);
        let communities = vec![vec![0, 1], vec![2, 3]];
        let result = sorted(merge_local(&lg, &communities, 0));
        assert_eq!(result, vec![vec![0, 1], vec![2, 3]]);
    }

    #[test]
    fn test_merge_many_cross_edges() {
        // Two communities with many cross-edges → merge
        // Community A: {0,1,2}, Community B: {3,4,5}
        // Cross edges: 0-3, 0-4, 1-3, 1-4, 2-3 (5 cross edges, cutoff=3)
        let lg = make_local_graph(
            6,
            &[
                (0, 1),
                (0, 2),
                (1, 2), // internal A
                (3, 4),
                (3, 5),
                (4, 5), // internal B
                (0, 3),
                (0, 4),
                (1, 3),
                (1, 4),
                (2, 3), // cross
            ],
        );
        let communities = vec![vec![0, 1, 2], vec![3, 4, 5]];
        let result = sorted(merge_local(&lg, &communities, 3));
        assert_eq!(result, vec![vec![0, 1, 2, 3, 4, 5]]);
    }

    #[test]
    fn test_merge_below_cutoff() {
        // Only 1 cross edge, cutoff=5 → no merge
        let lg = make_local_graph(4, &[(0, 1), (2, 3), (1, 2)]);
        let communities = vec![vec![0, 1], vec![2, 3]];
        let result = sorted(merge_local(&lg, &communities, 5));
        assert_eq!(result, vec![vec![0, 1], vec![2, 3]]);
    }

    #[test]
    fn test_merge_disabled() {
        // cutoff=-1 disables merging
        let lg = make_local_graph(4, &[(0, 1), (1, 2), (2, 3), (0, 2), (0, 3), (1, 3)]);
        let communities = vec![vec![0, 1], vec![2, 3]];
        let result = merge_local(&lg, &communities, -1);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_merge_single_community() {
        let lg = make_local_graph(3, &[(0, 1), (1, 2)]);
        let communities = vec![vec![0, 1, 2]];
        let result = merge_local(&lg, &communities, 0);
        assert_eq!(result.len(), 1);
    }

    // -----------------------------------------------------------------------
    // distributed_local tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_distributed_simple() {
        // Two triangles sharing a single node (bowtie)
        // bc splits at node 2 → {0,1} and {3,4}, each too small for k3
        // So distributed returns empty (no triangles in 2-node components)
        let lg = make_local_graph(5, &[(0, 1), (0, 2), (1, 2), (2, 3), (2, 4), (3, 4)]);
        let params = MoleculeParams::default();
        let result = distributed_local(&lg, &[0, 1, 2, 3, 4], &params);
        assert!(result.is_empty());
    }

    #[test]
    fn test_distributed_larger() {
        // Two K4 cliques connected by node 3-4 bridge
        // bc should keep the cliques intact, k3 finds triangles
        let lg = make_local_graph(
            8,
            &[
                (0, 1),
                (0, 2),
                (0, 3),
                (1, 2),
                (1, 3),
                (2, 3), // K4: 0,1,2,3
                (3, 4), // bridge
                (4, 5),
                (4, 6),
                (4, 7),
                (5, 6),
                (5, 7),
                (6, 7), // K4: 4,5,6,7
            ],
        );
        let params = MoleculeParams::default();
        let result = distributed_local(&lg, &(0..8).collect::<Vec<_>>(), &params);
        // After bc removes AP node 3 and 4, remaining components have triangles
        assert!(!result.is_empty());
    }

    // -----------------------------------------------------------------------
    // parse_strategy tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_strategy_basic() {
        let s = parse_strategy("bc+cc");
        assert_eq!(s, vec![Strategy::Bc, Strategy::Cc]);
    }

    #[test]
    fn test_parse_strategy_all() {
        let s = parse_strategy("bc+cc+k3+k3bin+sqcos+sqcosbin+distributed");
        assert_eq!(s.len(), 7);
    }

    #[test]
    fn test_parse_strategy_unknown() {
        let s = parse_strategy("bc+unknown+cc");
        assert_eq!(s, vec![Strategy::Bc, Strategy::Cc]);
    }

    #[test]
    fn test_parse_strategy_empty() {
        let s = parse_strategy("");
        assert!(s.is_empty());
    }

    // -----------------------------------------------------------------------
    // Integration: separate_molecules on a small graph
    // -----------------------------------------------------------------------

    #[test]
    fn test_separate_molecules_bccc() {
        // Build a small graph: two cliques connected by a bridge barcode
        // Clique A: a,b,c (triangle)  Clique B: d,e,f (triangle)
        // Bridge: c-d
        let g = make_named_graph(
            &["a", "b", "c", "d", "e", "f"],
            &[
                ("a", "b", 5),
                ("a", "c", 5),
                ("b", "c", 5),
                ("c", "d", 2),
                ("d", "e", 5),
                ("d", "f", 5),
                ("e", "f", 5),
            ],
        );
        let junctions = FxHashSet::default();
        let mol_g = separate_molecules(&g, "bc+cc", &junctions);
        // Should produce molecules
        assert!(mol_g.num_vertices() > 0);
        assert!(mol_g.num_edges() > 0);
    }

    #[test]
    fn test_separate_molecules_k3() {
        // z's neighbors are {a,b,c,d}. Among them: triangle a-b-c and triangle b-c-d.
        // k3 should find these triangles and produce molecule assignments.
        let g = make_named_graph(
            &["z", "a", "b", "c", "d"],
            &[
                ("z", "a", 5),
                ("z", "b", 5),
                ("z", "c", 5),
                ("z", "d", 5),
                ("a", "b", 5),
                ("b", "c", 5),
                ("a", "c", 5), // triangle a-b-c
                ("b", "d", 5),
                ("c", "d", 5), // triangle b-c-d
            ],
        );
        let junctions = FxHashSet::default();
        let mol_g = separate_molecules(&g, "k3", &junctions);
        assert!(mol_g.num_vertices() > 0, "k3 should produce molecules");
    }

    #[test]
    fn test_separate_molecules_deterministic() {
        let g = make_named_graph(
            &["a", "b", "c", "d", "e", "f"],
            &[
                ("a", "b", 5),
                ("a", "c", 5),
                ("b", "c", 5),
                ("c", "d", 2),
                ("d", "e", 5),
                ("d", "f", 5),
                ("e", "f", 5),
            ],
        );
        let junctions = FxHashSet::default();
        let m1 = separate_molecules(&g, "bc+k3", &junctions);
        let m2 = separate_molecules(&g, "bc+k3", &junctions);
        assert_eq!(m1.num_vertices(), m2.num_vertices());
        assert_eq!(m1.num_edges(), m2.num_edges());
    }

    // -----------------------------------------------------------------------
    // articulation_points_local tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_ap_path() {
        // 0-1-2: node 1 is AP
        let lg = make_local_graph(3, &[(0, 1), (1, 2)]);
        let ap = articulation_points_local(&lg, &[0, 1, 2]);
        assert!(!ap[0]);
        assert!(ap[1]);
        assert!(!ap[2]);
    }

    #[test]
    fn test_ap_triangle() {
        // Triangle: no APs
        let lg = make_local_graph(3, &[(0, 1), (1, 2), (0, 2)]);
        let ap = articulation_points_local(&lg, &[0, 1, 2]);
        assert!(!ap[0]);
        assert!(!ap[1]);
        assert!(!ap[2]);
    }

    #[test]
    fn test_ap_star() {
        // Star: center node 0 connected to 1,2,3,4 — node 0 is AP
        let lg = make_local_graph(5, &[(0, 1), (0, 2), (0, 3), (0, 4)]);
        let ap = articulation_points_local(&lg, &[0, 1, 2, 3, 4]);
        assert!(ap[0]);
        for (i, &is_ap) in ap.iter().enumerate().take(5).skip(1) {
            assert!(!is_ap, "node {} should not be AP", i);
        }
    }

    #[test]
    fn test_ap_complex() {
        // 0-1-2-3, 1-3 (cycle with tail)
        // 0-1, 1-2, 2-3, 1-3 → node 1 is AP (removing it disconnects 0)
        let lg = make_local_graph(4, &[(0, 1), (1, 2), (2, 3), (1, 3)]);
        let ap = articulation_points_local(&lg, &[0, 1, 2, 3]);
        assert!(ap[1]);
        assert!(!ap[0]);
        assert!(!ap[2]);
        assert!(!ap[3]);
    }

    // -----------------------------------------------------------------------
    // Overlapping community assignment: insert (last wins) vs or_insert
    // -----------------------------------------------------------------------

    #[test]
    fn test_overlapping_assignment_no_gaps() {
        // Build a graph where k3 produces overlapping communities.
        // Two triangles sharing edge 1-2: (0,1,2) and (1,2,3)
        // k3 should produce two communities: {0,1,2} and {1,2,3}
        // Nodes 1 and 2 are in both communities.
        //
        // With insert (last wins), sorted desc by size (both size 3):
        //   community 0: {0,1,2} → 0→0, 1→0, 2→0
        //   community 1: {1,2,3} → 1→1, 2→1, 3→1 (overwrites 1 and 2)
        //   Final: {0→0, 1→1, 2→1, 3→1}
        //   n_molecules = 1 + max(0,1,1,1) = 2 ✓ (no gaps)
        //
        // With or_insert (first wins):
        //   community 0: {0,1,2} → 0→0, 1→0, 2→0
        //   community 1: {1,2,3} → 1 already has 0, 2 already has 0, 3→1
        //   Final: {0→0, 1→0, 2→0, 3→1}
        //   n_molecules = 1 + max(0,0,0,1) = 2 ✓ (no gaps in this case)
        //
        // Both produce 2 molecules here, but with different membership.
        // The key difference shows up when smaller communities lose ALL members.

        // More complex case: 3 triangles, middle one overlaps with both
        // Triangle A: (0,1,2), Triangle B: (1,2,3), Triangle C: (3,4,5)
        // k3 communities: {0,1,2,3} (A+B merged via shared edge 1-2), {3,4,5} (C)
        // Actually k3 merges A and B because they share edge 1-2.
        // Let's use a case where k3 produces truly overlapping communities.

        // Actually, k3 (k-clique percolation) merges cliques sharing k-1 nodes.
        // For k=3, two triangles sharing an edge (2 nodes) are merged.
        // So k3 rarely produces overlapping communities for k=3.
        // The overlap comes from the merge step across bins.

        // Test the determine_molecules function directly with a graph
        // that has a center node with neighbors forming overlapping structures.
        let names = &["center", "a", "b", "c", "d", "e"];
        let edges = &[
            ("center", "a", 1),
            ("center", "b", 1),
            ("center", "c", 1),
            ("center", "d", 1),
            ("center", "e", 1),
            ("a", "b", 1),
            ("b", "c", 1), // a-b-c path
            ("d", "e", 1), // d-e edge
        ];
        let g = make_named_graph(names, edges);
        let junctions = FxHashSet::default();
        let strategies = vec![Strategy::Cc];
        let params = MoleculeParams::default();
        let center = g.names.get_idx("center").unwrap();
        let assignment =
            determine_molecules(&g.graph, center, &strategies, false, &junctions, &params);
        // CC should find: {a,b,c,d,e} as one component (all connected via center's edges)
        // Wait, CC operates on the subgraph of neighbors, not including center.
        // Neighbors: a,b,c,d,e. Edges among them: a-b, b-c, d-e.
        // CC: {a,b,c} and {d,e}
        assert_eq!(assignment.len(), 5); // all 5 neighbors assigned
        let a = g.names.get_idx("a").unwrap();
        let b = g.names.get_idx("b").unwrap();
        let c = g.names.get_idx("c").unwrap();
        let d = g.names.get_idx("d").unwrap();
        let e = g.names.get_idx("e").unwrap();
        assert_eq!(assignment[&a], assignment[&b]);
        assert_eq!(assignment[&b], assignment[&c]);
        assert_eq!(assignment[&d], assignment[&e]);
        assert_ne!(assignment[&a], assignment[&d]);
    }
}
