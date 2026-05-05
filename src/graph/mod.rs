/// Graph types and algorithms for the Physlr overlap graph.
///
/// Uses petgraph's `StableUnGraph` so node indices remain valid after removals.
use petgraph::stable_graph::{NodeIndex, StableUnGraph};
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use rustc_hash::{FxHashMap, FxHashSet};

/// Properties stored on each vertex.
#[derive(Debug, Clone)]
pub struct VertexProp {
    pub m: u32,
    pub mol: Option<u32>,
}

/// Properties stored on each edge.
#[derive(Debug, Clone, Copy)]
pub struct EdgeProp {
    pub m: u32,
}

/// The core overlap graph type.
pub type OverlapGraph = StableUnGraph<VertexProp, EdgeProp>;

/// Bidirectional mapping between vertex names and graph node indices.
#[derive(Debug, Clone)]
pub struct NameIndex {
    pub name_to_idx: FxHashMap<String, NodeIndex>,
    pub idx_to_name: FxHashMap<NodeIndex, String>,
}

impl NameIndex {
    pub fn new() -> Self {
        Self {
            name_to_idx: FxHashMap::default(),
            idx_to_name: FxHashMap::default(),
        }
    }

    pub fn insert(&mut self, name: String, idx: NodeIndex) {
        self.name_to_idx.insert(name.clone(), idx);
        self.idx_to_name.insert(idx, name);
    }

    pub fn get_idx(&self, name: &str) -> Option<NodeIndex> {
        self.name_to_idx.get(name).copied()
    }

    pub fn get_name(&self, idx: NodeIndex) -> Option<&str> {
        self.idx_to_name.get(&idx).map(|s| s.as_str())
    }

    pub fn len(&self) -> usize {
        self.name_to_idx.len()
    }

    pub fn is_empty(&self) -> bool {
        self.name_to_idx.is_empty()
    }
}

impl Default for NameIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// A graph with its name index.
#[derive(Debug, Clone)]
pub struct NamedGraph {
    pub graph: OverlapGraph,
    pub names: NameIndex,
}

impl NamedGraph {
    pub fn new() -> Self {
        Self {
            graph: OverlapGraph::default(),
            names: NameIndex::new(),
        }
    }

    /// Add a vertex, returning its index. If the name already exists, return the existing index.
    pub fn add_vertex(&mut self, name: &str, m: u32) -> NodeIndex {
        if let Some(idx) = self.names.get_idx(name) {
            return idx;
        }
        let idx = self.graph.add_node(VertexProp { m, mol: None });
        self.names.insert(name.to_string(), idx);
        idx
    }

    /// Add an edge between two named vertices.
    pub fn add_edge(&mut self, u: NodeIndex, v: NodeIndex, m: u32) {
        self.graph.add_edge(u, v, EdgeProp { m });
    }

    pub fn num_vertices(&self) -> usize {
        self.graph.node_count()
    }

    pub fn num_edges(&self) -> usize {
        self.graph.edge_count()
    }

    /// Remove isolated vertices (degree 0). Returns count removed.
    pub fn remove_singletons(&mut self) -> usize {
        let singletons: Vec<NodeIndex> = self
            .graph
            .node_indices()
            .filter(|&n| self.graph.edges(n).next().is_none())
            .collect();
        let count = singletons.len();
        for idx in &singletons {
            if let Some(name) = self.names.idx_to_name.remove(idx) {
                self.names.name_to_idx.remove(&name);
            }
            self.graph.remove_node(*idx);
        }
        count
    }

    /// Remove edges with weight < min_m. Returns count removed.
    pub fn filter_edges(&mut self, min_m: u32) -> usize {
        if min_m == 0 {
            return 0;
        }
        let to_remove: Vec<petgraph::stable_graph::EdgeIndex> = self
            .graph
            .edge_indices()
            .filter(|&e| self.graph[e].m < min_m)
            .collect();
        let count = to_remove.len();
        for e in to_remove {
            self.graph.remove_edge(e);
        }
        count
    }

    /// Remove components smaller than min_size. Returns (components_removed, vertices_removed).
    pub fn remove_small_components(&mut self, min_size: usize) -> (usize, usize) {
        if min_size < 2 {
            return (0, 0);
        }
        let components = connected_components(&self.graph);
        let mut to_remove = Vec::new();
        let mut comp_count = 0;
        for comp in &components {
            if comp.len() < min_size {
                to_remove.extend(comp);
                comp_count += 1;
            }
        }
        let vert_count = to_remove.len();
        for idx in &to_remove {
            if let Some(name) = self.names.idx_to_name.remove(idx) {
                self.names.name_to_idx.remove(&name);
            }
            self.graph.remove_node(*idx);
        }
        (comp_count, vert_count)
    }
}

impl Default for NamedGraph {
    fn default() -> Self {
        Self::new()
    }
}

/// Find connected components of an undirected graph.
pub fn connected_components(graph: &OverlapGraph) -> Vec<Vec<NodeIndex>> {
    let mut visited = FxHashSet::default();
    let mut components = Vec::new();

    for node in graph.node_indices() {
        if visited.contains(&node) {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![node];
        while let Some(n) = stack.pop() {
            if !visited.insert(n) {
                continue;
            }
            component.push(n);
            for neighbor in graph.neighbors(n) {
                if !visited.contains(&neighbor) {
                    stack.push(neighbor);
                }
            }
        }
        components.push(component);
    }
    components
}

/// Compute the maximum spanning tree using Kruskal's algorithm.
/// Returns a new graph containing only the MST edges, with a node mapping.
pub fn maximum_spanning_tree(
    graph: &OverlapGraph,
) -> (OverlapGraph, FxHashMap<NodeIndex, NodeIndex>) {
    let mut edges: Vec<(NodeIndex, NodeIndex, u32)> = graph
        .edge_references()
        .map(|e| (e.source(), e.target(), e.weight().m))
        .collect();
    edges.sort_by_key(|e| std::cmp::Reverse(e.2));

    // Assign contiguous IDs for union-find
    let node_list: Vec<NodeIndex> = graph.node_indices().collect();
    let mut node_to_uf: FxHashMap<NodeIndex, usize> = FxHashMap::default();
    for (i, &n) in node_list.iter().enumerate() {
        node_to_uf.insert(n, i);
    }
    let n = node_list.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank = vec![0u32; n];

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut [usize], rank: &mut [u32], x: usize, y: usize) -> bool {
        let rx = find(parent, x);
        let ry = find(parent, y);
        if rx == ry {
            return false;
        }
        if rank[rx] < rank[ry] {
            parent[rx] = ry;
        } else if rank[rx] > rank[ry] {
            parent[ry] = rx;
        } else {
            parent[ry] = rx;
            rank[rx] += 1;
        }
        true
    }

    let mut mst = OverlapGraph::default();
    let mut node_map: FxHashMap<NodeIndex, NodeIndex> = FxHashMap::default();
    for &node in &node_list {
        let new_idx = mst.add_node(graph[node].clone());
        node_map.insert(node, new_idx);
    }

    for (u, v, w) in edges {
        let uf_u = node_to_uf[&u];
        let uf_v = node_to_uf[&v];
        if union(&mut parent, &mut rank, uf_u, uf_v) {
            mst.add_edge(node_map[&u], node_map[&v], EdgeProp { m: w });
        }
    }

    (mst, node_map)
}

/// Measure branch lengths in a tree using BFS-based message passing.
pub fn measure_branch_lengths(
    tree: &OverlapGraph,
    component: &[NodeIndex],
) -> FxHashMap<(NodeIndex, NodeIndex), usize> {
    if component.is_empty() {
        return FxHashMap::default();
    }

    let mut messages: FxHashMap<(NodeIndex, NodeIndex), usize> = FxHashMap::default();
    let root = component[0];

    // BFS ordering
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(root);
    let mut bfs_visited = FxHashSet::default();
    bfs_visited.insert(root);
    let mut bfs_edges = Vec::new();

    while let Some(node) = queue.pop_front() {
        for neighbor in tree.neighbors(node) {
            if bfs_visited.insert(neighbor) {
                bfs_edges.push((node, neighbor));
                queue.push_back(neighbor);
            }
        }
    }

    // Gather: leaves to root
    for &(parent, child) in bfs_edges.iter().rev() {
        let degree = tree.neighbors(child).count();
        if degree == 1 {
            messages.insert((parent, child), 1);
        } else {
            let max_branch = tree
                .neighbors(child)
                .filter(|&n| n != parent)
                .map(|n| *messages.get(&(child, n)).unwrap_or(&0))
                .max()
                .unwrap_or(0);
            messages.insert((parent, child), 1 + max_branch);
        }
    }

    // Distribute: root to leaves
    for &(parent, child) in bfs_edges.iter() {
        let max_branch = tree
            .neighbors(parent)
            .filter(|&n| n != child)
            .map(|n| *messages.get(&(parent, n)).unwrap_or(&0))
            .max()
            .unwrap_or(0);
        if tree.neighbors(parent).count() == 1 {
            messages.insert((child, parent), 1);
        } else {
            messages.insert((child, parent), 1 + max_branch);
        }
    }

    messages
}

/// Find the diameter path of a tree using weighted distances.
///
/// Matches original: `diameter_of_tree(g, weight="m")` + `nx.shortest_path(g, u, v, weight="m")`
/// Uses Dijkstra-like BFS with edge weight "m" to find the two farthest nodes,
/// then traces the path between them.
pub fn diameter_path(tree: &OverlapGraph, component: &[NodeIndex]) -> Vec<NodeIndex> {
    if component.is_empty() {
        return Vec::new();
    }
    if component.len() == 1 {
        return vec![component[0]];
    }

    let comp_set: FxHashSet<NodeIndex> = component.iter().copied().collect();
    let start = component[0];
    let (far1, _) = weighted_farthest(tree, start, &comp_set);
    let (far2, _) = weighted_farthest(tree, far1, &comp_set);
    tree_path(tree, far1, far2, &comp_set)
}

/// Find the farthest node from `start` using weighted distances (edge weight = m).
/// In a tree, BFS with accumulated weights is equivalent to Dijkstra.
fn weighted_farthest(
    tree: &OverlapGraph,
    start: NodeIndex,
    valid: &FxHashSet<NodeIndex>,
) -> (NodeIndex, u64) {
    let mut visited = FxHashSet::default();
    let mut queue = std::collections::VecDeque::new();
    queue.push_back((start, 0u64));
    visited.insert(start);
    let mut farthest = (start, 0u64);

    while let Some((node, dist)) = queue.pop_front() {
        if dist > farthest.1 {
            farthest = (node, dist);
        }
        for edge in tree.edges(node) {
            let neighbor = if edge.source() == node {
                edge.target()
            } else {
                edge.source()
            };
            if valid.contains(&neighbor) && visited.insert(neighbor) {
                queue.push_back((neighbor, dist + edge.weight().m as u64));
            }
        }
    }
    farthest
}

/// Find the unique path between two nodes in a tree.
fn tree_path(
    tree: &OverlapGraph,
    start: NodeIndex,
    end: NodeIndex,
    valid: &FxHashSet<NodeIndex>,
) -> Vec<NodeIndex> {
    let mut visited = FxHashSet::default();
    let mut parent_map: FxHashMap<NodeIndex, NodeIndex> = FxHashMap::default();
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(start);
    visited.insert(start);

    while let Some(node) = queue.pop_front() {
        if node == end {
            break;
        }
        for neighbor in tree.neighbors(node) {
            if valid.contains(&neighbor) && visited.insert(neighbor) {
                parent_map.insert(neighbor, node);
                queue.push_back(neighbor);
            }
        }
    }

    let mut path = Vec::new();
    let mut current = end;
    path.push(current);
    while current != start {
        if let Some(&p) = parent_map.get(&current) {
            path.push(p);
            current = p;
        } else {
            return Vec::new();
        }
    }
    path.reverse();
    path
}

/// Prune short branches from a tree. Returns count of pruned vertices.
///
/// A branch is "short" if it hangs off a node with degree >= 3 and has
/// fewer than `min_branch_size` nodes.
pub fn prune_branches(tree: &mut OverlapGraph, min_branch_size: usize) -> usize {
    if min_branch_size == 0 {
        return 0;
    }

    let mut total_pruned = 0;
    loop {
        // Find leaves
        let leaves: Vec<NodeIndex> = tree
            .node_indices()
            .filter(|&n| tree.neighbors(n).count() == 1)
            .collect();

        let mut pruned_this_round = 0;

        for leaf in leaves {
            // Walk from leaf toward the interior, collecting the branch
            let mut branch = Vec::new();
            let mut current = leaf;
            loop {
                let degree = tree.neighbors(current).count();
                if degree == 0 {
                    // Already removed
                    branch.push(current);
                    break;
                }
                if degree >= 3 {
                    // Reached a junction — don't include it in the branch
                    break;
                }
                branch.push(current);
                if degree == 1 && current != leaf {
                    // Reached another leaf (the other end of a path component)
                    break;
                }
                // degree == 2 or (degree == 1 and current == leaf): continue walking
                let neighbors: Vec<NodeIndex> = tree.neighbors(current).collect();
                let next = if neighbors.len() == 1 {
                    neighbors[0]
                } else {
                    // degree 2: pick the neighbor not already in branch
                    *neighbors
                        .iter()
                        .find(|&&n| !branch.contains(&n))
                        .unwrap_or(&neighbors[0])
                };
                if branch.contains(&next) {
                    break;
                }
                current = next;
            }

            // Check if this branch is short enough to prune
            let junction_degree = tree.neighbors(current).count();
            if branch.len() < min_branch_size && junction_degree >= 3 {
                for &node in &branch {
                    tree.remove_node(node);
                    pruned_this_round += 1;
                }
            }
        }

        total_pruned += pruned_this_round;
        if pruned_this_round == 0 {
            break;
        }
    }
    total_pruned
}

/// Collect all nodes reachable from `start` via BFS.
pub fn bfs_collect(graph: &OverlapGraph, start: NodeIndex) -> Vec<NodeIndex> {
    let mut visited = Vec::new();
    let mut queue = std::collections::VecDeque::new();
    let mut seen = FxHashSet::default();
    queue.push_back(start);
    seen.insert(start);

    while let Some(node) = queue.pop_front() {
        visited.push(node);
        for neighbor in graph.neighbors(node) {
            if seen.insert(neighbor) {
                queue.push_back(neighbor);
            }
        }
    }
    visited
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_graph(names: &[&str], edges: &[(&str, &str, u32)]) -> NamedGraph {
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

    // -----------------------------------------------------------------------
    // NamedGraph basic operations
    // -----------------------------------------------------------------------

    #[test]
    fn test_add_vertex() {
        let mut g = NamedGraph::new();
        let idx = g.add_vertex("a", 10);
        assert_eq!(g.num_vertices(), 1);
        assert_eq!(g.names.get_name(idx), Some("a"));
    }

    #[test]
    fn test_add_vertex_duplicate() {
        let mut g = NamedGraph::new();
        let idx1 = g.add_vertex("a", 10);
        let idx2 = g.add_vertex("a", 20);
        assert_eq!(idx1, idx2);
        assert_eq!(g.num_vertices(), 1);
    }

    #[test]
    fn test_add_edge() {
        let mut g = NamedGraph::new();
        let a = g.add_vertex("a", 10);
        let b = g.add_vertex("b", 20);
        g.add_edge(a, b, 5);
        assert_eq!(g.num_edges(), 1);
    }

    #[test]
    fn test_remove_singletons() {
        let mut g = NamedGraph::new();
        let a = g.add_vertex("a", 10);
        let b = g.add_vertex("b", 20);
        g.add_vertex("c", 30); // singleton
        g.add_edge(a, b, 5);
        let removed = g.remove_singletons();
        assert_eq!(removed, 1);
        assert_eq!(g.num_vertices(), 2);
    }

    #[test]
    fn test_filter_edges() {
        let mut g = make_graph(
            &["a", "b", "c"],
            &[("a", "b", 3), ("b", "c", 7)],
        );
        let removed = g.filter_edges(5);
        assert_eq!(removed, 1);
        assert_eq!(g.num_edges(), 1);
    }

    #[test]
    fn test_filter_edges_zero() {
        let mut g = make_graph(&["a", "b"], &[("a", "b", 1)]);
        let removed = g.filter_edges(0);
        assert_eq!(removed, 0);
    }

    #[test]
    fn test_remove_small_components() {
        let mut g = make_graph(
            &["a", "b", "c", "d", "e"],
            &[("a", "b", 5), ("b", "c", 5), ("d", "e", 5)],
        );
        let (comp_removed, vert_removed) = g.remove_small_components(3);
        assert_eq!(comp_removed, 1);
        assert_eq!(vert_removed, 2);
        assert_eq!(g.num_vertices(), 3);
    }

    // -----------------------------------------------------------------------
    // Connected components
    // -----------------------------------------------------------------------

    #[test]
    fn test_connected_components_single() {
        let g = make_graph(&["a", "b", "c"], &[("a", "b", 5), ("b", "c", 5)]);
        let cc = connected_components(&g.graph);
        assert_eq!(cc.len(), 1);
        assert_eq!(cc[0].len(), 3);
    }

    #[test]
    fn test_connected_components_multiple() {
        let g = make_graph(
            &["a", "b", "c", "d"],
            &[("a", "b", 5), ("c", "d", 5)],
        );
        let cc = connected_components(&g.graph);
        assert_eq!(cc.len(), 2);
    }

    #[test]
    fn test_connected_components_isolated() {
        let mut g = NamedGraph::new();
        g.add_vertex("a", 10);
        g.add_vertex("b", 20);
        g.add_vertex("c", 30);
        let cc = connected_components(&g.graph);
        assert_eq!(cc.len(), 3);
    }

    // -----------------------------------------------------------------------
    // Maximum spanning tree
    // -----------------------------------------------------------------------

    #[test]
    fn test_mst_triangle() {
        // Triangle with weights 3, 5, 7 → MST picks edges 7 and 5
        let g = make_graph(
            &["a", "b", "c"],
            &[("a", "b", 3), ("b", "c", 5), ("a", "c", 7)],
        );
        let (mst, _) = maximum_spanning_tree(&g.graph);
        assert_eq!(mst.node_count(), 3);
        assert_eq!(mst.edge_count(), 2);
        // Total weight should be 12 (7 + 5)
        let total: u32 = mst.edge_references().map(|e| e.weight().m).sum();
        assert_eq!(total, 12);
    }

    #[test]
    fn test_mst_path() {
        // Already a tree: a-b-c
        let g = make_graph(&["a", "b", "c"], &[("a", "b", 5), ("b", "c", 3)]);
        let (mst, _) = maximum_spanning_tree(&g.graph);
        assert_eq!(mst.edge_count(), 2);
    }

    #[test]
    fn test_mst_disconnected() {
        // Two components → MST is a forest
        let g = make_graph(
            &["a", "b", "c", "d"],
            &[("a", "b", 5), ("c", "d", 3)],
        );
        let (mst, _) = maximum_spanning_tree(&g.graph);
        assert_eq!(mst.edge_count(), 2);
        assert_eq!(mst.node_count(), 4);
    }

    // -----------------------------------------------------------------------
    // Diameter path
    // -----------------------------------------------------------------------

    #[test]
    fn test_diameter_path_line() {
        // Path: a(w=10)-b(w=5)-c(w=3)-d
        let g = make_graph(
            &["a", "b", "c", "d"],
            &[("a", "b", 10), ("b", "c", 5), ("c", "d", 3)],
        );
        let nodes: Vec<NodeIndex> = g.graph.node_indices().collect();
        let path = diameter_path(&g.graph, &nodes);
        assert_eq!(path.len(), 4);
    }

    #[test]
    fn test_diameter_path_single() {
        let mut g = NamedGraph::new();
        let a = g.add_vertex("a", 10);
        let path = diameter_path(&g.graph, &[a]);
        assert_eq!(path, vec![a]);
    }

    #[test]
    fn test_diameter_path_star() {
        // Star: center b connected to a, c, d with different weights
        let g = make_graph(
            &["a", "b", "c", "d"],
            &[("b", "a", 10), ("b", "c", 5), ("b", "d", 3)],
        );
        let nodes: Vec<NodeIndex> = g.graph.node_indices().collect();
        let path = diameter_path(&g.graph, &nodes);
        // Diameter should be a-b-c (weight 15) or a-b-d (weight 13)
        // → a-b-c is the diameter
        assert_eq!(path.len(), 3);
    }

    // -----------------------------------------------------------------------
    // Branch pruning
    // -----------------------------------------------------------------------

    #[test]
    fn test_prune_no_branches() {
        // Simple path: nothing to prune (no degree-3 junctions)
        let g = make_graph(&["a", "b", "c"], &[("a", "b", 5), ("b", "c", 5)]);
        let mut tree = g.graph;
        let pruned = prune_branches(&mut tree, 5);
        assert_eq!(pruned, 0);
    }

    #[test]
    fn test_prune_short_branch() {
        // T-shape: a-b-c with d hanging off b
        // Branch d-b has length 1 < min_branch_size=2
        let g = make_graph(
            &["a", "b", "c", "d"],
            &[("a", "b", 5), ("b", "c", 5), ("b", "d", 5)],
        );
        let mut tree = g.graph;
        let pruned = prune_branches(&mut tree, 2);
        // d should be pruned (branch length 1 < 2, hanging off junction b)
        assert_eq!(pruned, 1);
        assert_eq!(tree.node_count(), 3);
    }

    #[test]
    fn test_prune_disabled() {
        let g = make_graph(
            &["a", "b", "c", "d"],
            &[("a", "b", 5), ("b", "c", 5), ("b", "d", 5)],
        );
        let mut tree = g.graph;
        let pruned = prune_branches(&mut tree, 0);
        assert_eq!(pruned, 0);
    }

    // -----------------------------------------------------------------------
    // Branch length measurement
    // -----------------------------------------------------------------------

    #[test]
    fn test_measure_branch_lengths_path() {
        let g = make_graph(&["a", "b", "c"], &[("a", "b", 5), ("b", "c", 5)]);
        let nodes: Vec<NodeIndex> = g.graph.node_indices().collect();
        let messages = measure_branch_lengths(&g.graph, &nodes);
        // Should have messages for all directed edges
        assert!(!messages.is_empty());
    }

    #[test]
    fn test_measure_branch_lengths_empty() {
        let g = NamedGraph::new();
        let messages = measure_branch_lengths(&g.graph, &[]);
        assert!(messages.is_empty());
    }

    // -----------------------------------------------------------------------
    // BFS collect
    // -----------------------------------------------------------------------

    #[test]
    fn test_bfs_collect() {
        let g = make_graph(
            &["a", "b", "c", "d"],
            &[("a", "b", 5), ("b", "c", 5), ("c", "d", 5)],
        );
        let start = g.names.get_idx("a").unwrap();
        let collected = bfs_collect(&g.graph, start);
        assert_eq!(collected.len(), 4);
    }

    #[test]
    fn test_bfs_collect_disconnected() {
        let g = make_graph(
            &["a", "b", "c", "d"],
            &[("a", "b", 5), ("c", "d", 5)],
        );
        let start = g.names.get_idx("a").unwrap();
        let collected = bfs_collect(&g.graph, start);
        assert_eq!(collected.len(), 2); // only a's component
    }

    // -----------------------------------------------------------------------
    // NameIndex
    // -----------------------------------------------------------------------

    #[test]
    fn test_name_index() {
        let mut ni = NameIndex::new();
        let idx = NodeIndex::new(0);
        ni.insert("foo".to_string(), idx);
        assert_eq!(ni.get_idx("foo"), Some(idx));
        assert_eq!(ni.get_name(idx), Some("foo"));
        assert_eq!(ni.len(), 1);
        assert!(!ni.is_empty());
    }

    #[test]
    fn test_name_index_missing() {
        let ni = NameIndex::new();
        assert_eq!(ni.get_idx("missing"), None);
        assert_eq!(ni.get_name(NodeIndex::new(99)), None);
        assert!(ni.is_empty());
    }
}
