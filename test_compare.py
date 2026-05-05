#!/usr/bin/env python3
"""
Compare original Physlr molecule separation with physlr-next on a small graph.
Tests the distributed+sqcosbin strategy.
"""
import sys
sys.path.insert(0, '/workspaces/physlr')

import networkx as nx
import numpy as np
import random

# Reproduce the original Physlr functions
def detect_communities_biconnected_components(g, node_set):
    cut_vertices = set(nx.articulation_points(g.subgraph(node_set)))
    components = list(nx.connected_components(g.subgraph(node_set - cut_vertices)))
    return components

def detect_communities_k_clique(g, node_set, k):
    return list(nx.algorithms.community.k_clique_communities(g.subgraph(node_set), k))

def partition_subgraph_into_bins_randomly(node_set, max_size=50):
    bins_count = 1 + len(node_set) // max_size
    node_list = list(node_set)
    random.shuffle(node_list)
    size, leftover = divmod(len(node_set), bins_count)
    bins = [node_list[0 + size * i: size * (i + 1)] for i in range(bins_count)]
    edge = size * bins_count
    for i in range(leftover):
        bins[i % bins_count].append(node_list[edge + i])
    return [set(x) for x in bins]

def merge_communities(g, communities, node_set=0, strategy=0, cutoff=20):
    mode = 1
    if cutoff == -1:
        return communities
    if len(communities) == 1 and (node_set == 0 or strategy != 1):
        return communities
    if strategy == 1:
        raise NotImplementedError("Louvain merge not implemented here")
    merge_network = nx.Graph()
    for i in range(len(communities)):
        merge_network.add_node(i)
    for i, com1 in enumerate(communities):
        for j, com2 in enumerate(communities):
            if i >= j:
                continue
            if mode == 1:
                cross = (nx.number_of_edges(g.subgraph(com1.union(com2))) -
                         nx.number_of_edges(g.subgraph(com1)) -
                         nx.number_of_edges(g.subgraph(com2)))
                if cross > cutoff:
                    merge_network.add_edge(i, j)
    return [{barcode for j in i for barcode in communities[j]}
            for i in nx.connected_components(merge_network)]

def determine_molecules_partition_split_merge(g, component):
    return [merged
            for bi_connected_component in
            detect_communities_biconnected_components(g, component)
            for merged in
            merge_communities(
                g, [cluster
                    for bin_set in
                    partition_subgraph_into_bins_randomly(bi_connected_component)
                    for bi_con2 in
                    detect_communities_biconnected_components(g, bin_set)
                    for cluster in
                    detect_communities_k_clique(g, bi_con2, k=3)],
                bi_connected_component, strategy=0)
            ]

def detect_communities_cosine_of_squared(g, node_set, squaring=True, threshold=0.75):
    from sklearn.metrics.pairwise import cosine_similarity
    import scipy as sp

    if len(node_set) < 10:  # skip_small default
        return [set(node_set)]
    if len(node_set) == 1:
        return [set(node_set)]

    communities = []
    if len(node_set) > 1:
        adj_array = nx.adjacency_matrix(g.subgraph(node_set)).toarray()
        if squaring:
            new_adj = np.multiply(
                cosine_similarity(
                    sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) >= threshold, adj_array)
        else:
            new_adj = np.multiply(cosine_similarity(adj_array) >= threshold, adj_array)
        edges_to_remove = np.argwhere(new_adj != adj_array)
        barcode_dict = dict(zip(range(len(node_set)), list(node_set)))
        edges_to_remove_barcode = [(barcode_dict[i], barcode_dict[j])
                                   for i, j in edges_to_remove]
        sub_graph_copy = nx.Graph(g.subgraph(node_set))
        sub_graph_copy.remove_edges_from(edges_to_remove_barcode)
        cos_components = list(nx.connected_components(sub_graph_copy))
        for com in cos_components:
            communities.append(com)
    return communities

def run_distributed_sqcosbin(g, u, sqcost=0.75):
    """Run the full distributed+sqcosbin pipeline for vertex u."""
    neighbors = set(g[u].keys())
    
    # Step 1: distributed
    communities = [neighbors]
    communities_temp = []
    for component in communities:
        communities_temp.extend(
            determine_molecules_partition_split_merge(g, component))
    communities = communities_temp
    
    # Step 2: sqcosbin
    communities_temp = []
    for component in communities:
        communities_temp.extend(
            merge_communities(
                g, [cluster
                    for bin_set in
                    partition_subgraph_into_bins_randomly(component)
                    for cluster in
                    detect_communities_cosine_of_squared(
                        g, bin_set, squaring=True, threshold=sqcost)
                    ]
            )
        )
    communities = communities_temp
    
    communities.sort(key=len, reverse=True)
    assignment = {}
    for i, vs in enumerate(communities):
        if len(vs) > 1:
            for v in vs:
                if v not in assignment:
                    assignment[v] = i
    return assignment

# Build a test graph
random.seed(42)
g = nx.Graph()

# Create a barcode with ~30 neighbors that form 2-3 molecules
center = "BX_CENTER"
g.add_node(center, m=100)

# Molecule 1: 15 nodes, well-connected
mol1 = [f"M1_{i}" for i in range(15)]
for n in mol1:
    g.add_node(n, m=50)
    g.add_edge(center, n, m=20)
for i in range(len(mol1)):
    for j in range(i+1, len(mol1)):
        if random.random() < 0.5:
            g.add_edge(mol1[i], mol1[j], m=10)

# Molecule 2: 15 nodes, well-connected
mol2 = [f"M2_{i}" for i in range(15)]
for n in mol2:
    g.add_node(n, m=50)
    g.add_edge(center, n, m=20)
for i in range(len(mol2)):
    for j in range(i+1, len(mol2)):
        if random.random() < 0.5:
            g.add_edge(mol2[i], mol2[j], m=10)

# A few cross-edges between molecules (chimeric)
for _ in range(3):
    g.add_edge(random.choice(mol1), random.choice(mol2), m=5)

print(f"Graph: {g.number_of_nodes()} nodes, {g.number_of_edges()} edges")
print(f"Center has {len(list(g.neighbors(center)))} neighbors")
print()

# Run original
assignment = run_distributed_sqcosbin(g, center)
n_mols = max(assignment.values()) + 1 if assignment else 0
print(f"Original distributed+sqcosbin: {n_mols} molecules")
for mol_id in range(n_mols):
    members = [v for v, m in assignment.items() if m == mol_id]
    print(f"  Molecule {mol_id}: {len(members)} members")
    m1_count = sum(1 for m in members if m.startswith("M1_"))
    m2_count = sum(1 for m in members if m.startswith("M2_"))
    print(f"    M1: {m1_count}, M2: {m2_count}")

unassigned = set(g.neighbors(center)) - set(assignment.keys())
print(f"  Unassigned: {len(unassigned)}")
