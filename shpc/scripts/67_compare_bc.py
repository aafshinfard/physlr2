#!/usr/bin/env python3
"""
Compare biconnected components between Python and Rust for specific barcodes.
Outputs the neighbor subgraph in a format that can be loaded by Rust for comparison.
"""
import sys
import os
import timeit

PHYSLR_DIR = "/home/afshinfa/miniconda3/envs/physlr-original/bin/share/physlr-1.0.4-8"
sys.path.insert(0, PHYSLR_DIR)

import networkx as nx
from physlr.physlr import Physlr


def read_overlap_graph(tsv_file):
    g = nx.Graph()
    header_skipped = False
    with open(tsv_file) as f:
        in_edges = False
        for line in f:
            line = line.strip()
            if not line:
                in_edges = True
                header_skipped = False
                continue
            if not header_skipped:
                header_skipped = True
                if line.startswith('U\t'):
                    continue
            parts = line.split('\t')
            if not in_edges:
                barcode = parts[0]
                g.add_node(barcode)
            else:
                if len(parts) >= 3:
                    u, v = parts[0], parts[1]
                    try:
                        m = int(parts[2])
                    except ValueError:
                        continue
                    g.add_edge(u, v, m=m)
    return g


def analyze_barcode(g, u):
    """Analyze biconnected components for a barcode's neighbor subgraph."""
    neighbors = set(g[u].keys())
    sub = g.subgraph(neighbors)

    print(f"\nBarcode: {u}")
    print(f"Degree: {len(neighbors)}")
    print(f"Edges in neighbor subgraph: {sub.number_of_edges()}")

    # Articulation points
    aps = set(nx.articulation_points(sub))
    print(f"Articulation points: {len(aps)}")

    # Connected components after removing APs
    remaining = neighbors - aps
    sub_no_ap = g.subgraph(remaining)
    components = list(nx.connected_components(sub_no_ap))
    sizes = sorted([len(c) for c in components], reverse=True)
    print(f"BC components: {len(components)}")
    print(f"BC sizes: {sizes[:20]}")
    print(f"Total nodes in BC: {sum(sizes)}")

    # For the largest BC component, do the inner pipeline
    if components:
        largest = max(components, key=len)
        print(f"\nLargest BC component (size {len(largest)}):")

        # Partition into bins
        import random
        random.seed(42)
        bins = Physlr.partition_subgraph_into_bins_randomly(largest)
        print(f"  Bins: {len(bins)}, sizes: {sorted([len(b) for b in bins], reverse=True)}")

        # Inner bc + k3 for each bin
        total_inner_bc = 0
        total_k3 = 0
        for bi, bin_set in enumerate(bins):
            inner_sub = g.subgraph(bin_set)
            inner_aps = set(nx.articulation_points(inner_sub))
            inner_remaining = bin_set - inner_aps
            inner_components = list(nx.connected_components(g.subgraph(inner_remaining)))

            k3_clusters = []
            for inner_bc in inner_components:
                k3_result = Physlr.detect_communities_k_clique(g, inner_bc, k=3)
                k3_clusters.extend(k3_result)

            total_inner_bc += len(inner_components)
            total_k3 += len(k3_clusters)

            if bi < 3:
                print(f"  Bin {bi} (size {len(bin_set)}): {len(inner_aps)} APs, "
                      f"{len(inner_components)} inner BCs, {len(k3_clusters)} k3 clusters")

        print(f"  Total inner BCs: {total_inner_bc}")
        print(f"  Total k3 clusters: {total_k3}")

    # Also output the neighbor subgraph as edge list for Rust comparison
    # Sort neighbors for deterministic output
    sorted_neighbors = sorted(neighbors)
    name_to_idx = {n: i for i, n in enumerate(sorted_neighbors)}

    print(f"\n# Edge list (local indices) for Rust comparison:")
    print(f"# n={len(sorted_neighbors)}")
    edge_count = 0
    for u_node in sorted_neighbors:
        for v_node in sub[u_node]:
            if name_to_idx[u_node] < name_to_idx[v_node]:
                edge_count += 1
    print(f"# edges={edge_count}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 67_compare_bc.py <overlap.filtered.tsv> [barcode1,barcode2,...]")
        sys.exit(1)

    tsv_file = sys.argv[1]
    barcodes = sys.argv[2].split(",") if len(sys.argv) > 2 else None

    class Args:
        pass
    args = Args()
    args.skip_small = 10
    Physlr.args = args

    print(f"Loading graph...", file=sys.stderr)
    t0 = timeit.default_timer()
    g = read_overlap_graph(tsv_file)
    print(f"Loaded in {timeit.default_timer() - t0:.1f}s: V={g.number_of_nodes()} E={g.number_of_edges()}", file=sys.stderr)

    if barcodes is None:
        degree_list = sorted(g.degree(), key=lambda x: x[1], reverse=True)
        barcodes = [n for n, d in degree_list[:5]]

    for barcode in barcodes:
        if barcode not in g:
            print(f"Barcode {barcode} not found")
            continue
        analyze_barcode(g, barcode)
        print("=" * 60)


if __name__ == "__main__":
    main()
