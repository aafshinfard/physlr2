#!/usr/bin/env python3
"""
Trace molecule separation for specific barcodes in the original Physlr.
Outputs detailed community structure at each step.

Usage: python3 61_trace_single_barcode.py <overlap.filtered.tsv> [barcode1,barcode2,...]
"""
import sys
import os
import random
import timeit
import collections

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
                # Skip header line (U\tm or U\tV\tm)
                if line.startswith('U\t'):
                    continue
            parts = line.split('\t')
            if not in_edges:
                barcode = parts[0]
                mxs = parts[1] if len(parts) > 1 else ""
                m = len(mxs.split()) if mxs else 0
                g.add_node(barcode, m=m)
            else:
                if len(parts) >= 3:
                    u, v = parts[0], parts[1]
                    try:
                        m = int(parts[2])
                    except ValueError:
                        continue  # skip header
                    g.add_edge(u, v, m=m)
    return g


def trace_determine_molecules(g, u, strategy):
    """Run determine_molecules with detailed tracing."""
    neighbors = list(g[u].keys())
    print(f"\nBarcode: {u}")
    print(f"Degree: {len(neighbors)}")

    if not neighbors:
        print("No neighbors, skipping")
        return

    communities = [set(g[u].keys())]
    print(f"Initial: 1 community with {len(communities[0])} nodes")

    alg_list = strategy.split("+")
    for algorithm in alg_list:
        communities_temp = []
        print(f"\n--- Strategy: {algorithm} ---")
        print(f"Input: {len(communities)} communities, sizes: {sorted([len(c) for c in communities], reverse=True)[:20]}")

        if algorithm == "distributed":
            for ci, component in enumerate(communities):
                result = Physlr.determine_molecules_partition_split_merge(g, component)
                communities_temp.extend(result)
                if ci < 5:  # Only trace first 5
                    print(f"  Component {ci} (size {len(component)}) → {len(result)} communities, sizes: {sorted([len(c) for c in result], reverse=True)[:10]}")
        elif algorithm == "sqcosbin":
            for ci, component in enumerate(communities):
                bins = Physlr.partition_subgraph_into_bins_randomly(component)
                clusters = []
                for bin_set in bins:
                    clusters.extend(
                        Physlr.detect_communities_cosine_of_squared(
                            g, bin_set, squaring=True, threshold=Physlr.args.sqcost))
                merged = Physlr.merge_communities(g, clusters)
                communities_temp.extend(merged)
                if ci < 5:
                    print(f"  Component {ci} (size {len(component)}) → {len(bins)} bins → {len(clusters)} sqcos clusters → {len(merged)} merged, sizes: {sorted([len(c) for c in merged], reverse=True)[:10]}")
        elif algorithm == "bc":
            for component in communities:
                communities_temp.extend(
                    Physlr.detect_communities_biconnected_components(g, component))
        elif algorithm == "cc":
            for component in communities:
                communities_temp.extend(
                    Physlr.detect_communities_connected_components(g, component))
        elif algorithm == "k3":
            for component in communities:
                communities_temp.extend(
                    Physlr.detect_communities_k_clique(g, component, k=3))

        communities = communities_temp
        sizes = sorted([len(c) for c in communities], reverse=True)
        print(f"Output: {len(communities)} communities, sizes: {sizes[:20]}")
        print(f"  Total nodes in communities: {sum(len(c) for c in communities)}")
        print(f"  Communities with >1 member: {sum(1 for c in communities if len(c) > 1)}")

    # Final assignment
    communities_sorted = sorted(communities, key=len, reverse=True)
    assignment = {v: i for i, vs in enumerate(communities_sorted) if len(vs) > 1 for v in vs}
    n_mols = 1 + max(assignment.values()) if assignment else 0
    n_assigned = len(assignment)
    n_unassigned = len(neighbors) - n_assigned

    print(f"\n--- Final ---")
    print(f"n_molecules: {n_mols}")
    print(f"n_assigned: {n_assigned}")
    print(f"n_unassigned: {n_unassigned}")
    print(f"Distinct molecule IDs: {len(set(assignment.values()))}")

    return n_mols


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 61_trace_single_barcode.py <overlap.filtered.tsv> [barcode1,barcode2,...]")
        sys.exit(1)

    tsv_file = sys.argv[1]
    barcodes = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    strategy = sys.argv[3] if len(sys.argv) > 3 else "distributed+sqcosbin"

    class Args:
        pass
    args = Args()
    args.strategy = strategy
    args.sqcost = 0.75
    args.skip_small = 10
    Physlr.args = args

    print(f"Loading graph from {tsv_file}...", file=sys.stderr)
    t0 = timeit.default_timer()
    g = read_overlap_graph(tsv_file)
    print(f"Loaded in {timeit.default_timer() - t0:.1f}s: V={g.number_of_nodes()} E={g.number_of_edges()}", file=sys.stderr)

    if barcodes is None:
        # Pick barcodes with highest degree
        degree_list = sorted(g.degree(), key=lambda x: x[1], reverse=True)
        barcodes = [n for n, d in degree_list[:5]]
        # Also pick some with moderate degree
        mid = [n for n, d in degree_list if 20 <= d <= 30]
        random.seed(42)
        random.shuffle(mid)
        barcodes.extend(mid[:5])

    print(f"\nTracing {len(barcodes)} barcodes with strategy={strategy}")
    print(f"{'='*60}")

    total_mols = 0
    for barcode in barcodes:
        if barcode not in g:
            print(f"\nBarcode {barcode} not found in graph")
            continue
        random.seed(42)  # Reset seed for reproducibility
        n_mols = trace_determine_molecules(g, barcode, strategy)
        if n_mols:
            total_mols += n_mols
        print(f"{'='*60}")

    print(f"\nTotal molecules across {len(barcodes)} barcodes: {total_mols}")


if __name__ == "__main__":
    main()
