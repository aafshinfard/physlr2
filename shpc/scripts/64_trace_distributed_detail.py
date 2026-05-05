#!/usr/bin/env python3
"""
Trace the distributed pipeline step by step for a single barcode.
"""
import sys
import os
import random
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
                mxs = parts[1] if len(parts) > 1 else ""
                m = len(mxs.split()) if mxs else 0
                g.add_node(barcode, m=m)
            else:
                if len(parts) >= 3:
                    u, v = parts[0], parts[1]
                    try:
                        m = int(parts[2])
                    except ValueError:
                        continue
                    g.add_edge(u, v, m=m)
    return g


def trace_distributed(g, u):
    """Trace the distributed pipeline step by step."""
    neighbors = set(g[u].keys())
    print(f"Barcode: {u}, degree: {len(neighbors)}")

    # Step 1: bc
    bc_components = Physlr.detect_communities_biconnected_components(g, neighbors)
    print(f"\nStep 1 (bc): {len(bc_components)} biconnected components")
    bc_sizes = sorted([len(c) for c in bc_components], reverse=True)
    print(f"  Sizes: {bc_sizes[:20]}")
    print(f"  Total nodes: {sum(len(c) for c in bc_components)}")
    print(f"  Nodes lost (articulation points): {len(neighbors) - sum(len(c) for c in bc_components)}")

    total_clusters = 0
    total_merged = 0
    all_merged = []

    for bci, bc_comp in enumerate(bc_components):
        if bci >= 3 and len(bc_components) > 5:
            continue  # Only trace first 3

        # Step 2: bin
        random.seed(42)
        bins = Physlr.partition_subgraph_into_bins_randomly(bc_comp)
        if bci < 3:
            print(f"\n  BC component {bci} (size {len(bc_comp)}):")
            print(f"    Step 2 (bin): {len(bins)} bins, sizes: {sorted([len(b) for b in bins], reverse=True)[:10]}")

        # Step 3+4: bc + k3 within each bin
        clusters = []
        for bin_set in bins:
            inner_bcs = Physlr.detect_communities_biconnected_components(g, bin_set)
            for inner_bc in inner_bcs:
                k3_result = Physlr.detect_communities_k_clique(g, inner_bc, k=3)
                clusters.extend(k3_result)

        if bci < 3:
            print(f"    Step 3+4 (bc+k3): {len(clusters)} clusters, sizes: {sorted([len(c) for c in clusters], reverse=True)[:10]}")
            print(f"    Total nodes in clusters: {sum(len(c) for c in clusters)}")

        total_clusters += len(clusters)

        # Step 5: merge
        merged = Physlr.merge_communities(g, clusters, bc_comp, strategy=0)
        if bci < 3:
            print(f"    Step 5 (merge): {len(merged)} merged communities, sizes: {sorted([len(c) for c in merged], reverse=True)[:10]}")
            print(f"    Total nodes in merged: {sum(len(c) for c in merged)}")

        total_merged += len(merged)
        all_merged.extend(merged)

    # Process remaining bc components
    for bci in range(min(3, len(bc_components)), len(bc_components)):
        bc_comp = bc_components[bci]
        random.seed(42)
        bins = Physlr.partition_subgraph_into_bins_randomly(bc_comp)
        clusters = []
        for bin_set in bins:
            inner_bcs = Physlr.detect_communities_biconnected_components(g, bin_set)
            for inner_bc in inner_bcs:
                k3_result = Physlr.detect_communities_k_clique(g, inner_bc, k=3)
                clusters.extend(k3_result)
        total_clusters += len(clusters)
        merged = Physlr.merge_communities(g, clusters, bc_comp, strategy=0)
        total_merged += len(merged)
        all_merged.extend(merged)

    print(f"\n  Summary:")
    print(f"    Total clusters (pre-merge): {total_clusters}")
    print(f"    Total merged communities: {total_merged}")
    print(f"    Total nodes in all communities: {sum(len(c) for c in all_merged)}")
    print(f"    Communities with >1 member: {sum(1 for c in all_merged if len(c) > 1)}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 64_trace_distributed_detail.py <overlap.filtered.tsv> [barcode1,barcode2,...]")
        sys.exit(1)

    tsv_file = sys.argv[1]
    barcodes = sys.argv[2].split(",") if len(sys.argv) > 2 else None

    class Args:
        pass
    args = Args()
    args.strategy = "distributed"
    args.sqcost = 0.75
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
        trace_distributed(g, barcode)
        print("=" * 60)


if __name__ == "__main__":
    main()
