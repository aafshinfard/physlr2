#!/usr/bin/env python3
"""
Compare molecule separation between original Physlr and physlr-next.

Extracts a subgraph, runs original Python molecule separation,
and outputs detailed diagnostics for comparison with Rust.

Usage: python3 59_compare_molecules.py <overlap.filtered.tsv> [num_barcodes] [strategy]
"""
import sys
import os
import random
import timeit
import collections

# Add physlr to path
PHYSLR_DIR = "/home/afshinfa/miniconda3/envs/physlr-original/bin/share/physlr-1.0.4-8"
sys.path.insert(0, PHYSLR_DIR)

import networkx as nx
from physlr.physlr import Physlr


def read_overlap_graph(tsv_file):
    """Read overlap graph from TSV file."""
    g = nx.Graph()
    with open(tsv_file) as f:
        in_edges = False
        for line in f:
            line = line.strip()
            if not line:
                in_edges = True
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
                    m = int(parts[2])
                    g.add_edge(u, v, m=m)
    return g


def analyze_barcode(g, u, strategy):
    """Run molecule separation on a single barcode and return detailed info."""
    neighbors = list(g[u].keys())
    if not neighbors:
        return None

    # Run determine_molecules
    _, assignment = Physlr.determine_molecules(g, u, [], strategy)

    n_mols = len(set(assignment.values())) if assignment else 0
    n_assigned = len(assignment)
    n_unassigned = len(neighbors) - n_assigned

    # Get community sizes
    mol_sizes = collections.Counter(assignment.values())

    return {
        'barcode': u,
        'degree': len(neighbors),
        'n_molecules': n_mols,
        'n_assigned': n_assigned,
        'n_unassigned': n_unassigned,
        'mol_sizes': sorted(mol_sizes.values(), reverse=True) if mol_sizes else [],
    }


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 59_compare_molecules.py <overlap.filtered.tsv> [num_barcodes] [strategy]")
        sys.exit(1)

    tsv_file = sys.argv[1]
    num_barcodes = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    strategy = sys.argv[3] if len(sys.argv) > 3 else "distributed+sqcosbin"

    # Set up args
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

    # Select barcodes with varying degrees
    degree_list = sorted(g.degree(), key=lambda x: x[1], reverse=True)

    # Pick barcodes from different degree ranges
    test_barcodes = []
    # Top 100 highest degree
    test_barcodes.extend([n for n, d in degree_list[:100]])
    # Random sample from degree >= 20
    high_degree = [n for n, d in degree_list if d >= 20]
    random.seed(42)
    random.shuffle(high_degree)
    test_barcodes.extend(high_degree[:min(300, len(high_degree))])
    # Random sample from degree >= 5
    mid_degree = [n for n, d in degree_list if 5 <= d < 20]
    random.shuffle(mid_degree)
    test_barcodes.extend(mid_degree[:min(300, len(mid_degree))])
    # Random sample from all
    all_nodes = list(g.nodes())
    random.shuffle(all_nodes)
    test_barcodes.extend(all_nodes[:min(300, len(all_nodes))])

    # Deduplicate while preserving order
    seen = set()
    unique_barcodes = []
    for b in test_barcodes:
        if b not in seen:
            seen.add(b)
            unique_barcodes.append(b)
    test_barcodes = unique_barcodes[:num_barcodes]

    print(f"Testing {len(test_barcodes)} barcodes with strategy={strategy}", file=sys.stderr)

    # Run molecule separation
    total_mols = 0
    total_assigned = 0
    total_unassigned = 0
    max_mols = 0
    gt10 = 0
    gt20 = 0
    gt50 = 0
    degree_hist = collections.Counter()
    mol_hist = collections.Counter()

    # Output per-barcode results for comparison
    print("barcode\tdegree\tn_molecules\tn_assigned\tn_unassigned\tmol_sizes")
    t0 = timeit.default_timer()
    for i, u in enumerate(test_barcodes):
        if i % 100 == 0:
            print(f"  Processing {i}/{len(test_barcodes)}...", file=sys.stderr)
        result = analyze_barcode(g, u, strategy)
        if result is None:
            continue

        total_mols += result['n_molecules']
        total_assigned += result['n_assigned']
        total_unassigned += result['n_unassigned']
        max_mols = max(max_mols, result['n_molecules'])
        if result['n_molecules'] > 10:
            gt10 += 1
        if result['n_molecules'] > 20:
            gt20 += 1
        if result['n_molecules'] > 50:
            gt50 += 1
        degree_hist[result['degree']] += 1
        mol_hist[result['n_molecules']] += 1

        sizes_str = ",".join(str(s) for s in result['mol_sizes'])
        print(f"{result['barcode']}\t{result['degree']}\t{result['n_molecules']}\t{result['n_assigned']}\t{result['n_unassigned']}\t{sizes_str}")

    elapsed = timeit.default_timer() - t0
    print(f"\nDone in {elapsed:.1f}s", file=sys.stderr)

    print(f"\n=== Summary ===", file=sys.stderr)
    print(f"Strategy: {strategy}", file=sys.stderr)
    print(f"Barcodes tested: {len(test_barcodes)}", file=sys.stderr)
    print(f"Total molecules: {total_mols}", file=sys.stderr)
    print(f"Total assigned: {total_assigned}", file=sys.stderr)
    print(f"Total unassigned: {total_unassigned}", file=sys.stderr)
    print(f"Max molecules per barcode: {max_mols}", file=sys.stderr)
    print(f"Barcodes with >10 molecules: {gt10}", file=sys.stderr)
    print(f"Barcodes with >20 molecules: {gt20}", file=sys.stderr)
    print(f"Barcodes with >50 molecules: {gt50}", file=sys.stderr)
    print(f"Mean molecules/barcode: {total_mols/max(len(test_barcodes),1):.3f}", file=sys.stderr)

    print(f"\n=== Molecule count distribution ===", file=sys.stderr)
    for n_mols in sorted(mol_hist.keys()):
        print(f"  {n_mols} molecules: {mol_hist[n_mols]} barcodes", file=sys.stderr)


if __name__ == "__main__":
    main()
