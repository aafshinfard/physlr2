#!/usr/bin/env python3
"""
Test the impact of random vs deterministic binning on molecule separation.

Runs the original Physlr's molecule separation on a subset of the overlap graph,
comparing random.shuffle binning (original) vs hash-based deterministic binning.

Usage: python3 58_test_binning_impact.py <overlap.filtered.tsv> [num_barcodes]
"""
import sys
import os
import random
import timeit

# Add physlr to path
PHYSLR_DIR = "/home/afshinfa/miniconda3/envs/physlr-original/bin/share/physlr-1.0.4-8"
sys.path.insert(0, PHYSLR_DIR)

import networkx as nx
from physlr.physlr import Physlr

def deterministic_partition(node_set, max_size=50):
    """Deterministic binning using hash-based sort (matching Rust implementation)."""
    bins_count = 1 + len(node_set) // max_size
    node_list = list(node_set)
    # Sort by hash (simulating Rust's wrapping_mul approach)
    # Use Python's hash for determinism
    node_list.sort(key=lambda x: hash(x) & 0xFFFFFFFFFFFFFFFF)
    size, leftover = divmod(len(node_set), bins_count)
    bins = [node_list[0 + size * i: size * (i + 1)] for i in range(bins_count)]
    edge = size * bins_count
    for i in range(leftover):
        bins[i % bins_count].append(node_list[edge + i])
    return [set(x) for x in bins]


def run_molecules_with_strategy(g, nodes, strategy_name, use_deterministic_binning=False):
    """Run molecule separation on a set of nodes with the given strategy."""
    # Monkey-patch the binning function if needed
    original_partition = Physlr.partition_subgraph_into_bins_randomly
    if use_deterministic_binning:
        Physlr.partition_subgraph_into_bins_randomly = staticmethod(deterministic_partition)

    total_molecules = 0
    total_singletons = 0
    molecule_counts = []

    for u in nodes:
        neighbors = list(g[u].keys())
        if not neighbors:
            continue
        _, assignment = Physlr.determine_molecules(g, u, None, strategy_name)
        n_mols = len(set(assignment.values())) if assignment else 1
        unassigned = len(neighbors) - len(assignment)
        molecule_counts.append(n_mols)
        total_molecules += n_mols
        total_singletons += unassigned

    # Restore original
    Physlr.partition_subgraph_into_bins_randomly = original_partition

    return total_molecules, total_singletons, molecule_counts


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 58_test_binning_impact.py <overlap.filtered.tsv> [num_barcodes]")
        sys.exit(1)

    tsv_file = sys.argv[1]
    num_barcodes = int(sys.argv[2]) if len(sys.argv) > 2 else 10000

    # Set up args
    class Args:
        strategy = "distributed+sqcosbin"
        sqcost = 0.75
        skip_small = 10
    Physlr.args = Args()

    print(f"Loading graph from {tsv_file}...", file=sys.stderr)
    t0 = timeit.default_timer()

    # Read the overlap graph
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
                # Vertex line: barcode\tmin1 min2 ...
                barcode = parts[0]
                g.add_node(barcode)
            else:
                # Edge line: u\tv\tm
                if len(parts) >= 3:
                    u, v = parts[0], parts[1]
                    m = int(parts[2])
                    g.add_edge(u, v, m=m)

    print(f"Loaded in {timeit.default_timer() - t0:.1f}s: V={g.number_of_nodes()} E={g.number_of_edges()}", file=sys.stderr)

    # Select a subset of barcodes with enough neighbors
    candidates = [n for n in g.nodes() if g.degree(n) >= 5]
    random.seed(42)
    random.shuffle(candidates)
    test_nodes = candidates[:num_barcodes]

    print(f"Testing on {len(test_nodes)} barcodes with degree >= 5", file=sys.stderr)

    # Test with random binning (original)
    print("\n=== Random binning (original) ===")
    t0 = timeit.default_timer()
    random.seed(42)  # Set seed for reproducibility
    mols_random, singles_random, counts_random = run_molecules_with_strategy(
        g, test_nodes, "distributed+sqcosbin", use_deterministic_binning=False)
    elapsed = timeit.default_timer() - t0
    print(f"Total molecules: {mols_random}")
    print(f"Total singletons: {singles_random}")
    print(f"Mean molecules/barcode: {mols_random/len(test_nodes):.3f}")
    print(f"Time: {elapsed:.1f}s")

    # Test with deterministic binning (Rust-style)
    print("\n=== Deterministic binning (Rust-style) ===")
    t0 = timeit.default_timer()
    mols_det, singles_det, counts_det = run_molecules_with_strategy(
        g, test_nodes, "distributed+sqcosbin", use_deterministic_binning=True)
    elapsed = timeit.default_timer() - t0
    print(f"Total molecules: {mols_det}")
    print(f"Total singletons: {singles_det}")
    print(f"Mean molecules/barcode: {mols_det/len(test_nodes):.3f}")
    print(f"Time: {elapsed:.1f}s")

    # Compare
    print(f"\n=== Comparison ===")
    print(f"Random: {mols_random} molecules, {singles_random} singletons")
    print(f"Deterministic: {mols_det} molecules, {singles_det} singletons")
    print(f"Ratio (det/random): {mols_det/max(mols_random,1):.3f}")

    # Also test with random binning but different seed
    print("\n=== Random binning (different seed) ===")
    t0 = timeit.default_timer()
    random.seed(123)
    mols_random2, singles_random2, counts_random2 = run_molecules_with_strategy(
        g, test_nodes, "distributed+sqcosbin", use_deterministic_binning=False)
    elapsed = timeit.default_timer() - t0
    print(f"Total molecules: {mols_random2}")
    print(f"Total singletons: {singles_random2}")
    print(f"Mean molecules/barcode: {mols_random2/len(test_nodes):.3f}")
    print(f"Ratio (seed123/seed42): {mols_random2/max(mols_random,1):.3f}")

    # Test distributed only (no sqcosbin)
    print("\n=== Distributed only (random binning) ===")
    t0 = timeit.default_timer()
    random.seed(42)
    mols_dist, singles_dist, counts_dist = run_molecules_with_strategy(
        g, test_nodes, "distributed", use_deterministic_binning=False)
    elapsed = timeit.default_timer() - t0
    print(f"Total molecules: {mols_dist}")
    print(f"Total singletons: {singles_dist}")
    print(f"Mean molecules/barcode: {mols_dist/len(test_nodes):.3f}")
    print(f"Time: {elapsed:.1f}s")

    print("\n=== Distributed only (deterministic binning) ===")
    t0 = timeit.default_timer()
    mols_dist_det, singles_dist_det, counts_dist_det = run_molecules_with_strategy(
        g, test_nodes, "distributed", use_deterministic_binning=True)
    elapsed = timeit.default_timer() - t0
    print(f"Total molecules: {mols_dist_det}")
    print(f"Total singletons: {singles_dist_det}")
    print(f"Mean molecules/barcode: {mols_dist_det/len(test_nodes):.3f}")
    print(f"Ratio (det/random): {mols_dist_det/max(mols_dist,1):.3f}")


if __name__ == "__main__":
    main()
