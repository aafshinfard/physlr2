#!/usr/bin/env python3
"""
Generate a synthetic linked-read dataset for testing the Physlr pipeline.

Creates:
  - A synthetic reference genome (multiple chromosomes)
  - Synthetic linked reads with barcodes (FASTQ with BX:Z: tags)
  - Two fragmented draft assemblies for scaffolding

The synthetic genome has known structure so we can verify the pipeline output.
"""

import argparse
import os
import random
import sys

BASES = "ACGT"


def random_seq(length: int) -> str:
    return "".join(random.choice(BASES) for _ in range(length))


def reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def generate_reference(n_chroms: int, chrom_sizes: list[int], outpath: str):
    """Generate a synthetic reference genome."""
    chroms = {}
    with open(outpath, "w") as f:
        for i in range(n_chroms):
            name = f"chr{i+1}"
            seq = random_seq(chrom_sizes[i % len(chrom_sizes)])
            chroms[name] = seq
            f.write(f">{name} LN:{len(seq)}\n")
            # Write in 80-char lines
            for j in range(0, len(seq), 80):
                f.write(seq[j : j + 80] + "\n")
    return chroms


def generate_linked_reads(
    chroms: dict[str, str],
    n_molecules_per_chrom: int,
    molecule_size: int,
    reads_per_molecule: int,
    read_length: int,
    outpath: str,
):
    """Generate synthetic linked reads with barcodes."""
    barcode_id = 0
    read_id = 0

    with open(outpath, "w") as f:
        for chrom_name, chrom_seq in chroms.items():
            chrom_len = len(chrom_seq)
            if chrom_len < molecule_size:
                molecule_size_actual = chrom_len // 2
            else:
                molecule_size_actual = molecule_size

            for mol_i in range(n_molecules_per_chrom):
                barcode_id += 1
                barcode = f"BARCODE{barcode_id:06d}-1"

                # Random molecule position
                start = random.randint(0, max(0, chrom_len - molecule_size_actual))
                end = min(start + molecule_size_actual, chrom_len)

                for _ in range(reads_per_molecule):
                    read_id += 1
                    # Random read position within molecule
                    rstart = random.randint(start, max(start, end - read_length))
                    rend = min(rstart + read_length, end)
                    seq = chrom_seq[rstart:rend]

                    # Randomly reverse complement
                    if random.random() < 0.5:
                        seq = reverse_complement(seq)

                    qual = "I" * len(seq)
                    f.write(f"@read{read_id} BX:Z:{barcode}\n")
                    f.write(seq + "\n")
                    f.write("+\n")
                    f.write(qual + "\n")

    return barcode_id, read_id


def generate_draft_assembly(
    chroms: dict[str, str], chunk_size: int, label: str, outpath: str
):
    """Generate a fragmented draft assembly by splitting chromosomes."""
    contig_id = 0
    with open(outpath, "w") as f:
        for chrom_name, chrom_seq in chroms.items():
            for i in range(0, len(chrom_seq), chunk_size):
                contig_id += 1
                seq = chrom_seq[i : i + chunk_size]
                if len(seq) < 100:
                    continue
                f.write(f">{label}_contig{contig_id}\n")
                for j in range(0, len(seq), 80):
                    f.write(seq[j : j + 80] + "\n")
    return contig_id


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o", "--outdir", default="data/synthetic", help="Output directory"
    )
    parser.add_argument(
        "--n-chroms", type=int, default=3, help="Number of chromosomes"
    )
    parser.add_argument(
        "--chrom-size",
        type=int,
        default=100_000,
        help="Size of each chromosome (bp)",
    )
    parser.add_argument(
        "--n-molecules",
        type=int,
        default=50,
        help="Molecules per chromosome",
    )
    parser.add_argument(
        "--molecule-size",
        type=int,
        default=50_000,
        help="Average molecule size (bp)",
    )
    parser.add_argument(
        "--reads-per-mol",
        type=int,
        default=20,
        help="Reads per molecule",
    )
    parser.add_argument(
        "--read-length", type=int, default=150, help="Read length (bp)"
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    random.seed(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    chrom_sizes = [args.chrom_size] * args.n_chroms

    # Generate reference
    ref_path = os.path.join(args.outdir, "reference.fa")
    print(f"Generating reference genome ({args.n_chroms} chroms, {args.chrom_size} bp each)...")
    chroms = generate_reference(args.n_chroms, chrom_sizes, ref_path)
    print(f"  Written to {ref_path}")

    # Generate linked reads
    reads_path = os.path.join(args.outdir, "linked_reads.fq")
    print(f"Generating linked reads...")
    n_barcodes, n_reads = generate_linked_reads(
        chroms,
        args.n_molecules,
        args.molecule_size,
        args.reads_per_mol,
        args.read_length,
        reads_path,
    )
    print(f"  {n_barcodes} barcodes, {n_reads} reads -> {reads_path}")

    # Generate draft assemblies
    draft1_path = os.path.join(args.outdir, "draft1.fa")
    print(f"Generating draft assembly 1 (10kb fragments)...")
    n1 = generate_draft_assembly(chroms, 10_000, "d1", draft1_path)
    print(f"  {n1} contigs -> {draft1_path}")

    draft2_path = os.path.join(args.outdir, "draft2.fa")
    print(f"Generating draft assembly 2 (5kb fragments)...")
    n2 = generate_draft_assembly(chroms, 5_000, "d2", draft2_path)
    print(f"  {n2} contigs -> {draft2_path}")

    # Generate config
    genome_size = sum(len(s) for s in chroms.values())
    config_path = os.path.join(args.outdir, "config.yaml")
    with open(config_path, "w") as f:
        f.write(f"""# Physlr test configuration (synthetic data)
physlr_bin: "target/release/physlr"
outdir: "{args.outdir}/output"
prefix: "synthetic"

linked_reads: "{reads_path}"
reference: "{ref_path}"
genome_size: {genome_size}

draft_assemblies:
  draft1: "{draft1_path}"
  draft2: "{draft2_path}"

# Parameters (tuned for small synthetic data)
k: 20
w: 10
min_bx_count: 2
max_bx_count: 5000
min_overlap: 2
edge_percentile: 0
prune_branches: 2
min_path_size: 3
min_map_score: 2
gap_size: 100
""")
    print(f"\nConfig written to {config_path}")
    print(f"Genome size: {genome_size} bp")
    print(f"\nTo run the pipeline:")
    print(f"  physlr physical-map {reads_path} -o {args.outdir}/output -p synthetic \\")
    print(f"    --min-bx-count 2 --min-overlap 2 --edge-percentile 0 \\")
    print(f"    --prune-branches 2 --min-path-size 3")


if __name__ == "__main__":
    main()
