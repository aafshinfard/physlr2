#!/usr/bin/env python3
"""
Simulate linked reads from a reference genome.

Generates FASTQ reads with BX:Z: barcode tags that mimic 10x Chromium
linked-read sequencing. Each barcode (GEM) contains multiple molecules
from different genomic regions, matching real 10x data characteristics.
"""

import argparse
import random
import sys

BASES = "ACGT"
BARCODE_CHARS = "ACGT"


def read_fasta(path):
    """Read a FASTA file into a list of (name, sequence) tuples."""
    seqs = []
    name = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(seq_parts).upper()))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        seqs.append((name, "".join(seq_parts).upper()))
    return seqs


def random_barcode():
    """Generate a random 16-mer barcode like 10x Chromium."""
    return "".join(random.choice(BARCODE_CHARS) for _ in range(16)) + "-1"


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def simulate(
    ref_seqs,
    n_barcodes,
    molecules_per_barcode_mean,
    molecule_size_mean,
    molecule_size_std,
    reads_per_molecule,
    read_length,
    error_rate,
    output_path,
):
    """Simulate linked reads with multiple molecules per barcode."""
    read_id = 0
    total_molecules = 0

    with open(output_path, "w") as out:
        for bc_i in range(n_barcodes):
            barcode = random_barcode()

            # Each barcode (GEM) contains multiple molecules
            n_mols = max(1, int(random.expovariate(1.0 / molecules_per_barcode_mean)))
            n_mols = min(n_mols, 30)  # Cap at 30 molecules per GEM

            for _ in range(n_mols):
                total_molecules += 1

                # Pick a random chromosome weighted by length
                weights = [len(s) for _, s in ref_seqs]
                total = sum(weights)
                r = random.random() * total
                cumulative = 0
                chrom_idx = 0
                for i, w in enumerate(weights):
                    cumulative += w
                    if r <= cumulative:
                        chrom_idx = i
                        break

                chrom_name, chrom_seq = ref_seqs[chrom_idx]
                chrom_len = len(chrom_seq)

                # Molecule size
                mol_size = max(
                    1000, int(random.gauss(molecule_size_mean, molecule_size_std))
                )
                mol_size = min(mol_size, chrom_len)

                # Molecule position
                mol_start = random.randint(0, max(0, chrom_len - mol_size))
                mol_end = min(mol_start + mol_size, chrom_len)

                for _ in range(reads_per_molecule):
                    read_id += 1

                    # Read position within molecule
                    rstart = random.randint(
                        mol_start, max(mol_start, mol_end - read_length)
                    )
                    rend = min(rstart + read_length, mol_end)
                    seq = chrom_seq[rstart:rend]

                    if len(seq) < 50:
                        continue

                    # Add errors
                    if error_rate > 0:
                        seq_list = list(seq)
                        for i in range(len(seq_list)):
                            if random.random() < error_rate:
                                seq_list[i] = random.choice(
                                    [b for b in BASES if b != seq_list[i]]
                                )
                        seq = "".join(seq_list)

                    # Random orientation
                    if random.random() < 0.5:
                        seq = reverse_complement(seq)

                    qual = "I" * len(seq)
                    out.write(f"@read{read_id} BX:Z:{barcode}\n")
                    out.write(seq + "\n")
                    out.write("+\n")
                    out.write(qual + "\n")

            if (bc_i + 1) % 500 == 0:
                print(
                    f"  Generated {bc_i + 1}/{n_barcodes} barcodes, "
                    f"{total_molecules} molecules, {read_id} reads",
                    file=sys.stderr,
                )

    return read_id, total_molecules


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", help="Reference FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTQ file")
    parser.add_argument(
        "-n",
        "--n-barcodes",
        type=int,
        default=500,
        help="Number of barcodes (GEMs) to simulate",
    )
    parser.add_argument(
        "--mols-per-barcode",
        type=int,
        default=5,
        help="Mean molecules per barcode (GEM)",
    )
    parser.add_argument(
        "--molecule-size",
        type=int,
        default=50000,
        help="Mean molecule size (bp)",
    )
    parser.add_argument(
        "--molecule-std",
        type=int,
        default=20000,
        help="Molecule size std dev (bp)",
    )
    parser.add_argument(
        "--reads-per-mol",
        type=int,
        default=15,
        help="Reads per molecule",
    )
    parser.add_argument(
        "--read-length", type=int, default=150, help="Read length (bp)"
    )
    parser.add_argument(
        "--error-rate",
        type=float,
        default=0.001,
        help="Per-base error rate",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    random.seed(args.seed)

    print(f"Reading reference from {args.reference}...", file=sys.stderr)
    ref_seqs = read_fasta(args.reference)
    total_bp = sum(len(s) for _, s in ref_seqs)
    print(
        f"  {len(ref_seqs)} sequences, {total_bp} bp total",
        file=sys.stderr,
    )

    print(
        f"Simulating {args.n_barcodes} barcodes with ~{args.mols_per_barcode} molecules each...",
        file=sys.stderr,
    )
    n_reads, n_molecules = simulate(
        ref_seqs,
        args.n_barcodes,
        args.mols_per_barcode,
        args.molecule_size,
        args.molecule_std,
        args.reads_per_mol,
        args.read_length,
        args.error_rate,
        args.output,
    )

    coverage = (n_reads * args.read_length) / total_bp
    print(
        f"Generated {n_reads} reads from {n_molecules} molecules "
        f"in {args.n_barcodes} barcodes ({coverage:.1f}x coverage)",
        file=sys.stderr,
    )
    print(f"Output: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
