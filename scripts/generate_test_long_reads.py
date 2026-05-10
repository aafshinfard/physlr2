#!/usr/bin/env python3
"""Generate synthetic long-read FASTQ for testing physlr2 long-read mode.

Simulates ONT-like reads from a synthetic genome:
- Reads are 5-50 kbp (ONT ultra-long distribution)
- Homopolymer runs inserted to test compression
- ~10x coverage of a 500 kbp synthetic genome
- Overlapping reads to produce meaningful overlap graphs
"""

import random
import sys
import os

def generate_genome(size, seed=42):
    """Generate a random genome with realistic homopolymer content."""
    random.seed(seed)
    bases = "ACGT"
    genome = []
    i = 0
    while i < size:
        base = random.choice(bases)
        # ~5% chance of a homopolymer run (3-15 bases)
        if random.random() < 0.05:
            run_len = random.randint(3, 15)
            genome.extend([base] * min(run_len, size - i))
            i += min(run_len, size - i)
        else:
            genome.append(base)
            i += 1
    return "".join(genome)


def simulate_reads(genome, coverage=10, min_len=5000, max_len=50000, seed=42):
    """Simulate long reads from a genome at given coverage."""
    random.seed(seed)
    genome_len = len(genome)
    total_bases = genome_len * coverage
    reads = []
    bases_generated = 0
    read_id = 0

    while bases_generated < total_bases:
        read_len = random.randint(min_len, max_len)
        if read_len > genome_len:
            read_len = genome_len
        start = random.randint(0, genome_len - read_len)
        seq = genome[start:start + read_len]

        # Simulate ~5% error rate (substitutions only, for simplicity)
        seq_list = list(seq)
        for j in range(len(seq_list)):
            if random.random() < 0.05:
                seq_list[j] = random.choice("ACGT")
        seq = "".join(seq_list)

        # Random strand
        if random.random() < 0.5:
            seq = reverse_complement(seq)

        reads.append((f"read_{read_id:06d}", seq, start, read_len))
        bases_generated += read_len
        read_id += 1

    return reads


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def write_fastq(reads, output_path):
    """Write reads as FASTQ."""
    with open(output_path, "w") as f:
        for name, seq, start, length in reads:
            qual = "I" * len(seq)  # Phred 40
            f.write(f"@{name} start={start} length={length}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{qual}\n")


def main():
    output = sys.argv[1] if len(sys.argv) > 1 else "tests/data/test_long_reads.fq"
    genome_size = 500_000  # 500 kbp
    coverage = 10

    print(f"Generating {genome_size/1000:.0f} kbp genome...", file=sys.stderr)
    genome = generate_genome(genome_size)
    print(f"Genome: {len(genome)} bp", file=sys.stderr)

    print(f"Simulating {coverage}x long reads...", file=sys.stderr)
    reads = simulate_reads(genome, coverage=coverage)
    print(f"Generated {len(reads)} reads", file=sys.stderr)

    total_bases = sum(len(seq) for _, seq, _, _ in reads)
    print(f"Total bases: {total_bases/1e6:.1f} Mbp ({total_bases/genome_size:.1f}x)", file=sys.stderr)

    os.makedirs(os.path.dirname(output) if os.path.dirname(output) else ".", exist_ok=True)
    write_fastq(reads, output)
    print(f"Wrote {output}", file=sys.stderr)


if __name__ == "__main__":
    main()
