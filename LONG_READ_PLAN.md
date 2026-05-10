# Long-Read Physical Map Support — Development Plan

## Overview

Extend physlr2 to construct physical maps from long-read sequencing data (ONT, PacBio).
Each long read is treated as its own "barcode" (read name = barcode ID), and the existing
pipeline (index → filter → overlap → molecules → backbone) is reused with long-read-specific
adaptations.

Based on: Afshinfard, A. (2024). *Data Mining Techniques for De Novo Genome Assembly and Analysis*. PhD thesis, UBC.

## Protocol Detection

A new `--protocol` flag with three modes:

| Value    | Behavior |
|----------|----------|
| `auto`   | Inspect first N reads: if BX:Z: or stLFR barcodes found → linked; otherwise → long |
| `linked` | Force linked-read mode (current defaults) |
| `long`   | Force long-read mode (adjusted defaults below) |

Auto-detection heuristic: sample the first 1000 reads. If ≥10% have a BX:Z: tag or
stLFR `#barcode` pattern, classify as linked reads. Otherwise, classify as long reads.

## Long-Read Adaptations

### 1. Read Name as Barcode
Already partially supported: when `extract_barcode()` returns `None`, the read name is
used. For long-read mode, skip barcode extraction entirely and always use the read name.

### 2. Homopolymer Compression (Reversible)
Collapse runs of identical bases before minimizer extraction to reduce indel noise in
homopolymer regions. The compression must be reversible — maintain a coordinate mapping
from compressed positions back to original positions so downstream tools can report
original coordinates.

Implementation: `homopolymer_compress(seq) -> (compressed_seq, coord_map)` where
`coord_map[compressed_pos] = original_pos`.

### 3. Ordered Minimizers
Long reads produce ordered minimizers (position matters), unlike linked reads where
minimizers from multiple reads under one barcode form an unordered set. The indexing
step stores `Vec<u64>` (ordered) instead of `FxHashSet<u64>` (unordered).

### 4. Parameter Defaults for Long Reads
- `min_count` (min minimizers per read): 50 (vs 100 for linked)
- `max_count` (max minimizers per read): 50000 (vs 5000 for linked)
- `min_overlap` (min shared minimizers for edge): 5 (vs 10 for linked)

These are starting points; will be tuned with CHM13 data.

### 5. Chimeric Read Detection (via Molecule Separation)
A chimeric long read spans two unrelated genomic regions. Its minimizers will cluster
into two groups with no overlap. The existing molecule separation machinery (bc+cc strategy)
naturally handles this: each read name is a "barcode", and if its neighborhood in the
overlap graph has two disconnected components, molecule separation splits it into two
sub-reads.

### 6. Overlap Correction via RANSAC (Optional)
For each edge in the overlap graph, fit a linear regression to the positions of shared
minimizers in both reads. True overlaps have slope ≈ ±1. Reject edges where the 95% CI
of the slope doesn't cover ±1 or the margin of error exceeds 0.1.

This step is **optional** (off by default) because it can be slow on large datasets.
Enable with `--refine-overlap`.

### 7. Double Bloom Filter
Two-pass repeat masking: first pass counts k-mer frequencies, second pass builds a
Bloom filter of repeats. Reduces false positives compared to single-pass. Already
partially implemented in the `repeat` module.

## Test Data

CHM13 human cell line, ONT ultra-long reads (~390 Gbp, ~126x coverage):

- **rel8** (Guppy 5.0.7): `s3://human-pangenomics/T2T/CHM13/nanopore/rel8-guppy-5.0.7/reads.fastq.gz`
- **Guppy 6.3.7 HAC**: NCBI SRA accession `SRX19306105`

For development/testing: extract reads mapping to chr1 (~8% of genome, ~10x subset).

## Expected Results (from thesis, rel3 data)

| Parameter | NG50 (Mb) | Misassemblies |
|-----------|-----------|---------------|
| m=70      | ~83       | 14            |
| m=85      | ~40.7     | 2             |

## Implementation Phases

1. **Protocol detection + CLI** — `--protocol auto/linked/long`, parameter defaults
2. **Indexing changes** — homopolymer compression, ordered minimizers, read-name-as-barcode
3. **Pipeline integration** — wire long-read mode through `physical-map` and `pipeline` commands
4. **Testing on CHM13 chr1** — validate first few steps work end-to-end
5. **Overlap refinement** — optional RANSAC step (later phase)
