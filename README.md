[![Published in DNA](https://img.shields.io/badge/Published%20in-DNA-blue.svg)](https://doi.org/10.3390/dna2020009)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<p align="center">
  <img src="https://raw.githubusercontent.com/bcgsc/physlr/master/physlr-logo-transparent.png" alt="Physlr logo" width="400">
</p>

# Physlr 2

**A ground-up Rust rewrite of [Physlr](https://github.com/bcgsc/physlr) for constructing physical maps from linked reads.**

Physlr 2 takes linked-read sequencing data (10x Genomics Chromium or MGI stLFR) and constructs *de novo* physical maps — ordered sets of molecules along each chromosome. These physical maps can then scaffold draft genome assemblies to chromosome-level contiguity.

## Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline](#pipeline)
- [Commands](#commands)
- [Parameters](#parameters)
- [Reproducing Results](#reproducing-results)
- [Comparison with Physlr v1](#comparison-with-physlr-v1)
- [Project Structure](#project-structure)
- [Citation](#citation)
- [License](#license)

## Overview

<p align="center">
  <img src="https://raw.githubusercontent.com/bcgsc/physlr/master/physlr-stages.png" alt="Physlr pipeline stages" height="200">
</p>

Physlr 2 performs two main tasks:

1. **Physical map construction** — builds an ordered map of molecules along each chromosome from linked-read barcodes.
2. **Assembly scaffolding** — uses the physical map to order and orient contigs from a draft genome assembly.

An optional **merge-paths** post-processing step can further improve physical map contiguity by identifying non-backbone "bridge" molecules that share minimizers with the endpoints of adjacent backbone paths, providing evidence to merge them. On human stLFR data (NA12878 + NA24143), merge-paths adds 9 true-positive merges with zero false positives using the default parameters.

### Physical Map Visualization

Backbone paths mapped to the GRCh38 reference genome for two human cell lines:

| | NA12878 (stLFR) | NA24143 (stLFR) |
|:---|:---:|:---:|
| **Backbone** | [<img src="results/na12878_backbone_v023.png" height="200">](results/na12878_backbone_v023.png) | [<img src="results/na24143_backbone_v023.png" height="200">](results/na24143_backbone_v023.png) |
| **Reference** | [<img src="results/na12878_reference_v023.png" height="200">](results/na12878_reference_v023.png) | [<img src="results/na24143_reference_v023.png" height="200">](results/na24143_reference_v023.png) |

*Backbone: paths colored by reference chromosome. Reference: chromosomes colored by backbone path. Click to enlarge.*

## Installation

### Install Rust

If you do not have Rust installed, install it via [rustup](https://rustup.rs/):

```bash
curl --proto "=https" --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

### Build from Source

```bash
git clone https://github.com/aafshinfard/physlr2.git
cd physlr2
cargo build --release
```

The compiled binary will be at `target/release/physlr`. You can either:

**Option A — Add to PATH (recommended):**

```bash
export PATH="$(pwd)/target/release:$PATH"
# Add the line above to ~/.bashrc or ~/.zshrc to make it permanent
```

**Option B — Install to ~/.cargo/bin:**

```bash
cargo install --path .
```

**Option C — Copy to a system directory:**

```bash
sudo cp target/release/physlr /usr/local/bin/
```

### Verify Installation

```bash
physlr --version
physlr --help
```

### Dependencies

**Build:**
- [Rust](https://www.rust-lang.org/tools/install) ≥ 1.70 (includes cargo)

**Runtime (optional):**
- [btllib](https://github.com/bcgsc/btllib) — for `indexlr` minimizer extraction backend (`conda install -c bioconda btllib`)
- Python ≥ 3.8 with matplotlib — for backbone-vs-reference visualization (`scripts/plotpaf.py`)
- [Snakemake](https://snakemake.readthedocs.io/) ≥ 7.0 — for automated workflow (`workflow/Snakefile`)
- [QUAST](https://github.com/ablab/quast) — for reference-based assembly evaluation

## Quick Start

### One-command physical map

Build a physical map directly from linked-read FASTQ files:

```bash
physlr physical-map reads.R1.fq.gz reads.R2.fq.gz -o output/ -p mygenome
```

This runs the full pipeline (indexing, filtering, overlap, molecule separation, backbone extraction) and produces the physical map in `output/`.

To also run the optional **merge-paths** step, which uses bridge molecules to merge adjacent backbone paths and improve contiguity:

```bash
physlr physical-map reads.R1.fq.gz reads.R2.fq.gz -o output/ -p mygenome --merge-paths
```

### Scaffolding

Scaffold a draft assembly using the physical map output:

```bash
physlr scaffolds output/mygenome.backbone.path output/mygenome.filtered.tsv draft.fa -o output/ -p mygenome
```

The scaffolds command takes three positional arguments:
1. The backbone path file produced by `physical-map`
2. The filtered minimizer TSV produced by `physical-map`
3. The draft assembly FASTA to scaffold

To include NG50 in the output metrics, add the expected genome size (optional):

```bash
physlr scaffolds output/mygenome.backbone.path output/mygenome.filtered.tsv draft.fa \
  -o output/ -p mygenome -g 3088269832
```

### Step-by-step CLI

```bash
# 1. Index minimizers from linked reads
physlr index reads.fq.gz -o reads.mxs.tsv -k 32 -w 32

# 2. Filter barcodes and minimizers
physlr filter-minimizers reads.mxs.tsv -o filtered.mxs.tsv -n 100 -N 5000

# 3. Compute barcode overlap graph
physlr overlap filtered.mxs.tsv -o overlap.tsv

# 4. Filter edges by percentile
physlr filter-overlap overlap.tsv -o overlap.filtered.tsv -p 85

# 5. Separate barcodes into molecules
physlr molecules overlap.filtered.tsv -o mol.tsv --strategy bc+cc

# 6. Extract backbone paths (physical map)
physlr backbone mol.tsv -o backbone.path

# 7. (Optional) Merge adjacent backbone paths
physlr split-minimizers mol.tsv filtered.mxs.tsv -o split.mxs.tsv
physlr merge-paths backbone.path split.mxs.tsv -o merged.path

# 8. (Optional) Scaffold a draft assembly using the physical map
physlr index-contigs draft.fa -o draft.mxs.tsv
physlr map backbone.path filtered.mxs.tsv draft.mxs.tsv -o map.bed
physlr bed-to-path map.bed -o scaffold.path
physlr path-to-fasta draft.fa scaffold.path -o scaffolds.fa
```

### Snakemake Workflow

```bash
cd workflow/
# Edit config.yaml with your file paths
snakemake -s Snakefile --configfile config.yaml -j 8
```

## Pipeline

```
Linked reads (FASTQ + barcodes)
  │
  ├── index              Extract (k,w)-minimizers per barcode
  ├── filter-minimizers  Remove low/high-count barcodes, singleton minimizers
  ├── overlap            Compute barcode overlap graph (shared minimizers)
  ├── filter-overlap     Remove low-weight edges by percentile
  ├── molecules          Separate barcodes into individual molecules
  ├── backbone           MST → prune branches → extract paths
  │
  └──► Physical map (backbone paths)
         │
         ├── merge-paths       (Optional) Merge adjacent paths via bridge molecules
         ├── map / map-paf     Map contigs or reference to the physical map
         ├── bed-to-path       Convert mappings to scaffold paths
         ├── path-to-fasta     Produce scaffolded FASTA
         │
         └──► Scaffolded assembly
```

## Commands

| Command | Description |
|---------|-------------|
| **Indexing** | |
| `index` | Extract (k,w)-minimizers from FASTA/FASTQ, grouped by barcode |
| `index-contigs` | Extract ordered minimizers from FASTA contigs or reference |
| `repeat-filter` | Detect repetitive k-mers and build a Bloom filter |
| **Graph Construction** | |
| `filter-minimizers` | Filter barcodes by count; remove singleton/repetitive minimizers |
| `overlap` | Compute barcode overlap graph from shared minimizers |
| `filter-overlap` | Remove low-weight edges by percentile |
| **Molecule Separation** | |
| `molecules` | Separate barcodes into individual molecules |
| `split-minimizers` | Assign barcode minimizers to individual molecules |
| `trace-molecules` | Diagnostic: trace molecule separation for specific barcodes |
| **Physical Map** | |
| `backbone` | Extract backbone paths from the molecule overlap graph |
| `merge-paths` | Merge adjacent backbone paths using bridge molecule evidence |
| **Scaffolding** | |
| `map` | Map query sequences to the physical map (BED output) |
| `map-paf` | Map sequences to the physical map (PAF output, for visualization) |
| `bed-to-path` | Convert BED mappings to scaffold paths |
| `path-to-fasta` | Produce scaffolded FASTA from scaffold paths |
| **Reporting** | |
| `metrics` | Compute assembly metrics (N50, NG50, etc.) |
| `path-metrics` | Compute physical map metrics |
| `backbone-dot` | Generate DOT visualization of backbone paths |
| **Pipelines** | |
| `physical-map` | Run the full physical map pipeline |
| `scaffolds` | Run the full scaffolding pipeline |

## Parameters

Most parameters have sensible defaults and do not need to be changed for typical use. The tables below are for advanced users who want to tune the pipeline.

### Core Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-k` | 32 | K-mer size for minimizer extraction |
| `-w` | 32 | Window size for minimizer extraction |
| `-t` | auto | Number of threads (auto-detected, capped at 16) |
| `-v` | 1 | Verbosity: 0 = silent, 1 = info, 2 = debug |

### Barcode and Edge Filtering

These control which barcodes and edges are kept in the overlap graph.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-count` / `-n` | 100 | Minimum minimizers per barcode (removes low-coverage barcodes) |
| `--max-count` / `-N` | 5000 | Maximum minimizers per barcode (removes chimeric/noisy barcodes) |
| `--min-shared` | 10 | Minimum shared minimizers to create an overlap edge |
| `--percentile` / `-p` | 90 | Remove the bottom N% of edges by weight. Use ~85 for stLFR, ~92.5 for 10x Chromium |

### Backbone Extraction

These control how the minimum spanning tree is pruned to extract backbone paths.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--prune-branches` | 10 | Remove branches shorter than this from MST junctions |
| `--prune-bridges` | 10 | Remove bridge edges connecting components smaller than this |
| `--prune-junctions` | 200 | Remove junction branches shorter than this |
| `--min-component-size` | 50 | Discard backbone paths shorter than this (in molecules) |

### Merge-Paths (optional)

These control the optional post-processing step that merges adjacent backbone paths. The defaults were optimized on human stLFR data (NA12878 + NA24143) for zero false positives.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--endpoint-depth` | 25 | Number of molecules from each path end used as endpoints |
| `--min-endpoint-hits` | 4 | A bridge molecule must connect to ≥ this many endpoint molecules on each side |
| `--min-bridges` | 2 | Minimum bridge molecules required to accept a merge |
| `--min-shared-mx` | 3 | Minimum shared minimizers between a bridge and an endpoint molecule |
| `--max-connections` | 2 | A bridge molecule connecting > this many paths is discarded (specificity filter) |
| `--max-links-per-endpoint` | 1 | Endpoints with > this many candidate links are discarded (ambiguity filter) |
| `--min-bridge-density` | 0.01 | Minimum ratio of bridges to the shorter path length |

### Scaffolding and Reporting

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-score` / `-n` | 10 | Minimum mapping score when mapping contigs to the physical map |
| `--gap-size` | 100 | Number of Ns inserted between scaffolded contigs |
| `-g` | — | (Optional) Expected genome size in bp. Only used for NG50 in metrics output |

## Reproducing Results

See [REPRODUCING.md](REPRODUCING.md) for step-by-step instructions to reproduce the NA12878 and NA24143 physical map results shown above, including data download links and the complete pipeline script.

## Comparison with Physlr v1

| Aspect | Physlr v1 | Physlr 2 |
|--------|-----------|----------|
| Language | Python + C++ | Rust |
| Graph library | NetworkX | petgraph |
| Pipeline driver | GNU Make | Snakemake + CLI subcommands |
| Minimizer extraction | External (indexlr/btllib) | Built-in + btllib support |
| Parallelism | Limited | rayon (data-parallel) |
| Memory | Python dicts + NetworkX | Compact hash maps (rustc-hash) |
| Molecule separation | Louvain communities | Biconnected components / k-clique / cosine similarity |
| Path merging | — | Bridge molecule evidence (merge-paths) |

## Project Structure

```
physlr2/
├── Cargo.toml                 # Rust package manifest
├── src/
│   ├── main.rs                # CLI entry point (clap)
│   ├── lib.rs                 # Library root
│   ├── minimizer/mod.rs       # Minimizer extraction and filtering
│   ├── overlap/mod.rs         # Barcode overlap computation
│   ├── molecules/mod.rs       # Molecule separation
│   ├── graph/mod.rs           # Graph algorithms (MST, pruning)
│   ├── backbone/mod.rs        # Backbone extraction + merge-paths
│   ├── map/mod.rs             # Mapping to physical map
│   ├── scaffold/mod.rs        # Assembly scaffolding
│   ├── repeat/mod.rs          # Repeat k-mer detection (Bloom filter)
│   ├── report/mod.rs          # Metrics and reporting
│   └── io/mod.rs              # File I/O (TSV, FASTA, BED, gzip)
├── scripts/
│   ├── plotpaf.py             # Backbone-vs-reference visualization
│   ├── find-ntcard-mode.py    # K-mer histogram mode finder
│   └── profile_pipeline.sh    # Pipeline profiling (time + memory)
├── workflow/
│   ├── Snakefile              # Snakemake pipeline
│   └── scripts/               # Workflow helper scripts
├── results/                   # Example result plots
├── REPRODUCING.md             # Reproducibility instructions
└── LICENSE                    # GPL-3.0
```

## Citation

If you use Physlr, please cite:

> Afshinfard, A., Jackman, S.D., Wong, J., Coombe, L., Nikolic, V., Chu, J., Mohamadi, H., & Birol, I. (2022).
> Physlr: Next-Generation Physical Maps. *DNA*, 2(2), 116–130.
> [https://doi.org/10.3390/dna2020009](https://doi.org/10.3390/dna2020009)

## Support

[Create a new issue on GitHub.](https://github.com/aafshinfard/physlr2/issues)

## License

[GNU General Public License v3.0](LICENSE)
