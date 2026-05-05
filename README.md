# Physlr v2

**Next-generation physical maps from linked reads.**

Physlr constructs de novo physical maps using linked reads (10x Genomics Chromium or MGI stLFR) and uses them to scaffold genome assemblies. This is a ground-up rewrite of the [original Physlr](https://github.com/bcgsc/physlr) in Rust for performance and maintainability.

## Overview

Physlr takes linked-read sequencing data and:

1. **Builds a physical map** — an ordered set of molecules along each chromosome
2. **Scaffolds draft assemblies** — orders and orients contigs using the physical map

### Pipeline stages

```
Linked reads (FASTQ + barcodes)
  │
  ├─ Index minimizers
  ├─ Filter barcodes & minimizers
  ├─ Compute barcode overlap graph
  ├─ Separate barcodes into molecules
  ├─ Maximum spanning tree → prune → backbone paths
  │
  └─► Physical map (backbone paths)
        │
        ├─ Map draft contigs to physical map
        ├─ Order & orient contigs
        │
        └─► Scaffolded assembly (FASTA)
```

## Installation

### From source (requires Rust 1.70+)

```bash
git clone <repo-url>
cd physlr-next
cargo build --release
# Binary at target/release/physlr
```

### Dependencies

- **Rust** ≥ 1.70 (build)
- **Python 3** ≥ 3.10 (Snakemake workflow, visualization scripts)
- **Snakemake** ≥ 7.0 (optional, for automated pipeline)
- **Quast** (optional, for reference-based assembly evaluation)
- **Graphviz** (optional, for rendering DOT visualizations)

## Quick start

### Step-by-step CLI

```bash
# 1. Index minimizers from linked reads
physlr index reads.fq.gz -o reads.mxs.tsv -k 32 -w 32

# 2. Filter barcodes and minimizers
physlr filter-minimizers reads.mxs.tsv -o filtered.mxs.tsv -n 100 -N 5000

# 3. Compute overlap graph
physlr overlap filtered.mxs.tsv -o overlap.tsv --min-shared 10

# 4. Separate molecules
physlr molecules overlap.tsv -o mol.tsv

# 5. Extract backbone paths (physical map)
physlr backbone mol.tsv -o backbone.path

# 6. Index draft assembly contigs
physlr index-contigs draft.fa -o draft.mxs.tsv -k 32 -w 32

# 7. Map contigs to physical map
physlr map backbone.path filtered.mxs.tsv draft.mxs.tsv -o map.bed

# 8. Convert to scaffold paths
physlr bed-to-path map.bed -o scaffold.path

# 9. Produce scaffolded FASTA
physlr path-to-fasta draft.fa scaffold.path -o scaffolds.fa
```

### One-command pipelines

```bash
# Physical map only
physlr physical-map reads.mxs.tsv -o output/ -p mygenome

# Full scaffolding
physlr scaffolds reads.mxs.tsv draft.fa draft.mxs.tsv -o output/ -p mygenome -g 300000000
```

### Snakemake workflow

```bash
# Edit workflow/config.yaml with your file paths and parameters
snakemake -s workflow/Snakefile --configfile config.yaml -j 8
```

## Commands

| Command | Description |
|---------|-------------|
| `index` | Extract (k,w)-minimizers from FASTA/FASTQ, grouping by barcode |
| `index-contigs` | Extract ordered minimizers from FASTA (for contigs/reference) |
| `filter-minimizers` | Filter barcodes by count, remove singleton/repetitive minimizers |
| `overlap` | Compute barcode overlap graph from shared minimizers |
| `filter-overlap` | Remove low-weight edges by percentile |
| `molecules` | Separate barcodes into individual molecules |
| `backbone` | Extract backbone paths (physical map) from overlap graph |
| `map` | Map query sequences to the physical map |
| `bed-to-path` | Convert BED mappings to scaffold paths |
| `path-to-fasta` | Produce scaffolded FASTA from paths |
| `metrics` | Compute assembly metrics (N50, NG50, etc.) |
| `path-metrics` | Compute physical map metrics |
| `backbone-dot` | Generate DOT visualization of backbone paths |
| `physical-map` | Run the full physical map pipeline |
| `scaffolds` | Run the full scaffolding pipeline |

## Parameters

Key parameters and their defaults:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-k` | 32 | K-mer size for minimizers |
| `-w` | 32 | Window size for minimizers |
| `--min-bx-count` | 100 | Minimum minimizers per barcode |
| `--max-bx-count` | 5000 | Maximum minimizers per barcode |
| `--min-overlap` / `--min-shared` | 10 | Minimum shared minimizers for an edge |
| `--edge-percentile` | 90 | Percentile of edges to remove |
| `--prune-branches` | 10 | Minimum branch size in MST |
| `--min-path-size` / `--min-component-size` | 50 | Minimum backbone path length |
| `--min-map-score` | 10 | Minimum mapping score |
| `--gap-size` | 100 | Gap size (Ns) between scaffolded contigs |
| `-g` | — | Expected genome size (for NG50 reporting) |

## Project structure

```
physlr-next/
├── Cargo.toml              # Rust package manifest
├── src/
│   ├── main.rs             # CLI entry point (clap)
│   ├── lib.rs              # Library root
│   ├── minimizer/mod.rs    # Minimizer extraction and filtering
│   ├── overlap/mod.rs      # Barcode overlap computation (parallel)
│   ├── molecules/mod.rs    # Molecule separation (biconnected components)
│   ├── graph/mod.rs        # Graph types and algorithms (MST, pruning, paths)
│   ├── backbone/mod.rs     # Backbone path extraction
│   ├── map/mod.rs          # Mapping sequences to physical map
│   ├── scaffold/mod.rs     # Scaffolding (ordering, orienting, joining)
│   ├── report/mod.rs       # Metrics computation and reporting
│   └── io/mod.rs           # File I/O (TSV, FASTA, BED, gzip)
├── workflow/
│   ├── Snakefile           # Snakemake pipeline
│   └── scripts/            # Helper scripts (data download, visualization)
├── data/                   # Test datasets
└── tests/                  # Integration tests
```

## Design decisions

| Aspect | Original Physlr | Physlr v2 |
|--------|----------------|-----------|
| Language | Python + C++ | Rust (single language) |
| Graph library | NetworkX (Python) | petgraph (Rust) |
| Pipeline | GNU Make | Snakemake + CLI subcommands |
| Minimizer extraction | External (indexlr/btllib) | Built-in |
| Parallelism | Limited | rayon (data-parallel) |
| Memory | Python dicts + NetworkX overhead | Compact hash maps (rustc-hash) |
| Index stability | N/A | StableUnGraph (safe node removal) |

## Testing

```bash
# Generate synthetic test data
python3 workflow/scripts/generate_synthetic_test.py -o data/synthetic

# Run pipeline on synthetic data
physlr physical-map data/synthetic/linked_reads.fq \
  -o data/synthetic/output -p test \
  --min-bx-count 2 --min-overlap 2 --edge-percentile 0 \
  --prune-branches 2 --min-path-size 3

# Compute metrics
physlr metrics data/synthetic/output/test.scaffolds.fa -g 300000
```

## Citation

If you use Physlr, please cite:

> Afshinfard, A., Jackman, S.D., Wong, J., Coombe, L., Nikolic, V., Chu, J., Mohamadi, H., & Birol, I. (2022).
> Physlr: Next-Generation Physical Maps. *DNA*, 2(2), 116–130.
> https://doi.org/10.3390/dna2020009

## License

GPL-3.0
