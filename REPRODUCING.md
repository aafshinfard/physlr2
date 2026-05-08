# Reproducing the Physical Map Results

This document provides step-by-step instructions to reproduce the Physlr 2 physical map results for two human cell lines: **NA12878** and **NA24143** using stLFR linked reads.

## Contents

- [Requirements](#requirements)
- [Data Sources](#data-sources)
- [Quick Start](#quick-start)
- [Step-by-Step Instructions](#step-by-step-instructions)
- [Expected Results](#expected-results)
- [Troubleshooting](#troubleshooting)

## Requirements

### Hardware

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| RAM | 120 GB | 128 GB (for a human-size genome) |
| CPU cores | 8 | 16 |
| Disk space | 600 GB | 1 TB |
| Wall time | 24 hours | 12 hours |

The overlap computation step is the most memory-intensive. Disk space is dominated by the input FASTQ files (~250 GB per sample).

### Software

- **Physlr 2** (this repo)
- **Python 3.8+** with matplotlib (for visualization)
- **btllib** with `indexlr` (for minimizer extraction; install via `conda install -c bioconda btllib`)
- **samtools** (optional, for reference indexing)
- **QUAST** (optional, for assembly evaluation; `conda install -c bioconda quast`)

## Data Sources

### Linked Reads (stLFR)

| Sample | File | URL | Size |
|--------|------|-----|------|
| NA12878 R1 | `stLFR.split_read.1.fq.gz` | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/stLFR.split_read.1.fq.gz) | ~120 GB |
| NA12878 R2 | `stLFR.split_read.2.fq.gz` | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/stLFR.split_read.2.fq.gz) | ~120 GB |
| NA24143 R1 | `stLFR_NA24143_split_read.1.fq.gz` | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/stLFR/stLFR_NA24143_split_read.1.fq.gz) | ~110 GB |
| NA24143 R2 | `stLFR_NA24143_split_read.2.fq.gz` | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/stLFR/stLFR_NA24143_split_read.2.fq.gz) | ~110 GB |

### Reference Genome

| File | URL |
|------|-----|
| GRCh38 (no alt) | [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) |

### Draft Assemblies (for scaffolding evaluation)

| Sample | Assembly | URL |
|--------|----------|-----|
| NA12878 | ONT Shasta | [BCGSC](https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na12878.ont.baseline.fa) |
| NA12878 | PE+MPET ABySS | [BCGSC](https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na12878.pempet.baseline.fa) |
| NA12878 | stLFR Supernova | [BCGSC](https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na12878.stlfr.baseline.fa) |
| NA12878 | PacBio Peregrine | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/JasonChin_Peregrine_PacBioCCS_assembly_05072019/NA12878.ccs.peregrine.fa.gz) |
| NA24143 | ONT Shasta | [BCGSC](https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na24143.ont.baseline.fa) |
| NA24143 | PE+MPET ABySS | [BCGSC](https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na24143.pempet.baseline.fa) |
| NA24143 | stLFR Supernova | [BCGSC](https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na24143.stlfr.baseline.fa) |
| NA24143 | PacBio Falcon | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/MtSinai_PacBio_Assembly_falcon_03282016/NA24143_hg004_falcon.fa) |

## Quick Start

The provided script automates the entire pipeline for both samples:

```bash
# 1. Build physlr
cargo build --release

# 2. Run the full reproducibility pipeline
#    This downloads data, runs the pipeline, and generates plots.
bash scripts/reproduce.sh /path/to/workdir
```

The script is resumable — it skips completed steps and uses `wget -c` for partial downloads. Re-run it if interrupted.

## Step-by-Step Instructions

### 1. Download Data

```bash
WORKDIR=/path/to/workdir
mkdir -p $WORKDIR/data

# Reference genome
wget -c -O $WORKDIR/data/grch38.fa.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
gunzip $WORKDIR/data/grch38.fa.gz

# NA12878 stLFR reads
wget -c -O $WORKDIR/data/na12878.R1.fq.gz \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/stLFR.split_read.1.fq.gz"
wget -c -O $WORKDIR/data/na12878.R2.fq.gz \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/stLFR.split_read.2.fq.gz"

# NA24143 stLFR reads
wget -c -O $WORKDIR/data/na24143.R1.fq.gz \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/stLFR/stLFR_NA24143_split_read.1.fq.gz"
wget -c -O $WORKDIR/data/na24143.R2.fq.gz \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/stLFR/stLFR_NA24143_split_read.2.fq.gz"
```

### 2. Build Physlr

```bash
cd /path/to/physlr2
cargo build --release
export PATH="$(pwd)/target/release:$PATH"
```

### 3. Run the Pipeline

For each sample (`na12878` or `na24143`):

```bash
SAMPLE=na12878
OUTDIR=$WORKDIR/output/$SAMPLE
mkdir -p $OUTDIR && cd $OUTDIR

# Index minimizers (k=32, w=32)
physlr index -k 32 -w 32 \
  $WORKDIR/data/${SAMPLE}.R1.fq.gz $WORKDIR/data/${SAMPLE}.R2.fq.gz \
  -o ${SAMPLE}.reads.tsv

# Filter barcodes
physlr filter-minimizers ${SAMPLE}.reads.tsv \
  -o ${SAMPLE}.filtered.tsv -n 100 -N 5000

# Compute overlap graph
physlr overlap ${SAMPLE}.filtered.tsv \
  -o ${SAMPLE}.overlap.tsv

# Filter edges (85th percentile for stLFR)
physlr filter-overlap ${SAMPLE}.overlap.tsv \
  -o ${SAMPLE}.overlap.filtered.tsv -p 85

# Separate molecules
physlr molecules ${SAMPLE}.overlap.filtered.tsv \
  -o ${SAMPLE}.molecules.tsv --strategy distributed+sqcosbin

# Extract backbone paths
physlr backbone ${SAMPLE}.molecules.tsv \
  -o ${SAMPLE}.backbone.path \
  --prune-branches 10 --prune-bridges 10 --prune-junctions 200

# Split minimizers (for merge-paths and visualization)
physlr split-minimizers ${SAMPLE}.molecules.tsv ${SAMPLE}.filtered.tsv \
  -o ${SAMPLE}.split.tsv

# Merge adjacent backbone paths
physlr merge-paths ${SAMPLE}.backbone.path ${SAMPLE}.split.tsv \
  -o ${SAMPLE}.merged.path

# Physical map metrics
physlr path-metrics ${SAMPLE}.merged.path
```

### 4. Visualize Against Reference

```bash
# Index reference minimizers (requires indexlr from btllib)
indexlr -t 5 -k 32 -w 32 -o ${SAMPLE}.ref.tsv $WORKDIR/data/grch38.fa

# Map backbone to reference
physlr map-paf ${SAMPLE}.merged.path ${SAMPLE}.split.tsv ${SAMPLE}.ref.tsv \
  -o ${SAMPLE}.paf -n 1 --mx-type split

# Generate plots
python3 /path/to/physlr2/scripts/plotpaf.py ${SAMPLE}.paf ${SAMPLE} 1
# Produces: ${SAMPLE}.backbone.png and ${SAMPLE}.reference.png
```

## Expected Results

### Physical Map Statistics

| Metric | NA12878 | NA24143 |
|--------|---------|---------|
| Backbone paths | ~180 | ~60 |
| Merged paths (after merge-paths) | ~170 | ~55 |
| Largest path (molecules) | ~6,000 | ~5,000 |

### Merge-Paths Performance (v0.23 defaults)

| Sample | True Positives | False Positives |
|--------|---------------|-----------------|
| NA12878 | 4 | 0 |
| NA24143 | 5 | 0 |

True/false positives are determined by mapping merged paths to the GRCh38 reference and checking whether merged path pairs map to the same chromosome in the correct order.

### Output Files

| File | Description |
|------|-------------|
| `*.backbone.path` | Backbone paths (physical map, before merging) |
| `*.merged.path` | Merged backbone paths (after merge-paths) |
| `*.backbone.png` | Backbone coverage plot |
| `*.reference.png` | Reference coverage plot |
| `*.paf` | PAF alignment of backbone to reference |

## Troubleshooting

**Out of memory during overlap computation:**
Increase `--min-count` in `filter-minimizers` to reduce the number of barcodes, or increase available RAM.

**Slow minimizer indexing:**
Use `--indexer btllib` with `-t` for multi-threaded indexing if btllib/indexlr is installed.

**Missing plots:**
Ensure Python 3 with matplotlib is installed: `pip install matplotlib`.

**Resuming interrupted runs:**
The `scripts/reproduce.sh` script checks for existing output files and skips completed steps. Re-run it to resume.

