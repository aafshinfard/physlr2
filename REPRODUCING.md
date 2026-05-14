# Reproducing the Physical Map Results

Step-by-step instructions to reproduce the Physlr 2 physical map results for **NA12878** and **NA24143** using stLFR linked reads.

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
- **btllib** with `indexlr`, `ntcard`, `nthits` (for repeat filtering and minimizer extraction; `conda install -c bioconda btllib`)
- **Python 3.8+** with matplotlib (for visualization)
- **samtools** (optional, for reference indexing)

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

## Quick Start

The provided script automates the entire pipeline for both samples:

```bash
# 1. Build physlr
cargo build --release

# 2. Run the full reproducibility pipeline
#    Downloads data, runs the pipeline, and generates plots.
bash scripts/reproduce.sh /path/to/workdir
```

The script is resumable — it skips completed steps and uses `wget -c` for partial downloads.

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

### 3. Run the Physical Map Pipeline

The `physical-map` command runs the full pipeline from FASTQ to backbone paths, using external tools (ntcard, nthits, indexlr from btllib) for repeat filtering and minimizer extraction:

```bash
SAMPLE=na12878
OUTDIR=$WORKDIR/output/$SAMPLE

physlr physical-map \
  -k 32 -w 32 -t 16 \
  -o $OUTDIR -p $SAMPLE \
  --reference $WORKDIR/data/grch38.fa \
  $WORKDIR/data/${SAMPLE}.R1.fq.gz $WORKDIR/data/${SAMPLE}.R2.fq.gz
```

This produces:
- `$SAMPLE.backbone.path` — backbone paths (physical map)
- `$SAMPLE.filtered.tsv` — filtered minimizers
- `$SAMPLE.mol.tsv` — molecule graph
- `$SAMPLE.backbone.paf` — PAF mapping to reference (when `--reference` is given)
- `$SAMPLE.positions.tsv` — position map for coordinate conversion
- Backbone and reference plots (when `--reference` is given and plotpaf.py is found)

### 4. Merge Adjacent Paths (Optional)

Merge-paths identifies non-backbone "bridge" molecules to merge adjacent backbone paths:

```bash
# Split minimizers by molecule (needed for merge-paths)
physlr split-minimizers \
  $OUTDIR/$SAMPLE.mol.tsv $OUTDIR/$SAMPLE.filtered.tsv \
  -o $OUTDIR/$SAMPLE.split.tsv -t 16

# Merge adjacent paths
physlr merge-paths \
  $OUTDIR/$SAMPLE.backbone.path $OUTDIR/$SAMPLE.split.tsv \
  -o $OUTDIR/$SAMPLE.merged.path

# Report metrics
physlr path-metrics $OUTDIR/$SAMPLE.merged.path
```

### 5. Visualize Merged Paths Against Reference

```bash
# Map merged paths to reference
physlr map-paf \
  $OUTDIR/$SAMPLE.merged.path $OUTDIR/$SAMPLE.split.tsv $OUTDIR/$SAMPLE.ref.tsv \
  -o $OUTDIR/$SAMPLE.merged.paf -n 1 --mx-type split

# Generate plots
python3 scripts/plotpaf.py $OUTDIR/$SAMPLE.merged.paf $OUTDIR/$SAMPLE.merged \
  --positions $OUTDIR/$SAMPLE.positions.tsv
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

### Output Files

| File | Description |
|------|-------------|
| `*.backbone.path` | Backbone paths (physical map, before merging) |
| `*.merged.path` | Merged backbone paths (after merge-paths) |
| `*.backbone.paf` | PAF alignment of backbone to reference |
| `*.backbone.png` | Backbone coverage plot |
| `*.reference.png` | Reference coverage plot |

## Troubleshooting

**Out of memory during overlap computation:**
Increase `--min-bx-count` to reduce the number of barcodes, or increase available RAM.

**Missing external tools (ntcard, nthits, indexlr):**
Install btllib: `conda install -c bioconda btllib`

**Missing plots:**
Ensure Python 3 with matplotlib is installed: `pip install matplotlib`.

**Resuming interrupted runs:**
The `scripts/reproduce.sh` script checks for existing output files and skips completed steps. Re-run it to resume.
