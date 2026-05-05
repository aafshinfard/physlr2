# Physlr Benchmarking on SHPC

Benchmark scripts for running the Physlr Rust rewrite on real human genome data using SLURM.

## Directory Structure

```
shpc/
├── README.md
├── scripts/
│   ├── common.sh                 # Shared functions and variables
│   ├── 00_download_data.sh       # Download all input data
│   ├── 01_build_physlr.sh        # Build the Physlr binary
│   ├── 02_preprocess_10x.sh      # Preprocess 10x Chromium reads (SLURM)
│   ├── 10_na12878_stlfr.sh       # NA12878 stLFR pipeline (SLURM)
│   ├── 11_na24143_stlfr.sh       # NA24143 stLFR pipeline (SLURM)
│   ├── 12_na24143_10x.sh         # NA24143 10x Chromium pipeline (SLURM)
│   └── 20_collect_results.sh     # Collect and summarize results
├── data/                         # Downloaded input data (created by 00_download_data.sh)
├── output/                       # Pipeline output (created by pipeline scripts)
├── logs/                         # SLURM job logs
└── physlr-next/                  # Physlr source code (copy from workspace)
```

## Setup

### 1. Copy files to SHPC

```bash
SHPC_DIR=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1

# Copy the benchmark scripts
scp -r shpc/* user@shpc:${SHPC_DIR}/

# Copy the Physlr source code
scp -r /workspaces/physlr-next user@shpc:${SHPC_DIR}/physlr-next
```

### 2. Build Physlr

```bash
cd ${SHPC_DIR}
sbatch scripts/01_build_physlr.sh
```

This requires Rust (install via `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh` if not available).
Runs as a SLURM job (~5 min). Check `logs/build_*.out` for status.

### 3. Download data

```bash
cd ${SHPC_DIR}
sbatch scripts/00_download_data.sh
```

Runs as a SLURM job (hours, network-dependent). Uses `wget -c` for resumable downloads — resubmit if it times out.

**Disk space required**: ~600 GB total
- Reference genome: ~3 GB
- NA12878 stLFR reads: ~250 GB
- NA24143 stLFR reads: ~220 GB
- NA24143 10x Chromium reads: ~80 GB
- Draft assemblies: ~10 GB

The download script uses `wget -c` for resumable downloads. Run it again if interrupted.

### 4. Preprocess 10x Chromium reads (NA24143 only)

The 10x Chromium data needs barcode demultiplexing before Physlr can use it. The script supports two methods:

**Option A: longranger basic** (recommended if available)
```bash
# Ensure longranger is in PATH
sbatch scripts/02_preprocess_10x.sh
```

**Option B: Manual barcode extraction** (fallback)
The script automatically falls back to manual extraction if longranger is not found. This extracts the 16bp barcode from the first 16bp of R1 and adds it as a `BX:Z:` tag.

**stLFR reads do NOT need preprocessing** — Physlr handles the `#barcode` format natively.

## Running the Benchmarks

### Test 1: NA12878 with stLFR reads

```bash
cd ${SHPC_DIR}
sbatch scripts/10_na12878_stlfr.sh
```

Scaffolds 4 draft assemblies:
- ONT Shasta
- PacBio Peregrine
- PE+MPET ABySS
- stLFR Supernova

### Test 2: NA24143 with stLFR reads

```bash
sbatch scripts/11_na24143_stlfr.sh
```

Scaffolds 4 draft assemblies:
- ONT Shasta
- PacBio Falcon
- PE+MPET ABySS
- stLFR Supernova

### Test 3: NA24143 with 10x Chromium reads

```bash
sbatch scripts/12_na24143_10x.sh
```

Requires preprocessing step (02) to be complete first. Scaffolds 5 draft assemblies:
- ONT Shasta
- PacBio Falcon
- PE+MPET ABySS
- stLFR Supernova
- Chromium Supernova

## Resource Requirements

| Job | Memory | CPUs | Wall time | Est. Runtime |
|-----|--------|------|-----------|--------------|
| Build Physlr | 8 GB | 4 | 1 hour | ~5 min |
| Download data | 4 GB | 1 | 24 hours | hours (network) |
| Preprocess 10x | 64 GB | 16 | 24 hours | 6-12 hours |
| NA12878 stLFR | 450 GB | 28 | 24 hours | 12-24 hours |
| NA24143 stLFR | 450 GB | 28 | 24 hours | 12-24 hours |
| NA24143 10x | 450 GB | 28 | 24 hours | 12-24 hours |

All scripts run as SLURM jobs, so they survive SSH disconnects.
All scripts are resumable — they skip completed steps and `wget -c` resumes partial downloads. Resubmit if a job hits the 24h wall time.

## Collecting Results

After all jobs complete:

```bash
bash scripts/20_collect_results.sh
```

This produces `results_summary.tsv` with QUAST metrics for all baseline and scaffolded assemblies, plus physical map statistics.

## Pipeline Parameters

| Parameter | stLFR | 10x Chromium | Notes |
|-----------|-------|-------------|-------|
| k (k-mer size) | 40 | 40 | From original Physlr paper |
| w (window size) | 32 | 32 | |
| Overlap percentile | 85 | 92.5 | stLFR has lower barcode diversity |
| Min barcode count | 100 | 100 | |
| Max barcode count | 5000 | 5000 | |
| Molecule strategy | bc+cc | bc+cc | Biconnected components + connected components |
| Prune branches | 10 | 10 | |
| Prune bridges | 10 | 10 | |

## Evaluation

Each pipeline script runs:

1. **Physical map construction**: Index → filter → overlap → molecules → backbone
2. **Physical map metrics**: Path count, molecule count, path length distribution
3. **Scaffolding**: Map contigs → scaffold paths → scaffolded FASTA
4. **QUAST evaluation**: NG50, NGA50, scaffold count, misassemblies (requires `quast` in PATH)

The comparison between baseline and scaffolded assemblies shows the improvement from Physlr scaffolding.

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# View job output in real time
tail -f logs/na12878_stlfr_*.out

# Check completed job
sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS
```

## Troubleshooting

- **Out of memory**: The overlap step is the most memory-intensive. If 450 GB is insufficient, try increasing `--min-count` in the filter step to reduce the number of barcodes.
- **Slow indexing**: The indexing step processes all reads sequentially. For very large datasets (>100x coverage), consider subsampling.
- **QUAST not found**: Install via `pip install quast` or `conda install -c bioconda quast`. QUAST evaluation is optional — the pipeline still produces scaffolds without it.
- **Pipeline resumes**: All scripts check for existing output files and skip completed steps. Safe to re-run after failures.
