#!/bin/bash
#SBATCH --job-name=next-dist-na12878
#SBATCH --output=./logs/next_distributed_na12878.%j.out
#SBATCH --error=./logs/next_distributed_na12878.%j.err
#SBATCH --mem=300G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
PHYSLR="${BASEDIR}/physlr-next/target/release/physlr"
PREV_OUTDIR="${BASEDIR}/output/na12878_stlfr"
OUTDIR="${BASEDIR}/output/na12878_stlfr_distributed"

source activate physlr-tools 2>/dev/null || conda activate physlr-tools 2>/dev/null || true

echo "============================================================"
echo "Physlr-next (distributed only) — NA12878 stLFR"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "============================================================"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

ln -sf "${PREV_OUTDIR}/sample.overlap.filtered.tsv" sample.overlap.filtered.tsv 2>/dev/null || true

echo "=== Molecules (distributed only) ==="
START=$(date +%s)
"${PHYSLR}" molecules sample.overlap.filtered.tsv \
    -o sample.molecules.tsv \
    --strategy "distributed"
echo "Done in $(($(date +%s) - START))s"

echo ""
echo "=== Backbone ==="
START=$(date +%s)
"${PHYSLR}" backbone sample.molecules.tsv \
    -o sample.backbone.tsv \
    --prune-branches 10 --prune-bridges 10
echo "Done in $(($(date +%s) - START))s"

echo ""
echo "=== Summary ==="
echo "Backbone paths: $(wc -l < sample.backbone.tsv)"
echo "Total backbone nodes: $(awk '{print NF}' sample.backbone.tsv | paste -sd+ | bc)"
echo "Top 10 path sizes:"
awk '{print NF}' sample.backbone.tsv | sort -rn | head -10

echo ""
echo "Job complete: $(date)"
