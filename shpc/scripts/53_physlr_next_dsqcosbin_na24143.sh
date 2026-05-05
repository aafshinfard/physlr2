#!/bin/bash
#SBATCH --job-name=next-dsqcos-na24143
#SBATCH --output=./logs/next_dsqcosbin_na24143.%j.out
#SBATCH --error=./logs/next_dsqcosbin_na24143.%j.err
#SBATCH --mem=300G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
PHYSLR="${BASEDIR}/physlr-next/target/release/physlr"
SCRIPT_DIR="${BASEDIR}/physlr-next/scripts"
PREV_OUTDIR="${BASEDIR}/output/na24143_stlfr"
OUTDIR="${BASEDIR}/output/na24143_stlfr_dsqcosbin"
export REFERENCE="${BASEDIR}/data/reference/grch38.fa"

source activate physlr-tools 2>/dev/null || conda activate physlr-tools 2>/dev/null || true

echo "============================================================"
echo "Physlr-next (distributed+sqcosbin) — NA24143 stLFR"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "============================================================"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

if [ ! -f sample.overlap.filtered.tsv ]; then
    ln -sf "${PREV_OUTDIR}/sample.overlap.filtered.tsv" sample.overlap.filtered.tsv
fi
if [ ! -f sample.filtered.tsv ]; then
    ln -sf "${PREV_OUTDIR}/sample.filtered.tsv" sample.filtered.tsv
fi

# Step 5: Separate molecules with distributed+sqcosbin
echo ""
echo "================================================================"
echo "  Step 5: Molecules (distributed+sqcosbin)"
echo "================================================================"
START=$(date +%s)
"${PHYSLR}" molecules sample.overlap.filtered.tsv \
    -o sample.molecules.tsv \
    --strategy "distributed+sqcosbin"
echo "Step 5 done in $(($(date +%s) - START))s"

# Step 6: Extract backbone
echo ""
echo "================================================================"
echo "  Step 6: Backbone"
echo "================================================================"
START=$(date +%s)
"${PHYSLR}" backbone sample.molecules.tsv \
    -o sample.backbone.tsv \
    --prune-branches 10 --prune-bridges 10
echo "Step 6 done in $(($(date +%s) - START))s"

# Summary
echo ""
echo "================================================================"
echo "  Summary"
echo "================================================================"
echo "Backbone paths: $(wc -l < sample.backbone.tsv)"
echo "Total backbone nodes: $(awk '{print NF}' sample.backbone.tsv | paste -sd+ | bc)"
echo "Top 10 path sizes:"
awk '{print NF}' sample.backbone.tsv | sort -rn | head -10

# Step 7: Visualization
if [ -n "${REFERENCE}" ] && [ -f "${REFERENCE}" ]; then
    echo ""
    echo "================================================================"
    echo "  Step 7: Visualization"
    echo "================================================================"

    if [ ! -f sample.ref.tsv ]; then
        ln -sf "${PREV_OUTDIR}/sample.ref.tsv" sample.ref.tsv
    fi

    START=$(date +%s)
    "${PHYSLR}" map-paf sample.backbone.tsv sample.filtered.tsv sample.ref.tsv \
        -o sample.backbone.paf -n 1
    echo "map-paf done in $(($(date +%s) - START))s"

    START=$(date +%s)
    python3 "${SCRIPT_DIR}/plotpaf.py" sample.backbone.paf sample.backbone 1
    echo "Plotting done in $(($(date +%s) - START))s"
fi

echo ""
echo "Job complete: $(date)"
