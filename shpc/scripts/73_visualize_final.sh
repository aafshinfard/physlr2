#!/bin/bash
#SBATCH --job-name=viz-final
#SBATCH --output=./logs/viz_final.%j.out
#SBATCH --error=./logs/viz_final.%j.err
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
PHYSLR="${BASEDIR}/physlr"
CONDA_BIN="/home/afshinfa/miniconda3/envs/physlr-tools/bin"
PLOTPAF="${BASEDIR}/physlr-next/scripts/plotpaf.py"

echo "============================================================"
echo "Visualization: final pipeline results"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "============================================================"

for CELL in na12878_stlfr na24143_stlfr; do
    OUTDIR="${BASEDIR}/output/${CELL}_final"
    PREV="${BASEDIR}/output/${CELL}"

    echo ""
    echo "=== ${CELL} ==="

    cd "${OUTDIR}"

    # Symlink reference index from previous run
    ln -sf "${PREV}/sample.ref.tsv" sample.ref.tsv 2>/dev/null || true

    if [ ! -f sample.ref.tsv ]; then
        echo "No reference index found for ${CELL}, skipping"
        continue
    fi

    # Map backbone to reference (barcode minimizers)
    echo "--- Map PAF (barcode minimizers) ---"
    START=$(date +%s)
    ${PHYSLR} -v2 map-paf \
        sample.backbone.path \
        sample.filtered.tsv \
        sample.ref.tsv \
        -o sample.backbone.paf \
        -n 1 2>&1
    echo "Done in $(($(date +%s) - START))s"

    # Map backbone to reference (split minimizers)
    echo "--- Map PAF (split minimizers) ---"
    START=$(date +%s)
    ${PHYSLR} -v2 map-paf \
        sample.backbone.path \
        sample.split.tsv \
        sample.ref.tsv \
        -o sample.backbone.split.paf \
        -n 1 \
        --mx-type split 2>&1
    echo "Done in $(($(date +%s) - START))s"

    # Generate plots
    echo "--- Plots ---"
    ${CONDA_BIN}/python3 ${PLOTPAF} sample.backbone.paf sample.backbone 1 2>&1
    ${CONDA_BIN}/python3 ${PLOTPAF} sample.backbone.split.paf sample.backbone.split 1 2>&1

    echo "Generated: sample.backbone.*.png and sample.backbone.split.*.png"
    ls -la sample.backbone*.png 2>/dev/null || echo "No PNG files found"
done

echo ""
echo "Job complete: $(date)"
