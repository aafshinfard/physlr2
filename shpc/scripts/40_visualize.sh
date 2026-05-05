#!/bin/bash
#SBATCH --job-name=physlr-visualize
#SBATCH --output=./logs/visualize.%j.out
#SBATCH --error=./logs/visualize.%j.err
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
PHYSLR="${BASEDIR}/physlr-next/target/release/physlr"
SCRIPT_DIR="${BASEDIR}/scripts"
REFERENCE="${BASEDIR}/data/reference/grch38.fa"
K=40
W=32

# Activate conda env
source activate physlr-tools 2>/dev/null || conda activate physlr-tools 2>/dev/null || true

echo "============================================================"
echo "Physlr Visualization — Backbone vs Reference"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "============================================================"

for DATASET in na12878_stlfr na24143_stlfr; do
    OUTDIR="${BASEDIR}/output/${DATASET}"
    echo ""
    echo "=== Processing ${DATASET} ==="

    if [ ! -f "${OUTDIR}/sample.backbone.tsv" ]; then
        echo "ERROR: ${OUTDIR}/sample.backbone.tsv not found. Skipping."
        continue
    fi

    if [ ! -f "${OUTDIR}/sample.filtered.tsv" ]; then
        echo "ERROR: ${OUTDIR}/sample.filtered.tsv not found. Skipping."
        continue
    fi

    cd "${OUTDIR}"

    # Step 1: Index reference with indexlr (if not already done)
    if [ ! -f "sample.ref.tsv" ]; then
        echo "Indexing reference..."
        indexlr -t5 -k "${K}" -w "${W}" -o "sample.ref.tsv" "${REFERENCE}"
        echo "Reference indexed: $(wc -l < sample.ref.tsv) contigs"
    else
        echo "Reference already indexed: $(wc -l < sample.ref.tsv) contigs"
    fi

    # Step 2: Map backbone to reference
    echo "Mapping backbone to reference..."
    "${PHYSLR}" map-paf \
        "sample.backbone.tsv" \
        "sample.filtered.tsv" \
        "sample.ref.tsv" \
        -o "sample.backbone.paf" \
        -n 1
    echo "PAF records: $(wc -l < sample.backbone.paf)"

    # Step 3: Generate plots
    echo "Generating plots..."
    python3 "${SCRIPT_DIR}/plotpaf.py" \
        "sample.backbone.paf" \
        "sample.backbone" \
        1

    echo "Done: ${DATASET}"
    ls -lh sample.backbone.*.png 2>/dev/null || echo "No PNG files generated"
done

echo ""
echo "============================================================"
echo "Visualization complete: $(date)"
echo "============================================================"
