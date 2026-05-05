#!/bin/bash
#SBATCH --job-name=full-na24143
#SBATCH --output=./logs/full_na24143.%j.out
#SBATCH --error=./logs/full_na24143.%j.err
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
PREV_OUTDIR="${BASEDIR}/output/na24143_stlfr"
OUTDIR="${BASEDIR}/output/na24143_stlfr_final"
PHYSLR="${BASEDIR}/physlr"
REF="${BASEDIR}/data/GRCh38.fa"
CONDA_BIN="/home/afshinfa/miniconda3/envs/physlr-tools/bin"
K=32
W=32

echo "============================================================"
echo "physlr-next full pipeline (k3 fix + split-minimizers) — NA24143"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "============================================================"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Symlink pre-computed files from previous run
ln -sf "${PREV_OUTDIR}/sample.overlap.filtered.tsv" sample.overlap.filtered.tsv 2>/dev/null || true
ln -sf "${PREV_OUTDIR}/sample.filtered.tsv" sample.filtered.tsv 2>/dev/null || true

echo "=== Step 5: Molecules (distributed+sqcosbin) ==="
START=$(date +%s)
${PHYSLR} -v2 -t16 molecules \
    --strategy distributed+sqcosbin \
    --sqcos-threshold 0.75 \
    --skip-small 10 \
    --bin-max-size 50 \
    --merge-cutoff 20 \
    sample.overlap.filtered.tsv \
    -o sample.molecules.tsv 2>&1
echo "Molecules done in $(($(date +%s) - START))s"

BLANK=$(grep -n '^$' sample.molecules.tsv | head -1 | cut -d: -f1)
TOTAL=$(wc -l < sample.molecules.tsv)
VERTS=$((BLANK - 2))
EDGES=$((TOTAL - BLANK - 1))
echo "Molecule graph: V=${VERTS} E=${EDGES}"

echo ""
echo "=== Step 6: Backbone ==="
START=$(date +%s)
${PHYSLR} -v2 backbone \
    --prune-branches 10 \
    --prune-bridges 10 \
    --prune-junctions 200 \
    --min-component-size 50 \
    sample.molecules.tsv \
    -o sample.backbone.path 2>&1
echo "Backbone done in $(($(date +%s) - START))s"
echo "Backbone paths: $(wc -l < sample.backbone.path)"
echo "Top 10 path sizes:"
awk '{print NF}' sample.backbone.path | sort -rn | head -10

echo ""
echo "=== Step 7: Split minimizers ==="
START=$(date +%s)
${PHYSLR} -v2 -t16 split-minimizers \
    sample.molecules.tsv \
    sample.filtered.tsv \
    -o sample.split.tsv 2>&1
echo "Split minimizers done in $(($(date +%s) - START))s"

echo ""
echo "=== Step 8: Visualization ==="
if [ -f "${REF}" ]; then
    if [ ! -f sample.ref.tsv ]; then
        START=$(date +%s)
        ${CONDA_BIN}/indexlr -t5 -k ${K} -w ${W} \
            -o sample.ref.tsv \
            "${REF}" 2>&1
        echo "Index ref done in $(($(date +%s) - START))s"
    else
        echo "Reference index already exists, skipping"
    fi

    echo ""
    echo "=== Step 9a: Map backbone (barcode minimizers) ==="
    START=$(date +%s)
    ${PHYSLR} -v2 map-paf \
        sample.backbone.path \
        sample.filtered.tsv \
        sample.ref.tsv \
        -o sample.backbone.paf \
        -n 1 2>&1
    echo "Map PAF done in $(($(date +%s) - START))s"

    echo ""
    echo "=== Step 9b: Map backbone (split minimizers) ==="
    START=$(date +%s)
    ${PHYSLR} -v2 map-paf \
        sample.backbone.path \
        sample.split.tsv \
        sample.ref.tsv \
        -o sample.backbone.split.paf \
        -n 1 \
        --mx-type split 2>&1
    echo "Map PAF split done in $(($(date +%s) - START))s"

    echo ""
    echo "=== Step 10: Generate plots ==="
    START=$(date +%s)
    ${CONDA_BIN}/python3 ${BASEDIR}/physlr-next/scripts/plotpaf.py \
        sample.backbone.paf \
        sample.backbone \
        1 2>&1
    ${CONDA_BIN}/python3 ${BASEDIR}/physlr-next/scripts/plotpaf.py \
        sample.backbone.split.paf \
        sample.backbone.split \
        1 2>&1
    echo "Plots done in $(($(date +%s) - START))s"
else
    echo "Reference not found at ${REF}, skipping visualization"
fi

echo ""
echo "=== Summary ==="
echo "Original (cpython): V=3,254,106 E=60,662,752 Paths=87"
echo "physlr-next (final): V=${VERTS} E=${EDGES} Paths=$(wc -l < sample.backbone.path)"
echo ""
echo "Job complete: $(date)"
