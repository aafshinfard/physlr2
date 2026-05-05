#!/bin/bash
#SBATCH --job-name=orig-dist-na12878
#SBATCH --output=./logs/orig_distributed_na12878.%j.out
#SBATCH --error=./logs/orig_distributed_na12878.%j.err
#SBATCH --mem=300G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
PREV_OUTDIR="${BASEDIR}/output/na12878_stlfr"
OUTDIR="${BASEDIR}/output/na12878_stlfr_original_distributed"

CONDA_PYTHON="/home/afshinfa/miniconda3/envs/physlr-cpython/bin/python3"
PHYSLR_DIR="/home/afshinfa/miniconda3/envs/physlr-original/bin/share/physlr-1.0.4-8"
export PYTHONPATH="${PHYSLR_DIR}"
PHYSLR_PY="${CONDA_PYTHON} ${PHYSLR_DIR}/bin/physlr"

echo "============================================================"
echo "Original Physlr (distributed only) — NA12878"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "============================================================"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

ln -sf "${PREV_OUTDIR}/sample.overlap.filtered.tsv" sample.overlap.filtered.tsv 2>/dev/null || true

echo "=== Molecules (distributed only) ==="
START=$(date +%s)
${PHYSLR_PY} molecules -t1 --separation-strategy=distributed sample.overlap.filtered.tsv > sample.molecules.tsv
echo "Done in $(($(date +%s) - START))s"

BLANK=$(grep -n '^$' sample.molecules.tsv | head -1 | cut -d: -f1)
TOTAL=$(wc -l < sample.molecules.tsv)
VERTS=$((BLANK - 2))
EDGES=$((TOTAL - BLANK - 1))
echo "Molecule graph: V=${VERTS} E=${EDGES}"

echo ""
echo "=== Backbone ==="
START=$(date +%s)
${PHYSLR_PY} backbone --prune-branches=10 --prune-bridges=10 --prune-junctions=200 sample.molecules.tsv > sample.backbone.path
echo "Done in $(($(date +%s) - START))s"
echo "Backbone paths: $(wc -l < sample.backbone.path)"
echo "Top 5 path sizes:"
awk '{print NF}' sample.backbone.path | sort -rn | head -5

echo ""
echo "Job complete: $(date)"
