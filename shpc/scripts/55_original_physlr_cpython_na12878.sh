#!/bin/bash
#SBATCH --job-name=orig-cpython-na12878
#SBATCH --output=./logs/orig_cpython_na12878.%j.out
#SBATCH --error=./logs/orig_cpython_na12878.%j.err
#SBATCH --mem=300G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
PREV_OUTDIR="${BASEDIR}/output/na12878_stlfr"
OUTDIR="${BASEDIR}/output/na12878_stlfr_original_cpython"
SCRIPT_DIR="${BASEDIR}/physlr-next/scripts"
PHYSLR_RUST="${BASEDIR}/physlr-next/target/release/physlr"
export REFERENCE="${BASEDIR}/data/reference/grch38.fa"

# Use cpython env for sklearn/scipy support — use explicit path
CONDA_PYTHON="/home/afshinfa/miniconda3/envs/physlr-cpython/bin/python3"

# Physlr Python code location (from physlr-original conda package)
PHYSLR_DIR="/home/afshinfa/miniconda3/envs/physlr-original/bin/share/physlr-1.0.4-8"
export PYTHONPATH="${PHYSLR_DIR}"
PHYSLR_PY="${CONDA_PYTHON} ${PHYSLR_DIR}/bin/physlr"

# Also need C++ tools
export PATH="${BASEDIR}/physlr-src:${HOME}/miniconda3/envs/physlr-tools/bin:${PATH}"

echo "============================================================"
echo "Original Physlr (cpython) — NA12878 stLFR"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "============================================================"
echo "python3: $(${CONDA_PYTHON} --version)"
echo "sklearn: $(${CONDA_PYTHON} -c 'import sklearn; print(sklearn.__version__)')"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Reuse the filtered overlap from physlr-next (identical to original)
ln -sf "${PREV_OUTDIR}/sample.overlap.filtered.tsv" sample.overlap.filtered.tsv 2>/dev/null || true
ln -sf "${PREV_OUTDIR}/sample.filtered.tsv" sample.filtered.tsv 2>/dev/null || true

# Step 1: Molecules (distributed+sqcosbin) — the key comparison
echo ""
echo "================================================================"
echo "  Step 1: Molecules (distributed+sqcosbin via cpython)"
echo "================================================================"
START=$(date +%s)
${PHYSLR_PY} molecules -t4 --separation-strategy=distributed+sqcosbin sample.overlap.filtered.tsv > sample.molecules.tsv
echo "Step 1 done in $(($(date +%s) - START))s"

# Step 2: Backbone
echo ""
echo "================================================================"
echo "  Step 2: Backbone"
echo "================================================================"
START=$(date +%s)
${PHYSLR_PY} backbone --prune-branches=10 --prune-bridges=10 --prune-junctions=200 sample.molecules.tsv > sample.backbone.path
echo "Step 2 done in $(($(date +%s) - START))s"

# Summary
echo ""
echo "================================================================"
echo "  Summary"
echo "================================================================"
echo "Backbone paths: $(wc -l < sample.backbone.path)"
echo "Total backbone nodes: $(awk '{print NF}' sample.backbone.path | paste -sd+ | bc)"
echo "Top 10 path sizes:"
awk '{print NF}' sample.backbone.path | sort -rn | head -10

# Molecule graph stats
BLANK=$(grep -n '^$' sample.molecules.tsv | head -1 | cut -d: -f1)
TOTAL=$(wc -l < sample.molecules.tsv)
VERTS=$((BLANK - 2))
EDGES=$((TOTAL - BLANK - 1))
echo ""
echo "Molecule graph: V=${VERTS} E=${EDGES}"

# Step 3: Visualization
if [ -n "${REFERENCE}" ] && [ -f "${REFERENCE}" ]; then
    echo ""
    echo "================================================================"
    echo "  Step 3: Visualization"
    echo "================================================================"

    if [ ! -f sample.ref.tsv ]; then
        ln -sf "${PREV_OUTDIR}/sample.ref.tsv" sample.ref.tsv
    fi

    START=$(date +%s)
    "${PHYSLR_RUST}" map-paf sample.backbone.path sample.filtered.tsv sample.ref.tsv \
        -o sample.backbone.paf -n 1
    echo "map-paf done in $(($(date +%s) - START))s"

    START=$(date +%s)
    ${CONDA_PYTHON} "${SCRIPT_DIR}/plotpaf.py" sample.backbone.paf sample.backbone 1
    echo "Plotting done in $(($(date +%s) - START))s"
fi

echo ""
echo "Job complete: $(date)"
