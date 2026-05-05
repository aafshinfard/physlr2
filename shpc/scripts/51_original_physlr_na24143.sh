#!/bin/bash
#SBATCH --job-name=orig-physlr-na24143
#SBATCH --output=./logs/orig_physlr_na24143.%j.out
#SBATCH --error=./logs/orig_physlr_na24143.%j.err
#SBATCH --mem=300G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
SRCDIR="${BASEDIR}/physlr-src"
OUTDIR="${BASEDIR}/output/na24143_stlfr_original"
RUST_OUTDIR="${BASEDIR}/output/na24143_stlfr"

# Activate physlr-original conda env (has pypy3 + physlr python + ntcard + nthits)
source activate physlr-original 2>/dev/null || conda activate physlr-original 2>/dev/null || true

# Also need indexlr from physlr-tools — add it to PATH
export PATH="${BASEDIR}/physlr-src:${HOME}/miniconda3/envs/physlr-tools/bin:${PATH}"

# Reference genome for visualization
export REFERENCE="${BASEDIR}/data/reference/grch38.fa"

echo "============================================================"
echo "Original Physlr — NA24143 stLFR"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}M"
echo "Start: $(date)"
echo "============================================================"

# Physlr Python location
PHYSLR_DIR="/home/afshinfa/miniconda3/envs/physlr-original/bin/share/physlr-1.0.4-8"
export PYTHONPATH="${PHYSLR_DIR}"
PHYSLR_PY="pypy3 ${PHYSLR_DIR}/bin/physlr"

# Verify tools
echo "physlr: $(${PHYSLR_PY} version 2>&1 || true)"
echo "physlr-overlap: $(which physlr-overlap)"
echo "physlr-filter-bxmx: $(which physlr-filter-bxmx)"
echo "indexlr: $(which indexlr)"
echo "ntcard: $(which ntcard)"
echo "pypy3: $(pypy3 --version 2>&1 | head -1)"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

K=40
W=32
M=85  # stLFR default overlap percentile
STRATEGY="distributed+sqcosbin"  # original default
THREADS=16

# ── Reuse indexlr output from physlr-next run ────────────────────────────────
READS_TSV="${RUST_OUTDIR}/sample.reads.tsv"
if [ ! -f "${READS_TSV}" ]; then
    echo "ERROR: Cannot find ${READS_TSV} — run physlr-next pipeline first"
    exit 1
fi
echo "Reusing indexlr output: ${READS_TSV} ($(wc -l < "${READS_TSV}") lines)"

# ── Step 5: Filter minimizers (physlr-filter-bxmx) ──────────────────────────
echo ""
echo "================================================================"
echo "  Step 5: Filter minimizers (physlr-filter-bxmx)"
echo "================================================================"
START=$(date +%s)
physlr-filter-bxmx -n100 -N5000 -o sample.filtered.tsv "${READS_TSV}"
echo "Step 5 done in $(($(date +%s) - START))s"
echo "Filtered barcodes: $(wc -l < sample.filtered.tsv)"

# ── Step 6: Compute overlaps (physlr-overlap) ───────────────────────────────
echo ""
echo "================================================================"
echo "  Step 6: Compute overlaps (physlr-overlap)"
echo "================================================================"
START=$(date +%s)
physlr-overlap -t${THREADS} -m10 sample.filtered.tsv > sample.overlap.tsv
echo "Step 6 done in $(($(date +%s) - START))s"

# ── Step 7: Filter overlaps (physlr filter-overlap) ─────────────────────────
echo ""
echo "================================================================"
echo "  Step 7: Filter overlaps (physlr filter-overlap)"
echo "================================================================"
START=$(date +%s)
${PHYSLR_PY} filter-overlap --minimizer-overlap ${M} sample.overlap.tsv > sample.overlap.filtered.tsv
echo "Step 7 done in $(($(date +%s) - START))s"

# ── Step 8: Separate molecules (physlr molecules) ───────────────────────────
echo ""
echo "================================================================"
echo "  Step 8: Separate molecules (strategy=${STRATEGY})"
echo "================================================================"
START=$(date +%s)
${PHYSLR_PY} molecules -t4 --separation-strategy=${STRATEGY} sample.overlap.filtered.tsv > sample.molecules.tsv
echo "Step 8 done in $(($(date +%s) - START))s"

# ── Step 9: Extract backbone (physlr backbone) ──────────────────────────────
echo ""
echo "================================================================"
echo "  Step 9: Extract backbone"
echo "================================================================"
START=$(date +%s)
${PHYSLR_PY} backbone --prune-branches=10 --prune-bridges=10 --prune-junctions=200 sample.molecules.tsv > sample.backbone.path
echo "Step 9 done in $(($(date +%s) - START))s"

# ── Summary ─────────────────────────────────────────────────────────────────
echo ""
echo "================================================================"
echo "  Summary"
echo "================================================================"
echo "Backbone paths: $(wc -l < sample.backbone.path)"
echo "Total backbone nodes: $(awk '{print NF}' sample.backbone.path | paste -sd+ | bc)"
echo "Top 10 path sizes:"
awk '{print NF}' sample.backbone.path | sort -rn | head -10
echo ""

# ── Step 10: Visualization (if reference available) ─────────────────────────
if [ -n "${REFERENCE}" ] && [ -f "${REFERENCE}" ]; then
    echo "================================================================"
    echo "  Step 10: Backbone-to-reference visualization"
    echo "================================================================"

    PHYSLR_RUST="${BASEDIR}/physlr-next/target/release/physlr"
    SCRIPT_DIR="${BASEDIR}/physlr-next/scripts"

    # Index reference
    if [ ! -f sample.ref.tsv ]; then
        START=$(date +%s)
        indexlr -t5 -k ${K} -w ${W} -o sample.ref.tsv "${REFERENCE}"
        echo "Reference indexing done in $(($(date +%s) - START))s"
    fi

    # Map backbone to reference using physlr-next's map-paf
    START=$(date +%s)
    "${PHYSLR_RUST}" map-paf sample.backbone.path sample.filtered.tsv sample.ref.tsv \
        -o sample.backbone.paf -n 1
    echo "map-paf done in $(($(date +%s) - START))s"

    # Generate plots
    START=$(date +%s)
    python3 "${SCRIPT_DIR}/plotpaf.py" sample.backbone.paf sample.backbone 1
    echo "Plotting done in $(($(date +%s) - START))s"
fi

echo ""
echo "Job complete: $(date)"
