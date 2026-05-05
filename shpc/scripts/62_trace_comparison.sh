#!/bin/bash
#SBATCH --job-name=trace-compare
#SBATCH --output=./logs/trace_comparison.%j.out
#SBATCH --error=./logs/trace_comparison.%j.err
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
OVERLAP="${BASEDIR}/output/na12878_stlfr/sample.overlap.filtered.tsv"
PHYSLR="${BASEDIR}/physlr"
CONDA_PYTHON="/home/afshinfa/miniconda3/envs/physlr-cpython/bin/python3"
PHYSLR_DIR="/home/afshinfa/miniconda3/envs/physlr-original/bin/share/physlr-1.0.4-8"
TRACE_PY="${BASEDIR}/physlr-next/shpc/scripts/61_trace_single_barcode.py"

echo "============================================================"
echo "Trace comparison: Original vs Rust molecule separation"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "============================================================"

# First, find the top 5 highest-degree barcodes
echo "=== Finding top barcodes ==="
TOP_BARCODES=$(${PHYSLR} -t1 trace-molecules ${OVERLAP} --barcodes top5 2>&1 | grep "^Barcode:" | head -5 | awk '{print $2}' | sed 's/(degree=.*//' | paste -sd,)
echo "Top barcodes: ${TOP_BARCODES}"

# If that didn't work, use a simpler approach
if [ -z "${TOP_BARCODES}" ]; then
    echo "Falling back to manual barcode selection..."
    # Get barcodes with highest degree from the TSV
    TOP_BARCODES=$(awk 'NR>1 && !/^\t/ && !/^$/ {print $1}' ${OVERLAP} | head -100000 | while read bc; do
        echo "$bc $(grep -c "^${bc}\t\|	${bc}\t" ${OVERLAP} 2>/dev/null || echo 0)"
    done | sort -k2 -rn | head -5 | awk '{print $1}' | paste -sd,)
    echo "Fallback barcodes: ${TOP_BARCODES}"
fi

echo ""
echo "============================================================"
echo "=== Rust trace ==="
echo "============================================================"
${PHYSLR} -t1 trace-molecules ${OVERLAP} \
    --barcodes "${TOP_BARCODES}" \
    --strategy "distributed+sqcosbin" 2>&1

echo ""
echo "============================================================"
echo "=== Python (original) trace ==="
echo "============================================================"
export PYTHONPATH="${PHYSLR_DIR}"
${CONDA_PYTHON} ${TRACE_PY} ${OVERLAP} "${TOP_BARCODES}" "distributed+sqcosbin" 2>&1

echo ""
echo "============================================================"
echo "=== Also test distributed only ==="
echo "============================================================"
echo "--- Rust distributed only ---"
${PHYSLR} -t1 trace-molecules ${OVERLAP} \
    --barcodes "${TOP_BARCODES}" \
    --strategy "distributed" 2>&1

echo ""
echo "--- Python distributed only ---"
${CONDA_PYTHON} ${TRACE_PY} ${OVERLAP} "${TOP_BARCODES}" "distributed" 2>&1

echo ""
echo "Job complete: $(date)"
