#!/bin/bash
#SBATCH --job-name=trace-dist-detail
#SBATCH --output=./logs/trace_dist_detail.%j.out
#SBATCH --error=./logs/trace_dist_detail.%j.err
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
OVERLAP="${BASEDIR}/output/na12878_stlfr/sample.overlap.filtered.tsv"
CONDA_PYTHON="/home/afshinfa/miniconda3/envs/physlr-cpython/bin/python3"
PHYSLR_DIR="/home/afshinfa/miniconda3/envs/physlr-original/bin/share/physlr-1.0.4-8"
TRACE_PY="${BASEDIR}/physlr-next/shpc/scripts/64_trace_distributed_detail.py"
BARCODES="195_574_231,602_1381_479,165_210_584,669_506_651,126_292_1400"

echo "============================================================"
echo "Detailed distributed pipeline trace (Python)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "============================================================"

export PYTHONPATH="${PHYSLR_DIR}"
${CONDA_PYTHON} ${TRACE_PY} ${OVERLAP} "${BARCODES}" 2>&1

echo ""
echo "Job complete: $(date)"
