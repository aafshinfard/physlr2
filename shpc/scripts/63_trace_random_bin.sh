#!/bin/bash
#SBATCH --job-name=trace-randbin
#SBATCH --output=./logs/trace_randbin.%j.out
#SBATCH --error=./logs/trace_randbin.%j.err
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
OVERLAP="${BASEDIR}/output/na12878_stlfr/sample.overlap.filtered.tsv"
PHYSLR="${BASEDIR}/physlr"
BARCODES="195_574_231,602_1381_479,165_210_584,669_506_651,126_292_1400"

echo "============================================================"
echo "Trace with random binning fix"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "============================================================"

echo "=== Rust trace (random binning) ==="
${PHYSLR} -t1 trace-molecules ${OVERLAP} \
    --barcodes "${BARCODES}" \
    --strategy "distributed+sqcosbin" 2>&1

echo ""
echo "=== Rust trace (distributed only) ==="
${PHYSLR} -t1 trace-molecules ${OVERLAP} \
    --barcodes "${BARCODES}" \
    --strategy "distributed" 2>&1

echo ""
echo "=== Reference: Original Python results ==="
echo "Barcode 195_574_231 (deg=1071): Python distributed → 26 communities, 222 nodes, final=28 mols"
echo "Barcode 602_1381_479 (deg=1017): Python distributed → 10 communities, 40 nodes, final=10 mols"
echo "Barcode 165_210_584 (deg=962): Python distributed → 14 communities, 65 nodes, final=14 mols"
echo "Barcode 669_506_651 (deg=870): Python distributed → 24 communities, 161 nodes, final=24 mols"
echo "Barcode 126_292_1400 (deg=858): Python distributed → 9 communities, 34 nodes, final=9 mols"

echo ""
echo "Job complete: $(date)"
