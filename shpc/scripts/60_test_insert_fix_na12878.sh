#!/bin/bash
#SBATCH --job-name=fix-insert-na12878
#SBATCH --output=./logs/fix_insert_na12878.%j.out
#SBATCH --error=./logs/fix_insert_na12878.%j.err
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
PREV_OUTDIR="${BASEDIR}/output/na12878_stlfr"
OUTDIR="${BASEDIR}/output/na12878_stlfr_insert_fix"
PHYSLR="${BASEDIR}/physlr"

echo "============================================================"
echo "physlr-next with insert fix (or_insert → insert) — NA12878"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "============================================================"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Symlink the filtered overlap graph from previous run
ln -sf "${PREV_OUTDIR}/sample.overlap.filtered.tsv" sample.overlap.filtered.tsv 2>/dev/null || true

echo "=== Step 6: Molecules (distributed+sqcosbin) ==="
START=$(date +%s)
${PHYSLR} -v2 -t16 molecules \
    --strategy distributed+sqcosbin \
    --sqcos-threshold 0.75 \
    --skip-small 10 \
    --bin-max-size 50 \
    --merge-cutoff 20 \
    sample.overlap.filtered.tsv \
    -o sample.molecules.tsv 2>&1 | tee molecules.log
echo "Done in $(($(date +%s) - START))s"

# Parse molecule graph stats
BLANK=$(grep -n '^$' sample.molecules.tsv | head -1 | cut -d: -f1)
TOTAL=$(wc -l < sample.molecules.tsv)
VERTS=$((BLANK - 2))
EDGES=$((TOTAL - BLANK - 1))
echo "Molecule graph: V=${VERTS} E=${EDGES}"

echo ""
echo "=== Step 7: Backbone ==="
START=$(date +%s)
${PHYSLR} -v2 backbone \
    --prune-branches 10 \
    --prune-bridges 10 \
    --prune-junctions 200 \
    --min-component-size 50 \
    sample.molecules.tsv \
    -o sample.backbone.path 2>&1 | tee backbone.log
echo "Done in $(($(date +%s) - START))s"
echo "Backbone paths: $(wc -l < sample.backbone.path)"
echo "Top 10 path sizes:"
awk '{print NF}' sample.backbone.path | sort -rn | head -10

echo ""
echo "=== Comparison ==="
echo "Original (cpython): V=4,282,315 E=47,923,428 Paths=183"
echo "physlr-next (old or_insert): V=8,071,572 E=51,563,282 Paths=1,294"
echo "physlr-next (insert fix): V=${VERTS} E=${EDGES} Paths=$(wc -l < sample.backbone.path)"

echo ""
echo "Job complete: $(date)"
