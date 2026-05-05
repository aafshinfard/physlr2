#!/bin/bash
#SBATCH --job-name=physlr-na12878-stlfr
#SBATCH --output=./logs/na12878_stlfr_profile.%j.out
#SBATCH --error=./logs/na12878_stlfr_profile.%j.err
#SBATCH --mem=500G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
export PHYSLR="${BASEDIR}/physlr-next/target/release/physlr"
export PATH="${HOME}/.cargo/bin:${PATH}"

# Activate conda env with ntcard, nthits, indexlr
source activate physlr-tools 2>/dev/null || conda activate physlr-tools 2>/dev/null || true
export MAKEBF="${BASEDIR}/physlr-makebf"
echo "ntcard: $(which ntcard 2>/dev/null || echo 'NOT FOUND')"
echo "nthits: $(which nthits 2>/dev/null || echo 'NOT FOUND')"
echo "indexlr: $(which indexlr 2>/dev/null || echo 'NOT FOUND')"
echo "physlr-makebf: ${MAKEBF}"

mkdir -p "${BASEDIR}/logs"

echo "============================================================"
echo "Physlr NA12878 stLFR — Profiled Physical Map"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}M"
echo "Start: $(date)"
echo "============================================================"

# Space-separated list of input files — profile_pipeline.sh word-splits this
INPUT_FQ="${BASEDIR}/data/na12878_stlfr/na12878.stlfr.R1.fq.gz ${BASEDIR}/data/na12878_stlfr/na12878.stlfr.R2.fq.gz"
OUTDIR="${BASEDIR}/output/na12878_stlfr"

# Reference genome for backbone visualization
export REFERENCE="${BASEDIR}/data/reference/grch38.fa"

"${BASEDIR}/scripts/profile_pipeline.sh" \
    "${INPUT_FQ}" \
    "${OUTDIR}" \
    40 32 85 bc+cc builtin

echo ""
echo "Job complete: $(date)"
