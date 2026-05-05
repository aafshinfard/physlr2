#!/bin/bash
#SBATCH --job-name=physlr-na24143-stlfr
#SBATCH --output=./logs/na24143_stlfr_profile.%j.out
#SBATCH --error=./logs/na24143_stlfr_profile.%j.err
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
echo "Physlr NA24143 stLFR — Profiled Physical Map"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}M"
echo "Start: $(date)"
echo "============================================================"

INPUT_FQ="${BASEDIR}/data/na24143_stlfr/na24143.stlfr.R1.fq.gz ${BASEDIR}/data/na24143_stlfr/na24143.stlfr.R2.fq.gz"
OUTDIR="${BASEDIR}/output/na24143_stlfr"

# Reference genome for backbone visualization
export REFERENCE="${BASEDIR}/data/reference/grch38.fa"

"${BASEDIR}/scripts/profile_pipeline.sh" \
    "${INPUT_FQ}" \
    "${OUTDIR}" \
    40 32 85 bc+cc builtin

echo ""
echo "Job complete: $(date)"
