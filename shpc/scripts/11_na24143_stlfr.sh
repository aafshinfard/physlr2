#!/bin/bash
#SBATCH --job-name=physlr-na24143-stlfr
#SBATCH --output=./logs/na24143_stlfr.%j.out
#SBATCH --error=./logs/na24143_stlfr.%j.err
#SBATCH --mem=450G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

# NA24143 stLFR: Build physical map and scaffold draft assemblies.
#
# Input: stLFR linked reads (barcode in read name after #)
# Drafts: ONT Shasta, PacBio Falcon, PE+MPET ABySS, stLFR Supernova
#
# Usage: sbatch scripts/11_na24143_stlfr.sh

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
source "${BASEDIR}/scripts/common.sh"
mkdir -p "${BASEDIR}/logs"

SAMPLE="na24143_stlfr"
OUTDIR="${BASEDIR}/output/${SAMPLE}"
INPUT_FQ="${DATADIR}/na24143_stlfr/na24143.stlfr.R1.fq.gz ${DATADIR}/na24143_stlfr/na24143.stlfr.R2.fq.gz"

# stLFR parameters
K=40
W=32
OVERLAP_PCT=85
MOL_STRATEGY="bc+cc"

echo "============================================================"
echo "Physlr: NA24143 stLFR"
echo "Start: $(date)"
echo "Threads: ${THREADS}"
echo "Memory: ${SLURM_MEM_PER_NODE:-450G}"
echo "============================================================"

# ─── Step 1: Build physical map ──────────────────────────────────────────────

run_physical_map "${SAMPLE}" "${INPUT_FQ}" "${OUTDIR}" \
    "${K}" "${W}" "${OVERLAP_PCT}" "${MOL_STRATEGY}"

# ─── Step 2: Scaffold draft assemblies ───────────────────────────────────────

DRAFTS=(
    "na24143.ont.baseline.fa:ont"
    "na24143.pacbio.baseline.fa:pacbio"
    "na24143.pempet.baseline.fa:pempet"
    "na24143.stlfr.baseline.fa:supernova"
)

for entry in "${DRAFTS[@]}"; do
    draft_file="${entry%%:*}"
    draft_label="${entry##*:}"
    draft_path="${ASMDIR}/${draft_file}"

    if [ ! -f "${draft_path}" ]; then
        echo "WARNING: Draft assembly not found: ${draft_path}, skipping ${draft_label}"
        continue
    fi

    scaffold_assembly "${SAMPLE}" "${draft_path}" "${draft_label}" \
        "${OUTDIR}" "${K}" "${W}"
done

# ─── Step 3: Evaluate with QUAST ─────────────────────────────────────────────

for entry in "${DRAFTS[@]}"; do
    draft_file="${entry%%:*}"
    draft_label="${entry##*:}"
    draft_path="${ASMDIR}/${draft_file}"

    if [ -f "${draft_path}" ]; then
        run_quast "${draft_path}" "${SAMPLE}.${draft_label}.baseline" "${OUTDIR}"
    fi
done

for entry in "${DRAFTS[@]}"; do
    draft_label="${entry##*:}"
    scaffold="${OUTDIR}/${SAMPLE}.${draft_label}.scaffold.fa"

    if [ -f "${scaffold}" ]; then
        run_quast "${scaffold}" "${SAMPLE}.${draft_label}.scaffold" "${OUTDIR}"
    fi
done

# ─── Step 4: Summary ─────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo "NA24143 stLFR: Complete"
echo "End: $(date)"
echo "Output: ${OUTDIR}"
echo "============================================================"

echo ""
echo "--- Baseline vs Scaffolded ---"
for entry in "${DRAFTS[@]}"; do
    draft_label="${entry##*:}"
    baseline_quast="${OUTDIR}/quast_${SAMPLE}.${draft_label}.baseline/transposed_report.tsv"
    scaffold_quast="${OUTDIR}/quast_${SAMPLE}.${draft_label}.scaffold/transposed_report.tsv"

    if [ -f "${baseline_quast}" ] && [ -f "${scaffold_quast}" ]; then
        echo ""
        echo "=== ${draft_label} ==="
        echo "Baseline:"
        grep -E "NG50|NGA50|# scaffolds" "${baseline_quast}" || true
        echo "Scaffolded:"
        grep -E "NG50|NGA50|# scaffolds" "${scaffold_quast}" || true
    fi
done
