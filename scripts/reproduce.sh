#!/usr/bin/env bash
# Reproduce the Physlr 2 physical map results for NA12878 and NA24143.
#
# Downloads stLFR linked reads and reference genome, runs the full
# physical map pipeline, and generates backbone-vs-reference plots.
#
# Usage:
#   bash scripts/reproduce.sh <workdir> [sample]
#
# Arguments:
#   workdir   Working directory for data and output (needs ~600 GB free)
#   sample    Optional: "na12878", "na24143", or "all" (default: all)
#
# The script is resumable — it skips completed steps and uses wget -c
# for partial downloads. Re-run it if interrupted.
#
# Requirements:
#   - physlr binary (cargo build --release)
#   - btllib (indexlr, ntcard, nthits): conda install -c bioconda btllib
#   - Python 3 with matplotlib (for plots)
#   - ~200 GB RAM, 16 CPU cores, 12-24 hours per sample

set -euo pipefail

# ── Arguments ────────────────────────────────────────────────────────────────

WORKDIR="${1:?Usage: $0 <workdir> [sample]}"
SAMPLE="${2:-all}"

# ── Locate physlr binary ────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [ -x "${REPO_DIR}/target/release/physlr" ]; then
    PHYSLR="${REPO_DIR}/target/release/physlr"
elif command -v physlr &>/dev/null; then
    PHYSLR="$(command -v physlr)"
else
    echo "ERROR: physlr binary not found."
    echo "Build it with: cd ${REPO_DIR} && cargo build --release"
    exit 1
fi
echo "Using physlr: ${PHYSLR}"

PLOTPAF="${REPO_DIR}/scripts/plotpaf.py"
if [ ! -f "${PLOTPAF}" ]; then
    echo "WARNING: plotpaf.py not found at ${PLOTPAF}, plots will be skipped"
fi

# ── Check external tool dependencies ────────────────────────────────────────

for tool in indexlr ntcard nthits; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "ERROR: ${tool} not found. Install btllib: conda install -c bioconda btllib"
        exit 1
    fi
done
echo "External tools: indexlr=$(command -v indexlr), ntcard=$(command -v ntcard), nthits=$(command -v nthits)"

# ── Parameters ───────────────────────────────────────────────────────────────

K=32
W=32
THREADS="${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 8)}"
THREADS=$((THREADS > 16 ? 16 : THREADS))

DATADIR="${WORKDIR}/data"
mkdir -p "${DATADIR}"

# ── Download helper ──────────────────────────────────────────────────────────

download() {
    local url="$1"
    local dest="$2"
    if [ -f "${dest}" ] && [ -s "${dest}" ]; then
        echo "  [skip] $(basename "${dest}") already exists"
        return 0
    fi
    echo "  [download] $(basename "${dest}")"
    mkdir -p "$(dirname "${dest}")"
    wget -c -q --show-progress -t 3 -T 120 "${url}" -O "${dest}"
}

# ── Step runner ──────────────────────────────────────────────────────────────

run_step() {
    local name="$1"
    local output="$2"
    shift 2
    if [ -f "${output}" ] && [ -s "${output}" ]; then
        echo "  [skip] ${name} ($(basename "${output}") exists)"
        return 0
    fi
    echo ""
    echo "  === ${name} ==="
    local start=$(date +%s)
    "$@"
    local elapsed=$(( $(date +%s) - start ))
    echo "  --- ${name}: ${elapsed}s ---"
}

# ── Download reference genome ────────────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Downloading reference genome"
echo "============================================================"

download \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" \
    "${DATADIR}/grch38.fa.gz"

if [ -f "${DATADIR}/grch38.fa.gz" ] && [ ! -f "${DATADIR}/grch38.fa" ]; then
    echo "  Decompressing reference..."
    gunzip -k "${DATADIR}/grch38.fa.gz"
fi

# ── Per-sample pipeline ──────────────────────────────────────────────────────

run_sample() {
    local sample="$1"
    local r1_url="$2"
    local r2_url="$3"

    echo ""
    echo "============================================================"
    echo "  ${sample}: Downloading reads"
    echo "============================================================"

    download "${r1_url}" "${DATADIR}/${sample}.R1.fq.gz"
    download "${r2_url}" "${DATADIR}/${sample}.R2.fq.gz"

    local OUTDIR="${WORKDIR}/output/${sample}"
    mkdir -p "${OUTDIR}"
    cd "${OUTDIR}"

    echo ""
    echo "============================================================"
    echo "  ${sample}: Building physical map"
    echo "============================================================"

    # Build physical map with reference mapping.
    # Uses external tools (ntcard → nthits → physlr-makebf → indexlr)
    # for repeat filtering and minimizer extraction.
    local REF_ARG=""
    if [ -f "${DATADIR}/grch38.fa" ]; then
        REF_ARG="--reference ${DATADIR}/grch38.fa"
    fi

    run_step "Physical map" "${OUTDIR}/${sample}.backbone.path" \
        "${PHYSLR}" physical-map \
        -k ${K} -w ${W} -t ${THREADS} \
        -o "${OUTDIR}" -p "${sample}" \
        ${REF_ARG} \
        "${DATADIR}/${sample}.R1.fq.gz" "${DATADIR}/${sample}.R2.fq.gz"

    # Merge adjacent backbone paths using bridge molecule evidence
    run_step "Split minimizers" "${OUTDIR}/${sample}.split.tsv" \
        "${PHYSLR}" split-minimizers \
        "${OUTDIR}/${sample}.mol.tsv" "${OUTDIR}/${sample}.filtered.tsv" \
        -o "${OUTDIR}/${sample}.split.tsv" -t ${THREADS}

    run_step "Merge paths" "${OUTDIR}/${sample}.merged.path" \
        "${PHYSLR}" merge-paths \
        "${OUTDIR}/${sample}.backbone.path" "${OUTDIR}/${sample}.split.tsv" \
        -o "${OUTDIR}/${sample}.merged.path"

    # Merged path metrics
    echo ""
    echo "  === Physical map metrics (after merge) ==="
    "${PHYSLR}" path-metrics "${OUTDIR}/${sample}.merged.path"

    # Re-map merged paths to reference for final visualization
    if [ -f "${DATADIR}/grch38.fa" ] && [ -f "${OUTDIR}/${sample}.ref.tsv" ]; then
        run_step "Map merged paths to reference" "${OUTDIR}/${sample}.merged.paf" \
            "${PHYSLR}" map-paf \
            "${OUTDIR}/${sample}.merged.path" "${OUTDIR}/${sample}.split.tsv" \
            "${OUTDIR}/${sample}.ref.tsv" \
            -o "${OUTDIR}/${sample}.merged.paf" -n 1 --mx-type split

        if [ -f "${PLOTPAF}" ]; then
            run_step "Generate merged plots" "${OUTDIR}/${sample}.merged.backbone.png" \
                python3 "${PLOTPAF}" "${OUTDIR}/${sample}.merged.paf" \
                "${OUTDIR}/${sample}.merged" \
                --positions "${OUTDIR}/${sample}.positions.tsv"
            echo "  Plots: ${OUTDIR}/${sample}.merged.backbone.png"
            echo "         ${OUTDIR}/${sample}.merged.reference.png"
        fi
    fi

    echo ""
    echo "============================================================"
    echo "  ${sample}: Complete"
    echo "  Output: ${OUTDIR}"
    echo "============================================================"
}

# ── Run selected samples ─────────────────────────────────────────────────────

if [ "${SAMPLE}" = "all" ] || [ "${SAMPLE}" = "na12878" ]; then
    run_sample "na12878" \
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/stLFR.split_read.1.fq.gz" \
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/stLFR.split_read.2.fq.gz"
fi

if [ "${SAMPLE}" = "all" ] || [ "${SAMPLE}" = "na24143" ]; then
    run_sample "na24143" \
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/stLFR/stLFR_NA24143_split_read.1.fq.gz" \
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/stLFR/stLFR_NA24143_split_read.2.fq.gz"
fi

echo ""
echo "============================================================"
echo "  All done."
echo "============================================================"
