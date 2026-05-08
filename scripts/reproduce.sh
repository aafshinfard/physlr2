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
#   - indexlr from btllib (conda install -c bioconda btllib)
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

# ── Locate indexlr ───────────────────────────────────────────────────────────

INDEXLR=""
if command -v indexlr &>/dev/null; then
    INDEXLR="$(command -v indexlr)"
    echo "Using indexlr: ${INDEXLR}"
else
    echo "WARNING: indexlr not found. Reference visualization will use built-in indexer."
    echo "  Install btllib for faster indexing: conda install -c bioconda btllib"
fi

# ── Parameters ───────────────────────────────────────────────────────────────

K=32
W=32
OVERLAP_PCT=85
MOL_STRATEGY="distributed+sqcosbin"
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
    echo "  ${sample}: Running pipeline"
    echo "============================================================"

    # Step 1: Index minimizers
    run_step "Index minimizers" "${sample}.reads.tsv" \
        "${PHYSLR}" index -k ${K} -w ${W} -t ${THREADS} \
        "${DATADIR}/${sample}.R1.fq.gz" "${DATADIR}/${sample}.R2.fq.gz" \
        -o "${sample}.reads.tsv"

    # Step 2: Filter barcodes
    run_step "Filter minimizers" "${sample}.filtered.tsv" \
        "${PHYSLR}" filter-minimizers "${sample}.reads.tsv" \
        -o "${sample}.filtered.tsv" -n 100 -N 5000

    # Step 3: Compute overlaps
    run_step "Compute overlaps" "${sample}.overlap.tsv" \
        "${PHYSLR}" overlap "${sample}.filtered.tsv" \
        -o "${sample}.overlap.tsv" -t ${THREADS}

    # Step 4: Filter edges
    run_step "Filter overlaps" "${sample}.overlap.filtered.tsv" \
        "${PHYSLR}" filter-overlap "${sample}.overlap.tsv" \
        -o "${sample}.overlap.filtered.tsv" -p ${OVERLAP_PCT}

    # Step 5: Separate molecules
    run_step "Separate molecules" "${sample}.molecules.tsv" \
        "${PHYSLR}" molecules "${sample}.overlap.filtered.tsv" \
        -o "${sample}.molecules.tsv" \
        --strategy "${MOL_STRATEGY}" -t ${THREADS}

    # Step 6: Extract backbone
    run_step "Extract backbone" "${sample}.backbone.path" \
        "${PHYSLR}" backbone "${sample}.molecules.tsv" \
        -o "${sample}.backbone.path" \
        --prune-branches 10 --prune-bridges 10 --prune-junctions 200

    # Step 7: Split minimizers
    run_step "Split minimizers" "${sample}.split.tsv" \
        "${PHYSLR}" split-minimizers "${sample}.molecules.tsv" "${sample}.filtered.tsv" \
        -o "${sample}.split.tsv" -t ${THREADS}

    # Step 8: Merge paths
    run_step "Merge paths" "${sample}.merged.path" \
        "${PHYSLR}" merge-paths "${sample}.backbone.path" "${sample}.split.tsv" \
        -o "${sample}.merged.path"

    # Step 9: Physical map metrics
    echo ""
    echo "  === Physical map metrics ==="
    "${PHYSLR}" path-metrics "${sample}.merged.path"

    # Step 10: Visualize against reference
    if [ -f "${DATADIR}/grch38.fa" ]; then
        echo ""
        echo "  === Reference visualization ==="

        # Index reference
        if [ ! -f "${sample}.ref.tsv" ]; then
            if [ -n "${INDEXLR}" ]; then
                run_step "Index reference (indexlr)" "${sample}.ref.tsv" \
                    "${INDEXLR}" -t 5 -k ${K} -w ${W} \
                    -o "${sample}.ref.tsv" "${DATADIR}/grch38.fa"
            else
                run_step "Index reference (built-in)" "${sample}.ref.tsv" \
                    "${PHYSLR}" index-contigs "${DATADIR}/grch38.fa" \
                    -o "${sample}.ref.tsv" -k ${K} -w ${W}
            fi
        fi

        # Map backbone to reference
        run_step "Map to reference" "${sample}.paf" \
            "${PHYSLR}" map-paf \
            "${sample}.merged.path" "${sample}.split.tsv" "${sample}.ref.tsv" \
            -o "${sample}.paf" -n 1 --mx-type split

        # Generate plots
        if [ -f "${PLOTPAF}" ]; then
            run_step "Generate plots" "${sample}.backbone.png" \
                python3 "${PLOTPAF}" "${sample}.paf" "${sample}" 1
            echo "  Plots: ${OUTDIR}/${sample}.backbone.png"
            echo "         ${OUTDIR}/${sample}.reference.png"
        fi
    else
        echo "  [skip] Reference not available, skipping visualization"
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
