#!/bin/bash
# Self-profiling wrapper for the Physlr pipeline.
#
# Runs the full physical map pipeline while recording wall-clock time
# and peak RSS (resident set size) for each step. Produces:
#   <outdir>/profile.tsv   — tab-separated profiling data
#   <outdir>/profile.png   — bar chart of time and RAM per step
#
# Usage:
#   ./scripts/profile_pipeline.sh <input_fq> <outdir> [k] [w] [overlap_pct] [strategy]
#
# Example:
#   ./scripts/profile_pipeline.sh reads.fq.gz output/ 40 32 85 bc+cc
#
# Requires: /usr/bin/time (GNU time, not shell builtin), gnuplot (for plotting)

set -euo pipefail

# ── Arguments ────────────────────────────────────────────────────────────────

INPUT_FQ="${1:?Usage: $0 <input_fq> <outdir> [k] [w] [overlap_pct] [strategy] [indexer]}"
OUTDIR="${2:?Usage: $0 <input_fq> <outdir> [k] [w] [overlap_pct] [strategy] [indexer]}"
K="${3:-40}"
W="${4:-32}"
OVERLAP_PCT="${5:-85}"
STRATEGY="${6:-bc+cc}"
INDEXER="${7:-btllib}"

# ── Locate physlr binary ────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Try to find physlr relative to this script (inside the repo)
if [ -x "${SCRIPT_DIR}/../target/release/physlr" ]; then
    PHYSLR="${SCRIPT_DIR}/../target/release/physlr"
elif [ -n "${PHYSLR:-}" ] && [ -x "${PHYSLR}" ]; then
    : # use PHYSLR from environment
else
    echo "ERROR: Cannot find physlr binary. Set PHYSLR env var or build with cargo build --release"
    exit 1
fi
PHYSLR="$(realpath "${PHYSLR}")"
echo "Using physlr: ${PHYSLR}"

# ── GNU time ─────────────────────────────────────────────────────────────────

GNU_TIME=""
if command -v /usr/bin/time &>/dev/null; then
    GNU_TIME="/usr/bin/time"
elif command -v gtime &>/dev/null; then
    GNU_TIME="gtime"
else
    echo "WARNING: GNU time not found. RAM tracking will be unavailable."
fi

# ── Setup ────────────────────────────────────────────────────────────────────

mkdir -p "${OUTDIR}"
PROFILE_TSV="${OUTDIR}/profile.tsv"
SAMPLE="sample"
TIME_TMP="${OUTDIR}/.time_tmp"

echo -e "step\twall_sec\tpeak_rss_mb\texit_code" > "${PROFILE_TSV}"

# ── Profiled step runner ─────────────────────────────────────────────────────

# Run a command, record wall time and peak RSS.
# Usage: profile_step "Step Name" cmd arg1 arg2 ...
profile_step() {
    local step_name="$1"
    shift

    echo ""
    echo "================================================================"
    echo "  ${step_name}"
    echo "  Command: $*"
    echo "================================================================"

    local start_time=$(date +%s)
    local exit_code=0
    local peak_rss_kb=0

    if [ -n "${GNU_TIME}" ]; then
        # GNU time -v outputs "Maximum resident set size (kbytes): NNNN"
        ${GNU_TIME} -v "$@" 2> "${TIME_TMP}" || exit_code=$?
        peak_rss_kb=$(grep "Maximum resident set size" "${TIME_TMP}" | awk '{print $NF}') || peak_rss_kb=0
        # Print stderr from time (contains useful info)
        cat "${TIME_TMP}" >&2
    else
        "$@" || exit_code=$?
    fi

    local end_time=$(date +%s)
    local wall_sec=$((end_time - start_time))
    local peak_rss_mb=$(awk "BEGIN {printf \"%.1f\", ${peak_rss_kb}/1024}")

    echo "--- ${step_name}: ${wall_sec}s, peak RSS: ${peak_rss_mb} MB ---"
    echo -e "${step_name}\t${wall_sec}\t${peak_rss_mb}\t${exit_code}" >> "${PROFILE_TSV}"

    if [ "${exit_code}" -ne 0 ]; then
        echo "ERROR: ${step_name} failed with exit code ${exit_code}"
        exit "${exit_code}"
    fi
}

# ── Pipeline ─────────────────────────────────────────────────────────────────

echo "============================================================"
echo "Physlr Profiled Pipeline"
echo "Start: $(date)"
echo "Input: ${INPUT_FQ}"
echo "Output: ${OUTDIR}"
echo "Parameters: k=${K} w=${W} overlap_pct=${OVERLAP_PCT} strategy=${STRATEGY} indexer=${INDEXER}"
echo "============================================================"

cd "${OUTDIR}"

# Number of threads (from SLURM or default 32)
THREADS="${SLURM_CPUS_PER_TASK:-32}"

# Bloom filter size in bytes (default 10GB, matching original)
BF_SIZE="${BF_SIZE:-10000000000}"

# ── Step 0a: Count k-mers with ntcard ────────────────────────────────────────
# ntcard produces a histogram of k-mer frequencies using streaming HyperLogLog.
# Multi-threaded, single pass, constant memory.
profile_step "0a_ntcard" \
    ntcard -t "${THREADS}" -k "${K}" -o "${SAMPLE}_k${K}.histogram" ${INPUT_FQ}

# ── Step 0b: Find histogram mode ─────────────────────────────────────────────
# Identical to: python physlr find-ntcard-mode <histogram>
MODE=$(python3 "${SCRIPT_DIR}/find-ntcard-mode.py" "${SAMPLE}_k${K}.histogram")
REPEAT_THRESHOLD=$((MODE * 3))
echo "K-mer histogram mode: ${MODE}"
echo "Repeat threshold: ${REPEAT_THRESHOLD} (mode=${MODE} x 3)"

# ── Step 0c: Find repetitive k-mers with nthits ─────────────────────────────
# nthits outputs a text file of k-mers above the threshold.
# Multi-threaded, uses ntHash rolling hash.
profile_step "0c_nthits" \
    nthits -t "${THREADS}" -k "${K}" -c "${REPEAT_THRESHOLD}" \
    -p "${SAMPLE}" ${INPUT_FQ}

# ── Step 0d: Build Bloom filter from repetitive k-mers ──────────────────────
# physlr-makebf converts the nthits text output to a btllib Bloom filter.
# MAKEBF env var must point to the physlr-makebf binary.
if [ -z "${MAKEBF}" ]; then
    echo "ERROR: MAKEBF env var not set. Set it to the path of physlr-makebf."
    exit 1
fi
profile_step "0d_makebf" \
    "${MAKEBF}" -k "${K}" -b "${BF_SIZE}" -t "${THREADS}" \
    -o "${SAMPLE}.k${K}.bf" \
    "${SAMPLE}_k${K}.rep"

# ── Step 1: Index minimizers with indexlr ────────────────────────────────────
# btllib indexlr extracts (k,w)-minimizers, filtering repetitive k-mers via BF.
# Multi-threaded (max 5 threads per input file), uses ntHash rolling hash.
# Accepts multiple input files directly.
profile_step "1_index" \
    indexlr --bx -t5 -k "${K}" -w "${W}" \
    -r "${SAMPLE}.k${K}.bf" \
    -o "${SAMPLE}.reads.tsv" \
    ${INPUT_FQ}

# Step 2: Filter minimizers
profile_step "2_filter_minimizers" \
    "${PHYSLR}" filter-minimizers "${SAMPLE}.reads.tsv" \
    -o "${SAMPLE}.filtered.tsv" \
    --min-count 100 --max-count 5000

# Step 3: Compute overlaps
profile_step "3_overlap" \
    "${PHYSLR}" overlap "${SAMPLE}.filtered.tsv" \
    -o "${SAMPLE}.overlap.tsv"

# Step 4: Filter overlaps
profile_step "4_filter_overlap" \
    "${PHYSLR}" filter-overlap "${SAMPLE}.overlap.tsv" \
    -o "${SAMPLE}.overlap.filtered.tsv" \
    -p "${OVERLAP_PCT}"

# Step 5: Separate molecules
profile_step "5_molecules" \
    "${PHYSLR}" molecules "${SAMPLE}.overlap.filtered.tsv" \
    -o "${SAMPLE}.molecules.tsv" \
    --strategy "${STRATEGY}"

# Step 6: Extract backbone
profile_step "6_backbone" \
    "${PHYSLR}" backbone "${SAMPLE}.molecules.tsv" \
    -o "${SAMPLE}.backbone.tsv" \
    --prune-branches 10 --prune-bridges 10

# Step 7: Metrics
profile_step "7_metrics" \
    "${PHYSLR}" path-metrics "${SAMPLE}.backbone.tsv"

# Step 7b: Split minimizers (molecule-level minimizer assignment)
profile_step "7b_split_minimizers" \
    "${PHYSLR}" split-minimizers \
    "${SAMPLE}.molecules.tsv" \
    "${SAMPLE}.filtered.tsv" \
    -o "${SAMPLE}.split.tsv"

# ── Step 8-10: Backbone-vs-reference visualization (optional) ────────────────
# Requires: REFERENCE env var pointing to a reference FASTA file.
# If not set, visualization steps are skipped.
if [ -n "${REFERENCE}" ] && [ -f "${REFERENCE}" ]; then
    echo ""
    echo "=== Visualization: mapping backbone to reference ==="

    # Step 8: Index reference with indexlr
    profile_step "8_index_ref" \
        indexlr -t5 -k "${K}" -w "${W}" \
        -o "${SAMPLE}.ref.tsv" \
        "${REFERENCE}"

    # Step 9a: Map backbone to reference using barcode minimizers
    profile_step "9a_map_paf" \
        "${PHYSLR}" map-paf \
        "${SAMPLE}.backbone.tsv" \
        "${SAMPLE}.filtered.tsv" \
        "${SAMPLE}.ref.tsv" \
        -o "${SAMPLE}.backbone.paf" \
        -n 1

    # Step 9b: Map backbone to reference using split minimizers
    profile_step "9b_map_paf_split" \
        "${PHYSLR}" map-paf \
        "${SAMPLE}.backbone.tsv" \
        "${SAMPLE}.split.tsv" \
        "${SAMPLE}.ref.tsv" \
        -o "${SAMPLE}.backbone.split.paf" \
        -n 1 \
        --mx-type split

    # Step 10a: Generate plots (barcode minimizers)
    profile_step "10a_plot" \
        python3 "${SCRIPT_DIR}/plotpaf.py" \
        "${SAMPLE}.backbone.paf" \
        "${SAMPLE}.backbone" \
        1

    # Step 10b: Generate plots (split minimizers)
    profile_step "10b_plot_split" \
        python3 "${SCRIPT_DIR}/plotpaf.py" \
        "${SAMPLE}.backbone.split.paf" \
        "${SAMPLE}.backbone.split" \
        1
else
    echo ""
    echo "Skipping visualization (REFERENCE env var not set or file not found)."
    echo "Set REFERENCE=/path/to/reference.fa to enable."
fi

# ── Summary ──────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo "Pipeline complete: $(date)"
echo "============================================================"
echo ""
echo "Profile data: ${PROFILE_TSV}"
cat "${PROFILE_TSV}"

# ── Generate plot ────────────────────────────────────────────────────────────

generate_plot() {
    local tsv="$1"
    local png="${tsv%.tsv}.png"

    # Check for gnuplot
    if ! command -v gnuplot &>/dev/null; then
        echo ""
        echo "WARNING: gnuplot not found. Skipping plot generation."
        echo "Install gnuplot and re-run, or use the TSV data directly."
        # Generate a simple ASCII bar chart as fallback
        generate_ascii_chart "$tsv"
        return
    fi

    local gp_file="${OUTDIR}/.profile_plot.gp"
    cat > "${gp_file}" << 'GNUPLOT_EOF'
set terminal pngcairo size 1200,600 enhanced font "Arial,11"
set output outfile

set style data histogram
set style histogram clustered gap 1
set style fill solid 0.8 border -1
set boxwidth 0.9

set xlabel "Pipeline Step" font ",12"
set xtics rotate by -30 font ",10"
set grid y

set multiplot layout 1,2 title "Physlr Pipeline Profile" font ",14"

# Left panel: Wall time
set ylabel "Wall Time (seconds)" font ",12"
set title "Wall Clock Time" font ",12"
plot datafile using 2:xtic(1) title "Time (s)" lc rgb "#4472C4" notitle

# Right panel: Peak RSS
set ylabel "Peak RSS (MB)" font ",12"
set title "Peak Memory (RSS)" font ",12"
plot datafile using 3:xtic(1) title "RSS (MB)" lc rgb "#ED7D31" notitle

unset multiplot
GNUPLOT_EOF

    gnuplot -e "datafile='${tsv}'; outfile='${png}'" "${gp_file}" 2>/dev/null

    if [ -f "${png}" ]; then
        echo "Plot saved: ${png}"
    else
        echo "WARNING: gnuplot failed. Falling back to ASCII chart."
        generate_ascii_chart "$tsv"
    fi

    rm -f "${gp_file}"
}

generate_ascii_chart() {
    local tsv="$1"
    echo ""
    echo "═══════════════════════════════════════════════════════════════"
    echo "  Pipeline Profile Summary"
    echo "═══════════════════════════════════════════════════════════════"
    echo ""

    # Find max values for scaling
    local max_time=$(awk -F'\t' 'NR>1 {if($2>m)m=$2} END{print m+0}' "$tsv")
    local max_rss=$(awk -F'\t' 'NR>1 {if($3>m)m=$3} END{print m+0}' "$tsv")
    local bar_width=40

    printf "  %-25s %10s %12s\n" "Step" "Time (s)" "Peak RSS (MB)"
    printf "  %-25s %10s %12s\n" "-------------------------" "----------" "-------------"

    local total_time=0
    local max_rss_overall=0

    while IFS=$'\t' read -r step wall rss exit_code; do
        [ "$step" = "step" ] && continue  # skip header

        total_time=$((total_time + wall))
        if awk "BEGIN{exit ($rss > $max_rss_overall) ? 0 : 1}" 2>/dev/null; then
            max_rss_overall="$rss"
        fi

        # Time bar
        if [ "$max_time" -gt 0 ] 2>/dev/null; then
            local time_bar_len=$(awk "BEGIN{printf \"%d\", ($wall/$max_time)*$bar_width}")
        else
            local time_bar_len=0
        fi
        local time_bar=$(printf '%*s' "$time_bar_len" '' | tr ' ' '█')

        # RSS bar
        if awk "BEGIN{exit ($max_rss > 0) ? 0 : 1}" 2>/dev/null; then
            local rss_bar_len=$(awk "BEGIN{printf \"%d\", ($rss/$max_rss)*$bar_width}")
        else
            local rss_bar_len=0
        fi
        local rss_bar=$(printf '%*s' "$rss_bar_len" '' | tr ' ' '▓')

        printf "  %-25s %8ss  %10s MB\n" "$step" "$wall" "$rss"
        printf "    Time: %-${bar_width}s\n" "$time_bar"
        printf "    RSS:  %-${bar_width}s\n" "$rss_bar"
    done < "$tsv"

    echo ""
    printf "  %-25s %8ss  %10s MB\n" "TOTAL" "$total_time" "$max_rss_overall"
    echo "═══════════════════════════════════════════════════════════════"
}

generate_plot "${PROFILE_TSV}"

# Cleanup
rm -f "${TIME_TMP}"

echo ""
echo "Done. Results in: ${OUTDIR}/"
