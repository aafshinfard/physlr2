#!/bin/bash
set -euo pipefail

# Collect and summarize results from all Physlr test runs.
#
# Usage: bash scripts/20_collect_results.sh

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
OUTDIR="${BASEDIR}/output"
REPORT="${BASEDIR}/results_summary.tsv"

echo "=== Collecting Physlr Benchmark Results ==="
echo ""

# Header
printf "Sample\tDraft\tStage\tSequences\tTotal_Length\tN50\tNG50\tNGA50\tL50\n" > "${REPORT}"

collect_quast() {
    local sample="$1"
    local draft="$2"
    local stage="$3"
    local quast_dir="$4"

    local tsv="${quast_dir}/transposed_report.tsv"
    if [ ! -f "${tsv}" ]; then
        return
    fi

    local seqs=$(awk -F'\t' 'NR==2{print $2}' "${tsv}" 2>/dev/null || echo "NA")
    local total=$(awk -F'\t' '/Total length\t/{print $2}' "${tsv}" 2>/dev/null | head -1 || echo "NA")
    local n50=$(awk -F'\t' '/^N50\t/{print $2}' "${tsv}" 2>/dev/null || echo "NA")
    local ng50=$(awk -F'\t' '/^NG50\t/{print $2}' "${tsv}" 2>/dev/null || echo "NA")
    local nga50=$(awk -F'\t' '/^NGA50\t/{print $2}' "${tsv}" 2>/dev/null || echo "NA")
    local l50=$(awk -F'\t' '/^L50\t/{print $2}' "${tsv}" 2>/dev/null || echo "NA")

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${sample}" "${draft}" "${stage}" \
        "${seqs}" "${total}" "${n50}" "${ng50}" "${nga50}" "${l50}" >> "${REPORT}"
}

# Collect from all samples
for sample_dir in "${OUTDIR}"/*/; do
    sample=$(basename "${sample_dir}")

    # Find all QUAST directories
    for quast_dir in "${sample_dir}"/quast_*/; do
        [ -d "${quast_dir}" ] || continue
        label=$(basename "${quast_dir}" | sed 's/^quast_//')

        # Parse sample.draft.stage from label
        if [[ "${label}" == *".baseline" ]]; then
            draft=$(echo "${label}" | sed "s/^${sample}\.\(.*\)\.baseline$/\1/")
            collect_quast "${sample}" "${draft}" "baseline" "${quast_dir}"
        elif [[ "${label}" == *".scaffold" ]]; then
            draft=$(echo "${label}" | sed "s/^${sample}\.\(.*\)\.scaffold$/\1/")
            collect_quast "${sample}" "${draft}" "scaffold" "${quast_dir}"
        fi
    done
done

echo "Results written to: ${REPORT}"
echo ""
echo "--- Summary ---"
column -t -s$'\t' "${REPORT}"

# Also collect physical map metrics
echo ""
echo "--- Physical Map Metrics ---"
for sample_dir in "${OUTDIR}"/*/; do
    sample=$(basename "${sample_dir}")
    backbone="${sample_dir}/${sample}.backbone.tsv"
    if [ -f "${backbone}" ]; then
        echo ""
        echo "=== ${sample} ==="
        wc -l "${backbone}" | awk '{print "Backbone paths: " $1}'
        awk '{print NF}' "${backbone}" | sort -rn | head -5 | awk 'NR==1{print "Longest path: " $1 " molecules"}'
        awk '{total+=NF} END{print "Total molecules: " total}' "${backbone}"
    fi
done

echo ""
echo "=== Done ==="
