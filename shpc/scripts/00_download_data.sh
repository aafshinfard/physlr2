#!/bin/bash
#SBATCH --job-name=physlr-download
#SBATCH --output=./logs/download.%j.out
#SBATCH --error=./logs/download.%j.err
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

# Download all test data for Physlr benchmarking.
# Uses parallel wget (8 concurrent downloads) to maximize throughput.
# All downloads are resumable — resubmit if the job hits the wall time.
#
# Usage: sbatch scripts/00_download_data.sh

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
DATADIR="${BASEDIR}/data"
JOBS=8

echo "=== Downloading Physlr test data to ${DATADIR} ==="
echo "Using ${JOBS} parallel downloads"
echo "Start: $(date)"

# Download helper: skip if file exists and is non-empty
dl() {
    local url="$1"
    local dest="$2"
    if [ -f "${dest}" ] && [ -s "${dest}" ]; then
        echo "[skip] ${dest} already exists"
        return 0
    fi
    echo "[download] ${dest}"
    mkdir -p "$(dirname "${dest}")"
    wget -c -q --show-progress -t 3 -T 60 "${url}" -O "${dest}" || {
        echo "[WARN] Failed: ${dest}, will retry on resubmit"
        return 0
    }
}
export -f dl

# ─── Build download manifest ─────────────────────────────────────────────────

MANIFEST=$(mktemp)

# Reference genome
cat >> "${MANIFEST}" << 'EOF'
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz|reference/grch38.fa.gz
EOF

# NA12878 stLFR reads
cat >> "${MANIFEST}" << 'EOF'
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/stLFR.split_read.1.fq.gz|na12878_stlfr/na12878.stlfr.R1.fq.gz
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/stLFR.split_read.2.fq.gz|na12878_stlfr/na12878.stlfr.R2.fq.gz
EOF

# NA24143 stLFR reads
cat >> "${MANIFEST}" << 'EOF'
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/stLFR/stLFR_NA24143_split_read.1.fq.gz|na24143_stlfr/na24143.stlfr.R1.fq.gz
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/stLFR/stLFR_NA24143_split_read.2.fq.gz|na24143_stlfr/na24143.stlfr.R2.fq.gz
EOF

# Draft assemblies — NA12878
cat >> "${MANIFEST}" << 'EOF'
https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na12878.ont.baseline.fa|assemblies/na12878.ont.baseline.fa
https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na12878.pempet.baseline.fa|assemblies/na12878.pempet.baseline.fa
https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na12878.stlfr.baseline.fa|assemblies/na12878.stlfr.baseline.fa
https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/JasonChin_Peregrine_PacBioCCS_assembly_05072019/NA12878.ccs.peregrine.fa.gz|assemblies/na12878.pacbio.baseline.fa.gz
EOF

# Draft assemblies — NA24143
cat >> "${MANIFEST}" << 'EOF'
https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na24143.ont.baseline.fa|assemblies/na24143.ont.baseline.fa
https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na24143.pempet.baseline.fa|assemblies/na24143.pempet.baseline.fa
https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na24143.stlfr.baseline.fa|assemblies/na24143.stlfr.baseline.fa
https://www.bcgsc.ca/downloads/btl/physlr/assemblies/hsapiens/baseline/na24143.10x.baseline.fa|assemblies/na24143.10x.baseline.fa
https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/MtSinai_PacBio_Assembly_falcon_03282016/NA24143_hg004_falcon.fa|assemblies/na24143.pacbio.baseline.fa
EOF

# ─── Run parallel downloads ──────────────────────────────────────────────────

echo ""
echo "--- Downloading $(wc -l < "${MANIFEST}") files with ${JOBS} parallel jobs ---"

cat "${MANIFEST}" | xargs -P "${JOBS}" -I {} bash -c '
    url="${1%%|*}"
    dest="'"${DATADIR}"'/${1##*|}"
    dl "${url}" "${dest}"
' _ {}

rm -f "${MANIFEST}"

# ─── Post-processing ─────────────────────────────────────────────────────────

echo ""
echo "--- Post-processing ---"

# Decompress reference if needed
if [ -f "${DATADIR}/reference/grch38.fa.gz" ] && [ ! -f "${DATADIR}/reference/grch38.fa" ]; then
    echo "Decompressing reference genome..."
    gunzip "${DATADIR}/reference/grch38.fa.gz"
fi

# Index reference if samtools available
if [ -f "${DATADIR}/reference/grch38.fa" ] && [ ! -f "${DATADIR}/reference/grch38.fa.fai" ]; then
    if command -v samtools &>/dev/null; then
        echo "Indexing reference genome..."
        samtools faidx "${DATADIR}/reference/grch38.fa"
    fi
fi

# Decompress PacBio assembly if needed
if [ -f "${DATADIR}/assemblies/na12878.pacbio.baseline.fa.gz" ] && [ ! -f "${DATADIR}/assemblies/na12878.pacbio.baseline.fa" ]; then
    echo "Decompressing NA12878 PacBio assembly..."
    gunzip "${DATADIR}/assemblies/na12878.pacbio.baseline.fa.gz"
fi

# ─── Download 10x Chromium reads ─────────────────────────────────────────────

echo ""
echo "--- Downloading NA24143 10x Chromium reads ---"

CHROMIUM_DIR="${DATADIR}/na24143_10x"
mkdir -p "${CHROMIUM_DIR}"
BASE_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/10Xgenomics_ChromiumGenome/NA24143.fastqs"

# Build list of 10x files and download in parallel
CHROMIUM_MANIFEST=$(mktemp)
for si in ACCGGCTC CGTCCTAG GAGTTAGT TTAAAGCA; do
    for lane in 001 002 003 004 005 006 007 008; do
        for read_type in RA I1; do
            fname=$(wget -q -O - "${BASE_URL}/" 2>/dev/null \
                | grep -oP "read-${read_type}_si-${si}_lane-${lane}-chunk-\d+\.fastq\.gz" \
                | head -1 || true)
            if [ -n "${fname}" ]; then
                echo "${BASE_URL}/${fname}|na24143_10x/${fname}" >> "${CHROMIUM_MANIFEST}"
            fi
        done
    done
done

if [ -s "${CHROMIUM_MANIFEST}" ]; then
    echo "Downloading $(wc -l < "${CHROMIUM_MANIFEST}") 10x Chromium files..."
    cat "${CHROMIUM_MANIFEST}" | xargs -P "${JOBS}" -I {} bash -c '
        url="${1%%|*}"
        dest="'"${DATADIR}"'/${1##*|}"
        dl "${url}" "${dest}"
    ' _ {}
fi
rm -f "${CHROMIUM_MANIFEST}"

# ─── Summary ─────────────────────────────────────────────────────────────────

echo ""
echo "=== Download complete ==="
echo "End: $(date)"
echo ""
echo "--- Disk usage ---"
du -sh "${DATADIR}"/*/ 2>/dev/null || true
echo ""
du -sh "${DATADIR}" 2>/dev/null || true
