#!/bin/bash
#SBATCH --job-name=physlr-preprocess-10x
#SBATCH --output=./logs/preprocess_10x.%j.out
#SBATCH --error=./logs/preprocess_10x.%j.err
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos 1d
#SBATCH --chdir=/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1
set -euo pipefail

# Preprocess 10x Chromium reads using longranger basic.
# Attaches BX:Z: barcode tags to FASTQ reads.
#
# Prerequisites: longranger must be in PATH.
# Install from: https://support.10xgenomics.com/genome-exome/software/downloads/latest
#
# Usage: sbatch scripts/02_preprocess_10x.sh

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
DATADIR="${BASEDIR}/data"
OUTDIR="${DATADIR}/na24143_10x"

mkdir -p "${BASEDIR}/logs"
cd "${OUTDIR}"

echo "=== Preprocessing NA24143 10x Chromium reads ==="

if [ -f na24143.10x.barcoded.fq.gz ]; then
    echo "Barcoded FASTQ already exists, skipping."
    exit 0
fi

# Check for longranger
if ! command -v longranger &>/dev/null; then
    echo "ERROR: longranger not found in PATH."
    echo "Install from https://support.10xgenomics.com/genome-exome/software/downloads/latest"
    echo ""
    echo "Alternative: use the manual barcode extraction below."
    echo ""

    # Manual barcode extraction fallback:
    # 10x Chromium stores the 16bp barcode in the first 16bp of R1.
    # We can extract it and add as BX:Z: tag.
    echo "--- Using manual barcode extraction ---"

    # Interleave R1 and R2, extract barcode from first 16bp of R1
    # The read-RA files contain interleaved R1/R2 pairs
    echo "Concatenating and processing read-RA files..."

    # Process all RA files: extract barcode from first 16bp, trim R1, add BX:Z: tag
    zcat read-RA_si-*_lane-*-chunk-*.fastq.gz \
    | awk '
    BEGIN { OFS="\n" }
    NR%8==1 { # R1 header
        header=$0
        getline seq    # R1 sequence
        getline plus   # +
        getline qual   # R1 quality
        barcode=substr(seq,1,16)
        seq=substr(seq,17)
        qual=substr(qual,17)
        # Add BX:Z: tag to header
        sub(/ .*/, "", header)
        print header " BX:Z:" barcode "-1"
        print seq
        print "+"
        print qual
    }
    NR%8==5 { # R2 header
        header=$0
        getline seq
        getline plus
        getline qual
        sub(/ .*/, "", header)
        print header " BX:Z:" barcode "-1"
        print seq
        print "+"
        print qual
    }
    ' | gzip > na24143.10x.barcoded.fq.gz

    echo "Done. Output: na24143.10x.barcoded.fq.gz"
    exit 0
fi

# Use longranger basic
echo "Running longranger basic..."
longranger basic \
    --id=na24143_10x \
    --fastqs="${OUTDIR}" \
    --localcores=${SLURM_CPUS_PER_TASK:-16}

# Extract the barcoded FASTQ
if [ -f na24143_10x/outs/barcoded.fastq.gz ]; then
    mv na24143_10x/outs/barcoded.fastq.gz na24143.10x.barcoded.fq.gz
    echo "Done. Output: na24143.10x.barcoded.fq.gz"
else
    echo "ERROR: longranger output not found"
    exit 1
fi

echo "=== Preprocessing complete ==="
