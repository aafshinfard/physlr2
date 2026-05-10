#!/bin/bash
# Download CHM13 ONT data for physlr2 long-read testing.
#
# Usage:
#   bash scripts/download_chm13_ont.sh [output_dir]
#
# Downloads:
#   1. rel8 (Guppy 5.0.7) — full dataset re-basecalled with Guppy 5
#   2. Guppy 6.3.7 HAC — from NCBI SRA (requires sra-tools)
#
# After download, extracts chr1 reads for development testing.

set -euo pipefail

OUTDIR="${1:-/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1/data/chm13_ont}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== Downloading CHM13 ONT data to $OUTDIR ==="

# --- rel8: Guppy 5.0.7 ---
REL8_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/nanopore/rel8-guppy-5.0.7/reads.fastq.gz"
REL8_MD5="39262716285ad6efb1f39d374f57dd4e"
REL8_FILE="rel8-guppy5.fastq.gz"

if [ ! -f "$REL8_FILE" ]; then
    echo "Downloading rel8 (Guppy 5.0.7)..."
    wget -c "$REL8_URL" -O "$REL8_FILE"
    echo "Verifying checksum..."
    echo "$REL8_MD5  $REL8_FILE" | md5sum -c -
else
    echo "rel8 already downloaded: $REL8_FILE"
fi

# --- Guppy 6.3.7 HAC from SRA ---
GUPPY6_ACC="SRX19306105"
GUPPY6_FILE="guppy6.fastq.gz"

if [ ! -f "$GUPPY6_FILE" ]; then
    if command -v fasterq-dump &> /dev/null; then
        echo "Downloading Guppy 6.3.7 HAC from SRA ($GUPPY6_ACC)..."
        # Get the SRR accession from SRX
        SRR=$(esearch -db sra -query "$GUPPY6_ACC" | efetch -format runinfo | grep "^SRR" | head -1 | cut -d',' -f1)
        if [ -z "$SRR" ]; then
            echo "Warning: Could not resolve SRR accession from $GUPPY6_ACC"
            echo "Try manually: fasterq-dump SRR22585867 && gzip SRR22585867.fastq"
            SRR="SRR22585867"  # Known accession
        fi
        echo "Downloading $SRR..."
        fasterq-dump "$SRR" -O . -t /tmp -e 8
        gzip "${SRR}.fastq"
        mv "${SRR}.fastq.gz" "$GUPPY6_FILE"
    else
        echo "Warning: sra-tools (fasterq-dump) not found. Skipping Guppy 6.3.7 download."
        echo "Install with: conda install -c bioconda sra-tools"
        echo "Then run: fasterq-dump SRR22585867 -O $OUTDIR && gzip $OUTDIR/SRR22585867.fastq"
    fi
else
    echo "Guppy 6.3.7 already downloaded: $GUPPY6_FILE"
fi

# --- Extract chr1 subset for development testing ---
# Requires: minimap2, samtools, seqtk
CHM13_REF="chm13v2.0.fa"
CHR1_READS="chr1_reads.fastq.gz"

if [ ! -f "$CHR1_READS" ] && [ -f "$REL8_FILE" ]; then
    if command -v minimap2 &> /dev/null && command -v samtools &> /dev/null; then
        echo "Extracting chr1 reads from rel8..."

        # Download CHM13 v2.0 reference if needed
        if [ ! -f "$CHM13_REF" ]; then
            echo "Downloading CHM13 v2.0 reference..."
            wget -c "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz" -O "${CHM13_REF}.gz"
            gunzip "${CHM13_REF}.gz"
        fi

        # Map reads to reference, extract chr1
        echo "Mapping reads to reference (this will take a while)..."
        minimap2 -t 16 -a -x map-ont "$CHM13_REF" "$REL8_FILE" \
            | samtools view -b -F 4 - \
            | samtools sort -@ 4 -o rel8_sorted.bam -
        samtools index rel8_sorted.bam

        echo "Extracting chr1 read names..."
        samtools view rel8_sorted.bam chr1 | cut -f1 | sort -u > chr1_readnames.txt
        echo "Found $(wc -l < chr1_readnames.txt) chr1 reads"

        echo "Extracting chr1 reads from FASTQ..."
        seqtk subseq "$REL8_FILE" chr1_readnames.txt | gzip > "$CHR1_READS"

        echo "chr1 subset: $CHR1_READS"
        echo "Cleaning up intermediate files..."
        rm -f rel8_sorted.bam rel8_sorted.bam.bai chr1_readnames.txt
    else
        echo "Warning: minimap2/samtools not found. Cannot extract chr1 subset."
        echo "Install with: conda install -c bioconda minimap2 samtools seqtk"
    fi
else
    echo "chr1 reads already extracted or source not available"
fi

echo "=== Done ==="
ls -lh "$OUTDIR"/*.fastq.gz 2>/dev/null || echo "No FASTQ files found"
