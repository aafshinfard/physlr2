#!/usr/bin/env bash
# Download test datasets for Physlr.
#
# Datasets:
#   fly_chr4: D. melanogaster chromosome 4 (~1.3 Mbp)
#     - 10x Genomics linked reads (subset mapped to chr4)
#     - Reference genome (chr4 only)
#     - Supernova assembly (chr4 scaffold)
#     - ABySS draft assembly (for scaffolding comparison)
#
# Usage: ./download_test_data.sh <dataset> <outdir>
#   dataset: fly_chr4 | fly_chr2R
#   outdir: output directory (default: data/)

set -euo pipefail

DATASET="${1:-fly_chr4}"
OUTDIR="${2:-data}"

mkdir -p "$OUTDIR"

case "$DATASET" in
  fly_chr4)
    echo "=== Downloading D. melanogaster chr4 test data ==="

    # Reference genome (BDGP Release 6 / dm6)
    REF="$OUTDIR/dm6_chr4.fa"
    if [ ! -f "$REF" ]; then
      echo "Downloading reference genome (chr4)..."
      curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/chromosomes/chr4.fa.gz" \
        | gunzip > "$REF"
      echo "Reference: $(grep -c '^>' "$REF") sequences"
    fi

    # Supernova assembly
    SUPERNOVA="$OUTDIR/fly_supernova.fa"
    if [ ! -f "$SUPERNOVA" ]; then
      echo "Downloading Supernova assembly..."
      curl -sL "http://cf.10xgenomics.com/samples/assembly/2.1.0/fly/fly_pseudohap.fasta.gz" \
        | gunzip > "$SUPERNOVA"
      echo "Supernova: $(grep -c '^>' "$SUPERNOVA") sequences"
    fi

    # Extract chr4 scaffold from Supernova (scaffold 49 maps to chr4)
    SUPERNOVA_CHR4="$OUTDIR/fly_supernova_chr4.fa"
    if [ ! -f "$SUPERNOVA_CHR4" ]; then
      echo "Extracting chr4 scaffold from Supernova assembly..."
      if command -v samtools &>/dev/null; then
        samtools faidx "$SUPERNOVA" 49 > "$SUPERNOVA_CHR4"
      else
        # Fallback: extract scaffold 49 with awk
        awk '/^>49 /{p=1} /^>/{if(p && !/^>49 /)exit} p' "$SUPERNOVA" > "$SUPERNOVA_CHR4"
      fi
      echo "Supernova chr4: $(grep -c '^>' "$SUPERNOVA_CHR4") sequences"
    fi

    # Create a fragmented draft assembly by splitting the Supernova chr4 scaffold
    # into smaller contigs (simulating a draft that needs scaffolding)
    DRAFT1="$OUTDIR/fly_draft1_chr4.fa"
    if [ ! -f "$DRAFT1" ]; then
      echo "Creating fragmented draft assembly (split at 10kb)..."
      python3 -c "
import sys
name = None
seq = []
for line in open('$SUPERNOVA_CHR4'):
    line = line.strip()
    if line.startswith('>'):
        if name and seq:
            full = ''.join(seq)
            chunk = 10000
            for i in range(0, len(full), chunk):
                print(f'>{name}_chunk{i//chunk}')
                print(full[i:i+chunk])
        name = line[1:].split()[0]
        seq = []
    else:
        seq.append(line)
if name and seq:
    full = ''.join(seq)
    chunk = 10000
    for i in range(0, len(full), chunk):
        print(f'>{name}_chunk{i//chunk}')
        print(full[i:i+chunk])
" > "$DRAFT1"
      echo "Draft1 (10kb chunks): $(grep -c '^>' "$DRAFT1") contigs"
    fi

    # Create a second draft with different fragmentation (5kb)
    DRAFT2="$OUTDIR/fly_draft2_chr4.fa"
    if [ ! -f "$DRAFT2" ]; then
      echo "Creating second fragmented draft assembly (split at 5kb)..."
      python3 -c "
import sys
name = None
seq = []
for line in open('$SUPERNOVA_CHR4'):
    line = line.strip()
    if line.startswith('>'):
        if name and seq:
            full = ''.join(seq)
            chunk = 5000
            for i in range(0, len(full), chunk):
                print(f'>{name}_frag{i//chunk}')
                print(full[i:i+chunk])
        name = line[1:].split()[0]
        seq = []
    else:
        seq.append(line)
if name and seq:
    full = ''.join(seq)
    chunk = 5000
    for i in range(0, len(full), chunk):
        print(f'>{name}_frag{i//chunk}')
        print(full[i:i+chunk])
" > "$DRAFT2"
      echo "Draft2 (5kb chunks): $(grep -c '^>' "$DRAFT2") contigs"
    fi

    echo ""
    echo "=== Test data ready in $OUTDIR ==="
    echo "Reference:  $REF"
    echo "Supernova:  $SUPERNOVA_CHR4"
    echo "Draft 1:    $DRAFT1 (10kb fragments)"
    echo "Draft 2:    $DRAFT2 (5kb fragments)"
    echo ""
    echo "NOTE: Linked reads must be downloaded separately (large file)."
    echo "  Full fly reads: http://s3-us-west-2.amazonaws.com/10x.files/samples/assembly/2.1.0/fly/fly_fastqs.tar"
    echo "  Then extract chr4 reads by mapping and filtering."
    ;;

  *)
    echo "Unknown dataset: $DATASET"
    echo "Available: fly_chr4"
    exit 1
    ;;
esac
