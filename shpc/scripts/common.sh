#!/bin/bash
# Common functions and variables for Physlr SLURM scripts.
# Source this file from other scripts: source scripts/common.sh

BASEDIR="/projects/site/dia/sbx/workspace/afshinfa/flexdev/physlr/rust/test1"
PHYSLR="${BASEDIR}/physlr-next/target/release/physlr"
DATADIR="${BASEDIR}/data"
REFDIR="${DATADIR}/reference"
ASMDIR="${DATADIR}/assemblies"
REFERENCE="${REFDIR}/grch38.fa"
GENOME_SIZE=3088269832  # GRCh38 primary assembly size

THREADS="${SLURM_CPUS_PER_TASK:-28}"

# Verify physlr binary exists
check_physlr() {
    if [ ! -x "${PHYSLR}" ]; then
        echo "ERROR: Physlr binary not found at ${PHYSLR}"
        echo "Run scripts/01_build_physlr.sh first."
        exit 1
    fi
}

# Run a pipeline step with timing
run_step() {
    local step_name="$1"
    shift
    echo ""
    echo "=== ${step_name} ==="
    echo "Command: $*"
    local start_time=$(date +%s)
    "$@"
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    echo "--- ${step_name}: ${elapsed}s ---"
}

# Run the full Physlr physical map pipeline.
# Arguments:
#   $1 - sample name (e.g., na12878_stlfr)
#   $2 - input FASTQ file(s), space-separated
#   $3 - output directory
#   $4 - k-mer size
#   $5 - window size
#   $6 - overlap percentile (e.g., 85 for stLFR, 92.5 for 10x)
#   $7 - molecule strategy (e.g., bc+cc)
run_physical_map() {
    local sample="$1"
    local input_fq="$2"
    local outdir="$3"
    local k="${4:-40}"
    local w="${5:-32}"
    local overlap_pct="${6:-85}"
    local mol_strategy="${7:-bc+cc}"

    mkdir -p "${outdir}"
    cd "${outdir}"

    check_physlr

    # Step 1: Index minimizers
    if [ ! -f "${sample}.reads.tsv" ]; then
        run_step "Index minimizers" \
            "${PHYSLR}" index -k "${k}" -w "${w}" ${input_fq} \
            -o "${sample}.reads.tsv"
    else
        echo "--- Skipping index (${sample}.reads.tsv exists) ---"
    fi

    # Step 2: Filter minimizers
    if [ ! -f "${sample}.filtered.tsv" ]; then
        run_step "Filter minimizers" \
            "${PHYSLR}" filter-minimizers "${sample}.reads.tsv" \
            -o "${sample}.filtered.tsv" \
            --min-count 100 --max-count 5000
    else
        echo "--- Skipping filter (${sample}.filtered.tsv exists) ---"
    fi

    # Step 3: Compute overlaps
    if [ ! -f "${sample}.overlap.tsv" ]; then
        run_step "Compute overlaps" \
            "${PHYSLR}" overlap "${sample}.filtered.tsv" \
            -o "${sample}.overlap.tsv"
    else
        echo "--- Skipping overlap (${sample}.overlap.tsv exists) ---"
    fi

    # Step 4: Filter overlaps by percentile
    if [ ! -f "${sample}.overlap.filtered.tsv" ]; then
        run_step "Filter overlaps (p=${overlap_pct})" \
            "${PHYSLR}" filter-overlap "${sample}.overlap.tsv" \
            -o "${sample}.overlap.filtered.tsv" \
            -p "${overlap_pct}"
    else
        echo "--- Skipping filter-overlap (${sample}.overlap.filtered.tsv exists) ---"
    fi

    # Step 5: Separate molecules
    if [ ! -f "${sample}.molecules.tsv" ]; then
        run_step "Separate molecules (strategy=${mol_strategy})" \
            "${PHYSLR}" molecules "${sample}.overlap.filtered.tsv" \
            -o "${sample}.molecules.tsv" \
            --strategy "${mol_strategy}"
    else
        echo "--- Skipping molecules (${sample}.molecules.tsv exists) ---"
    fi

    # Step 6: Extract backbone paths
    if [ ! -f "${sample}.backbone.tsv" ]; then
        run_step "Extract backbone paths" \
            "${PHYSLR}" backbone "${sample}.molecules.tsv" \
            -o "${sample}.backbone.tsv" \
            --prune-branches 10 --prune-bridges 10
    else
        echo "--- Skipping backbone (${sample}.backbone.tsv exists) ---"
    fi

    # Step 7: Physical map metrics
    run_step "Physical map metrics" \
        "${PHYSLR}" path-metrics "${sample}.backbone.tsv"

    echo ""
    echo "=== Physical map complete: ${sample} ==="
    echo "Backbone: ${outdir}/${sample}.backbone.tsv"
}

# Scaffold a draft assembly using the physical map.
# Arguments:
#   $1 - sample name
#   $2 - draft assembly FASTA
#   $3 - draft label (e.g., ont, pacbio, abyss)
#   $4 - output directory (where physical map files are)
#   $5 - k-mer size
#   $6 - window size
scaffold_assembly() {
    local sample="$1"
    local draft_fa="$2"
    local draft_label="$3"
    local outdir="$4"
    local k="${5:-40}"
    local w="${6:-32}"

    cd "${outdir}"
    check_physlr

    local prefix="${sample}.${draft_label}"

    # Index contigs
    if [ ! -f "${prefix}.contigs.tsv" ]; then
        run_step "Index contigs (${draft_label})" \
            "${PHYSLR}" index-contigs -k "${k}" -w "${w}" "${draft_fa}" \
            -o "${prefix}.contigs.tsv"
    fi

    # Map contigs to backbone
    if [ ! -f "${prefix}.map.bed" ]; then
        run_step "Map contigs to backbone (${draft_label})" \
            "${PHYSLR}" map \
            "${sample}.backbone.tsv" \
            "${sample}.filtered.tsv" \
            "${prefix}.contigs.tsv" \
            -o "${prefix}.map.bed"
    fi

    # Convert BED to scaffold paths
    if [ ! -f "${prefix}.path.tsv" ]; then
        run_step "BED to scaffold paths (${draft_label})" \
            "${PHYSLR}" bed-to-path "${prefix}.map.bed" \
            -o "${prefix}.path.tsv"
    fi

    # Produce scaffolded FASTA
    if [ ! -f "${prefix}.scaffold.fa" ]; then
        run_step "Scaffold FASTA (${draft_label})" \
            "${PHYSLR}" path-to-fasta "${draft_fa}" "${prefix}.path.tsv" \
            -o "${prefix}.scaffold.fa"
    fi

    # Metrics
    run_step "Scaffold metrics (${draft_label})" \
        "${PHYSLR}" metrics "${prefix}.scaffold.fa" \
        -g "${GENOME_SIZE}" -l "${prefix}"

    echo ""
    echo "=== Scaffolding complete: ${prefix} ==="
    echo "Scaffold: ${outdir}/${prefix}.scaffold.fa"
}

# Run QUAST evaluation.
# Arguments:
#   $1 - assembly FASTA
#   $2 - label
#   $3 - output directory
run_quast() {
    local assembly="$1"
    local label="$2"
    local outdir="$3"

    if ! command -v quast &>/dev/null; then
        echo "WARNING: quast not found, skipping QUAST evaluation for ${label}"
        return 0
    fi

    mkdir -p "${outdir}"

    run_step "QUAST (${label})" \
        quast-lg -t "${THREADS}" \
        --fast --large --scaffold-gap-max-size 100000 --min-identity 95 \
        -R "${REFERENCE}" \
        -o "${outdir}/quast_${label}" \
        "${assembly}"

    if [ -f "${outdir}/quast_${label}/transposed_report.tsv" ]; then
        echo "--- QUAST summary (${label}) ---"
        cat "${outdir}/quast_${label}/transposed_report.tsv"
    fi
}
