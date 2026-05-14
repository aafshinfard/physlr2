# Long-Read Support: Design & Development Plan

## Overview

Physlr2 supports linked reads (10X, stLFR) on `main`. This document describes the
design for adding long-read support (ONT, PacBio) on the `long-reads-v2` branch.

The two protocols share core algorithms (overlap graph, molecule separation, backbone
extraction) but differ in preprocessing, indexing, and coordinate handling.

---

## Architecture: Separate Pipeline Functions

Linked-read and long-read workflows are implemented as **separate pipeline functions**
called from a shared CLI handler:

```rust
Commands::PhysicalMap { ... } => {
    if proto == Protocol::Long {
        run_physical_map_long(&input, &outdir, ...)?;
    } else {
        run_physical_map_linked(&input, &outdir, ...)?;
    }
}
```

Each function has its own indexing, filtering, and mapping logic. They share:
- Overlap computation (`overlap::compute_overlap`)
- Edge filtering (`overlap::filter_edges_by_percentile`)
- Molecule separation (`molecules::separate_molecules`)
- Backbone extraction (`backbone::extract_named_backbones`)
- Report generation (`report::*`)

This keeps protocol-specific logic isolated while avoiding duplication of core algorithms.

---

## Homopolymer Compression (HPC) Coordinate Strategy

### Problem

Long reads benefit from homopolymer compression (HPC) to reduce noise from ONT/PacBio
systematic errors in homopolymer runs. However, the final physical map and downstream
outputs (PAF mapping, scaffolding) must use real genomic coordinates.

### Solution: Work in HPC space, translate back on small sequences

The key insight: **reads are huge (330 GB), reference and draft assembly are small (3 GB)**.
So we do everything in HPC space and only translate coordinates for the small sequences.

#### Pipeline flow for long reads:

```
1. HPC-compress reads (330 GB) → .hpc.fq
2. ntcard → nthits → makebf → indexlr on HPC reads
   → bx_to_mxs (minimizer hashes in HPC space)
3. Overlap → molecules → backbone (all in HPC space)
   → backbone paths (topology only, no coordinates)

4. Reference mapping (when --reference provided):
   a. HPC-compress reference (3 GB, seconds)
   b. Store per-chromosome coord_map: hpc_pos → original_pos
   c. Index HPC reference → extract ordered minimizers
   d. Generate position map from HPC reference (hpc_mx_index → hpc_bp)
   e. Map HPC reference to backbone → PAF in HPC coordinates
   f. Translate PAF coordinates back via coord_map → real bp coordinates
   g. Plot using real bp coordinates

5. Scaffolding (when --draft provided):
   a. HPC-compress draft assembly
   b. Store per-contig coord_map
   c. Index HPC contigs → map to backbone in HPC space
   d. Determine scaffold ordering (which contigs, what order/orientation)
   e. Apply ordering to original (non-HPC) contig sequences
   f. Output scaffolded FASTA with original sequences
```

#### Coordinate translation details:

`homopolymer_compress(seq)` returns `(compressed_seq, coord_map)` where
`coord_map[compressed_pos] = original_pos`.

For PAF records, translate:
- `Qstart`, `Qend` (reference query positions): look up in reference coord_map
- Minimizer indices in position map: generate position map from HPC reference,
  then translate the bp values through coord_map

For scaffolding:
- Scaffold ordering is determined by contig-to-backbone mapping (HPC space)
- Final FASTA uses original contig sequences — no coordinate translation needed,
  just the ordering and orientation

#### Data structures needed:

```rust
/// Per-sequence coordinate mapping from HPC to original space.
struct HpcCoordMap {
    name: String,
    coord_map: Vec<usize>,  // hpc_pos → original_pos
    original_len: usize,
    compressed_len: usize,
}

/// Collection of coord maps for a FASTA file (reference or assembly).
struct HpcIndex {
    sequences: Vec<HpcCoordMap>,
    /// Path to the HPC-compressed FASTA written to disk.
    compressed_path: PathBuf,
}
```

#### Functions to implement:

```
minimizer::hpc_compress_fasta(input_path, output_path) → HpcIndex
    Reads FASTA, writes HPC-compressed FASTA, returns coord maps.

map::translate_paf_coordinates(paf_records, hpc_index) → paf_records
    Translates PAF Qstart/Qend from HPC to original coordinates.

minimizer::translate_position_map(pos_map, hpc_index) → pos_map
    Translates position map bp values from HPC to original coordinates.
```

---

## Protocol Differences Summary

| Aspect | Linked Reads | Long Reads |
|--------|-------------|------------|
| Input format | FASTQ with BX:Z: or #barcode tags | FASTQ, read name = molecule ID |
| Preprocessing | None | Homopolymer compression |
| indexlr flags | `--bx` | `--long` (no `--bx`) |
| Minimizers per unit | ~100-5000 per barcode | ~50-50000 per read |
| Default min_count | 100 | 50 |
| Default max_count | 5000 | 50000 |
| Default min_overlap | 10 | 5 |
| Molecule separation | Needed (multiple molecules per barcode) | Each read is one molecule |
| Reference mapping | Direct (no HPC) | HPC + coordinate translation |
| Scaffolding | Direct (no HPC) | HPC mapping + original sequences |
| Memory | ~200 GB for human | ~600 GB for human (estimate) |

---

## Development Roadmap

### Phase 1: Basic long-read pipeline (current)
- [x] Protocol detection (auto/linked/long)
- [x] Homopolymer compression function
- [x] External tools pipeline with HPC preprocessing
- [x] indexlr `--long` flag support
- [x] Protocol-specific defaults (min_count, max_count, min_overlap)
- [ ] Separate pipeline functions (`run_physical_map_long`, `run_physical_map_linked`)
- [ ] HPC-aware reference mapping with coordinate translation
- [ ] HPC-aware scaffolding with coordinate translation
- [ ] End-to-end test on CHM13 ONT data

### Phase 2: Evaluation & tuning
- [ ] Evaluate physical map quality on CHM13 (path count, contiguity, misassemblies)
- [ ] Compare with/without HPC
- [ ] Tune parameters (min_count, max_count, min_overlap, edge_percentile)
- [ ] Evaluate molecule separation strategies for long reads
  - Each read is already a molecule — does bc+cc separation help or hurt?
  - May need a simplified or skipped molecule separation step
- [ ] Memory profiling and optimization

### Phase 3: Scaffolding & integration
- [ ] Test scaffolding with long-read physical maps
- [ ] Evaluate on multiple datasets (CHM13, HG002, etc.)
- [ ] Update documentation and README
- [ ] Merge to main

### Phase 4: Optimizations (future)
- [ ] Streaming HPC compression (avoid writing full .hpc.fq to disk)
- [ ] Parallel HPC compression
- [ ] Investigate whether molecule separation can be simplified for long reads
- [ ] Mixed-protocol support (linked + long reads together)

---

## Test Data

| Dataset | Species | Technology | Size | Location (ec-hub) |
|---------|---------|-----------|------|-------------------|
| CHM13 rel8 Guppy v5 | Human | ONT | 330 GB | `/projects/.../data/chm13_ont/rel8-guppy5.fastq.gz` |
| NA24143 stLFR | Human | stLFR linked | ~220 GB | `/projects/.../data/na24143_stlfr/` |
| NA12878 stLFR | Human | stLFR linked | ~240 GB | `/projects/.../data/na12878_stlfr/` |
| GRCh38 reference | Human | — | 3 GB | `/projects/.../data/reference/grch38.fa` |

---

## Branch Strategy

- `main`: Linked-read support only. Tagged releases (v0.29, etc.)
- `long-reads-v2`: Long-read development. Rebased on main periodically.
- Merge to main when Phase 1 is complete and results are validated.
