# Physlr: Physical Map Comparison

Comparison of physical maps produced by the original Physlr (Python) and Physlr 2 (Rust) on two human cell lines using stLFR linked reads.

## Versions Compared

| Version | Description |
|---------|-------------|
| **Original Physlr** | Python + C++ implementation ([bcgsc/physlr](https://github.com/bcgsc/physlr)) |
| **Physlr 2 v0.10** | Rust rewrite, backbone extraction only (no merge-paths) |
| **Physlr 2 v0.23** | Rust rewrite with merge-paths enabled (ed=25, eh=4, ml=1) |

## Datasets

| Sample | Technology | Source |
|--------|-----------|--------|
| NA12878 | stLFR | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/stLFR/) |
| NA24143 | stLFR | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/stLFR/) |

Reference: GRCh38 (no alt analysis set).

---

## NA12878

### Backbone View

Backbone paths colored by reference chromosome. Each horizontal bar is a backbone path; colors indicate which chromosome it maps to.

| Original Physlr | Physlr 2 v0.10 | Physlr 2 v0.23 |
|:---:|:---:|:---:|
| [<img src="results/comparison/original_na12878_backbone_v2.png" height="200">](results/comparison/original_na12878_backbone_v2.png) | [<img src="results/comparison/v010_na12878_backbone_v2.png" height="200">](results/comparison/v010_na12878_backbone_v2.png) | [<img src="results/comparison/v023_na12878_backbone_v2.png" height="200">](results/comparison/v023_na12878_backbone_v2.png) |

### Reference View

Reference chromosomes colored by backbone path. Shows how well the physical map covers each chromosome.

| Original Physlr | Physlr 2 v0.10 | Physlr 2 v0.23 |
|:---:|:---:|:---:|
| [<img src="results/comparison/original_na12878_reference_v2.png" height="200">](results/comparison/original_na12878_reference_v2.png) | [<img src="results/comparison/v010_na12878_reference_v2.png" height="200">](results/comparison/v010_na12878_reference_v2.png) | [<img src="results/comparison/v023_na12878_reference_v2.png" height="200">](results/comparison/v023_na12878_reference_v2.png) |

---

## NA24143

### Backbone View

| Original Physlr | Physlr 2 v0.10 | Physlr 2 v0.23 |
|:---:|:---:|:---:|
| [<img src="results/comparison/original_na24143_backbone_v2.png" height="200">](results/comparison/original_na24143_backbone_v2.png) | [<img src="results/comparison/v010_na24143_backbone_v2.png" height="200">](results/comparison/v010_na24143_backbone_v2.png) | [<img src="results/comparison/v023_na24143_backbone_v2.png" height="200">](results/comparison/v023_na24143_backbone_v2.png) |

### Reference View

| Original Physlr | Physlr 2 v0.10 | Physlr 2 v0.23 |
|:---:|:---:|:---:|
| [<img src="results/comparison/original_na24143_reference_v2.png" height="200">](results/comparison/original_na24143_reference_v2.png) | [<img src="results/comparison/v010_na24143_reference_v2.png" height="200">](results/comparison/v010_na24143_reference_v2.png) | [<img src="results/comparison/v023_na24143_reference_v2.png" height="200">](results/comparison/v023_na24143_reference_v2.png) |

---

## Merge-Paths Results (v0.23)

The merge-paths step identifies bridge molecules that connect adjacent backbone paths and merges them. Using the default parameters (endpoint-depth=25, min-endpoint-hits=4, min-bridges=1):

| Sample | True Positives | False Positives |
|--------|:-:|:-:|
| NA12878 | 4 | 0 |
| NA24143 | 5 | 0 |
| **Combined** | **9** | **0** |

<sub>Click any thumbnail above to view the full-resolution image.</sub>
