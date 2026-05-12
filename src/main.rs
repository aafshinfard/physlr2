use anyhow::Result;
use clap::{Parser, Subcommand};
use std::io::Write;

use physlr::backbone::BackboneConfig;
use physlr::scaffold::ScaffoldConfig;
use physlr::Timer;

#[derive(Parser)]
#[command(
    name = "physlr",
    version = "2.0.0",
    about = "Physlr: Next-generation physical maps from linked reads",
    long_about = "Physlr constructs de novo physical maps using linked reads (10X Genomics or \
                   MGI stLFR) and uses them to scaffold genome assemblies."
)]
struct Cli {
    /// Verbosity level (0=silent, 1=info, 2=debug)
    #[arg(short = 'v', long, default_value_t = 1)]
    verbose: u8,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = num_cpus())]
    threads: usize,

    #[command(subcommand)]
    command: Commands,
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get().min(16))
        .unwrap_or(1)
}

#[derive(Subcommand)]
enum Commands {
    /// Extract (k,w)-minimizers from FASTA/FASTQ, grouping by barcode
    Index {
        /// Input FASTA/FASTQ file(s) (may be gzipped)
        #[arg(required = true)]
        input: Vec<String>,
        /// Output minimizer TSV file
        #[arg(short, long, default_value = "-")]
        output: String,
        /// K-mer size
        #[arg(short, long, default_value_t = 32)]
        k: usize,
        /// Window size
        #[arg(short, long, default_value_t = 32)]
        w: usize,
        /// Minimizer extraction backend: "btllib" (default) or "builtin"
        #[arg(long, default_value = "btllib")]
        indexer: String,
        /// Number of threads for btllib indexlr (ignored for builtin)
        #[arg(short, long, default_value_t = 1)]
        threads: usize,
        /// Repeat Bloom filter file (from repeat-filter command). Skips repetitive k-mers.
        #[arg(long)]
        repeat_bf: Option<String>,
    },

    /// Detect repetitive k-mers and build a Bloom filter.
    /// Equivalent to ntcard + find-ntcard-mode + nthits + physlr-makebf.
    RepeatFilter {
        /// Input FASTA/FASTQ file(s) (may be gzipped)
        #[arg(required = true)]
        input: Vec<String>,
        /// Output Bloom filter file
        #[arg(short, long)]
        output: String,
        /// K-mer size (must match the k used for indexing)
        #[arg(short, long, default_value_t = 32)]
        k: usize,
        /// Repeat threshold multiplier: threshold = mode × multiplier [3]
        #[arg(long, default_value_t = 3)]
        multiplier: u32,
        /// Bloom filter size in bytes [10000000000] (10 GB)
        #[arg(long, default_value_t = 10_000_000_000)]
        bf_size: u64,
        /// Count-Min Sketch size in bytes [4000000000] (4 GB)
        #[arg(long, default_value_t = 4_000_000_000)]
        cms_size: u64,
        /// Also write the k-mer histogram to this file
        #[arg(long)]
        histogram: Option<String>,
    },

    /// Extract ordered minimizers from FASTA (for contigs/reference)
    IndexContigs {
        /// Input FASTA file
        input: String,
        /// Output minimizer TSV file
        #[arg(short, long, default_value = "-")]
        output: String,
        /// K-mer size
        #[arg(short, long, default_value_t = 32)]
        k: usize,
        /// Window size
        #[arg(short, long, default_value_t = 32)]
        w: usize,
    },

    /// Filter barcodes by minimizer count and remove singleton/repetitive minimizers
    FilterMinimizers {
        /// Input minimizer TSV file
        input: String,
        /// Output filtered TSV file
        #[arg(short, long, default_value = "-")]
        output: String,
        /// Minimum minimizers per barcode
        #[arg(short = 'n', long, default_value_t = 100)]
        min_count: usize,
        /// Maximum minimizers per barcode
        #[arg(short = 'N', long, default_value_t = 5000)]
        max_count: usize,
        /// Maximum minimizer frequency (0 = auto via IQR)
        #[arg(short = 'C', long, default_value_t = 0)]
        max_freq: u32,
    },

    /// Compute the barcode overlap graph
    Overlap {
        /// Input minimizer TSV file (filtered)
        input: String,
        /// Output overlap graph TSV file
        #[arg(short, long, default_value = "-")]
        output: String,
        /// Minimum shared minimizers for an edge
        #[arg(short, long, default_value_t = 10)]
        min_shared: u32,
    },

    /// Filter edges from an overlap graph by percentile
    FilterOverlap {
        /// Input overlap graph TSV
        input: String,
        /// Output filtered graph TSV
        #[arg(short, long, default_value = "-")]
        output: String,
        /// Percentile of edges to remove (0-100)
        #[arg(short, long, default_value_t = 90.0)]
        percentile: f64,
    },

    /// Separate barcodes into molecules
    Molecules {
        /// Input overlap graph TSV
        input: String,
        /// Output molecule graph TSV
        #[arg(short, long, default_value = "-")]
        output: String,
        /// Separation strategy: bc, cc, k3, k3bin, sqcos, sqcosbin, distributed joined with +
        #[arg(long, default_value = "bc+cc")]
        strategy: String,
        /// Cosine similarity threshold for sqcos/sqcosbin [0.75]
        #[arg(long, default_value_t = 0.75)]
        sqcos_threshold: f64,
        /// Skip sqcos splitting for components smaller than this [10]
        #[arg(long, default_value_t = 10)]
        skip_small: usize,
        /// Max bin size for random binning in k3bin/sqcosbin [50]
        #[arg(long, default_value_t = 50)]
        bin_max_size: usize,
        /// Merge cutoff: min cross-edges to merge communities (-1 = no merge) [20]
        #[arg(long, default_value_t = 20)]
        merge_cutoff: i64,
        /// Random seed for deterministic binning [42]
        #[arg(long, default_value_t = 42)]
        seed: u64,
    },

    /// Split minimizers by molecule: assign each barcode's minimizers to its
    /// individual molecules based on neighbor overlap in the molecule graph.
    /// Produces cleaner backbone-vs-reference visualizations.
    SplitMinimizers {
        /// Input molecule graph TSV
        mol_graph: String,
        /// Input barcode-to-minimizer TSV
        minimizers: String,
        /// Output molecule-to-minimizer TSV
        #[arg(short, long, default_value = "-")]
        output: String,
    },

    /// Trace molecule separation for specific barcodes (diagnostic)
    TraceMolecules {
        /// Input overlap graph TSV
        input: String,
        /// Comma-separated list of barcodes to trace (or "top5" for highest degree)
        #[arg(long, default_value = "top5")]
        barcodes: String,
        /// Separation strategy
        #[arg(long, default_value = "distributed+sqcosbin")]
        strategy: String,
        /// Cosine similarity threshold for sqcos/sqcosbin [0.75]
        #[arg(long, default_value_t = 0.75)]
        sqcos_threshold: f64,
        /// Skip sqcos splitting for components smaller than this [10]
        #[arg(long, default_value_t = 10)]
        skip_small: usize,
        /// Max bin size for random binning in k3bin/sqcosbin [50]
        #[arg(long, default_value_t = 50)]
        bin_max_size: usize,
        /// Merge cutoff: min cross-edges to merge communities (-1 = no merge) [20]
        #[arg(long, default_value_t = 20)]
        merge_cutoff: i64,
    },

    /// Extract backbone paths (physical map) from the overlap graph
    Backbone {
        /// Input overlap graph TSV (molecule-level)
        input: String,
        /// Output path file
        #[arg(short, long, default_value = "-")]
        output: String,
        /// Minimum branch size for pruning
        #[arg(long, default_value_t = 10)]
        prune_branches: usize,
        /// Minimum bridge size
        #[arg(long, default_value_t = 10)]
        prune_bridges: usize,
        /// Minimum junction branch size
        #[arg(long, default_value_t = 200)]
        prune_junctions: usize,
        /// Minimum path length to output
        #[arg(long, default_value_t = 50)]
        min_component_size: usize,
    },

    /// Merge adjacent backbone paths using split-minimizer bridge evidence.
    /// Optional post-processing step that finds non-backbone molecules
    /// bridging path endpoints to merge adjacent paths.
    MergePaths {
        /// Input backbone path file
        path_file: String,
        /// Split minimizer TSV (molecule-level minimizers)
        split_mxs: String,
        /// Output merged path file
        #[arg(short, long, default_value = "-")]
        output: String,
        /// Number of molecules from each path end to use as endpoints
        #[arg(long, default_value_t = 25)]
        endpoint_depth: usize,
        /// Minimum shared minimizers for a bridge molecule
        #[arg(long, default_value_t = 3)]
        min_shared_mx: usize,
        /// Minimum bridge molecules to accept a link
        #[arg(long, default_value_t = 2)]
        min_bridges: usize,
        /// Maximum paths a bridge molecule can connect (specificity filter)
        #[arg(long, default_value_t = 2)]
        max_connections: usize,
        /// Minimum path length to include in merging
        #[arg(long, default_value_t = 50)]
        min_path_size: usize,
        /// Maximum candidate links per endpoint (promiscuous endpoint filter)
        #[arg(long, default_value_t = 1)]
        max_links_per_endpoint: usize,
        /// Minimum bridge density (bridges / min_path_len). Set to 0 to disable.
        #[arg(long, default_value_t = 0.01)]
        min_bridge_density: f64,
        /// Minimum endpoint molecules a bridge must connect to on each side.
        /// Higher values require stronger neighborhood evidence per bridge.
        #[arg(long, default_value_t = 4)]
        min_endpoint_hits: usize,
    },

    /// Map sequences to the physical map
    Map {
        /// Backbone path file
        path_file: String,
        /// Target minimizer TSV (backbone molecules)
        target_mxs: String,
        /// Query minimizer TSV (contigs or reference)
        query_mxs: String,
        /// Output BED file
        #[arg(short, long, default_value = "-")]
        output: String,
        /// Minimum score for a mapping
        #[arg(short = 'n', long, default_value_t = 10)]
        min_score: u32,
    },

    /// Map sequences to the physical map and output PAF format.
    /// Used for backbone-vs-reference visualization.
    MapPaf {
        /// Backbone path file (.backbone.tsv)
        path_file: String,
        /// Target minimizer TSV (backbone molecules or split minimizers)
        target_mxs: String,
        /// Query minimizer TSV (reference, from indexlr)
        query_mxs: String,
        /// Output PAF file
        #[arg(short, long, default_value = "-")]
        output: String,
        /// Minimum score for a mapping
        #[arg(short = 'n', long, default_value_t = 10)]
        min_score: u32,
        /// IQR coefficient for whisker bounds [1.5]
        #[arg(long, default_value_t = 1.5)]
        coef: f64,
        /// Only output top N backbone paths by total score (0 = all)
        #[arg(long, default_value_t = 0)]
        top_n: usize,
        /// Minimizer type: "barcode" (default) or "split" (molecule-level)
        #[arg(long, default_value = "barcode")]
        mx_type: String,
    },

    /// Convert BED mappings to scaffold paths
    BedToPath {
        /// Input BED file
        input: String,
        /// Output path file
        #[arg(short, long, default_value = "-")]
        output: String,
        /// Minimum mapping score
        #[arg(short = 'n', long, default_value_t = 10)]
        min_score: u32,
    },

    /// Produce scaffolded FASTA from paths
    PathToFasta {
        /// Input FASTA file (draft assembly)
        fasta: String,
        /// Input path file
        path_file: String,
        /// Output FASTA file
        #[arg(short, long, default_value = "-")]
        output: String,
        /// Gap size (Ns between contigs)
        #[arg(long, default_value_t = 100)]
        gap_size: usize,
        /// Minimum scaffold length
        #[arg(long, default_value_t = 0)]
        min_length: usize,
    },

    /// Compute and report assembly metrics
    Metrics {
        /// Input FASTA file
        input: String,
        /// Expected genome size (for NG50)
        #[arg(short = 'g', long)]
        genome_size: Option<u64>,
        /// Label for the assembly
        #[arg(short, long, default_value = "assembly")]
        label: String,
    },

    /// Compute physical map metrics
    PathMetrics {
        /// Input path file
        input: String,
        /// Expected number of molecules (for NG50)
        #[arg(short = 'g', long)]
        expected_molecules: Option<usize>,
        /// Minimum path size to include
        #[arg(long, default_value_t = 1)]
        min_component_size: usize,
    },

    /// Generate a DOT visualization of backbone paths
    BackboneDot {
        /// Input path file
        input: String,
        /// Output DOT file
        #[arg(short, long, default_value = "-")]
        output: String,
    },

    /// End-to-end pipeline: build physical map from FASTQ and scaffold a draft assembly
    Pipeline {
        /// Input linked-read FASTQ file(s) (may be gzipped)
        #[arg(required = true)]
        input: Vec<String>,
        /// Draft assembly FASTA to scaffold
        #[arg(long)]
        draft: String,
        /// Output directory
        #[arg(short, long, default_value = ".")]
        outdir: String,
        /// Output prefix
        #[arg(short, long, default_value = "physlr")]
        prefix: String,
        /// K-mer size for minimizer extraction
        #[arg(short, long, default_value_t = 32)]
        k: usize,
        /// Window size for minimizer extraction
        #[arg(short, long, default_value_t = 32)]
        w: usize,
        /// Minimum minimizers per barcode
        #[arg(long, default_value_t = 100)]
        min_bx_count: usize,
        /// Maximum minimizers per barcode
        #[arg(long, default_value_t = 5000)]
        max_bx_count: usize,
        /// Minimum shared minimizers for overlap
        #[arg(long, default_value_t = 10)]
        min_overlap: u32,
        /// Edge filter percentile
        #[arg(long, default_value_t = 90.0)]
        edge_percentile: f64,
        /// Minimum branch size for pruning
        #[arg(long, default_value_t = 10)]
        prune_branches: usize,
        /// Minimum path length
        #[arg(long, default_value_t = 50)]
        min_path_size: usize,
        /// Minimum mapping score
        #[arg(long, default_value_t = 10)]
        min_map_score: u32,
        /// Gap size in scaffolds (Ns between contigs)
        #[arg(long, default_value_t = 100)]
        gap_size: usize,
        /// Expected genome size in bp (optional, for NG50 reporting only)
        #[arg(short = 'g', long)]
        genome_size: Option<u64>,
    },

    /// Run the full physical map pipeline from linked-read FASTQ files
    PhysicalMap {
        /// Input linked-read FASTQ file(s) (may be gzipped)
        #[arg(required = true)]
        input: Vec<String>,
        /// Output directory
        #[arg(short, long, default_value = ".")]
        outdir: String,
        /// Output prefix
        #[arg(short, long, default_value = "physlr")]
        prefix: String,
        /// K-mer size for minimizer extraction
        #[arg(short, long, default_value_t = 32)]
        k: usize,
        /// Window size for minimizer extraction
        #[arg(short, long, default_value_t = 32)]
        w: usize,
        /// Minimum minimizers per barcode
        #[arg(long, default_value_t = 100)]
        min_bx_count: usize,
        /// Maximum minimizers per barcode
        #[arg(long, default_value_t = 5000)]
        max_bx_count: usize,
        /// Minimum shared minimizers for overlap
        #[arg(long, default_value_t = 10)]
        min_overlap: u32,
        /// Edge filter percentile
        #[arg(long, default_value_t = 90.0)]
        edge_percentile: f64,
        /// Minimum branch size for pruning
        #[arg(long, default_value_t = 10)]
        prune_branches: usize,
        /// Minimum path length
        #[arg(long, default_value_t = 50)]
        min_path_size: usize,
    },

    /// Scaffold a draft assembly using a physical map produced by `physical-map`
    Scaffolds {
        /// Backbone path file (from physical-map output, e.g. physlr.backbone.path)
        path_file: String,
        /// Filtered minimizer TSV (from physical-map output, e.g. physlr.filtered.tsv)
        filtered_mxs: String,
        /// Draft assembly FASTA to scaffold
        draft: String,
        /// Output directory
        #[arg(short, long, default_value = ".")]
        outdir: String,
        /// Output prefix
        #[arg(short, long, default_value = "physlr")]
        prefix: String,
        /// K-mer size (must match the value used in physical-map)
        #[arg(short, long, default_value_t = 32)]
        k: usize,
        /// Window size (must match the value used in physical-map)
        #[arg(short, long, default_value_t = 32)]
        w: usize,
        /// Minimum mapping score
        #[arg(long, default_value_t = 10)]
        min_map_score: u32,
        /// Gap size in scaffolds (Ns between contigs)
        #[arg(long, default_value_t = 100)]
        gap_size: usize,
        /// Expected genome size in bp (optional, for NG50 reporting only)
        #[arg(short = 'g', long)]
        genome_size: Option<u64>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Configure logging
    let log_level = match cli.verbose {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        _ => log::LevelFilter::Debug,
    };
    env_logger::Builder::new()
        .filter_level(log_level)
        .format(|buf, record| writeln!(buf, "[physlr] {}", record.args()))
        .init();

    // Configure thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .ok();

    let timer = Timer::new();

    match cli.command {
        Commands::Index {
            input,
            output,
            k,
            w,
            indexer,
            threads,
            repeat_bf,
        } => {
            // Load repeat Bloom filter if provided
            let bf = if let Some(ref bf_path) = repeat_bf {
                timer.log(&format!("Loading repeat Bloom filter from {}", bf_path));
                let bf = physlr::repeat::BloomFilter::load(bf_path)?;
                timer.log(&format!(
                    "Loaded BF: {} bytes, FPR={:.4}%",
                    bf.size_bytes(),
                    bf.fpr() * 100.0
                ));
                Some(bf)
            } else {
                None
            };

            timer.log(&format!(
                "Indexing minimizers from {} file(s) (k={}, w={}, indexer={}, repeat_filter={})",
                input.len(),
                k,
                w,
                indexer,
                repeat_bf.is_some()
            ));

            let mut bx_to_mxs = rustc_hash::FxHashMap::default();
            for file in &input {
                timer.log(&format!("Processing {}", file));
                let file_mxs = if let Some(ref bf) = bf {
                    physlr::minimizer::index_file_with_repeat_filter(file, k, w, bf)?
                } else {
                    match indexer.as_str() {
                        "builtin" => physlr::minimizer::index_file(file, k, w)?,
                        "btllib" => physlr::minimizer::index_file_btllib(file, k, w, threads)?,
                        other => {
                            anyhow::bail!("Unknown indexer '{}'. Use 'builtin' or 'btllib'.", other)
                        }
                    }
                };
                // Merge: union minimizer sets per barcode
                for (barcode, mxs) in file_mxs {
                    bx_to_mxs
                        .entry(barcode)
                        .or_insert_with(rustc_hash::FxHashSet::default)
                        .extend(mxs);
                }
            }
            timer.log(&format!(
                "Indexed {} barcodes from {} file(s)",
                bx_to_mxs.len(),
                input.len()
            ));

            let mut writer = physlr::io::open_writer(&output)?;
            physlr::minimizer::write_minimizer_tsv(&bx_to_mxs, &mut *writer)?;
            timer.log("Done indexing");
        }

        Commands::RepeatFilter {
            input,
            output,
            k,
            multiplier,
            bf_size,
            cms_size,
            histogram,
        } => {
            timer.log(&format!(
                "Detecting repetitive k-mers from {} file(s) (k={}, multiplier={}, CMS={:.1}GB)",
                input.len(),
                k,
                multiplier,
                cms_size as f64 / 1_073_741_824.0
            ));

            let paths: Vec<&str> = input.iter().map(|s| s.as_str()).collect();
            let (bf, hist) =
                physlr::repeat::detect_repeats(&paths, k, multiplier, bf_size, cms_size)?;

            timer.log(&format!(
                "Repeat filter: BF size={} bytes, popcount={}, FPR={:.4}%",
                bf.size_bytes(),
                bf.popcount(),
                bf.fpr() * 100.0
            ));

            bf.save(&output)?;
            timer.log(&format!("Saved Bloom filter to {}", output));

            if let Some(hist_path) = histogram {
                physlr::repeat::write_histogram(&hist, &hist_path)?;
                timer.log(&format!("Saved histogram to {}", hist_path));
            }
        }

        Commands::IndexContigs {
            input,
            output,
            k,
            w,
        } => {
            timer.log(&format!(
                "Indexing contigs from {} (k={}, w={})",
                input, k, w
            ));
            let contig_mxs = physlr::minimizer::index_file_ordered(&input, k, w)?;
            timer.log(&format!("Indexed {} contigs", contig_mxs.len()));

            let mut writer = physlr::io::open_writer(&output)?;
            physlr::minimizer::write_minimizer_list_tsv(&contig_mxs, &mut *writer)?;
            timer.log("Done indexing contigs");
        }

        Commands::FilterMinimizers {
            input,
            output,
            min_count,
            max_count,
            max_freq,
        } => {
            timer.log(&format!("Reading minimizers from {}", input));
            let mut bx_to_mxs = physlr::io::read_minimizers(&input)?;
            timer.log(&format!("Read {} barcodes", bx_to_mxs.len()));

            physlr::minimizer::remove_singletons(&mut bx_to_mxs);
            physlr::minimizer::filter_barcodes(&mut bx_to_mxs, min_count, max_count);
            physlr::minimizer::remove_singletons(&mut bx_to_mxs);

            let max_freq_opt = if max_freq == 0 { None } else { Some(max_freq) };
            physlr::minimizer::remove_repetitive(&mut bx_to_mxs, max_freq_opt);

            let mut writer = physlr::io::open_writer(&output)?;
            physlr::io::write_minimizers(&bx_to_mxs, &mut *writer)?;
            timer.log("Done filtering minimizers");
        }

        Commands::Overlap {
            input,
            output,
            min_shared,
        } => {
            timer.log(&format!("Reading minimizers from {}", input));
            let bx_to_mxs = physlr::io::read_minimizers(&input)?;
            timer.log(&format!("Read {} barcodes", bx_to_mxs.len()));

            let g = physlr::overlap::compute_overlap(&bx_to_mxs, min_shared);

            let mut writer = physlr::io::open_writer(&output)?;
            physlr::io::write_graph_tsv(&g, &mut *writer)?;
            timer.log("Done computing overlaps");
        }

        Commands::FilterOverlap {
            input,
            output,
            percentile,
        } => {
            timer.log(&format!("Reading graph from {}", input));
            let mut g = physlr::io::read_graph_tsv(&input)?;
            timer.log(&format!(
                "Read graph: V={} E={}",
                g.num_vertices(),
                g.num_edges()
            ));

            physlr::overlap::filter_edges_by_percentile(&mut g, percentile);

            let mut writer = physlr::io::open_writer(&output)?;
            physlr::io::write_graph_tsv(&g, &mut *writer)?;
            timer.log("Done filtering overlap graph");
        }

        Commands::Molecules {
            input,
            output,
            strategy,
            sqcos_threshold,
            skip_small,
            bin_max_size,
            merge_cutoff,
            seed,
        } => {
            timer.log(&format!("Reading graph from {}", input));
            let g = physlr::io::read_graph_tsv(&input)?;
            timer.log(&format!(
                "Read graph: V={} E={}",
                g.num_vertices(),
                g.num_edges()
            ));

            let params = physlr::molecules::MoleculeParams {
                sqcos_threshold,
                skip_small,
                bin_max_size,
                merge_cutoff,
                seed,
            };
            let junctions = rustc_hash::FxHashSet::default();
            let mol_g = physlr::molecules::separate_molecules_with_params(
                &g, &strategy, &junctions, &params,
            );

            let mut writer = physlr::io::open_writer(&output)?;
            physlr::io::write_graph_tsv(&mol_g, &mut *writer)?;
            timer.log("Done separating molecules");
        }

        Commands::SplitMinimizers {
            mol_graph,
            minimizers,
            output,
        } => {
            timer.log(&format!("Reading molecule graph from {}", mol_graph));
            let g = physlr::io::read_graph_tsv(&mol_graph)?;
            timer.log(&format!(
                "Read molecule graph: V={} E={}",
                g.num_vertices(),
                g.num_edges()
            ));

            timer.log(&format!("Reading minimizers from {}", minimizers));
            let bx_mxs = physlr::io::read_minimizers_list(&minimizers)?;
            timer.log(&format!("Read minimizers for {} barcodes", bx_mxs.len()));

            let split = physlr::molecules::split_minimizers(&g, &bx_mxs);
            timer.log(&format!("Split minimizers for {} molecules", split.len()));

            let mut writer = physlr::io::open_writer(&output)?;
            for (mol_name, mxs) in &split {
                write!(writer, "{}\t", mol_name)?;
                for (i, mx) in mxs.iter().enumerate() {
                    if i > 0 {
                        write!(writer, " ")?;
                    }
                    write!(writer, "{}", mx)?;
                }
                writeln!(writer)?;
            }
            timer.log("Done splitting minimizers");
        }

        Commands::TraceMolecules {
            input,
            barcodes,
            strategy,
            sqcos_threshold,
            skip_small,
            bin_max_size,
            merge_cutoff,
        } => {
            timer.log(&format!("Reading graph from {}", input));
            let g = physlr::io::read_graph_tsv(&input)?;
            timer.log(&format!(
                "Read graph: V={} E={}",
                g.num_vertices(),
                g.num_edges()
            ));

            let params = physlr::molecules::MoleculeParams {
                sqcos_threshold,
                skip_small,
                bin_max_size,
                merge_cutoff,
                seed: 42,
            };

            let target_barcodes: Vec<String> = if barcodes == "top5" {
                // Pick 5 highest-degree barcodes
                let mut degrees: Vec<(String, usize)> = g
                    .graph
                    .node_indices()
                    .map(|ni| {
                        let name = g.names.get_name(ni).unwrap().to_string();
                        let deg = g.graph.neighbors(ni).count();
                        (name, deg)
                    })
                    .collect();
                degrees.sort_by_key(|(_, d)| std::cmp::Reverse(*d));
                degrees.into_iter().take(5).map(|(n, _)| n).collect()
            } else {
                barcodes.split(',').map(|s| s.trim().to_string()).collect()
            };

            physlr::molecules::trace_molecules(&g, &strategy, &params, &target_barcodes);
            timer.log("Done tracing molecules");
        }

        Commands::Backbone {
            input,
            output,
            prune_branches,
            prune_bridges,
            prune_junctions,
            min_component_size,
        } => {
            timer.log(&format!("Reading graph from {}", input));
            let g = physlr::io::read_graph_tsv(&input)?;
            timer.log(&format!(
                "Read graph: V={} E={}",
                g.num_vertices(),
                g.num_edges()
            ));

            let config = BackboneConfig {
                prune_branch_size: prune_branches,
                prune_bridge_size: prune_bridges,
                prune_junction_size: prune_junctions,
                min_path_size: min_component_size,
            };
            let paths = physlr::backbone::extract_named_backbones(&g, &config);

            let mut writer = physlr::io::open_writer(&output)?;
            physlr::io::write_paths(&paths, &mut *writer)?;
            timer.log("Done extracting backbones");
        }

        Commands::MergePaths {
            path_file,
            split_mxs,
            output,
            endpoint_depth,
            min_shared_mx,
            min_bridges,
            max_connections,
            min_path_size,
            max_links_per_endpoint,
            min_bridge_density,
            min_endpoint_hits,
        } => {
            timer.log("Loading backbone paths and split minimizers...");
            let paths = physlr::io::read_paths(&path_file)?;
            timer.log(&format!("Read {} backbone paths", paths.len()));

            let mxs = physlr::io::read_minimizers(&split_mxs)?;
            timer.log(&format!(
                "Read split minimizers for {} molecules",
                mxs.len()
            ));

            let config = physlr::backbone::MergePathsConfig {
                endpoint_depth,
                min_shared_mx,
                min_bridges,
                max_path_connections: max_connections,
                min_path_size,
                max_links_per_endpoint,
                min_bridge_density,
                min_endpoint_hits,
            };

            let merged = physlr::backbone::merge_paths(&paths, &mxs, &config);

            let mut writer = physlr::io::open_writer(&output)?;
            physlr::io::write_paths(&merged, &mut *writer)?;
            timer.log("Done merging paths");
        }

        Commands::Map {
            path_file,
            target_mxs,
            query_mxs,
            output,
            min_score,
        } => {
            timer.log("Loading data for mapping...");
            let backbones = physlr::io::read_paths(&path_file)?;
            let mol_mxs = physlr::io::read_minimizers(&target_mxs)?;
            let query = physlr::io::read_minimizers_list(&query_mxs)?;

            let mx_to_pos = physlr::map::index_backbone_minimizers(&backbones, &mol_mxs);
            let mappings = physlr::map::map_to_backbone(&query, &mx_to_pos, min_score);

            let mut writer = physlr::io::open_writer(&output)?;
            for m in &mappings {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    m.tid,
                    m.tpos,
                    m.tpos + 1,
                    m.qname,
                    m.score,
                    m.orientation
                )?;
            }
            timer.log("Done mapping");
        }

        Commands::MapPaf {
            path_file,
            target_mxs,
            query_mxs,
            output,
            min_score,
            coef,
            top_n,
            mx_type,
        } => {
            timer.log("Loading data for PAF mapping...");
            let backbones = physlr::map::read_backbone_paths(&path_file)?;
            let mol_mxs = if mx_type == "split" {
                // Split minimizers: molecule names are keys directly
                physlr::io::read_minimizers(&target_mxs)?
            } else {
                // Barcode minimizers: need to match molecule names to barcodes
                physlr::io::read_minimizers(&target_mxs)?
            };
            let query = physlr::io::read_minimizers_list(&query_mxs)?;

            timer.log(&format!(
                "Loaded {} backbone paths, {} molecules (mx_type={}), {} query sequences",
                backbones.len(),
                mol_mxs.len(),
                mx_type,
                query.len()
            ));

            let mx_to_pos = physlr::map::index_backbone_minimizers(&backbones, &mol_mxs);
            let records = physlr::map::map_paf(&query, &mx_to_pos, &backbones, min_score, coef);

            let mut writer = physlr::io::open_writer(&output)?;
            if top_n > 0 {
                let indices = physlr::map::filter_paf_top_n(&records, top_n);
                for i in indices {
                    records[i].write_to(&mut *writer)?;
                }
            } else {
                for r in &records {
                    r.write_to(&mut *writer)?;
                }
            }
            timer.log(&format!("Wrote {} PAF records", records.len()));
        }

        Commands::BedToPath {
            input,
            output,
            min_score,
        } => {
            timer.log(&format!("Reading BED from {}", input));
            let bed = physlr::io::read_bed(&input)?;
            let mappings: Vec<physlr::map::Mapping> = bed
                .into_iter()
                .map(|r| physlr::map::Mapping {
                    tid: r.tname.parse().unwrap_or(0),
                    tpos: r.tstart as usize,
                    qname: r.qname,
                    score: r.score,
                    orientation: r.orientation,
                })
                .collect();

            let paths = physlr::map::bed_to_scaffold_paths(&mappings, min_score);
            let named_paths: Vec<Vec<String>> = paths
                .into_iter()
                .map(|p| {
                    p.into_iter()
                        .map(|(name, ori)| format!("{}{}", name, ori))
                        .collect()
                })
                .collect();

            let mut writer = physlr::io::open_writer(&output)?;
            physlr::io::write_paths(&named_paths, &mut *writer)?;
            timer.log("Done converting BED to paths");
        }

        Commands::PathToFasta {
            fasta,
            path_file,
            output,
            gap_size,
            min_length,
        } => {
            timer.log("Loading data for scaffolding...");
            let seqs = physlr::io::read_fasta(&fasta)?;
            let paths = physlr::io::read_paths(&path_file)?;

            let scaffold_paths: Vec<Vec<(String, char)>> = paths
                .into_iter()
                .map(|p| {
                    p.into_iter()
                        .map(|s| {
                            let last = s.chars().last().unwrap_or('.');
                            if last == '+' || last == '-' || last == '.' {
                                (s[..s.len() - 1].to_string(), last)
                            } else {
                                (s, '+')
                            }
                        })
                        .collect()
                })
                .collect();

            let config = ScaffoldConfig {
                gap_size,
                min_score: 0,
                min_length,
            };
            let scaffolds = physlr::scaffold::scaffold_assembly(&seqs, &scaffold_paths, &config);

            let mut writer = physlr::io::open_writer(&output)?;
            physlr::io::write_fasta(&scaffolds, &mut *writer)?;
            timer.log("Done producing scaffolded FASTA");
        }

        Commands::Metrics {
            input,
            genome_size,
            label,
        } => {
            let seqs = physlr::io::read_fasta_ordered(&input)?;
            let metrics = physlr::report::compute_assembly_metrics(&seqs, genome_size);
            let mut stdout = std::io::stdout();
            physlr::report::write_metrics_tsv(&metrics, &label, &mut stdout)?;
        }

        Commands::PathMetrics {
            input,
            expected_molecules,
            min_component_size,
        } => {
            let paths = physlr::io::read_paths(&input)?;
            let filtered: Vec<Vec<String>> = paths
                .into_iter()
                .filter(|p| p.len() >= min_component_size)
                .collect();
            let metrics =
                physlr::report::compute_physical_map_metrics(&filtered, expected_molecules);
            let mut stdout = std::io::stdout();
            physlr::report::write_physical_map_metrics_tsv(&metrics, &mut stdout)?;
        }

        Commands::BackboneDot { input, output } => {
            let paths = physlr::io::read_paths(&input)?;
            let mut writer = physlr::io::open_writer(&output)?;
            physlr::report::backbone_to_dot(&paths, &mut *writer)?;
        }

        Commands::Pipeline {
            input,
            draft,
            outdir,
            prefix,
            k,
            w,
            min_bx_count,
            max_bx_count,
            min_overlap,
            edge_percentile,
            prune_branches,
            min_path_size,
            min_map_score,
            gap_size,
            genome_size,
        } => {
            std::fs::create_dir_all(&outdir)?;

            // ── Phase 1: Physical map ───────────────────────────────────
            timer.log("=== Phase 1: Building physical map ===");

            timer.log(&format!(
                "Step 1: Indexing minimizers from {} file(s) (k={}, w={})...",
                input.len(),
                k,
                w
            ));
            let mut bx_to_mxs = rustc_hash::FxHashMap::default();
            for file in &input {
                timer.log(&format!("  Processing {}", file));
                let file_mxs = physlr::minimizer::index_file(file, k, w)?;
                for (bx, mxs) in file_mxs {
                    bx_to_mxs
                        .entry(bx)
                        .or_insert_with(rustc_hash::FxHashSet::default)
                        .extend(mxs);
                }
            }
            timer.log(&format!("Indexed {} barcodes", bx_to_mxs.len()));

            timer.log("Step 2: Filtering minimizers...");
            physlr::minimizer::remove_singletons(&mut bx_to_mxs);
            physlr::minimizer::filter_barcodes(&mut bx_to_mxs, min_bx_count, max_bx_count);
            physlr::minimizer::remove_singletons(&mut bx_to_mxs);
            physlr::minimizer::remove_repetitive(&mut bx_to_mxs, None);

            let filtered_path = format!("{}/{}.filtered.tsv", outdir, prefix);
            let mut writer = physlr::io::open_writer(&filtered_path)?;
            physlr::io::write_minimizers(&bx_to_mxs, &mut *writer)?;
            drop(writer);

            timer.log("Step 3: Computing overlaps...");
            let mut g = physlr::overlap::compute_overlap(&bx_to_mxs, min_overlap);

            timer.log("Step 4: Filtering edges...");
            physlr::overlap::filter_edges_by_percentile(&mut g, edge_percentile);

            timer.log("Step 5: Separating molecules...");
            let mol_g = physlr::molecules::separate_molecules(
                &g,
                "bc+cc",
                &rustc_hash::FxHashSet::default(),
            );

            timer.log("Step 6: Extracting backbone paths...");
            let config = BackboneConfig {
                prune_branch_size: prune_branches,
                prune_bridge_size: 10,
                prune_junction_size: 200,
                min_path_size,
            };
            let backbones = physlr::backbone::extract_named_backbones(&mol_g, &config);

            let backbone_path = format!("{}/{}.backbone.path", outdir, prefix);
            let mut writer = physlr::io::open_writer(&backbone_path)?;
            physlr::io::write_paths(&backbones, &mut *writer)?;
            drop(writer);

            let pm_metrics = physlr::report::compute_physical_map_metrics(&backbones, None);
            timer.log(&format!(
                "Physical map: {} paths, {} total molecules",
                pm_metrics.num_paths, pm_metrics.total_molecules
            ));

            // ── Phase 2: Scaffolding ────────────────────────────────────
            timer.log("=== Phase 2: Scaffolding draft assembly ===");

            timer.log("Indexing draft assembly contigs...");
            let query_mxs = physlr::minimizer::index_file_ordered(&draft, k, w)?;
            timer.log(&format!("Indexed {} contigs", query_mxs.len()));

            timer.log("Mapping draft assembly to physical map...");
            let mx_to_pos = physlr::map::index_backbone_minimizers(&backbones, &bx_to_mxs);
            let mappings = physlr::map::map_to_backbone(&query_mxs, &mx_to_pos, min_map_score);
            let scaffold_paths = physlr::map::bed_to_scaffold_paths(&mappings, min_map_score);

            timer.log("Producing scaffolded assembly...");
            let draft_seqs = physlr::io::read_fasta(&draft)?;
            let scaffold_config = ScaffoldConfig {
                gap_size,
                min_score: min_map_score,
                min_length: 0,
            };
            let scaffolds =
                physlr::scaffold::scaffold_assembly(&draft_seqs, &scaffold_paths, &scaffold_config);

            let scaffolds_path = format!("{}/{}.scaffolds.fa", outdir, prefix);
            let mut writer = physlr::io::open_writer(&scaffolds_path)?;
            physlr::io::write_fasta(&scaffolds, &mut *writer)?;
            drop(writer);

            let draft_ordered = physlr::io::read_fasta_ordered(&draft)?;
            let before_metrics =
                physlr::report::compute_assembly_metrics(&draft_ordered, genome_size);
            let after_metrics = physlr::report::compute_assembly_metrics(&scaffolds, genome_size);

            let report_path = format!("{}/{}.report.json", outdir, prefix);
            let mut writer = physlr::io::open_writer(&report_path)?;
            physlr::report::write_json_report(
                &pm_metrics,
                Some(&before_metrics),
                Some(&after_metrics),
                &mut *writer,
            )?;
            drop(writer);

            let mut stdout = std::io::stdout();
            writeln!(stdout, "\n=== Before Scaffolding ===")?;
            physlr::report::write_metrics_tsv(&before_metrics, "draft", &mut stdout)?;
            writeln!(stdout, "\n=== After Scaffolding ===")?;
            physlr::report::write_metrics_tsv(&after_metrics, "physlr", &mut stdout)?;

            timer.log("Pipeline complete");
        }

        Commands::PhysicalMap {
            input,
            outdir,
            prefix,
            k,
            w,
            min_bx_count,
            max_bx_count,
            min_overlap,
            edge_percentile,
            prune_branches,
            min_path_size,
        } => {
            std::fs::create_dir_all(&outdir)?;

            timer.log(&format!(
                "Step 1: Indexing minimizers from {} file(s) (k={}, w={})...",
                input.len(),
                k,
                w
            ));
            let mut bx_to_mxs = rustc_hash::FxHashMap::default();
            for file in &input {
                timer.log(&format!("  Processing {}", file));
                let file_mxs = physlr::minimizer::index_file(file, k, w)?;
                for (bx, mxs) in file_mxs {
                    bx_to_mxs
                        .entry(bx)
                        .or_insert_with(rustc_hash::FxHashSet::default)
                        .extend(mxs);
                }
            }
            timer.log(&format!("Indexed {} barcodes", bx_to_mxs.len()));

            timer.log("Step 2: Filtering minimizers...");

            physlr::minimizer::remove_singletons(&mut bx_to_mxs);
            physlr::minimizer::filter_barcodes(&mut bx_to_mxs, min_bx_count, max_bx_count);
            physlr::minimizer::remove_singletons(&mut bx_to_mxs);
            physlr::minimizer::remove_repetitive(&mut bx_to_mxs, None);

            let filtered_path = format!("{}/{}.filtered.tsv", outdir, prefix);
            let mut writer = physlr::io::open_writer(&filtered_path)?;
            physlr::io::write_minimizers(&bx_to_mxs, &mut *writer)?;
            drop(writer);
            timer.log(&format!("Wrote filtered minimizers to {}", filtered_path));

            timer.log("Step 2: Computing overlaps...");
            let mut g = physlr::overlap::compute_overlap(&bx_to_mxs, min_overlap);

            timer.log("Step 3: Filtering edges...");
            physlr::overlap::filter_edges_by_percentile(&mut g, edge_percentile);

            let overlap_path = format!("{}/{}.overlap.tsv", outdir, prefix);
            let mut writer = physlr::io::open_writer(&overlap_path)?;
            physlr::io::write_graph_tsv(&g, &mut *writer)?;
            drop(writer);
            timer.log(&format!("Wrote overlap graph to {}", overlap_path));

            timer.log("Step 4: Separating molecules...");
            let mol_g = physlr::molecules::separate_molecules(
                &g,
                "bc+cc",
                &rustc_hash::FxHashSet::default(),
            );

            let mol_path = format!("{}/{}.mol.tsv", outdir, prefix);
            let mut writer = physlr::io::open_writer(&mol_path)?;
            physlr::io::write_graph_tsv(&mol_g, &mut *writer)?;
            drop(writer);

            timer.log("Step 5: Extracting backbone paths...");
            let config = BackboneConfig {
                prune_branch_size: prune_branches,
                prune_bridge_size: 10,
                prune_junction_size: 200,
                min_path_size,
            };
            let paths = physlr::backbone::extract_named_backbones(&mol_g, &config);

            let backbone_path = format!("{}/{}.backbone.path", outdir, prefix);
            let mut writer = physlr::io::open_writer(&backbone_path)?;
            physlr::io::write_paths(&paths, &mut *writer)?;
            drop(writer);

            let metrics = physlr::report::compute_physical_map_metrics(&paths, None);
            let metrics_path = format!("{}/{}.backbone.metrics.tsv", outdir, prefix);
            let mut writer = physlr::io::open_writer(&metrics_path)?;
            physlr::report::write_physical_map_metrics_tsv(&metrics, &mut *writer)?;
            drop(writer);

            let dot_path = format!("{}/{}.backbone.dot", outdir, prefix);
            let mut writer = physlr::io::open_writer(&dot_path)?;
            physlr::report::backbone_to_dot(&paths, &mut *writer)?;
            drop(writer);

            timer.log(&format!(
                "Physical map complete: {} paths, {} total molecules",
                metrics.num_paths, metrics.total_molecules
            ));
        }

        Commands::Scaffolds {
            path_file,
            filtered_mxs,
            draft,
            outdir,
            prefix,
            k,
            w,
            min_map_score,
            gap_size,
            genome_size,
        } => {
            std::fs::create_dir_all(&outdir)?;

            timer.log("Loading physical map...");
            let backbones = physlr::io::read_paths(&path_file)?;
            timer.log(&format!("Read {} backbone paths", backbones.len()));

            let bx_to_mxs = physlr::io::read_minimizers(&filtered_mxs)?;
            timer.log(&format!("Read {} filtered barcodes", bx_to_mxs.len()));

            timer.log("Indexing draft assembly contigs...");
            let query_mxs = physlr::minimizer::index_file_ordered(&draft, k, w)?;
            timer.log(&format!("Indexed {} contigs", query_mxs.len()));

            timer.log("Mapping draft assembly to physical map...");
            let mx_to_pos = physlr::map::index_backbone_minimizers(&backbones, &bx_to_mxs);
            let mappings = physlr::map::map_to_backbone(&query_mxs, &mx_to_pos, min_map_score);

            let scaffold_paths = physlr::map::bed_to_scaffold_paths(&mappings, min_map_score);

            timer.log("Producing scaffolded assembly...");
            let draft_seqs = physlr::io::read_fasta(&draft)?;
            let scaffold_config = ScaffoldConfig {
                gap_size,
                min_score: min_map_score,
                min_length: 0,
            };
            let scaffolds =
                physlr::scaffold::scaffold_assembly(&draft_seqs, &scaffold_paths, &scaffold_config);

            let scaffolds_path = format!("{}/{}.scaffolds.fa", outdir, prefix);
            let mut writer = physlr::io::open_writer(&scaffolds_path)?;
            physlr::io::write_fasta(&scaffolds, &mut *writer)?;
            drop(writer);

            let draft_ordered = physlr::io::read_fasta_ordered(&draft)?;
            let before_metrics =
                physlr::report::compute_assembly_metrics(&draft_ordered, genome_size);
            let after_metrics = physlr::report::compute_assembly_metrics(&scaffolds, genome_size);

            let pm_metrics = physlr::report::compute_physical_map_metrics(&backbones, None);

            let report_path = format!("{}/{}.report.json", outdir, prefix);
            let mut writer = physlr::io::open_writer(&report_path)?;
            physlr::report::write_json_report(
                &pm_metrics,
                Some(&before_metrics),
                Some(&after_metrics),
                &mut *writer,
            )?;
            drop(writer);

            let mut stdout = std::io::stdout();
            writeln!(stdout, "\n=== Before Scaffolding ===")?;
            physlr::report::write_metrics_tsv(&before_metrics, "draft", &mut stdout)?;
            writeln!(stdout, "\n=== After Scaffolding ===")?;
            physlr::report::write_metrics_tsv(&after_metrics, "physlr", &mut stdout)?;

            timer.log("Scaffolding complete");
        }
    }

    Ok(())
}
