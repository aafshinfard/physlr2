/// High-level pipeline functions for linked-read and long-read physical maps.
///
/// These are called from the CLI handlers in main.rs. Each protocol has its own
/// pipeline function that shares core algorithms (overlap, backbone) but differs
/// in preprocessing and coordinate handling.
use anyhow::Result;
use rustc_hash::FxHashSet;

use crate::backbone::BackboneConfig;
use crate::Timer;

/// Shared parameters for the physical map pipeline.
pub struct PhysicalMapParams {
    pub outdir: String,
    pub prefix: String,
    pub k: usize,
    pub w: usize,
    pub threads: usize,
    pub bf_size: u64,
    pub makebf_path: Option<String>,
    pub min_bx_count: usize,
    pub max_bx_count: usize,
    pub min_overlap: u32,
    pub edge_percentile: f64,
    pub prune_branches: usize,
    pub min_path_size: usize,
    pub reference: Option<String>,
}

/// Result of a physical map pipeline run.
pub struct PhysicalMapResult {
    pub backbone_path: String,
    pub num_paths: usize,
    pub total_molecules: usize,
}

/// Run the physical map pipeline for linked reads.
///
/// Steps: repeat-filter + indexlr → filter → overlap → molecules → backbone
/// Optional: reference mapping + plots (direct, no HPC)
pub fn run_physical_map_linked(
    input: &[String],
    params: &PhysicalMapParams,
    timer: &Timer,
) -> Result<PhysicalMapResult> {
    std::fs::create_dir_all(&params.outdir)?;

    // Steps 1-2: Repeat filtering + indexing
    timer.log("Steps 1-2: Repeat filtering and indexing via ntcard/nthits/indexlr...");
    let work_dir = std::path::Path::new(&params.outdir);
    let mut bx_to_mxs = crate::external::repeat_filter_and_index(
        input,
        params.k,
        params.w,
        params.threads,
        work_dir,
        &params.prefix,
        params.bf_size,
        3,
        false, // not long reads
        params.makebf_path.as_deref(),
    )?;
    timer.log(&format!("Indexed {} barcodes", bx_to_mxs.len()));

    // Steps 3-7: shared core pipeline
    let (paths, bx_to_mxs) = run_core_pipeline(&mut bx_to_mxs, params, timer)?;

    // Optional: reference mapping (direct, no HPC)
    if let Some(ref ref_path) = params.reference {
        map_reference_linked(ref_path, &paths, &bx_to_mxs, params, timer)?;
    }

    let metrics = crate::report::compute_physical_map_metrics(&paths, None);
    Ok(PhysicalMapResult {
        backbone_path: format!("{}/{}.backbone.path", params.outdir, params.prefix),
        num_paths: metrics.num_paths,
        total_molecules: metrics.total_molecules,
    })
}

/// Run the physical map pipeline for long reads.
///
/// Steps: HPC + repeat-filter + indexlr --long → filter → overlap → backbone
/// Optional: reference mapping with HPC coordinate translation
pub fn run_physical_map_long(
    input: &[String],
    params: &PhysicalMapParams,
    timer: &Timer,
) -> Result<PhysicalMapResult> {
    std::fs::create_dir_all(&params.outdir)?;

    // Steps 1-2: HPC + repeat filtering + indexing
    timer.log("Steps 1-2: HPC compression, repeat filtering, and indexing...");
    let work_dir = std::path::Path::new(&params.outdir);
    let mut bx_to_mxs = crate::external::repeat_filter_and_index(
        input,
        params.k,
        params.w,
        params.threads,
        work_dir,
        &params.prefix,
        params.bf_size,
        3,
        true, // long reads
        params.makebf_path.as_deref(),
    )?;
    timer.log(&format!("Indexed {} reads", bx_to_mxs.len()));

    // Steps 3-7: shared core pipeline
    let (paths, bx_to_mxs) = run_core_pipeline(&mut bx_to_mxs, params, timer)?;

    // Optional: reference mapping with HPC coordinate translation
    if let Some(ref ref_path) = params.reference {
        map_reference_long(ref_path, &paths, &bx_to_mxs, params, timer)?;
    }

    let metrics = crate::report::compute_physical_map_metrics(&paths, None);
    Ok(PhysicalMapResult {
        backbone_path: format!("{}/{}.backbone.path", params.outdir, params.prefix),
        num_paths: metrics.num_paths,
        total_molecules: metrics.total_molecules,
    })
}

/// Core pipeline shared between linked and long reads:
/// filter → overlap → edge filter → molecules → backbone
///
/// Returns (backbone_paths, filtered_bx_to_mxs).
fn run_core_pipeline(
    bx_to_mxs: &mut rustc_hash::FxHashMap<String, rustc_hash::FxHashSet<u64>>,
    params: &PhysicalMapParams,
    timer: &Timer,
) -> Result<(
    Vec<Vec<String>>,
    rustc_hash::FxHashMap<String, rustc_hash::FxHashSet<u64>>,
)> {
    timer.log("Step 3: Filtering minimizers...");
    crate::minimizer::remove_singletons(bx_to_mxs);
    crate::minimizer::filter_barcodes(bx_to_mxs, params.min_bx_count, params.max_bx_count);
    crate::minimizer::remove_singletons(bx_to_mxs);
    crate::minimizer::remove_repetitive(bx_to_mxs, None);

    let filtered_path = format!("{}/{}.filtered.tsv", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&filtered_path)?;
    crate::io::write_minimizers(bx_to_mxs, &mut *writer)?;
    drop(writer);
    timer.log(&format!("Wrote filtered minimizers to {}", filtered_path));

    timer.log("Step 4: Computing overlaps...");
    let mut g = crate::overlap::compute_overlap(bx_to_mxs, params.min_overlap);

    timer.log("Step 5: Filtering edges...");
    crate::overlap::filter_edges_by_percentile(&mut g, params.edge_percentile);

    let overlap_path = format!("{}/{}.overlap.tsv", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&overlap_path)?;
    crate::io::write_graph_tsv(&g, &mut *writer)?;
    drop(writer);
    timer.log(&format!("Wrote overlap graph to {}", overlap_path));

    timer.log("Step 6: Separating molecules...");
    let mol_g = crate::molecules::separate_molecules(&g, "bc+cc", &FxHashSet::default());

    let mol_path = format!("{}/{}.mol.tsv", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&mol_path)?;
    crate::io::write_graph_tsv(&mol_g, &mut *writer)?;
    drop(writer);

    timer.log("Step 7: Extracting backbone paths...");
    let config = BackboneConfig {
        prune_branch_size: params.prune_branches,
        prune_bridge_size: 10,
        prune_junction_size: 200,
        min_path_size: params.min_path_size,
    };
    let paths = crate::backbone::extract_named_backbones(&mol_g, &config);

    let backbone_path = format!("{}/{}.backbone.path", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&backbone_path)?;
    crate::io::write_paths(&paths, &mut *writer)?;
    drop(writer);

    let metrics = crate::report::compute_physical_map_metrics(&paths, None);
    let metrics_path = format!("{}/{}.backbone.metrics.tsv", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&metrics_path)?;
    crate::report::write_physical_map_metrics_tsv(&metrics, &mut *writer)?;
    drop(writer);

    let dot_path = format!("{}/{}.backbone.dot", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&dot_path)?;
    crate::report::backbone_to_dot(&paths, &mut *writer)?;
    drop(writer);

    timer.log(&format!(
        "Physical map complete: {} paths, {} total molecules",
        metrics.num_paths, metrics.total_molecules
    ));

    Ok((paths, bx_to_mxs.clone()))
}

/// Map backbones to reference for linked reads (direct, no HPC).
///
/// Uses indexlr to index the reference so minimizer hashes match the reads
/// (both use ntHash via indexlr).
pub fn map_reference_linked(
    ref_path: &str,
    paths: &[Vec<String>],
    bx_to_mxs: &rustc_hash::FxHashMap<String, rustc_hash::FxHashSet<u64>>,
    params: &PhysicalMapParams,
    timer: &Timer,
) -> Result<()> {
    timer.log(&format!(
        "Step 8: Mapping backbones to reference {}...",
        ref_path
    ));

    // Index reference with indexlr (ntHash) to match read minimizers
    let ref_mxs = crate::external::run_indexlr_reference_ordered(
        ref_path, params.k, params.w, params.threads,
    )?;
    timer.log(&format!("Indexed {} reference sequences via indexlr", ref_mxs.len()));

    let ref_tsv_path = format!("{}/{}.ref.tsv", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&ref_tsv_path)?;
    crate::external::write_minimizer_list_tsv_ordered(&ref_mxs, &mut *writer)?;
    drop(writer);

    // Position map via indexlr --pos
    let pos_path = format!("{}/{}.positions.tsv", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&pos_path)?;
    crate::external::run_indexlr_reference_positions(
        ref_path, params.k, params.w, params.threads, 100, None, &mut *writer,
    )?;
    drop(writer);
    timer.log(&format!("Wrote position map to {}", pos_path));

    let mx_to_pos = crate::map::index_backbone_minimizers(paths, bx_to_mxs);
    let paf_records = crate::map::map_paf(&ref_mxs, &mx_to_pos, paths, 10, 1.5);

    let paf_path = format!("{}/{}.backbone.paf", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&paf_path)?;
    for r in &paf_records {
        r.write_to(&mut *writer)?;
    }
    drop(writer);
    timer.log(&format!(
        "Wrote {} PAF records to {}",
        paf_records.len(),
        paf_path
    ));

    run_plotting(&paf_path, &params.outdir, &params.prefix, &pos_path, timer);
    Ok(())
}

/// Map backbones to reference for long reads (HPC + coordinate translation).
///
/// 1. HPC-compress the reference
/// 2. Index HPC reference via indexlr (ntHash, matching reads)
/// 3. Generate position map with translated coordinates
/// 4. Map HPC reference to backbone → PAF in HPC space
/// 5. Translate PAF coordinates back to original bp
/// 6. Plot using original coordinates
pub fn map_reference_long(
    ref_path: &str,
    paths: &[Vec<String>],
    bx_to_mxs: &rustc_hash::FxHashMap<String, rustc_hash::FxHashSet<u64>>,
    params: &PhysicalMapParams,
    timer: &Timer,
) -> Result<()> {
    timer.log(&format!(
        "Step 8: Mapping backbones to reference {} (HPC mode)...",
        ref_path
    ));

    // HPC-compress the reference
    let hpc_ref_path =
        std::path::Path::new(&params.outdir).join(format!("{}.ref.hpc.fa", params.prefix));
    timer.log("HPC-compressing reference...");
    let hpc_index = crate::minimizer::hpc_compress_fasta(ref_path, &hpc_ref_path)?;
    timer.log(&format!(
        "HPC-compressed {} reference sequences",
        hpc_index.sequences.len()
    ));

    // Index HPC reference with indexlr (ntHash) to match read minimizers
    let hpc_ref_str = hpc_ref_path.to_str().unwrap();
    let ref_mxs = crate::external::run_indexlr_reference_ordered(
        hpc_ref_str, params.k, params.w, params.threads,
    )?;
    timer.log(&format!(
        "Indexed {} HPC reference sequences via indexlr",
        ref_mxs.len()
    ));

    // Write reference minimizer TSV
    let ref_tsv_path = format!("{}/{}.ref.tsv", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&ref_tsv_path)?;
    crate::external::write_minimizer_list_tsv_ordered(&ref_mxs, &mut *writer)?;
    drop(writer);

    // Position map via indexlr --pos, with HPC→original coordinate translation
    let pos_path = format!("{}/{}.positions.tsv", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&pos_path)?;
    crate::external::run_indexlr_reference_positions(
        hpc_ref_str, params.k, params.w, params.threads, 100,
        Some(&hpc_index), &mut *writer,
    )?;
    drop(writer);
    timer.log(&format!(
        "Wrote position map (original coords) to {}",
        pos_path
    ));

    // Map HPC reference to backbone paths → PAF in HPC space
    let mx_to_pos = crate::map::index_backbone_minimizers(paths, bx_to_mxs);
    let mut paf_records = crate::map::map_paf(&ref_mxs, &mx_to_pos, paths, 10, 1.5);

    // Translate PAF coordinates from HPC to original space
    crate::map::translate_paf_from_hpc(&mut paf_records, &hpc_index);
    timer.log("Translated PAF coordinates to original space");

    let paf_path = format!("{}/{}.backbone.paf", params.outdir, params.prefix);
    let mut writer = crate::io::open_writer(&paf_path)?;
    for r in &paf_records {
        r.write_to(&mut *writer)?;
    }
    drop(writer);
    timer.log(&format!(
        "Wrote {} PAF records to {}",
        paf_records.len(),
        paf_path
    ));

    run_plotting(&paf_path, &params.outdir, &params.prefix, &pos_path, timer);
    Ok(())
}

/// Try to run plotpaf.py for visualization.
fn run_plotting(paf_path: &str, outdir: &str, prefix: &str, pos_path: &str, timer: &Timer) {
    timer.log("Generating plots...");
    let plot_prefix = format!("{}/{}", outdir, prefix);
    let exe_dir = std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(|d| d.to_path_buf()));
    let script_candidates = [
        exe_dir
            .as_ref()
            .map(|d| d.join("plotpaf.py"))
            .unwrap_or_default(),
        exe_dir
            .as_ref()
            .map(|d| d.join("../scripts/plotpaf.py"))
            .unwrap_or_default(),
        std::path::PathBuf::from("scripts/plotpaf.py"),
        std::path::PathBuf::from("plotpaf.py"),
    ];
    let plot_script = script_candidates.iter().find(|p| p.exists());
    if let Some(script) = plot_script {
        let status = std::process::Command::new("python3")
            .arg(script)
            .arg(paf_path)
            .arg(&plot_prefix)
            .arg("--positions")
            .arg(pos_path)
            .status();
        match status {
            Ok(s) if s.success() => {
                timer.log("Plots generated successfully");
            }
            Ok(s) => {
                eprintln!(
                    "Warning: plotting script exited with {}",
                    s.code().unwrap_or(-1)
                );
            }
            Err(e) => {
                eprintln!(
                    "Warning: could not run plotting script: {}. \
                     You can run it manually:\n  \
                     python3 {} {} {} --positions {}",
                    e,
                    script.display(),
                    paf_path,
                    plot_prefix,
                    pos_path
                );
            }
        }
    } else {
        eprintln!(
            "Note: plotpaf.py not found. To generate plots, run:\n  \
             python3 scripts/plotpaf.py {} {} --positions {}",
            paf_path, plot_prefix, pos_path
        );
    }
}
