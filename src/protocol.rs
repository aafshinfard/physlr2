/// Sequencing protocol detection and configuration.
///
/// Physlr supports two sequencing protocols:
/// - **Linked reads** (10x Genomics, stLFR): multiple short reads share a barcode tag.
///   Minimizers per barcode form an unordered set.
/// - **Long reads** (ONT, PacBio): each read is its own molecule.
///   The read name serves as the barcode. Minimizers are ordered.
///
/// The `--protocol` flag controls which mode is used:
/// - `auto`: inspect the first N reads to detect barcodes
/// - `linked`: force linked-read mode
/// - `long`: force long-read mode

use std::fmt;

/// Sequencing protocol.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Protocol {
    Linked,
    Long,
}

impl fmt::Display for Protocol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Protocol::Linked => write!(f, "linked"),
            Protocol::Long => write!(f, "long"),
        }
    }
}

/// Default parameters that differ between protocols.
#[derive(Debug, Clone)]
pub struct ProtocolDefaults {
    pub min_count: usize,
    pub max_count: usize,
    pub min_overlap: u32,
    pub homopolymer_compress: bool,
}

impl ProtocolDefaults {
    pub fn for_protocol(protocol: Protocol) -> Self {
        match protocol {
            Protocol::Linked => Self {
                min_count: 100,
                max_count: 5000,
                min_overlap: 10,
                homopolymer_compress: false,
            },
            Protocol::Long => Self {
                min_count: 50,
                max_count: 50000,
                min_overlap: 5,
                homopolymer_compress: true,
            },
        }
    }
}

/// Detect the sequencing protocol by inspecting the first `sample_size` reads.
///
/// Heuristic: if ≥10% of sampled reads have a BX:Z: tag or stLFR barcode pattern,
/// classify as linked reads. Otherwise, classify as long reads.
pub fn detect_protocol(path: &str, sample_size: usize) -> anyhow::Result<Protocol> {
    let mut reader = needletail::parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot open {}: {}", path, e))?;

    let mut total = 0usize;
    let mut with_barcode = 0usize;

    while let Some(record) = reader.next() {
        let record = record?;
        let header = std::str::from_utf8(record.id()).unwrap_or("");

        if has_barcode_tag(header) {
            with_barcode += 1;
        }

        total += 1;
        if total >= sample_size {
            break;
        }
    }

    if total == 0 {
        anyhow::bail!("No reads found in {}", path);
    }

    let fraction = with_barcode as f64 / total as f64;
    let protocol = if fraction >= 0.10 {
        Protocol::Linked
    } else {
        Protocol::Long
    };

    log::info!(
        "Auto-detected protocol: {} ({}/{} reads with barcodes, {:.1}%)",
        protocol,
        with_barcode,
        total,
        fraction * 100.0
    );

    Ok(protocol)
}

/// Check if a FASTQ header contains a barcode tag (BX:Z: or stLFR #barcode).
fn has_barcode_tag(header: &str) -> bool {
    // 10x Genomics: BX:Z:ACGTACGT-1
    for part in header.split_whitespace() {
        if part.starts_with("BX:Z:") {
            return true;
        }
    }

    // stLFR: readname#barcode (not #0 or #0_0_0)
    if let Some(pos) = header.find('#') {
        let bx = &header[pos + 1..];
        let bx = bx.split_whitespace().next().unwrap_or(bx);
        if !bx.is_empty() && bx != "0" && bx != "0_0_0" {
            return true;
        }
    }

    false
}

/// Resolve the protocol from the CLI `--protocol` value and input files.
pub fn resolve_protocol(protocol_str: &str, input_files: &[String]) -> anyhow::Result<Protocol> {
    match protocol_str {
        "linked" => Ok(Protocol::Linked),
        "long" => Ok(Protocol::Long),
        "auto" => {
            if input_files.is_empty() {
                anyhow::bail!("No input files provided for auto-detection");
            }
            detect_protocol(&input_files[0], 1000)
        }
        other => anyhow::bail!(
            "Unknown protocol '{}'. Use 'auto', 'linked', or 'long'.",
            other
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_has_barcode_10x() {
        assert!(has_barcode_tag("read1 BX:Z:ACGTACGT-1 RG:Z:foo"));
    }

    #[test]
    fn test_has_barcode_stlfr() {
        assert!(has_barcode_tag("read1#123_456_789"));
    }

    #[test]
    fn test_no_barcode_plain() {
        assert!(!has_barcode_tag("read1 length=5000"));
    }

    #[test]
    fn test_no_barcode_stlfr_unassigned() {
        assert!(!has_barcode_tag("read1#0"));
        assert!(!has_barcode_tag("read1#0_0_0"));
    }

    #[test]
    fn test_protocol_defaults_linked() {
        let d = ProtocolDefaults::for_protocol(Protocol::Linked);
        assert_eq!(d.min_count, 100);
        assert!(!d.homopolymer_compress);
    }

    #[test]
    fn test_protocol_defaults_long() {
        let d = ProtocolDefaults::for_protocol(Protocol::Long);
        assert_eq!(d.min_count, 50);
        assert!(d.homopolymer_compress);
    }

    #[test]
    fn test_resolve_protocol_explicit() {
        assert_eq!(
            resolve_protocol("linked", &[]).unwrap(),
            Protocol::Linked
        );
        assert_eq!(
            resolve_protocol("long", &[]).unwrap(),
            Protocol::Long
        );
    }

    #[test]
    fn test_resolve_protocol_invalid() {
        assert!(resolve_protocol("invalid", &[]).is_err());
    }
}
