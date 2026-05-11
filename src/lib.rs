//! # nanoget-rs
//!
//! A Rust library for extracting metrics from Oxford Nanopore sequencing data and alignments.
//!
//! This library provides functionality to extract useful metrics from:
//! - BAM/SAM/CRAM files (aligned reads)
//! - FASTQ files (with or without metadata)
//! - FASTA files
//! - Sequencing summary files
//!
//! ## Example
//!
//! ```rust,no_run
//! use nanoget_rs::{extract_metrics, FileType, ExtractArgs};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let args = ExtractArgs {
//!     files: vec!["reads.fastq".into()],
//!     file_type: FileType::Fastq,
//!     threads: 4,
//!     output_format: "json".to_string(),
//!     output: None,
//!     read_type: "1D".to_string(),
//!     barcoded: false,
//!     keep_supplementary: true,
//!     combine: "simple".to_string(),
//!     names: None,
//! };
//!
//! let metrics = extract_metrics(&args)?;
//! # Ok(())
//! # }
//! ```

pub mod cli;
pub mod error;
pub mod extract;
pub mod formats;
pub mod metrics;
pub mod utils;

pub use cli::{Cli, Commands, ExtractArgs};
pub use error::NanogetError;
pub use extract::extract_metrics;
pub use formats::FileType;
pub use metrics::{MetricsCollection, MetricsSummary, ReadMetrics, StatsSummary};

/// Convenience functions for common use cases
pub mod convenience {
    use super::*;
    use std::path::Path;

    // Default values as constants to avoid repeated allocations
    const DEFAULT_OUTPUT_FORMAT: &str = "json";
    const DEFAULT_READ_TYPE: &str = "1D";
    const DEFAULT_COMBINE: &str = "simple";
    const DEFAULT_THREADS: usize = 4;

    /// Create default ExtractArgs with the given files and file type
    fn default_args(files: Vec<std::path::PathBuf>, file_type: FileType) -> ExtractArgs {
        ExtractArgs {
            files,
            file_type,
            threads: DEFAULT_THREADS,
            output_format: DEFAULT_OUTPUT_FORMAT.to_string(),
            output: None,
            read_type: DEFAULT_READ_TYPE.to_string(),
            barcoded: false,
            keep_supplementary: true,
            combine: DEFAULT_COMBINE.to_string(),
            names: None,
        }
    }

    /// Extract metrics from a single FASTQ file with default settings
    pub fn extract_from_fastq<P: AsRef<Path>>(file: P) -> Result<MetricsCollection, NanogetError> {
        let args = default_args(vec![file.as_ref().to_path_buf()], FileType::Fastq);
        extract_metrics(&args)
    }

    /// Extract metrics from a single BAM file with default settings
    pub fn extract_from_bam<P: AsRef<Path>>(file: P) -> Result<MetricsCollection, NanogetError> {
        let args = default_args(vec![file.as_ref().to_path_buf()], FileType::Bam);
        extract_metrics(&args)
    }

    /// Extract metrics from a single FASTA file with default settings
    pub fn extract_from_fasta<P: AsRef<Path>>(file: P) -> Result<MetricsCollection, NanogetError> {
        let args = default_args(vec![file.as_ref().to_path_buf()], FileType::Fasta);
        extract_metrics(&args)
    }

    /// Extract metrics from multiple files with automatic format detection
    pub fn extract_from_files<P: AsRef<Path>>(
        files: Vec<P>,
        file_type: FileType,
        threads: Option<usize>,
    ) -> Result<MetricsCollection, NanogetError> {
        let mut args = default_args(
            files.into_iter().map(|p| p.as_ref().to_path_buf()).collect(),
            file_type,
        );
        if let Some(t) = threads {
            args.threads = t;
        }
        extract_metrics(&args)
    }
}
