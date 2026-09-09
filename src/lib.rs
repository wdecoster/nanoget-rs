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
//! ## Storage
//!
//! Metrics are stored columnar: one array per field rather than one struct per read. A
//! column an input format never populates is never allocated — a plain FASTQ sets 3 of
//! the 13 fields — and consumers that want a field get a slice rather than a copy.
//! [`ReadMetrics`] remains the row type, built transiently while parsing and never
//! stored; [`ReadView`] borrows one read back out. See [`columns`] for the layout.
//!
//! ## Example
//!
//! The formats of the input files are detected from their content, so `file_type` is
//! usually `None`:
//!
//! ```rust,no_run
//! use nanoget_rs::convenience::extract_auto;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let metrics = extract_auto(vec!["reads.fastq.gz", "more_reads.bam"])?;
//! println!("{} reads", metrics.len());
//!
//! // Columns are borrowed, not projected out of a sequence of structs.
//! let lengths: &[u32] = metrics.reads.lengths();
//! println!("longest {:?}", lengths.iter().max());
//!
//! // Row-wise access where it reads better.
//! for read in metrics.iter().take(5) {
//!     println!("{:?} {} bp", read.read_id(), read.length());
//! }
//! # Ok(())
//! # }
//! ```
//!
//! For full control, build [`ExtractArgs`] directly. Setting `file_type` to `Some(..)`
//! overrides detection — needed only to force a format or to select
//! [`FileType::FastqMinimal`], which is a processing mode rather than a format:
//!
//! ```rust,no_run
//! use nanoget_rs::{extract_metrics, CombineMethod, ExtractArgs, OutputFormat, ReadType};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let args = ExtractArgs {
//!     files: vec!["reads.fastq".into()],
//!     file_type: None,
//!     threads: 4,
//!     output_format: OutputFormat::Json,
//!     output: None,
//!     read_type: ReadType::OneD,
//!     barcoded: false,
//!     keep_supplementary: true,
//!     combine: CombineMethod::Simple,
//!     names: None,
//! };
//!
//! let metrics = extract_metrics(&args)?;
//! # Ok(())
//! # }
//! ```

pub mod cli;
pub mod columns;
pub mod error;
pub mod extract;
pub mod formats;
pub mod metrics;
pub mod utils;

pub use cli::{Cli, CombineMethod, Commands, ExtractArgs, OutputFormat, ReadType};
pub use columns::{ReadColumns, ReadColumnsBuilder, ReadView};
pub use error::NanogetError;
pub use extract::extract_metrics;
pub use formats::FileType;
pub use metrics::{MetricsCollection, MetricsSummary, ReadMetrics, StatsSummary};

/// Convenience functions for common use cases
pub mod convenience {
    use super::*;
    use std::path::Path;

    const DEFAULT_THREADS: usize = 4;

    /// Create default ExtractArgs with the given files and file type
    fn default_args(files: Vec<std::path::PathBuf>, file_type: Option<FileType>) -> ExtractArgs {
        ExtractArgs {
            files,
            file_type,
            threads: DEFAULT_THREADS,
            output_format: OutputFormat::default(),
            output: None,
            read_type: ReadType::default(),
            barcoded: false,
            keep_supplementary: true,
            combine: CombineMethod::default(),
            names: None,
        }
    }

    /// Extract metrics from files whose formats are detected from their content.
    pub fn extract_auto<P: AsRef<Path>>(files: Vec<P>) -> Result<MetricsCollection, NanogetError> {
        let args = default_args(
            files
                .into_iter()
                .map(|p| p.as_ref().to_path_buf())
                .collect(),
            None,
        );
        extract_metrics(&args)
    }

    /// Extract metrics from a single FASTQ file with default settings
    pub fn extract_from_fastq<P: AsRef<Path>>(file: P) -> Result<MetricsCollection, NanogetError> {
        let args = default_args(vec![file.as_ref().to_path_buf()], Some(FileType::Fastq));
        extract_metrics(&args)
    }

    /// Extract metrics from a single BAM file with default settings
    pub fn extract_from_bam<P: AsRef<Path>>(file: P) -> Result<MetricsCollection, NanogetError> {
        let args = default_args(vec![file.as_ref().to_path_buf()], Some(FileType::Bam));
        extract_metrics(&args)
    }

    /// Extract metrics from a single FASTA file with default settings
    pub fn extract_from_fasta<P: AsRef<Path>>(file: P) -> Result<MetricsCollection, NanogetError> {
        let args = default_args(vec![file.as_ref().to_path_buf()], Some(FileType::Fasta));
        extract_metrics(&args)
    }

    /// Extract metrics from multiple files with automatic format detection
    pub fn extract_from_files<P: AsRef<Path>>(
        files: Vec<P>,
        file_type: FileType,
        threads: Option<usize>,
    ) -> Result<MetricsCollection, NanogetError> {
        let mut args = default_args(
            files
                .into_iter()
                .map(|p| p.as_ref().to_path_buf())
                .collect(),
            Some(file_type),
        );
        if let Some(t) = threads {
            args.threads = t;
        }
        extract_metrics(&args)
    }
}
