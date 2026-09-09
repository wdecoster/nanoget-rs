use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(name = "nanoget")]
#[command(about = "Extract metrics from Oxford Nanopore sequencing data")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Extract metrics from sequencing files
    Extract(ExtractArgs),
}

/// Serialisation format for the extracted metrics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum, Serialize, Deserialize)]
pub enum OutputFormat {
    /// Pretty-printed JSON
    #[default]
    Json,
    /// Tab-separated values, with summary statistics as trailing comments
    Tsv,
}

/// How metrics from multiple input files are merged.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum, Serialize, Deserialize)]
pub enum CombineMethod {
    /// Concatenate all reads into one collection
    #[default]
    Simple,
    /// Concatenate, tagging each read with the name of the file it came from
    Track,
}

/// Which columns to read from a sequencing summary file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum, Serialize, Deserialize)]
pub enum ReadType {
    /// 1D sequencing: `sequence_length_template` / `mean_qscore_template`
    #[default]
    #[value(name = "1D", alias = "1d")]
    OneD,
    /// 2D sequencing: `sequence_length_2d` / `mean_qscore_2d`
    #[value(name = "2D", alias = "2d")]
    TwoD,
    /// 1D^2 sequencing: same columns as 2D
    #[value(name = "1D2", alias = "1d2")]
    OneD2,
}

impl ReadType {
    /// The (length, quality) column names this read type reads from a summary file.
    pub fn summary_columns(self) -> (&'static str, &'static str) {
        match self {
            Self::OneD => ("sequence_length_template", "mean_qscore_template"),
            Self::TwoD | Self::OneD2 => ("sequence_length_2d", "mean_qscore_2d"),
        }
    }
}

#[derive(Args, Debug, Clone)]
pub struct ExtractArgs {
    /// Input files to process, or `-` to read from stdin
    #[arg(required = true)]
    pub files: Vec<PathBuf>,

    /// Type of input files. Detected from the content of each file when omitted; give it
    /// explicitly to override detection, or to select `fastq-minimal`, which is a
    /// processing mode rather than a format and so cannot be detected.
    #[arg(short = 't', long, value_enum)]
    pub file_type: Option<crate::formats::FileType>,

    /// Number of threads to use for processing
    #[arg(short = 'j', long, default_value_t = 4)]
    pub threads: usize,

    /// Output format
    #[arg(short = 'f', long, value_enum, default_value_t = OutputFormat::Json)]
    pub output_format: OutputFormat,

    /// Output file (optional, defaults to stdout)
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,

    /// For summary files: which read type's columns to read
    #[arg(long, value_enum, default_value_t = ReadType::OneD)]
    pub read_type: ReadType,

    /// Include barcoded reads analysis
    #[arg(long)]
    pub barcoded: bool,

    /// Exclude supplementary alignments from BAM/CRAM. They are hard-clipped fragments
    /// of a read, so including them inflates read counts and yield. Kept by default.
    #[arg(
        long = "drop-supplementary",
        action = ArgAction::SetFalse,
        default_value_t = true
    )]
    pub keep_supplementary: bool,

    /// How to combine metrics from multiple files
    #[arg(long, value_enum, default_value_t = CombineMethod::Simple)]
    pub combine: CombineMethod,

    /// Names for datasets when using `--combine track`
    #[arg(long)]
    pub names: Option<Vec<String>>,
}
