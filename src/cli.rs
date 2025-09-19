use clap::{Parser, Subcommand, Args};
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

#[derive(Args)]
pub struct ExtractArgs {
    /// Input files to process
    #[arg(required = true)]
    pub files: Vec<PathBuf>,
    
    /// Type of input files
    #[arg(short = 't', long, value_enum)]
    pub file_type: crate::formats::FileType,
    
    /// Number of threads to use for processing
    #[arg(short = 'j', long, default_value = "4")]
    pub threads: usize,
    
    /// Output format (json, csv, tsv)
    #[arg(short = 'f', long, default_value = "json")]
    pub output_format: String,
    
    /// Output file (optional, defaults to stdout)
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,
    
    /// For summary files: read type (1D, 2D, 1D2)
    #[arg(long, default_value = "1D")]
    pub read_type: String,
    
    /// Include barcoded reads analysis
    #[arg(long)]
    pub barcoded: bool,
    
    /// Keep supplementary alignments (for BAM files)
    #[arg(long, default_value = "true")]
    pub keep_supplementary: bool,
    
    /// Process huge files without parallelization
    #[arg(long)]
    pub huge: bool,
    
    /// Combine multiple files: simple or track
    #[arg(long, default_value = "simple")]
    pub combine: String,
    
    /// Names for datasets when using track mode
    #[arg(long)]
    pub names: Option<Vec<String>>,
}