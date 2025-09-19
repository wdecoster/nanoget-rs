// Example of using nanoget-rs as a library

use nanoget_rs::{extract_metrics, ExtractArgs, FileType, ReadMetrics, MetricsCollection};
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Example 1: Extract metrics from a single FASTQ file
    let args = ExtractArgs {
        files: vec![PathBuf::from("reads.fastq")],
        file_type: FileType::Fastq,
        threads: 4,
        output_format: "json".to_string(),
        output: None,
        read_type: "1D".to_string(),
        barcoded: false,
        keep_supplementary: true,
        huge: false,
        combine: "simple".to_string(),
        names: None,
    };
    
    let metrics: MetricsCollection = extract_metrics(&args)?;
    
    // Access individual read metrics
    for read in &metrics.reads {
        println!("Read {}: {} bp, quality: {:.2}", 
                 read.read_id.as_deref().unwrap_or("unknown"),
                 read.length,
                 read.quality.unwrap_or(0.0));
    }
    
    // Access summary statistics
    println!("Total reads: {}", metrics.summary.read_count);
    println!("Mean length: {:.2}", metrics.summary.length_stats.mean);
    println!("Mean quality: {:.2}", 
             metrics.summary.quality_stats.as_ref()
                 .map(|q| q.mean)
                 .unwrap_or(0.0));
    
    // Example 2: Process multiple files with tracking
    let multi_args = ExtractArgs {
        files: vec![
            PathBuf::from("sample1.fastq"),
            PathBuf::from("sample2.fastq"),
        ],
        file_type: FileType::Fastq,
        threads: 8,
        output_format: "json".to_string(),
        output: None,
        read_type: "1D".to_string(),
        barcoded: false,
        keep_supplementary: true,
        huge: false,
        combine: "track".to_string(),
        names: Some(vec!["Sample1".to_string(), "Sample2".to_string()]),
    };
    
    let multi_metrics = extract_metrics(&multi_args)?;
    
    // Analyze per-dataset
    for dataset_name in ["Sample1", "Sample2"] {
        let dataset_reads: Vec<_> = multi_metrics.reads.iter()
            .filter(|r| r.dataset.as_deref() == Some(dataset_name))
            .collect();
        
        println!("{}: {} reads", dataset_name, dataset_reads.len());
    }
    
    Ok(())
}

// Example of creating custom analysis functions
fn analyze_read_length_distribution(metrics: &MetricsCollection) -> (f64, f64, f64) {
    let lengths: Vec<f64> = metrics.reads.iter()
        .map(|r| r.length as f64)
        .collect();
    
    let mean = lengths.iter().sum::<f64>() / lengths.len() as f64;
    let min = lengths.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max = lengths.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    
    (mean, min, max)
}

// Example of custom filtering
fn filter_high_quality_reads(metrics: &MetricsCollection, min_quality: f64) -> Vec<&ReadMetrics> {
    metrics.reads.iter()
        .filter(|read| read.quality.map(|q| q >= min_quality).unwrap_or(false))
        .collect()
}