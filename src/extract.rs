use crate::cli::ExtractArgs;
use crate::error::NanogetError;
use crate::formats::FileType;
use crate::metrics::{MetricsCollection, ReadMetrics};
use crate::utils;

use log::info;
use rayon::prelude::*;
use std::path::Path;

/// Main entry point for extracting metrics from files
pub fn extract_metrics(args: &ExtractArgs) -> Result<MetricsCollection, NanogetError> {
    info!(
        "Starting nanoget extraction with {} files",
        args.files.len()
    );

    // Validate input files
    for file in &args.files {
        utils::check_file_exists(file)?;
    }

    // Determine processing strategy based on file type and options
    // Process files in parallel
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .map_err(|e| NanogetError::ProcessingError(e.to_string()))?;

    let collections = thread_pool.install(|| {
        args.files
            .par_iter()
            .map(|file| process_single_file(file, &args.file_type, args))
            .collect::<Result<Vec<_>, _>>()
    })?;

    // Combine results
    let combined = MetricsCollection::combine(collections, &args.combine, args.names.clone());

    info!(
        "Extraction complete: {} reads processed",
        combined.summary.read_count
    );

    if combined.summary.read_count == 0 {
        return Err(NanogetError::ProcessingError(
            "No reads found in input files".to_string(),
        ));
    }

    Ok(combined)
}

/// Process a single file and return metrics
fn process_single_file(
    file: &Path,
    file_type: &FileType,
    args: &ExtractArgs,
) -> Result<MetricsCollection, NanogetError> {
    info!("Processing file: {}", file.display());

    let reads = match file_type {
        FileType::Fastq => process_fastq(file, false)?,
        FileType::FastqRich => process_fastq(file, true)?,
        FileType::FastqMinimal => process_fastq_minimal(file)?,
        FileType::Fasta => process_fasta(file)?,
        FileType::Bam => process_bam(file, args.keep_supplementary)?,
        FileType::Cram => process_cram(file, args.keep_supplementary)?,
        FileType::Ubam => process_ubam(file)?,
        FileType::Summary => process_summary(file, &args.read_type, args.barcoded)?,
    };

    Ok(MetricsCollection::new(reads))
}

/// Process FASTQ files
fn process_fastq(file: &Path, rich: bool) -> Result<Vec<ReadMetrics>, NanogetError> {
    use bio::io::fastq;

    let reader = utils::open_file(file)?;
    let fastq_reader = fastq::Reader::new(reader);
    let mut metrics = Vec::new();

    for (i, result) in fastq_reader.records().enumerate() {
        let record = result.map_err(|e| NanogetError::ParseError(e.to_string()))?;

        let read_id = record.id().to_string();
        let length = record.seq().len() as u32;
        let quality = utils::average_quality(record.qual());

        let mut read_metrics = ReadMetrics::new(Some(read_id), length);

        if let Some(q) = quality {
            read_metrics = read_metrics.with_quality(q);
        }

        // For rich FASTQ, try to extract additional metadata from the description
        if rich {
            let desc = record.desc().unwrap_or("");
            if let Some(metadata) = parse_rich_fastq_metadata(desc) {
                read_metrics = read_metrics.with_sequencing_metadata(
                    metadata.channel_id,
                    metadata.start_time,
                    metadata.duration,
                );
                read_metrics.run_id = metadata.run_id;
            }
        }

        metrics.push(read_metrics);

        if i % 10000 == 0 && i > 0 {
            info!("Processed {} reads from {}", i, file.display());
        }
    }

    info!(
        "Finished processing {} reads from {}",
        metrics.len(),
        file.display()
    );
    Ok(metrics)
}

/// Process FASTQ files with minimal information (length only)
fn process_fastq_minimal(file: &Path) -> Result<Vec<ReadMetrics>, NanogetError> {
    use bio::io::fastq;

    let reader = utils::open_file(file)?;
    let fastq_reader = fastq::Reader::new(reader);
    let mut metrics = Vec::new();

    for result in fastq_reader.records() {
        let record = result.map_err(|e| NanogetError::ParseError(e.to_string()))?;

        let length = record.seq().len() as u32;
        let read_metrics = ReadMetrics::new(None, length); // No read ID for minimal processing

        metrics.push(read_metrics);
    }

    Ok(metrics)
}

/// Process FASTA files
fn process_fasta(file: &Path) -> Result<Vec<ReadMetrics>, NanogetError> {
    use bio::io::fasta;

    let reader = utils::open_file(file)?;
    let fasta_reader = fasta::Reader::new(reader);
    let mut metrics = Vec::new();

    for result in fasta_reader.records() {
        let record = result.map_err(|e| NanogetError::ParseError(e.to_string()))?;

        let read_id = record.id().to_string();
        let length = record.seq().len() as u32;

        let read_metrics = ReadMetrics::new(Some(read_id), length);
        metrics.push(read_metrics);
    }

    Ok(metrics)
}

/// Process BAM files
fn process_bam(file: &Path, keep_supplementary: bool) -> Result<Vec<ReadMetrics>, NanogetError> {
    use rust_htslib::{bam, bam::Read};

    let mut bam_reader = bam::Reader::from_path(file)?;
    let mut metrics = Vec::new();

    for result in bam_reader.records() {
        let record = result?;

        // Skip unmapped reads
        if record.is_unmapped() {
            continue;
        }

        // Skip supplementary alignments if requested
        if !keep_supplementary && record.is_supplementary() {
            continue;
        }

        let read_id = String::from_utf8_lossy(record.qname()).to_string();
        let length = record.seq().len() as u32;
        let aligned_length = record.seq().len() as u32; // TODO: Calculate actual aligned length from CIGAR
        let mapping_quality = if record.mapq() == 255 {
            None
        } else {
            Some(record.mapq())
        };

        // Calculate quality scores
        let quality = record
            .qual()
            .iter()
            .any(|&q| q != 255)
            .then(|| utils::average_quality(record.qual()).unwrap_or(0.0));

        let aligned_quality = quality; // Same as overall quality for now

        // Calculate percent identity (simplified - would need CIGAR parsing for accuracy)
        let percent_identity = Some(95.0); // Placeholder - would calculate from CIGAR

        let read_metrics = ReadMetrics::new(Some(read_id), length)
            .with_quality(quality.unwrap_or(0.0))
            .with_alignment(
                aligned_length,
                aligned_quality,
                mapping_quality,
                percent_identity,
            );

        metrics.push(read_metrics);
    }

    Ok(metrics)
}

/// Process CRAM files (similar to BAM)
fn process_cram(file: &Path, keep_supplementary: bool) -> Result<Vec<ReadMetrics>, NanogetError> {
    // CRAM processing would be similar to BAM but with rust-htslib's CRAM support
    // For now, we'll use the same logic as BAM
    process_bam(file, keep_supplementary)
}

/// Process unaligned BAM files
fn process_ubam(file: &Path) -> Result<Vec<ReadMetrics>, NanogetError> {
    use rust_htslib::{bam, bam::Read};

    let mut bam_reader = bam::Reader::from_path(file)?;
    let mut metrics = Vec::new();

    for result in bam_reader.records() {
        let record = result?;

        let read_id = String::from_utf8_lossy(record.qname()).to_string();
        let length = record.seq().len() as u32;

        // Calculate quality scores
        let quality = record
            .qual()
            .iter()
            .any(|&q| q != 255)
            .then(|| utils::average_quality(record.qual()).unwrap_or(0.0));

        let mut read_metrics = ReadMetrics::new(Some(read_id), length);

        if let Some(q) = quality {
            read_metrics = read_metrics.with_quality(q);
        }

        metrics.push(read_metrics);
    }

    Ok(metrics)
}

/// Process sequencing summary files
fn process_summary(
    file: &Path,
    read_type: &str,
    barcoded: bool,
) -> Result<Vec<ReadMetrics>, NanogetError> {
    use csv::ReaderBuilder;
    use std::collections::HashMap;

    let reader = utils::open_file(file)?;
    let mut csv_reader = ReaderBuilder::new().delimiter(b'\t').from_reader(reader);

    // Get headers
    let headers = csv_reader.headers()?.clone();
    let mut metrics = Vec::new();

    for result in csv_reader.records() {
        let record = result?;
        let row: HashMap<&str, &str> = headers.iter().zip(record.iter()).collect();

        // Extract fields based on read type
        let (length_field, quality_field) = match read_type {
            "1D" => ("sequence_length_template", "mean_qscore_template"),
            "2D" | "1D2" => ("sequence_length_2d", "mean_qscore_2d"),
            _ => {
                return Err(NanogetError::InvalidInput(format!(
                    "Unsupported read type: {}",
                    read_type
                )))
            }
        };

        let length: u32 = row
            .get(length_field)
            .ok_or_else(|| NanogetError::ParseError(format!("Missing column: {}", length_field)))?
            .parse()
            .map_err(|e| NanogetError::ParseError(format!("Invalid length: {}", e)))?;

        let quality: f64 = row
            .get(quality_field)
            .ok_or_else(|| NanogetError::ParseError(format!("Missing column: {}", quality_field)))?
            .parse()
            .map_err(|e| NanogetError::ParseError(format!("Invalid quality: {}", e)))?;

        let channel_id: Option<u16> = row.get("channel").and_then(|s| s.parse().ok());

        let start_time = row
            .get("start_time")
            .and_then(|s| s.parse::<f64>().ok())
            .and_then(|timestamp| {
                use chrono::{TimeZone, Utc};
                Utc.timestamp_opt(timestamp as i64, (timestamp.fract() * 1e9) as u32)
                    .single()
            });

        let duration: Option<f64> = row.get("duration").and_then(|s| s.parse().ok());

        let barcode = if barcoded {
            row.get("barcode_arrangement").map(|s| s.to_string())
        } else {
            None
        };

        let mut read_metrics = ReadMetrics::new(None, length)
            .with_quality(quality)
            .with_sequencing_metadata(channel_id, start_time, duration);

        read_metrics.barcode = barcode;

        metrics.push(read_metrics);
    }

    Ok(metrics)
}

/// Metadata extracted from rich FASTQ descriptions
#[derive(Debug)]
struct RichFastqMetadata {
    channel_id: Option<u16>,
    start_time: Option<chrono::DateTime<chrono::Utc>>,
    duration: Option<f64>,
    run_id: Option<String>,
}

/// Parse metadata from rich FASTQ description lines
fn parse_rich_fastq_metadata(desc: &str) -> Option<RichFastqMetadata> {
    // Parse key=value pairs from the description
    let mut metadata = RichFastqMetadata {
        channel_id: None,
        start_time: None,
        duration: None,
        run_id: None,
    };

    for pair in desc.split_whitespace() {
        if let Some((key, value)) = pair.split_once('=') {
            match key {
                "ch" => {
                    metadata.channel_id = value.parse().ok();
                }
                "start_time" => {
                    if let Ok(timestamp) = value.parse::<f64>() {
                        use chrono::{TimeZone, Utc};
                        metadata.start_time = Utc
                            .timestamp_opt(timestamp as i64, (timestamp.fract() * 1e9) as u32)
                            .single();
                    }
                }
                "duration" => {
                    metadata.duration = value.parse().ok();
                }
                "runid" => {
                    metadata.run_id = Some(value.to_string());
                }
                _ => {} // Ignore unknown keys
            }
        }
    }

    // Return Some only if we found at least one piece of metadata
    if metadata.channel_id.is_some()
        || metadata.start_time.is_some()
        || metadata.duration.is_some()
        || metadata.run_id.is_some()
    {
        Some(metadata)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rich_fastq_metadata_parsing() {
        let desc = "ch=100 start_time=1234567890.5 duration=2.5 runid=test_run";
        let metadata = parse_rich_fastq_metadata(desc).unwrap();

        assert_eq!(metadata.channel_id, Some(100));
        assert_eq!(metadata.duration, Some(2.5));
        assert_eq!(metadata.run_id, Some("test_run".to_string()));
    }
}
