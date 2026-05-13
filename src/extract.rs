use crate::cli::ExtractArgs;
use crate::error::NanogetError;
use crate::formats::FileType;
use crate::metrics::{MetricsCollection, ReadMetrics};
use crate::utils;

use chrono::{DateTime, TimeZone, Utc};
use log::info;
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::Read as BamRead;
use std::io::Read;
use std::path::Path;

/// Safely parse a timestamp (seconds since epoch) to DateTime<Utc>
/// Handles nanosecond overflow by clamping to valid range
fn parse_timestamp(timestamp: f64) -> Option<DateTime<Utc>> {
    // Validate timestamp range
    if timestamp < 0.0 || timestamp > i64::MAX as f64 {
        return None;
    }

    let seconds = timestamp as i64;
    // Clamp nanoseconds to valid u32 range (0 to 999,999,999)
    let nanos = ((timestamp.fract().abs() * 1e9) as u32).min(999_999_999);

    Utc.timestamp_opt(seconds, nanos).single()
}

/// Main entry point for extracting metrics from files
pub fn extract_metrics(args: &ExtractArgs) -> Result<MetricsCollection, NanogetError> {
    // Stdin shortcut: single "-" path handled entirely here.
    if args.files.len() == 1 && args.files[0].as_os_str() == "-" {
        return extract_metrics_stdin(args);
    }

    info!(
        "Starting nanoget extraction with {} files",
        args.files.len()
    );

    // Validate input files
    for file in &args.files {
        utils::check_file_exists(file)?;
    }

    let collections = args
        .files
        .par_iter()
        .map(|file| process_single_file(file, &args.file_type, args))
        .collect::<Result<Vec<_>, _>>()?;

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
        FileType::Bam => process_bam(file, args.keep_supplementary, args.threads)?,
        FileType::Cram => process_bam(file, args.keep_supplementary, args.threads)?,
        FileType::Ubam => process_ubam(file)?,
        FileType::Summary => process_summary(file, &args.read_type, args.barcoded)?,
    };

    Ok(MetricsCollection::new(reads))
}

/// Process FASTQ files
fn process_fastq(file: &Path, rich: bool) -> Result<Vec<ReadMetrics>, NanogetError> {
    let reader = utils::open_file(file)?;
    process_fastq_from_reader(reader, rich)
}

fn process_fastq_from_reader<R: Read>(
    reader: R,
    rich: bool,
) -> Result<Vec<ReadMetrics>, NanogetError> {
    use bio::io::fastq;

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
            info!("Processed {} reads", i);
        }
    }

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
        metrics.push(ReadMetrics::new(None, record.seq().len() as u32));
    }

    Ok(metrics)
}

/// Process FASTA files
fn process_fasta(file: &Path) -> Result<Vec<ReadMetrics>, NanogetError> {
    let reader = utils::open_file(file)?;
    process_fasta_from_reader(reader)
}

fn process_fasta_from_reader<R: Read>(reader: R) -> Result<Vec<ReadMetrics>, NanogetError> {
    use bio::io::fasta;

    let fasta_reader = fasta::Reader::new(reader);
    let mut metrics = Vec::new();

    for result in fasta_reader.records() {
        let record = result.map_err(|e| NanogetError::ParseError(e.to_string()))?;
        metrics.push(ReadMetrics::new(
            Some(record.id().to_string()),
            record.seq().len() as u32,
        ));
    }

    Ok(metrics)
}

/// Get the NM (edit distance) tag from a BAM record
fn get_nm_tag(record: &rust_htslib::bam::Record) -> Option<u32> {
    match record.aux(b"NM") {
        Ok(value) => match value {
            Aux::U8(v) => Some(u32::from(v)),
            Aux::U16(v) => Some(u32::from(v)),
            Aux::U32(v) => Some(v),
            Aux::I8(v) => u32::try_from(v).ok(),
            Aux::I16(v) => u32::try_from(v).ok(),
            Aux::I32(v) => u32::try_from(v).ok(),
            _ => None,
        },
        Err(_) => None,
    }
}

/// Get the de (gap-compressed divergence) tag from a BAM record
/// This is provided by recent minimap2 versions
fn get_de_tag(record: &rust_htslib::bam::Record) -> Option<f64> {
    match record.aux(b"de") {
        Ok(value) => match value {
            Aux::Float(v) => Some(100.0 * (1.0 - v as f64)),
            _ => None,
        },
        Err(_) => None,
    }
}

/// Extract aligned length and gap-compressed identity with at most one CIGAR pass.
///
/// When the minimap2 `de` tag is present: one minimal CIGAR pass for aligned length only.
/// When absent: one combined CIGAR pass computing both values simultaneously.
fn alignment_stats(record: &rust_htslib::bam::Record) -> (u32, Option<f64>) {
    let mut aligned_len: u32 = 0;

    if let Some(identity) = get_de_tag(record) {
        // Minimal pass: aligned length only, no identity bookkeeping needed
        for entry in record.cigar().iter() {
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Ins(len) => {
                    aligned_len += len;
                }
                _ => {}
            }
        }
        return (aligned_len, Some(identity));
    }

    // No de tag: compute both in one pass
    let nm = get_nm_tag(record);
    let mut matches: u32 = 0;
    let mut gap_size: u32 = 0;
    let mut gap_count: u32 = 0;

    for entry in record.cigar().iter() {
        match entry {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                aligned_len += len;
                matches += len;
            }
            Cigar::Ins(len) => {
                aligned_len += len;
                gap_size += len;
                gap_count += 1;
            }
            Cigar::Del(len) => {
                gap_size += len;
                gap_count += 1;
            }
            _ => {}
        }
    }

    let identity = nm.and_then(|nm| {
        let denominator = matches + gap_count;
        if denominator == 0 {
            return None;
        }
        let numerator = nm.saturating_sub(gap_size) + gap_count;
        Some(100.0 * (1.0 - (numerator as f64 / denominator as f64)))
    });

    (aligned_len, identity)
}

/// Process BAM or CRAM files using sequential streaming with BGZF multi-threading.
///
/// htslib's BGZF threading pre-decompresses blocks on background threads while the
/// main thread processes records — much faster than chromosome-level parallelism,
/// which forces random seeks that break sequential BGZF streaming.
fn process_bam(
    file: &Path,
    keep_supplementary: bool,
    threads: usize,
) -> Result<Vec<ReadMetrics>, NanogetError> {
    let mut reader = if file.as_os_str() == "-" {
        rust_htslib::bam::Reader::from_stdin()?
    } else {
        rust_htslib::bam::Reader::from_path(file)?
    };
    // Use all-but-one thread for BGZF decompression; htslib manages the pool.
    let bgzf_threads = threads.saturating_sub(1);
    if bgzf_threads > 0 {
        reader
            .set_threads(bgzf_threads)
            .map_err(|e| NanogetError::ProcessingError(e.to_string()))?;
    }
    info!(
        "Processing {} with {} BGZF threads",
        file.display(),
        bgzf_threads
    );
    extract_bam_records(&mut reader, keep_supplementary)
}

/// Extract ReadMetrics from any type implementing bam::Read.
fn extract_bam_records<R: BamRead>(
    reader: &mut R,
    keep_supplementary: bool,
) -> Result<Vec<ReadMetrics>, NanogetError> {
    let mut metrics = Vec::new();

    for result in reader.records() {
        let record = result?;

        if record.is_unmapped() {
            continue;
        }
        if !keep_supplementary && record.is_supplementary() {
            continue;
        }

        let read_id = String::from_utf8_lossy(record.qname()).to_string();
        let length = record.seq().len() as u32;
        let (aligned_length, percent_identity) = alignment_stats(&record);
        let mapping_quality = if record.mapq() == 255 {
            None
        } else {
            Some(record.mapq())
        };
        let quality = utils::average_quality(record.qual());

        metrics.push(
            ReadMetrics::new(Some(read_id), length)
                .with_quality(quality.unwrap_or(0.0))
                .with_alignment(aligned_length, quality, mapping_quality, percent_identity),
        );
    }

    Ok(metrics)
}

/// Process unaligned BAM files
fn process_ubam(file: &Path) -> Result<Vec<ReadMetrics>, NanogetError> {
    use rust_htslib::{bam, bam::Read};

    let mut bam_reader = if file.as_os_str() == "-" {
        bam::Reader::from_stdin()?
    } else {
        bam::Reader::from_path(file)?
    };
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
    let reader = utils::open_file(file)?;
    process_summary_from_reader(reader, read_type, barcoded)
}

fn process_summary_from_reader<R: Read>(
    reader: R,
    read_type: &str,
    barcoded: bool,
) -> Result<Vec<ReadMetrics>, NanogetError> {
    use csv::ReaderBuilder;
    use std::collections::HashMap;

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
            .and_then(parse_timestamp);

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

/// Read from stdin: peek with fill_buf() to detect format, then route to the appropriate parser.
///
/// For text formats (FASTQ, FASTA, summary TSV): the BufReader is passed directly to the parser.
/// `fill_buf()` does not advance the BufReader's read position, so no bytes are lost.
///
/// For binary formats (BAM/CRAM): htslib reads from OS fd 0 directly, bypassing the BufReader.
/// We reconstruct stdin at the OS level by prepending the peeked bytes via a pipe + background thread.
fn extract_metrics_stdin(args: &ExtractArgs) -> Result<MetricsCollection, NanogetError> {
    use std::io::BufRead;

    let mut stdin_reader = std::io::BufReader::new(std::io::stdin());

    // Peek without consuming (BufReader internal buffer is filled, read position stays at 0).
    let file_type = {
        let peek = stdin_reader
            .fill_buf()
            .map_err(|e| NanogetError::ParseError(format!("Failed to read stdin: {}", e)))?;
        FileType::sniff_stdin_bytes(peek)?
    };

    info!("Detected stdin format: {:?}", file_type);

    let reads = match &file_type {
        FileType::Bam | FileType::Cram | FileType::Ubam => {
            // htslib reads from OS fd 0 directly, bypassing the BufReader.
            // Extract the peeked bytes and reconstruct fd 0 via a pipe so htslib
            // sees a complete, untruncated stream.
            let sniffed = stdin_reader.buffer().to_vec();
            drop(stdin_reader);
            reconstruct_stdin_prefix(sniffed)?;
            match file_type {
                FileType::Ubam => process_ubam(Path::new("-"))?,
                _ => process_bam(Path::new("-"), args.keep_supplementary, args.threads)?,
            }
        }
        _ => {
            // Text formats (FASTQ, FASTA, summary TSV) — may be gzip-compressed.
            // The BufReader still has all peeked bytes at position 0, so we can wrap
            // it in a GzDecoder if the stream is gzip-encoded.
            let is_plain_gzip = {
                let buf = stdin_reader.buffer();
                buf.len() >= 2
                    && buf[0] == 0x1f
                    && buf[1] == 0x8b
                    && !(buf.len() >= 16
                        && buf[3] & 0x04 != 0
                        && buf[12..16] == [0x42, 0x43, 0x02, 0x00])
            };
            let reader: Box<dyn Read> = if is_plain_gzip {
                Box::new(flate2::bufread::GzDecoder::new(stdin_reader))
            } else {
                Box::new(stdin_reader)
            };
            match file_type {
                FileType::Fastq | FileType::FastqRich => process_fastq_from_reader(reader, false)?,
                FileType::Fasta => process_fasta_from_reader(reader)?,
                FileType::Summary => {
                    process_summary_from_reader(reader, &args.read_type, args.barcoded)?
                }
                other => {
                    return Err(NanogetError::ParseError(format!(
                        "Format {:?} is not supported for stdin input",
                        other
                    )))
                }
            }
        }
    };

    Ok(MetricsCollection::new(reads))
}

/// Prepend `prefix` bytes to stdin by replacing fd 0 with a pipe whose write end is fed by a
/// background thread (prefix bytes first, then the rest of the original stdin).
///
/// This allows htslib — which reads from fd 0 directly — to see a complete, untruncated stream
/// even after we have consumed `prefix.len()` bytes from the OS stdin for format detection.
#[cfg(unix)]
fn reconstruct_stdin_prefix(prefix: Vec<u8>) -> Result<(), NanogetError> {
    use std::os::unix::io::FromRawFd;

    unsafe {
        // Save a dup of the current stdin before we replace it.
        let saved_stdin = libc::dup(0);
        if saved_stdin < 0 {
            return Err(NanogetError::ProcessingError(
                "Failed to dup stdin fd".into(),
            ));
        }

        // Create an anonymous pipe.
        let mut pipe_fds: [libc::c_int; 2] = [0; 2];
        if libc::pipe(pipe_fds.as_mut_ptr()) != 0 {
            libc::close(saved_stdin);
            return Err(NanogetError::ProcessingError(
                "Failed to create pipe for stdin reconstruction".into(),
            ));
        }
        let (read_fd, write_fd) = (pipe_fds[0], pipe_fds[1]);

        // Replace stdin (fd 0) with the read end of the pipe.
        if libc::dup2(read_fd, 0) < 0 {
            libc::close(read_fd);
            libc::close(write_fd);
            libc::close(saved_stdin);
            return Err(NanogetError::ProcessingError(
                "Failed to redirect stdin to pipe".into(),
            ));
        }
        libc::close(read_fd); // fd 0 is now the only reference to the read end.

        // Background thread: write prefix, then drain the original stdin into the write end.
        std::thread::spawn(move || {
            use std::io::Write;
            let mut writer = std::fs::File::from_raw_fd(write_fd);
            let mut orig = std::fs::File::from_raw_fd(saved_stdin);
            let _ = writer.write_all(&prefix);
            let _ = std::io::copy(&mut orig, &mut writer);
            // Both fds are closed when writer/orig drop, signalling EOF to the reader.
        });
    }

    Ok(())
}

#[cfg(not(unix))]
fn reconstruct_stdin_prefix(_prefix: Vec<u8>) -> Result<(), NanogetError> {
    Err(NanogetError::ProcessingError(
        "BAM/CRAM from stdin is only supported on Unix".into(),
    ))
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
                        metadata.start_time = parse_timestamp(timestamp);
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
