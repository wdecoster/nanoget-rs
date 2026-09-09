use crate::cli::{ExtractArgs, ReadType};
use crate::columns::{ReadColumns, ReadColumnsBuilder};
use crate::error::NanogetError;
use crate::formats::FileType;
use crate::metrics::{MetricsCollection, ReadMetrics};
use crate::utils;

use chrono::{DateTime, TimeZone, Utc};
use log::info;
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::Read as BamRead;
use rust_htslib::htslib::{
    hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS, sam_fields_SAM_AUX, sam_fields_SAM_CIGAR,
    sam_fields_SAM_FLAG, sam_fields_SAM_MAPQ, sam_fields_SAM_QNAME, sam_fields_SAM_QUAL,
    sam_fields_SAM_SEQ,
};
use std::io::Read;
use std::path::Path;

/// Safely parse a timestamp (seconds since epoch) to DateTime<Utc>
/// Handles nanosecond overflow by clamping to valid range
fn parse_timestamp(timestamp: f64) -> Option<DateTime<Utc>> {
    // Validate timestamp range. NaN fails every comparison, so it has to be rejected
    // explicitly — otherwise `NaN as i64` saturates to 0 and silently becomes the epoch.
    if !timestamp.is_finite() || timestamp < 0.0 || timestamp > i64::MAX as f64 {
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
        .map(|file| {
            // An explicit --file-type overrides detection; otherwise each file is
            // classified from its own content, so a mixed set still works.
            let file_type = match args.file_type {
                Some(ref t) => t.clone(),
                None => FileType::sniff(file)?,
            };
            process_single_file(file, &file_type, args)
        })
        .collect::<Result<Vec<ReadColumns>, _>>()?;

    // Combine results
    let combined = MetricsCollection::combine(collections, args.combine, args.names.clone());

    info!(
        "Extraction complete: {} reads processed",
        combined.summary().read_count
    );

    if combined.summary().read_count == 0 {
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
) -> Result<ReadColumns, NanogetError> {
    info!("Processing file: {} as {:?}", file.display(), file_type);

    let reads = match file_type {
        FileType::Fastq => process_fastq(file, false)?,
        FileType::FastqRich => process_fastq(file, true)?,
        FileType::FastqMinimal => process_fastq_minimal(file)?,
        FileType::Fasta => process_fasta(file)?,
        FileType::Bam => process_bam(file, args.keep_supplementary, args.threads)?,
        FileType::Cram => process_bam(file, args.keep_supplementary, args.threads)?,
        FileType::Ubam => process_ubam(file, args.threads)?,
        FileType::Summary => process_summary(file, args.read_type, args.barcoded)?,
    };

    Ok(reads)
}

/// Reject a malformed FASTQ record loudly.
///
/// `bio`'s reader does not validate what it parses, and corrupt or truncated FASTQ is
/// common enough in practice — a half-written file, two files concatenated mid-record,
/// an interrupted transfer — that silently accepting it and reporting plausible but
/// wrong metrics is worse than stopping. A sequence and quality line of different
/// lengths is the signature of nearly all of it.
///
/// Only O(1) checks are made: scanning every base for non-ASCII would cost a second
/// pass over the read for a corruption that a length mismatch almost always catches too.
fn check_fastq_record(record: &bio::io::fastq::Record, index: usize) -> Result<(), NanogetError> {
    let seq_len = record.seq().len();
    let qual_len = record.qual().len();
    if seq_len != qual_len {
        return Err(NanogetError::ParseError(format!(
            "Malformed FASTQ record {} ('{}'): sequence is {} bases \
             but the quality line is {} characters",
            index + 1,
            record.id(),
            seq_len,
            qual_len
        )));
    }
    if record.id().is_empty() {
        return Err(NanogetError::ParseError(format!(
            "Malformed FASTQ record {}: missing read id",
            index + 1
        )));
    }
    Ok(())
}

/// Process FASTQ files
fn process_fastq(file: &Path, rich: bool) -> Result<ReadColumns, NanogetError> {
    let reader = utils::open_file(file)?;
    process_fastq_from_reader(reader, rich)
}

fn process_fastq_from_reader<R: Read>(reader: R, rich: bool) -> Result<ReadColumns, NanogetError> {
    use bio::io::fastq;

    let fastq_reader = fastq::Reader::new(reader);
    let mut metrics = ReadColumnsBuilder::with_capacity(INITIAL_READ_CAPACITY);
    let mut skipped_empty = 0usize;

    for (i, result) in fastq_reader.records().enumerate() {
        let record = result.map_err(|e| NanogetError::ParseError(e.to_string()))?;
        check_fastq_record(&record, i)?;

        let length = record.seq().len() as u32;
        // Zero-length reads carry no usable metrics and drag the length distribution
        // down; python nanoget drops them too.
        if length == 0 {
            skipped_empty += 1;
            continue;
        }

        let read_id = record.id().to_string();
        // bio returns the raw FASTQ quality line, so it needs Phred+33 decoding.
        let quality = utils::average_quality_phred33(record.qual());

        let mut read_metrics = ReadMetrics::new(Some(read_id), length);

        if let Some(q) = quality {
            read_metrics = read_metrics.with_quality(q);
        }

        if rich {
            let desc = record.desc().unwrap_or("");
            // A record in a rich FASTQ with no recognisable metadata means something is
            // wrong with the file — a mixed concatenation, a stripped header, the wrong
            // --file-type — and silently emitting a read with no channel or start time
            // just makes the missing plots downstream unexplainable.
            let metadata = parse_rich_fastq_metadata(desc).ok_or_else(|| {
                NanogetError::ParseError(format!(
                    "FASTQ record {} ('{}') has no MinKNOW/albacore metadata in its \
                     header: expected key=value fields (ch=, start_time=, duration=, \
                     runid=) or SAM-style tags (ch:i:, st:Z:, du:f:, RG:Z:), found '{}'.\n\
                     The format was detected from the first record; if this file mixes \
                     rich and plain headers, pass --file-type fastq to read it without \
                     metadata.",
                    i + 1,
                    record.id(),
                    desc
                ))
            })?;
            read_metrics = read_metrics.with_sequencing_metadata(
                metadata.channel_id,
                metadata.start_time,
                metadata.duration,
            );
            read_metrics.run_id = metadata.run_id;
        }

        metrics.push(read_metrics);

        if i % 10000 == 0 && i > 0 {
            info!("Processed {} reads", i);
        }
    }

    log_skipped_empty(skipped_empty);
    Ok(metrics.finish())
}

/// Initial capacity for the per-file read accumulators.
///
/// A read count cannot be estimated from file size without being wrong by orders of
/// magnitude — record length varies from ~100 bases to ~100 kb — so this just skips the
/// first dozen doublings, which are the cheap ones. The expensive reallocations near the
/// end are avoided in `MetricsCollection::combine`, which knows the exact total.
const INITIAL_READ_CAPACITY: usize = 4096;

/// Report zero-length reads dropped during extraction, so a shrinking read count is
/// traceable rather than mysterious.
fn log_skipped_empty(skipped: usize) {
    if skipped > 0 {
        info!("Skipped {} zero-length read(s)", skipped);
    }
}

/// Process FASTQ files with minimal information (length only)
fn process_fastq_minimal(file: &Path) -> Result<ReadColumns, NanogetError> {
    let reader = utils::open_file(file)?;
    process_fastq_minimal_from_reader(reader)
}

fn process_fastq_minimal_from_reader<R: Read>(reader: R) -> Result<ReadColumns, NanogetError> {
    use bio::io::fastq;

    let fastq_reader = fastq::Reader::new(reader);
    let mut metrics = ReadColumnsBuilder::with_capacity(INITIAL_READ_CAPACITY);

    let mut skipped_empty = 0usize;

    for (i, result) in fastq_reader.records().enumerate() {
        let record = result.map_err(|e| NanogetError::ParseError(e.to_string()))?;
        check_fastq_record(&record, i)?;

        let length = record.seq().len() as u32;
        if length == 0 {
            skipped_empty += 1;
            continue;
        }
        metrics.push(ReadMetrics::new(None, length));
    }

    log_skipped_empty(skipped_empty);
    Ok(metrics.finish())
}

/// Process FASTA files
fn process_fasta(file: &Path) -> Result<ReadColumns, NanogetError> {
    let reader = utils::open_file(file)?;
    process_fasta_from_reader(reader)
}

fn process_fasta_from_reader<R: Read>(reader: R) -> Result<ReadColumns, NanogetError> {
    use bio::io::fasta;

    let fasta_reader = fasta::Reader::new(reader);
    let mut metrics = ReadColumnsBuilder::with_capacity(INITIAL_READ_CAPACITY);
    let mut skipped_empty = 0usize;

    for result in fasta_reader.records() {
        let record = result.map_err(|e| NanogetError::ParseError(e.to_string()))?;

        let length = record.seq().len() as u32;
        if length == 0 {
            skipped_empty += 1;
            continue;
        }
        metrics.push(ReadMetrics::new(Some(record.id().to_string()), length));
    }

    log_skipped_empty(skipped_empty);
    Ok(metrics.finish())
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
            // A non-finite divergence is not a measurement; treat it as absent so it
            // cannot propagate into the summary statistics.
            Aux::Float(v) => Some(100.0 * (1.0 - v as f64)).filter(|d| d.is_finite()),
            _ => None,
        },
        Err(_) => None,
    }
}

/// Extract aligned length and gap-compressed identity with at most one CIGAR pass.
///
/// When the minimap2 `de` tag is present: one minimal CIGAR pass for aligned length only.
/// When absent: one combined CIGAR pass computing both values simultaneously.
///
/// # Identity is gap-compressed, not BLAST-style
///
/// Identity is `1 - (NM - gap_bases + gap_count) / (matches + gap_count)`, which counts
/// each indel once regardless of its length. This deliberately differs from python
/// nanoget, which reports BLAST-style identity, `1 - NM / (matches + insertions +
/// deletions)`, charging every base of every indel. Gap-compressed identity is the more
/// meaningful measure for nanopore reads, whose errors are dominated by homopolymer
/// indels, and it is the definition behind the minimap2 `de` tag preferred above — so
/// the two sources agree rather than mixing conventions.
///
/// The consequence is that these numbers are **not comparable to python NanoPlot's**.
/// On the nanotest alignment (1115 reads) this crate reports a mean of 89.96 and a
/// minimum of 73.84 where python nanoget reports 86.35 and 49.09.
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
) -> Result<ReadColumns, NanogetError> {
    let mut reader = open_alignment_reader(file, threads)?;

    // For CRAM: tell htslib which fields we actually need so it can skip
    // decompressing the quality and mate-pair streams entirely.
    let is_cram =
        file.extension().and_then(|e| e.to_str()) == Some("cram") || file.as_os_str() == "-"; // stdin CRAM is handled safely — no-op on BAM
    if is_cram {
        #[allow(clippy::arithmetic_side_effects)]
        let fields = sam_fields_SAM_QNAME
            | sam_fields_SAM_FLAG
            | sam_fields_SAM_MAPQ
            | sam_fields_SAM_CIGAR
            | sam_fields_SAM_SEQ
            | sam_fields_SAM_AUX;
        reader
            .set_cram_options(hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS, fields)
            .map_err(|e| NanogetError::ProcessingError(e.to_string()))?;
    }

    extract_bam_records(&mut reader, keep_supplementary)
}

/// Open an alignment file (or stdin) and configure BGZF decompression threading.
///
/// htslib's BGZF threading pre-decompresses blocks on background threads while the main
/// thread processes records. Used by every alignment path — aligned and unaligned alike,
/// since a uBAM is just as BGZF-compressed as a BAM.
fn open_alignment_reader(
    file: &Path,
    threads: usize,
) -> Result<rust_htslib::bam::Reader, NanogetError> {
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
        "Reading {} with {} BGZF thread(s)",
        file.display(),
        bgzf_threads
    );
    Ok(reader)
}

/// Extract ReadMetrics from any type implementing bam::Read.
fn extract_bam_records<R: BamRead>(
    reader: &mut R,
    keep_supplementary: bool,
) -> Result<ReadColumns, NanogetError> {
    let mut metrics = ReadColumnsBuilder::with_capacity(INITIAL_READ_CAPACITY);
    let mut skipped_empty = 0usize;

    for result in reader.records() {
        let record = result?;

        // Secondary alignments are always excluded: they carry no full read
        // sequence (SEQ is '*' or hard-clipped) and would double-count reads.
        if record.is_unmapped() || record.is_secondary() {
            continue;
        }
        // Supplementary alignments are hard-clipped fragments of a read; including
        // them inflates read counts and yield, so they are excluded unless asked for.
        if !keep_supplementary && record.is_supplementary() {
            continue;
        }

        let length = record.seq().len() as u32;
        // A record with no stored SEQ ('*') carries no read to measure.
        if length == 0 {
            skipped_empty += 1;
            continue;
        }

        let read_id = String::from_utf8_lossy(record.qname()).to_string();
        let (aligned_length, percent_identity) = alignment_stats(&record);
        let mapping_quality = if record.mapq() == 255 {
            None
        } else {
            Some(record.mapq())
        };

        metrics.push(ReadMetrics::new(Some(read_id), length).with_alignment(
            aligned_length,
            None,
            mapping_quality,
            percent_identity,
        ));
    }

    log_skipped_empty(skipped_empty);
    Ok(metrics.finish())
}

/// Process unaligned BAM files
fn process_ubam(file: &Path, threads: usize) -> Result<ReadColumns, NanogetError> {
    let mut reader = open_alignment_reader(file, threads)?;
    extract_ubam_records(&mut reader)
}

/// Extract ReadMetrics from an unaligned BAM stream: no alignment fields, but
/// unlike the aligned path the per-read quality is available and is recorded.
fn extract_ubam_records<R: BamRead>(reader: &mut R) -> Result<ReadColumns, NanogetError> {
    let mut metrics = ReadColumnsBuilder::with_capacity(INITIAL_READ_CAPACITY);
    let mut skipped_empty = 0usize;

    for result in reader.records() {
        let record = result?;

        let length = record.seq().len() as u32;
        if length == 0 {
            skipped_empty += 1;
            continue;
        }

        let read_id = String::from_utf8_lossy(record.qname()).to_string();

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

    log_skipped_empty(skipped_empty);
    Ok(metrics.finish())
}

/// Process sequencing summary files
fn process_summary(
    file: &Path,
    read_type: ReadType,
    barcoded: bool,
) -> Result<ReadColumns, NanogetError> {
    let reader = utils::open_file(file)?;
    process_summary_from_reader(reader, read_type, barcoded)
}

fn process_summary_from_reader<R: Read>(
    reader: R,
    read_type: ReadType,
    barcoded: bool,
) -> Result<ReadColumns, NanogetError> {
    use csv::ReaderBuilder;

    let mut csv_reader = ReaderBuilder::new().delimiter(b'\t').from_reader(reader);
    let columns = SummaryColumns::resolve(csv_reader.headers()?, read_type, barcoded)?;

    let mut metrics = ReadColumnsBuilder::with_capacity(INITIAL_READ_CAPACITY);
    let mut skipped_empty = 0usize;

    for result in csv_reader.records() {
        let record = result?;

        let length: u32 = required_field(&record, columns.length, columns.length_name)?
            .parse()
            .map_err(|e| NanogetError::ParseError(format!("Invalid length: {}", e)))?;

        // Only reads with a >0 length are returned, matching python nanoget.
        if length == 0 {
            skipped_empty += 1;
            continue;
        }

        // The column must exist, but a `nan` in it means "not measured" rather than
        // "corrupt file" — some basecallers write one for a read that failed QC — so it
        // is recorded as an absent quality instead of a hard error.
        let quality: Option<f64> = required_field(&record, columns.quality, columns.quality_name)?
            .parse::<f64>()
            .map_err(|e| NanogetError::ParseError(format!("Invalid quality: {}", e)))
            .map(|q| if q.is_finite() { Some(q) } else { None })?;

        // Optional metadata columns. Absence is normal — not every summary file carries
        // every column — but a column that is present and unparseable is a real problem
        // and is reported rather than silently dropped.
        let channel_id: Option<u16> =
            optional_field(&record, columns.channel, "channel", |s| s.parse().ok())?;

        // Shares `parse_start_time` with the rich-FASTQ path, so both accept RFC3339
        // timestamps as well as seconds. Note that ONT writes seconds since the start of
        // the run here, which this (like python nanoget) reads as seconds since the Unix
        // epoch: the value is only meaningful relative to the other reads in the file.
        let start_time =
            optional_field(&record, columns.start_time, "start_time", parse_start_time)?;

        let duration: Option<f64> = optional_field(&record, columns.duration, "duration", |s| {
            s.parse().ok().filter(|d: &f64| d.is_finite())
        })?;

        // Routed through the same helper as the other optional columns, so a blank cell
        // is absence rather than a literal "" barcode in the distribution.
        let barcode: Option<String> =
            optional_field(&record, columns.barcode, "barcode_arrangement", |s| {
                Some(s.to_string())
            })?;

        let mut read_metrics = ReadMetrics::new(None, length)
            .with_sequencing_metadata(channel_id, start_time, duration);
        if let Some(q) = quality {
            read_metrics = read_metrics.with_quality(q);
        }
        read_metrics.barcode = barcode;

        metrics.push(read_metrics);
    }

    log_skipped_empty(skipped_empty);
    Ok(metrics.finish())
}

/// Column positions in a sequencing summary, resolved once from the header row.
///
/// Building a `HashMap<&str, &str>` per row — hashing and inserting every column, then
/// doing six string-keyed lookups — cost about 0.32 s per million rows, roughly a quarter
/// of the time to read a large summary file. Resolving positions once reduces each row to
/// a handful of slice indexes.
///
/// It also moves detection of a missing required column to the header, so a file with the
/// wrong columns fails immediately rather than on its first data row.
struct SummaryColumns {
    length: usize,
    length_name: &'static str,
    quality: usize,
    quality_name: &'static str,
    channel: Option<usize>,
    start_time: Option<usize>,
    duration: Option<usize>,
    barcode: Option<usize>,
}

impl SummaryColumns {
    fn resolve(
        headers: &csv::StringRecord,
        read_type: ReadType,
        barcoded: bool,
    ) -> Result<Self, NanogetError> {
        let (length_name, quality_name) = read_type.summary_columns();
        let index = |name: &str| headers.iter().position(|header| header == name);
        let required = |name: &str| {
            index(name).ok_or_else(|| NanogetError::ParseError(format!("Missing column: {}", name)))
        };

        Ok(Self {
            length: required(length_name)?,
            length_name,
            quality: required(quality_name)?,
            quality_name,
            channel: index("channel"),
            start_time: index("start_time"),
            duration: index("duration"),
            barcode: barcoded.then(|| index("barcode_arrangement")).flatten(),
        })
    }
}

/// Read a required field, which the header guaranteed exists. A row too short to reach it
/// is a malformed file rather than a missing column.
fn required_field<'a>(
    record: &'a csv::StringRecord,
    column: usize,
    name: &str,
) -> Result<&'a str, NanogetError> {
    record.get(column).ok_or_else(|| {
        NanogetError::ParseError(format!("Row is missing a value for column '{}'", name))
    })
}

/// Parse an optional summary field, distinguishing "not present" from "present but
/// unparseable".
///
/// A missing column yields `None`, as does a blank cell — real summary files leave
/// optional cells empty — and so does a non-finite number: `nan` and `inf` mean "not
/// measured", which some basecallers write for a read that failed QC, and which must not
/// reach the summary statistics either way.
///
/// A column that is present with a value `parse` cannot make sense of is an error naming
/// the column and the offending value, so a file that is not what it claims to be says so
/// instead of quietly producing fewer metrics.
fn optional_field<T>(
    record: &csv::StringRecord,
    column: Option<usize>,
    name: &str,
    parse: impl Fn(&str) -> Option<T>,
) -> Result<Option<T>, NanogetError> {
    match column.and_then(|i| record.get(i)) {
        None => Ok(None),
        Some(raw) if raw.trim().is_empty() => Ok(None),
        Some(raw) if is_non_finite_number(raw) => Ok(None),
        Some(raw) => parse(raw).map(Some).ok_or_else(|| {
            NanogetError::ParseError(format!(
                "Invalid value in summary column '{}': '{}'",
                name, raw
            ))
        }),
    }
}

/// True for a value that parses as a number but is not a finite one (`nan`, `inf`,
/// `-inf` and their spellings), i.e. a recorded non-measurement rather than corruption.
fn is_non_finite_number(raw: &str) -> bool {
    raw.trim()
        .parse::<f64>()
        .is_ok_and(|value| !value.is_finite())
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
        // An explicit --file-type overrides detection here too, which is the only way to
        // ask for `fastq-minimal` on a stream.
        match args.file_type {
            Some(ref t) => t.clone(),
            None => FileType::sniff_stdin_bytes(peek)?,
        }
    };

    info!("Reading stdin as {:?}", file_type);

    let reads = match &file_type {
        FileType::Bam | FileType::Cram | FileType::Ubam => {
            // htslib reads from OS fd 0 directly, bypassing the BufReader.
            // Extract the peeked bytes and reconstruct fd 0 via a pipe so htslib
            // sees a complete, untruncated stream.
            let sniffed = stdin_reader.buffer().to_vec();
            drop(stdin_reader);
            reconstruct_stdin_prefix(sniffed)?;
            process_stdin_alignments(args)?
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
                FileType::Fastq => process_fastq_from_reader(reader, false)?,
                FileType::FastqRich => process_fastq_from_reader(reader, true)?,
                FileType::FastqMinimal => process_fastq_minimal_from_reader(reader)?,
                FileType::Fasta => process_fasta_from_reader(reader)?,
                FileType::Summary => {
                    process_summary_from_reader(reader, args.read_type, args.barcoded)?
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

    let collection = MetricsCollection::new(reads);

    info!(
        "Extraction complete: {} reads processed",
        collection.summary().read_count
    );

    if collection.summary().read_count == 0 {
        return Err(NanogetError::ProcessingError(
            "No reads found in input files".to_string(),
        ));
    }

    Ok(collection)
}

/// Read an alignment stream from stdin, choosing the aligned or unaligned extractor
/// from the SAM header.
///
/// Magic bytes alone cannot tell an aligned BAM from an unaligned one, so
/// `sniff_stdin_bytes` reports every BGZF stream as `Bam`. The header settles it once
/// htslib can read it: a uBAM carries no @SQ lines, and every one of its records would
/// otherwise be silently discarded by the unmapped filter in `extract_bam_records`.
fn process_stdin_alignments(args: &ExtractArgs) -> Result<ReadColumns, NanogetError> {
    let mut reader = open_alignment_reader(Path::new("-"), args.threads)?;

    let unaligned = reader.header().target_count() == 0;

    // For CRAM: tell htslib which fields we actually need so it can skip decompressing
    // the streams we never touch. A no-op on BAM. The unaligned path reads the quality
    // scores, the aligned path does not.
    #[allow(clippy::arithmetic_side_effects)]
    let mut fields = sam_fields_SAM_QNAME
        | sam_fields_SAM_FLAG
        | sam_fields_SAM_MAPQ
        | sam_fields_SAM_CIGAR
        | sam_fields_SAM_SEQ
        | sam_fields_SAM_AUX;
    if unaligned {
        #[allow(clippy::arithmetic_side_effects)]
        {
            fields |= sam_fields_SAM_QUAL;
        }
    }
    reader
        .set_cram_options(hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS, fields)
        .map_err(|e| NanogetError::ProcessingError(e.to_string()))?;

    info!(
        "Reading stdin as {}",
        if unaligned {
            "unaligned BAM"
        } else {
            "aligned BAM/CRAM"
        }
    );

    if unaligned {
        extract_ubam_records(&mut reader)
    } else {
        extract_bam_records(&mut reader, args.keep_supplementary)
    }
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

/// Parse a read start time, accepting either an RFC3339 timestamp string
/// (real ONT output, e.g. "2019-12-23T13:44:31Z" and the SAM-style "st" tag)
/// or seconds since the Unix epoch as a float.
fn parse_start_time(value: &str) -> Option<DateTime<Utc>> {
    if let Ok(dt) = DateTime::parse_from_rfc3339(value) {
        return Some(dt.with_timezone(&Utc));
    }
    value.parse::<f64>().ok().and_then(parse_timestamp)
}

/// Parse metadata from rich FASTQ description lines.
///
/// Supports both the legacy "key=value" format (e.g. "ch=123") and the SAM-style
/// "tag:type:value" format introduced in MinKNOW 26.01 (e.g. "ch:i:123"), which also
/// renames fields: start_time -> "st", runid -> "RG" (with the runid as the first
/// underscore-separated component of the read group).
fn parse_rich_fastq_metadata(desc: &str) -> Option<RichFastqMetadata> {
    let mut metadata = RichFastqMetadata {
        channel_id: None,
        start_time: None,
        duration: None,
        run_id: None,
    };

    for field in desc.split_whitespace() {
        if let Some((key, value)) = field.split_once('=') {
            // Legacy albacore/MinKNOW format: key=value
            match key {
                "ch" => {
                    metadata.channel_id = value.parse().ok();
                }
                "start_time" => {
                    metadata.start_time = parse_start_time(value);
                }
                "duration" => {
                    metadata.duration = value.parse().ok().filter(|d: &f64| d.is_finite());
                }
                "runid" => {
                    metadata.run_id = Some(value.to_string());
                }
                _ => {} // Ignore unknown keys
            }
        } else {
            // SAM-style format (MinKNOW >= 26.01): tag:type:value
            let mut parts = field.splitn(3, ':');
            if let (Some(tag), Some(_ty), Some(value)) = (parts.next(), parts.next(), parts.next())
            {
                match tag {
                    "ch" => {
                        metadata.channel_id = value.parse().ok();
                    }
                    "st" => {
                        metadata.start_time = parse_start_time(value);
                    }
                    "du" => {
                        metadata.duration = value.parse().ok().filter(|d: &f64| d.is_finite());
                    }
                    "RG" => {
                        // RG holds "<runid>_<model>@<version>_<barcode>"
                        let runid = value.split('_').next().unwrap_or(value);
                        metadata.run_id = Some(runid.to_string());
                    }
                    _ => {} // Ignore unknown tags
                }
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

    #[test]
    fn test_rich_fastq_metadata_legacy_rfc3339_start_time() {
        let desc = "runid=ff83cfa read=19343 ch=53 start_time=2019-12-23T13:44:31Z";
        let metadata = parse_rich_fastq_metadata(desc).unwrap();

        assert_eq!(metadata.channel_id, Some(53));
        assert_eq!(metadata.run_id, Some("ff83cfa".to_string()));
        assert!(metadata.start_time.is_some());
    }

    #[test]
    fn test_rich_fastq_metadata_sam_format() {
        // MinKNOW >= 26.01 SAM-style "tag:type:value" header
        let desc = "ch:i:123 du:f:1.23 st:Z:2025-01-06T10:06:36.778368+00:00 \
                    RG:Z:e4994c62-93f9-439a-bc8f-d20c95a137a5_rna004_130bps_fast@v5.1.0_barcode02";
        let metadata = parse_rich_fastq_metadata(desc).unwrap();

        assert_eq!(metadata.channel_id, Some(123));
        assert_eq!(metadata.duration, Some(1.23));
        assert_eq!(
            metadata.run_id,
            Some("e4994c62-93f9-439a-bc8f-d20c95a137a5".to_string())
        );
        assert!(metadata.start_time.is_some());
    }
}
