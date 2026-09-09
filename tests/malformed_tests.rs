//! Corrupt FASTQ must fail loudly, and zero-length reads must be dropped.

use nanoget_rs::{
    extract_metrics, CombineMethod, ExtractArgs, FileType, NanogetError, OutputFormat, ReadType,
};
use std::io::Write;
use tempfile::NamedTempFile;

fn args(file: &std::path::Path, file_type: FileType) -> ExtractArgs {
    ExtractArgs {
        files: vec![file.to_path_buf()],
        file_type: Some(file_type),
        threads: 1,
        output_format: OutputFormat::Json,
        output: None,
        read_type: ReadType::OneD,
        barcoded: false,
        keep_supplementary: true,
        combine: CombineMethod::Simple,
        names: None,
    }
}

fn temp_with(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("temp file");
    file.write_all(content.as_bytes()).expect("write");
    file.flush().expect("flush");
    file
}

/// A sequence and quality line of different lengths is corruption, not something to
/// silently average over the shorter of the two.
#[test]
fn test_mismatched_quality_line_is_an_error() {
    let file = temp_with("@r1\nACGTACGTAC\n+\n55555\n");

    for file_type in [FileType::Fastq, FileType::FastqRich, FileType::FastqMinimal] {
        let err = extract_metrics(&args(file.path(), file_type.clone()))
            .expect_err("expected a parse error");
        let message = err.to_string();
        assert!(
            matches!(err, NanogetError::ParseError(_)),
            "{:?}: wrong variant: {}",
            file_type,
            message
        );
        // The message must name the offending record and both lengths.
        assert!(message.contains("record 1"), "{}", message);
        assert!(message.contains("'r1'"), "{}", message);
        assert!(message.contains("10 bases"), "{}", message);
        assert!(message.contains("5 characters"), "{}", message);
    }
}

/// The offending record must be identified by position, not just the first one.
#[test]
fn test_error_names_the_offending_record() {
    let file = temp_with("@good\nACGT\n+\nIIII\n@bad\nACGTACGT\n+\nII\n");
    let err = extract_metrics(&args(file.path(), FileType::Fastq)).expect_err("expected an error");
    let message = err.to_string();
    assert!(message.contains("record 2"), "{}", message);
    assert!(message.contains("'bad'"), "{}", message);
}

/// Well-formed input must still pass, including trailing CRLF, a description on the
/// header line, and the last record having no trailing newline.
#[test]
fn test_valid_fastq_is_not_rejected() {
    for (label, content) in [
        ("lf", "@r1\nACGT\n+\nIIII\n"),
        ("crlf", "@r1\r\nACGT\r\n+\r\nIIII\r\n"),
        ("with description", "@r1 ch=42 runid=abc\nACGT\n+\nIIII\n"),
        ("no trailing newline", "@r1\nACGT\n+\nIIII"),
    ] {
        let file = temp_with(content);
        let result = extract_metrics(&args(file.path(), FileType::Fastq));
        assert!(result.is_ok(), "{}: {:?}", label, result.err());
        assert_eq!(result.unwrap().summary().read_count, 1, "{}", label);
    }
}

/// Zero-length reads carry no usable metrics and would drag the length distribution to
/// zero; python nanoget drops them and so do we.
#[test]
fn test_zero_length_reads_are_dropped_from_fasta() {
    let file = temp_with(">empty\n\n>real\nACGTACGTAC\n");
    let result = extract_metrics(&args(file.path(), FileType::Fasta)).expect("extract");

    assert_eq!(result.summary().read_count, 1);
    assert_eq!(result.reads.get(0).unwrap().read_id(), Some("real"));
    assert_eq!(result.summary().length_stats.min, 10.0);
}

#[test]
fn test_zero_length_reads_are_dropped_from_summary() {
    let file = temp_with(
        "read_id\tchannel\tstart_time\tduration\tsequence_length_template\tmean_qscore_template\n\
         r1\t5\t100.0\t1.5\t1000\t12.5\n\
         r2\t6\t200.0\t2.0\t0\t0.0\n\
         r3\t7\t300.0\t2.5\t2000\t11.0\n",
    );
    let result = extract_metrics(&args(file.path(), FileType::Summary)).expect("extract");

    assert_eq!(result.summary().read_count, 2);
    assert_eq!(result.summary().length_stats.min, 1000.0);
    assert_eq!(result.summary().length_stats.max, 2000.0);
    // The dropped row's channel must not survive in the distribution either.
    let channels = result
        .summary()
        .channel_distribution
        .as_ref()
        .expect("channels");
    assert!(!channels.contains_key(&6), "channel 6 was dropped with r2");
}

/// A file of nothing but zero-length reads is empty input, not a successful extraction.
#[test]
fn test_all_zero_length_is_an_error() {
    let file = temp_with(">a\n\n>b\n\n");
    let err = extract_metrics(&args(file.path(), FileType::Fasta)).expect_err("expected an error");
    assert!(err.to_string().contains("No reads found"), "{}", err);
}

/// A summary column that is present but unparseable is a real problem: silently dropping
/// it leaves the user with missing plots and no reason given.
#[test]
fn test_unparseable_summary_column_is_an_error() {
    let file = temp_with(
        "read_id\tchannel\tstart_time\tsequence_length_template\tmean_qscore_template\n\
         r1\tnotanumber\t100.0\t1000\t12.5\n",
    );
    let err =
        extract_metrics(&args(file.path(), FileType::Summary)).expect_err("expected an error");
    let message = err.to_string();
    assert!(message.contains("channel"), "{}", message);
    assert!(message.contains("notanumber"), "{}", message);
}

/// A blank cell is absence, not corruption — real summary files leave optional cells
/// empty, and those reads must still be extracted.
#[test]
fn test_blank_summary_cells_are_treated_as_absent() {
    let file = temp_with(
        "read_id\tchannel\tstart_time\tsequence_length_template\tmean_qscore_template\n\
         r1\t\t\t1000\t12.5\n",
    );
    let result = extract_metrics(&args(file.path(), FileType::Summary)).expect("extract");

    assert_eq!(result.summary().read_count, 1);
    assert_eq!(result.reads.get(0).unwrap().channel_id(), None);
    assert_eq!(result.reads.get(0).unwrap().start_time(), None);
}

/// The summary parser shares `parse_start_time` with the rich-FASTQ path, so it accepts
/// RFC3339 timestamps as well as seconds.
#[test]
fn test_summary_accepts_rfc3339_start_time() {
    let file = temp_with(
        "read_id\tchannel\tstart_time\tsequence_length_template\tmean_qscore_template\n\
         r1\t5\t2019-12-23T13:44:31Z\t1000\t12.5\n",
    );
    let result = extract_metrics(&args(file.path(), FileType::Summary)).expect("extract");

    let start = result
        .reads
        .get(0)
        .unwrap()
        .start_time()
        .expect("start_time");
    assert_eq!(start.to_rfc3339(), "2019-12-23T13:44:31+00:00");
}

/// A record in a rich FASTQ with no metadata at all means something is wrong with the
/// file; emitting a read with no channel or start time makes the missing downstream
/// plots unexplainable.
#[test]
fn test_rich_fastq_without_metadata_is_an_error() {
    let file = temp_with(
        "@r1 ch=42 start_time=2019-12-23T13:44:31Z\nACGT\n+\nIIII\n\
         @r2\nACGT\n+\nIIII\n",
    );
    let err =
        extract_metrics(&args(file.path(), FileType::FastqRich)).expect_err("expected an error");
    let message = err.to_string();
    assert!(message.contains("record 2"), "{}", message);
    assert!(message.contains("'r2'"), "{}", message);

    // The same file read as plain FASTQ is fine: the metadata simply is not looked at.
    let ok = extract_metrics(&args(file.path(), FileType::Fastq)).expect("plain fastq");
    assert_eq!(ok.summary().read_count, 2);
}

/// `nan` in a quality column means "not measured", not "corrupt file" — some basecallers
/// write one for a read that failed QC. The read is kept with no quality rather than the
/// run being aborted, and the NaN must not reach the summary statistics.
#[test]
fn test_nan_summary_quality_is_treated_as_absent() {
    let file = temp_with(
        "read_id\tchannel\tstart_time\tsequence_length_template\tmean_qscore_template\n\
         r1\t5\t100.0\t1000\t12.5\n\
         r2\t6\t200.0\t2000\tnan\n",
    );
    let result = extract_metrics(&args(file.path(), FileType::Summary)).expect("extract");

    assert_eq!(result.summary().read_count, 2, "both reads are kept");
    assert_eq!(result.reads.get(1).unwrap().quality(), None);

    let quality = result
        .summary()
        .quality_stats
        .as_ref()
        .expect("quality stats");
    assert_eq!(quality.count, 1);
    assert!(quality.mean.is_finite());
}

/// A `nan` start_time must not silently become the Unix epoch (`NaN as i64` saturates
/// to 0), which is what happened before the finiteness guard.
#[test]
fn test_nan_summary_start_time_is_absent_not_epoch() {
    let file = temp_with(
        "read_id\tchannel\tstart_time\tsequence_length_template\tmean_qscore_template\n\
         r1\t5\tnan\t1000\t12.5\n",
    );
    let result = extract_metrics(&args(file.path(), FileType::Summary)).expect("extract");
    assert_eq!(result.reads.get(0).unwrap().start_time(), None);
}
