//! End-to-end tests for the stdin path, which the library API cannot reach.

use assert_cmd::Command;
use rust_htslib::bam;
use std::path::Path;

/// Write a small unaligned BAM: no @SQ lines, every record flagged unmapped.
fn write_ubam(path: &Path) {
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "unknown"),
    );
    let mut writer =
        bam::Writer::from_path(path, &header, bam::Format::Bam).expect("open ubam writer");

    let seq = b"ACGT".repeat(10);
    let qual = vec![20u8; seq.len()];
    for i in 0..5 {
        let mut record = bam::Record::new();
        record.set(format!("read{}", i).as_bytes(), None, &seq, &qual);
        record.set_flags(4); // unmapped
        record.set_tid(-1);
        record.set_pos(-1);
        writer.write(&record).expect("write ubam record");
    }
}

/// An unaligned BAM piped on stdin must be recognised from its header.
///
/// BGZF magic alone cannot distinguish it from an aligned BAM, and routing it to the
/// aligned extractor discarded every record via the unmapped filter — previously this
/// reported "Total reads: 0" and still exited 0.
#[test]
fn test_ubam_from_stdin_is_detected() {
    let dir = tempfile::tempdir().expect("tempdir");
    let ubam = dir.path().join("unaligned.bam");
    write_ubam(&ubam);
    let bytes = std::fs::read(&ubam).expect("read ubam");

    let output = Command::cargo_bin("nanoget")
        .expect("binary")
        .args(["extract", "-t", "ubam", "-f", "tsv", "-"])
        .write_stdin(bytes)
        .output()
        .expect("run nanoget");

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        output.status.success(),
        "exit {:?}, stderr: {}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        stdout.contains("# Total reads: 5"),
        "stdout was:\n{}",
        stdout
    );
    // The unaligned path records per-read quality; the aligned path does not.
    assert!(
        stdout.contains("\t20.000"),
        "quality missing from:\n{}",
        stdout
    );
}

/// The same uBAM by path must agree with the stdin route.
#[test]
fn test_ubam_from_path_matches_stdin() {
    let dir = tempfile::tempdir().expect("tempdir");
    let ubam = dir.path().join("unaligned.bam");
    write_ubam(&ubam);

    let metrics = nanoget_rs::convenience::extract_from_files(
        vec![ubam],
        nanoget_rs::FileType::Ubam,
        Some(1),
    )
    .expect("extract ubam");

    assert_eq!(metrics.summary().read_count, 5);
    assert!((metrics.reads.get(0).unwrap().quality().expect("quality") - 20.0).abs() < 0.01);
}

/// An empty collection from stdin must be an error, not a silent success — the file
/// path already errors this way.
#[test]
fn test_empty_stdin_errors() {
    let output = Command::cargo_bin("nanoget")
        .expect("binary")
        .args(["extract", "-t", "fastq", "-"])
        .write_stdin("")
        .output()
        .expect("run nanoget");

    assert!(!output.status.success(), "expected a non-zero exit");
}
