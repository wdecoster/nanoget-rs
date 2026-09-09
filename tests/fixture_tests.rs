//! End-to-end tests against the real ONT fixtures in the `nanotest` submodule.
//!
//! These cover the formats that had no end-to-end test at all: BAM, sequencing summary,
//! and compressed input through `utils::open_file`. Expected values are cross-checked
//! against python nanoget where it computes the same quantity.

use nanoget_rs::{
    extract_metrics, CombineMethod, ExtractArgs, FileType, MetricsCollection, NanogetError,
    OutputFormat, ReadType,
};
use std::path::{Path, PathBuf};

/// A fixture from the `nanotest` submodule.
///
/// Fails loudly rather than skipping: a silent skip made these tests pass on any machine
/// where the fixtures were missing, which is every CI runner that forgets the submodule.
fn nanotest(name: &str) -> PathBuf {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("nanotest")
        .join(name);
    assert!(
        path.exists(),
        "missing test fixture {}\n\
         The nanotest fixtures are a git submodule; check them out with:\n\
         \x20   git submodule update --init --depth 1",
        path.display()
    );
    path
}

fn extract(name: &str, file_type: Option<FileType>) -> Result<MetricsCollection, NanogetError> {
    extract_metrics(&ExtractArgs {
        files: vec![nanotest(name)],
        file_type,
        threads: 2,
        output_format: OutputFormat::Json,
        output: None,
        read_type: ReadType::OneD,
        barcoded: false,
        keep_supplementary: true,
        combine: CombineMethod::Simple,
        names: None,
    })
}

/// Formats are detected from content, including through gzip and BGZF.
#[test]
fn test_sniff_recognises_every_fixture() {
    for (name, expected) in [
        ("reads.fastq.gz", FileType::FastqRich),
        ("reads.fa.gz", FileType::Fasta),
        ("reads-mixed-timestamp.fastq", FileType::FastqRich),
        ("sequencing_summary.txt", FileType::Summary),
        ("alignment.bam", FileType::Bam),
        ("alignment_fasta.bam", FileType::Bam),
    ] {
        let got = FileType::sniff(&nanotest(name)).unwrap_or_else(|e| panic!("{}: {}", name, e));
        assert_eq!(got, expected, "{}", name);
    }
}

/// Gzipped rich FASTQ: read count, and the quality values that were wrong until the
/// Phred+33 fix. python nanoget on this exact file reports
/// count=371 mean=10.09 median=10.29 min=7.02 max=13.29.
#[test]
fn test_gzipped_rich_fastq() {
    let metrics = extract("reads.fastq.gz", None).expect("extract");

    assert_eq!(metrics.summary().read_count, 371);

    let quality = metrics.summary().quality_stats.as_ref().expect("quality");
    assert_eq!(quality.count, 371);
    assert!((quality.mean - 10.09).abs() < 0.01, "mean {}", quality.mean);
    assert!((quality.median - 10.29).abs() < 0.01);
    assert!((quality.min - 7.02).abs() < 0.01);
    assert!((quality.max - 13.29).abs() < 0.01);

    // Rich metadata must survive the gzip layer.
    assert!(metrics.reads.iter().all(|r| r.channel_id().is_some()));
    assert!(metrics.reads.iter().all(|r| r.start_time().is_some()));
    assert!(metrics.reads.iter().all(|r| r.run_id().is_some()));
    assert!(metrics.summary().channel_distribution.is_some());
}

/// Gzipped FASTA: same reads, lengths only.
#[test]
fn test_gzipped_fasta() {
    let metrics = extract("reads.fa.gz", None).expect("extract");

    assert_eq!(metrics.summary().read_count, 371);
    assert!(
        metrics.summary().quality_stats.is_none(),
        "FASTA has no quality"
    );
    assert!(metrics.reads.iter().all(|r| r.length() > 0));
}

/// The FASTA and FASTQ fixtures are the same reads, so their length distributions must
/// agree exactly — a good check that neither parser is dropping or miscounting bases.
#[test]
fn test_fasta_and_fastq_fixtures_agree_on_lengths() {
    let fastq = extract("reads.fastq.gz", None).expect("fastq");
    let fasta = extract("reads.fa.gz", None).expect("fasta");

    let mut a: Vec<u32> = fastq.reads.iter().map(|r| r.length()).collect();
    let mut b: Vec<u32> = fasta.reads.iter().map(|r| r.length()).collect();
    a.sort_unstable();
    b.sort_unstable();
    assert_eq!(a, b);
}

/// Aligned BAM: alignment fields are populated, quality deliberately is not (see the
/// "Known differences" section of the README).
#[test]
fn test_aligned_bam() {
    let metrics = extract("alignment.bam", None).expect("extract");

    assert_eq!(metrics.summary().read_count, 1115);
    assert!(metrics.reads.iter().all(|r| r.aligned_length().is_some()));
    assert!(metrics.reads.iter().all(|r| r.quality().is_none()));
    assert!(metrics.summary().quality_stats.is_none());

    // Gap-compressed identity, not python nanoget's BLAST-style identity, which reports
    // mean 86.35 / min 49.09 for the same reads.
    let identity = metrics
        .summary()
        .percent_identity_stats
        .as_ref()
        .expect("identity");
    assert_eq!(identity.count, 1115);
    assert!(
        (identity.mean - 89.96).abs() < 0.01,
        "mean {}",
        identity.mean
    );
    assert!((identity.min - 73.84).abs() < 0.01);
    assert!(identity.max <= 100.0);

    let mapq = metrics
        .summary()
        .mapping_quality_stats
        .as_ref()
        .expect("mapq");
    assert_eq!(mapq.count, 1115);
    assert!(mapq.max <= 60.0);
}

/// Supplementary alignments are hard-clipped fragments; excluding them lowers the count.
#[test]
fn test_bam_supplementary_toggle() {
    let mut args = ExtractArgs {
        files: vec![nanotest("alignment.bam")],
        file_type: Some(FileType::Bam),
        threads: 2,
        output_format: OutputFormat::Json,
        output: None,
        read_type: ReadType::OneD,
        barcoded: false,
        keep_supplementary: true,
        combine: CombineMethod::Simple,
        names: None,
    };
    let kept = extract_metrics(&args).expect("kept").summary().read_count;

    args.keep_supplementary = false;
    let dropped = extract_metrics(&args)
        .expect("dropped")
        .summary()
        .read_count;

    assert_eq!(kept, 1115);
    assert_eq!(dropped, 960);
    // Every remaining read must be a primary alignment, so yield can only fall.
    assert!(dropped < kept);
}

/// Sequencing summary: lengths, qualities, channels and times all come through.
#[test]
fn test_sequencing_summary() {
    let metrics = extract("sequencing_summary.txt", None).expect("extract");

    assert_eq!(metrics.summary().read_count, 371);
    assert!(metrics.reads.iter().all(|r| r.quality().is_some()));
    assert!(metrics.reads.iter().all(|r| r.channel_id().is_some()));
    assert!(metrics.reads.iter().all(|r| r.start_time().is_some()));
    assert!(metrics.reads.iter().all(|r| r.duration().is_some()));
    // Summary files carry no read ids in the columns we read.
    assert!(metrics.reads.iter().all(|r| r.read_id().is_none()));

    let quality = metrics.summary().quality_stats.as_ref().expect("quality");
    assert!((quality.mean - 10.10).abs() < 0.01, "mean {}", quality.mean);
}

/// The summary and the FASTQ describe the same run, so the read counts must match.
#[test]
fn test_summary_and_fastq_fixtures_agree_on_read_count() {
    let summary = extract("sequencing_summary.txt", None).expect("summary");
    let fastq = extract("reads.fastq.gz", None).expect("fastq");
    assert_eq!(summary.summary().read_count, fastq.summary().read_count);
}

/// A rich FASTQ with mixed timestamp spellings (RFC3339 with and without sub-second
/// precision) must parse both.
#[test]
fn test_mixed_timestamp_fastq() {
    let metrics = extract("reads-mixed-timestamp.fastq", None).expect("extract");

    assert_eq!(metrics.summary().read_count, 2);
    assert!(metrics.reads.iter().all(|r| r.start_time().is_some()));
    let times: Vec<_> = metrics
        .reads
        .iter()
        .filter_map(|r| r.start_time())
        .collect();
    assert_ne!(times[0], times[1]);
}

/// Multiple real files of different formats in one run, with dataset tracking.
#[test]
fn test_mixed_formats_with_tracking() {
    let metrics = extract_metrics(&ExtractArgs {
        files: vec![
            nanotest("reads.fastq.gz"),
            nanotest("alignment.bam"),
            nanotest("sequencing_summary.txt"),
        ],
        file_type: None,
        threads: 2,
        output_format: OutputFormat::Json,
        output: None,
        read_type: ReadType::OneD,
        barcoded: false,
        keep_supplementary: true,
        combine: CombineMethod::Track,
        names: Some(vec!["fq".into(), "bam".into(), "sum".into()]),
    })
    .expect("extract");

    assert_eq!(metrics.summary().read_count, 371 + 1115 + 371);
    assert_eq!(metrics.dataset_names(), vec!["bam", "fq", "sum"]);
    assert_eq!(metrics.reads_for_dataset("bam").len(), 1115);
    assert_eq!(metrics.reads_for_dataset("fq").len(), 371);
}

/// JSON output must be reproducible: `BTreeMap` distributions give a stable key order
/// where `HashMap` produced byte-different output on every run.
#[test]
fn test_json_output_is_deterministic() {
    let a = extract("sequencing_summary.txt", None)
        .expect("a")
        .to_json()
        .unwrap();
    let b = extract("sequencing_summary.txt", None)
        .expect("b")
        .to_json()
        .unwrap();
    assert_eq!(a, b);

    // And it round-trips through the crate's own Deserialize.
    let parsed: MetricsCollection = serde_json::from_str(&a).expect("round-trip");
    assert_eq!(parsed.summary().read_count, 371);
}
