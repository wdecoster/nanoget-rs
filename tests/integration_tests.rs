use nanoget_rs::{extract_metrics, CombineMethod, ExtractArgs, FileType, OutputFormat, ReadType};
use std::io::Write;
use tempfile::NamedTempFile;

fn fastq_args(file: &std::path::Path, file_type: FileType) -> ExtractArgs {
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

fn create_test_fastq() -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    writeln!(file, "@read1").unwrap();
    writeln!(file, "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG").unwrap();
    writeln!(file, "+").unwrap();
    writeln!(file, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII").unwrap();
    writeln!(file, "@read2").unwrap();
    writeln!(file, "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT").unwrap();
    writeln!(file, "+").unwrap();
    writeln!(file, "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ").unwrap();
    file
}

#[test]
fn test_fastq_extraction() {
    let temp_file = create_test_fastq();

    let args = ExtractArgs {
        files: vec![temp_file.path().to_path_buf()],
        file_type: Some(FileType::Fastq),
        threads: 1,
        output_format: OutputFormat::Json,
        output: None,
        read_type: ReadType::OneD,
        barcoded: false,
        keep_supplementary: true,
        combine: CombineMethod::Simple,
        names: None,
    };

    let result = extract_metrics(&args).expect("Failed to extract metrics");

    assert_eq!(result.summary().read_count, 2);
    assert_eq!(result.reads.len(), 2);

    // Check first read
    assert_eq!(result.reads.get(0).unwrap().read_id(), Some("read1"));
    assert_eq!(result.reads.get(0).unwrap().length(), 100);
    assert!(result.reads.get(0).unwrap().quality().is_some());

    // Check second read
    assert_eq!(result.reads.get(1).unwrap().read_id(), Some("read2"));
    assert_eq!(result.reads.get(1).unwrap().length(), 99);
    assert!(result.reads.get(1).unwrap().quality().is_some());

    // Check summary stats
    assert_eq!(result.summary().length_stats.count, 2);
    assert!(result.summary().length_stats.mean > 90.0);
    assert!(result.summary().quality_stats.is_some());
}

#[test]
fn test_fastq_minimal() {
    let temp_file = create_test_fastq();

    let args = ExtractArgs {
        files: vec![temp_file.path().to_path_buf()],
        file_type: Some(FileType::FastqMinimal),
        threads: 1,
        output_format: OutputFormat::Json,
        output: None,
        read_type: ReadType::OneD,
        barcoded: false,
        keep_supplementary: true,
        combine: CombineMethod::Simple,
        names: None,
    };

    let result = extract_metrics(&args).expect("Failed to extract metrics");

    assert_eq!(result.summary().read_count, 2);

    // In minimal mode, read IDs should be None
    assert_eq!(result.reads.get(0).unwrap().read_id(), None);
    assert_eq!(result.reads.get(1).unwrap().read_id(), None);
}

fn create_test_fasta() -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    writeln!(file, ">sequence1").unwrap();
    writeln!(file, "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG").unwrap();
    writeln!(file, ">sequence2").unwrap();
    writeln!(file, "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT").unwrap();
    file
}

#[test]
fn test_fasta_extraction() {
    let temp_file = create_test_fasta();

    let args = ExtractArgs {
        files: vec![temp_file.path().to_path_buf()],
        file_type: Some(FileType::Fasta),
        threads: 1,
        output_format: OutputFormat::Json,
        output: None,
        read_type: ReadType::OneD,
        barcoded: false,
        keep_supplementary: true,
        combine: CombineMethod::Simple,
        names: None,
    };

    let result = extract_metrics(&args).expect("Failed to extract metrics");

    assert_eq!(result.summary().read_count, 2);
    assert_eq!(result.reads.len(), 2);

    // Check first sequence
    assert_eq!(result.reads.get(0).unwrap().read_id(), Some("sequence1"));
    assert_eq!(result.reads.get(0).unwrap().length(), 100);
    assert!(result.reads.get(0).unwrap().quality().is_none()); // FASTA has no quality scores

    // Check second sequence
    assert_eq!(result.reads.get(1).unwrap().read_id(), Some("sequence2"));
    assert_eq!(result.reads.get(1).unwrap().length(), 99);
    assert!(result.reads.get(1).unwrap().quality().is_none());
}

#[test]
fn test_multiple_files_combination() {
    let temp_file1 = create_test_fastq();
    let temp_file2 = create_test_fastq();

    let args = ExtractArgs {
        files: vec![
            temp_file1.path().to_path_buf(),
            temp_file2.path().to_path_buf(),
        ],
        file_type: Some(FileType::Fastq),
        threads: 2,
        output_format: OutputFormat::Json,
        output: None,
        read_type: ReadType::OneD,
        barcoded: false,
        keep_supplementary: true,
        combine: CombineMethod::Simple,
        names: None,
    };

    let result = extract_metrics(&args).expect("Failed to extract metrics");

    // Should have reads from both files
    assert_eq!(result.summary().read_count, 4);
    assert_eq!(result.reads.len(), 4);
}

#[test]
fn test_track_combination() {
    let temp_file1 = create_test_fastq();
    let temp_file2 = create_test_fastq();

    let args = ExtractArgs {
        files: vec![
            temp_file1.path().to_path_buf(),
            temp_file2.path().to_path_buf(),
        ],
        file_type: Some(FileType::Fastq),
        threads: 2,
        output_format: OutputFormat::Json,
        output: None,
        read_type: ReadType::OneD,
        barcoded: false,
        keep_supplementary: true,
        combine: CombineMethod::Track,
        names: Some(vec!["sample1".to_string(), "sample2".to_string()]),
    };

    let result = extract_metrics(&args).expect("Failed to extract metrics");

    // Should have reads from both files with dataset tracking
    assert_eq!(result.summary().read_count, 4);
    assert_eq!(result.reads.len(), 4);

    // Check that dataset names are assigned
    let sample1_reads: Vec<_> = result
        .reads
        .iter()
        .filter(|r| r.dataset() == Some("sample1"))
        .collect();
    let sample2_reads: Vec<_> = result
        .reads
        .iter()
        .filter(|r| r.dataset() == Some("sample2"))
        .collect();

    assert_eq!(sample1_reads.len(), 2);
    assert_eq!(sample2_reads.len(), 2);
}

#[test]
fn test_tsv_output_format() {
    let temp_file = create_test_fastq();

    let args = ExtractArgs {
        files: vec![temp_file.path().to_path_buf()],
        file_type: Some(FileType::Fastq),
        threads: 1,
        output_format: OutputFormat::Tsv,
        output: None,
        read_type: ReadType::OneD,
        barcoded: false,
        keep_supplementary: true,
        combine: CombineMethod::Simple,
        names: None,
    };

    let metrics = extract_metrics(&args).expect("Failed to extract metrics");
    let tsv_output = metrics.to_tsv().expect("Failed to generate TSV output");

    // Check TSV format
    assert!(tsv_output.contains("read_id\tlength\tquality")); // Header with tabs
    assert!(tsv_output.contains("read1\t100\t")); // Data with tabs
    assert!(tsv_output.contains("read2\t99\t")); // Data with tabs
    assert!(tsv_output.contains("# Summary Statistics")); // Summary section
    assert!(tsv_output.contains("# Total reads: 2")); // Read count
    assert!(tsv_output.contains("# Length stats")); // Stats header
    assert!(tsv_output.contains("# Quality stats")); // Quality stats since FASTQ has quality
}

/// FASTQ quality lines are ASCII Phred+33; the extracted values must be the decoded
/// scores. Before this was fixed every value came back 33 too high and saturated at
/// the 60 cap, so '#'/'+'/'5' reported 35/43/53 instead of 2/10/20.
#[test]
fn test_fastq_quality_is_phred33_decoded() {
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    for (name, qual_char) in [("q2", '#'), ("q10", '+'), ("q20", '5')] {
        writeln!(file, "@{}", name).unwrap();
        writeln!(file, "{}", "A".repeat(10)).unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "{}", qual_char.to_string().repeat(10)).unwrap();
    }
    file.flush().unwrap();

    let result = extract_metrics(&fastq_args(file.path(), FileType::Fastq))
        .expect("Failed to extract metrics");

    let expected = [("q2", 2.0), ("q10", 10.0), ("q20", 20.0)];
    assert_eq!(result.reads.len(), expected.len());
    for (read, (name, want)) in result.reads.iter().zip(expected) {
        assert_eq!(read.read_id(), Some(name));
        let got = read.quality().expect("quality missing");
        assert!(
            (got - want).abs() < 0.01,
            "{}: got {}, want {}",
            name,
            got,
            want
        );
    }
}

/// A mixed-quality read must average in the probability domain, which is dominated by
/// the low-quality bases: a plain arithmetic mean of the decoded scores would give 18.2.
#[test]
fn test_fastq_quality_averages_in_probability_domain() {
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    writeln!(file, "@mixed").unwrap();
    writeln!(file, "AAAAAAAAAA").unwrap();
    writeln!(file, "+").unwrap();
    writeln!(file, "5555#55555").unwrap(); // nine Q20 bases and one Q2
    file.flush().unwrap();

    let result = extract_metrics(&fastq_args(file.path(), FileType::Fastq))
        .expect("Failed to extract metrics");

    let got = result
        .reads
        .get(0)
        .unwrap()
        .quality()
        .expect("quality missing");
    assert!((got - 11.42).abs() < 0.01, "got {}", got);
}
