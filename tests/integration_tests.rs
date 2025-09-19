use nanoget_rs::{extract_metrics, ExtractArgs, FileType};
use tempfile::NamedTempFile;
use std::io::Write;

fn create_test_fastq() -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    writeln!(file, "@read1").unwrap();
    writeln!(file, "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG").unwrap();
    writeln!(file, "+").unwrap();
    writeln!(file, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII").unwrap();
    writeln!(file, "@read2").unwrap();
    writeln!(file, "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT").unwrap();
    writeln!(file, "+").unwrap();
    writeln!(file, "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ").unwrap();
    file
}

#[test]
fn test_fastq_extraction() {
    let temp_file = create_test_fastq();
    
    let args = ExtractArgs {
        files: vec![temp_file.path().to_path_buf()],
        file_type: FileType::Fastq,
        threads: 1,
        output_format: "json".to_string(),
        output: None,
        read_type: "1D".to_string(),
        barcoded: false,
        keep_supplementary: true,
        huge: false,
        combine: "simple".to_string(),
        names: None,
    };
    
    let result = extract_metrics(&args).expect("Failed to extract metrics");
    
    assert_eq!(result.summary.read_count, 2);
    assert_eq!(result.reads.len(), 2);
    
    // Check first read
    assert_eq!(result.reads[0].read_id, Some("read1".to_string()));
    assert_eq!(result.reads[0].length, 100);
    assert!(result.reads[0].quality.is_some());
    
    // Check second read  
    assert_eq!(result.reads[1].read_id, Some("read2".to_string()));
    assert_eq!(result.reads[1].length, 99);
    assert!(result.reads[1].quality.is_some());
    
    // Check summary stats
    assert_eq!(result.summary.length_stats.count, 2);
    assert!(result.summary.length_stats.mean > 90.0);
    assert!(result.summary.quality_stats.is_some());
}

#[test]
fn test_fastq_minimal() {
    let temp_file = create_test_fastq();
    
    let args = ExtractArgs {
        files: vec![temp_file.path().to_path_buf()],
        file_type: FileType::FastqMinimal,
        threads: 1,
        output_format: "json".to_string(),
        output: None,
        read_type: "1D".to_string(),
        barcoded: false,
        keep_supplementary: true,
        huge: false,
        combine: "simple".to_string(),
        names: None,
    };
    
    let result = extract_metrics(&args).expect("Failed to extract metrics");
    
    assert_eq!(result.summary.read_count, 2);
    
    // In minimal mode, read IDs should be None
    assert_eq!(result.reads[0].read_id, None);
    assert_eq!(result.reads[1].read_id, None);
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
        file_type: FileType::Fasta,
        threads: 1,
        output_format: "json".to_string(),
        output: None,
        read_type: "1D".to_string(),
        barcoded: false,
        keep_supplementary: true,
        huge: false,
        combine: "simple".to_string(),
        names: None,
    };
    
    let result = extract_metrics(&args).expect("Failed to extract metrics");
    
    assert_eq!(result.summary.read_count, 2);
    assert_eq!(result.reads.len(), 2);
    
    // Check first sequence
    assert_eq!(result.reads[0].read_id, Some("sequence1".to_string()));
    assert_eq!(result.reads[0].length, 100);
    assert!(result.reads[0].quality.is_none()); // FASTA has no quality scores
    
    // Check second sequence
    assert_eq!(result.reads[1].read_id, Some("sequence2".to_string()));
    assert_eq!(result.reads[1].length, 99);
    assert!(result.reads[1].quality.is_none());
}

#[test]
fn test_multiple_files_combination() {
    let temp_file1 = create_test_fastq();
    let temp_file2 = create_test_fastq();
    
    let args = ExtractArgs {
        files: vec![temp_file1.path().to_path_buf(), temp_file2.path().to_path_buf()],
        file_type: FileType::Fastq,
        threads: 2,
        output_format: "json".to_string(),
        output: None,
        read_type: "1D".to_string(),
        barcoded: false,
        keep_supplementary: true,
        huge: false,
        combine: "simple".to_string(),
        names: None,
    };
    
    let result = extract_metrics(&args).expect("Failed to extract metrics");
    
    // Should have reads from both files
    assert_eq!(result.summary.read_count, 4);
    assert_eq!(result.reads.len(), 4);
}

#[test]
fn test_track_combination() {
    let temp_file1 = create_test_fastq();
    let temp_file2 = create_test_fastq();
    
    let args = ExtractArgs {
        files: vec![temp_file1.path().to_path_buf(), temp_file2.path().to_path_buf()],
        file_type: FileType::Fastq,
        threads: 2,
        output_format: "json".to_string(),
        output: None,
        read_type: "1D".to_string(),
        barcoded: false,
        keep_supplementary: true,
        huge: false,
        combine: "track".to_string(),
        names: Some(vec!["sample1".to_string(), "sample2".to_string()]),
    };
    
    let result = extract_metrics(&args).expect("Failed to extract metrics");
    
    // Should have reads from both files with dataset tracking
    assert_eq!(result.summary.read_count, 4);
    assert_eq!(result.reads.len(), 4);
    
    // Check that dataset names are assigned
    let sample1_reads: Vec<_> = result.reads.iter()
        .filter(|r| r.dataset == Some("sample1".to_string()))
        .collect();
    let sample2_reads: Vec<_> = result.reads.iter()
        .filter(|r| r.dataset == Some("sample2".to_string()))
        .collect();
    
    assert_eq!(sample1_reads.len(), 2);
    assert_eq!(sample2_reads.len(), 2);
}