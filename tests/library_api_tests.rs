// Integration test demonstrating library usage

use nanoget_rs::{convenience::*, MetricsCollection, ReadMetrics};
use tempfile::NamedTempFile;
use std::io::Write;

#[test]
fn test_convenience_api() {
    // Create a test FASTQ file
    let mut temp_file = NamedTempFile::new().expect("Failed to create temp file");
    writeln!(temp_file, "@read1").unwrap();
    writeln!(temp_file, "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG").unwrap();
    writeln!(temp_file, "+").unwrap();
    writeln!(temp_file, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII").unwrap();
    writeln!(temp_file, "@read2").unwrap();
    writeln!(temp_file, "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT").unwrap();
    writeln!(temp_file, "+").unwrap();
    writeln!(temp_file, "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ").unwrap();
    
    // Test convenience function
    let metrics = extract_from_fastq(temp_file.path()).expect("Failed to extract metrics");
    
    assert_eq!(metrics.summary.read_count, 2);
    assert_eq!(metrics.reads.len(), 2);
    
    // Test new collection methods
    let high_quality_reads = metrics.filter_by_quality(70.0);
    assert!(high_quality_reads.summary.read_count <= 2);
    
    let long_reads = metrics.filter_by_length(50);
    assert_eq!(long_reads.summary.read_count, 2); // Both reads are longer than 50bp
    
    // Test JSON export
    let json_output = metrics.to_json().expect("Failed to export JSON");
    assert!(json_output.contains("read_count"));
    assert!(json_output.contains("length_stats"));
    
    let compact_json = metrics.to_json_compact().expect("Failed to export compact JSON");
    assert!(compact_json.len() < json_output.len()); // Compact should be shorter
}

#[test]
fn test_dataset_functionality() {
    // Test the dataset methods even though we can't easily create tracked datasets in this test
    let reads = vec![
        ReadMetrics::new(Some("read1".to_string()), 100),
        ReadMetrics::new(Some("read2".to_string()), 200),
    ];
    
    let mut collection = MetricsCollection::new(reads);
    
    // Manually add dataset names to test the functionality
    collection.reads[0].dataset = Some("Sample1".to_string());
    collection.reads[1].dataset = Some("Sample2".to_string());
    
    let dataset_names = collection.dataset_names();
    assert_eq!(dataset_names, vec!["Sample1", "Sample2"]);
    
    let sample1_reads = collection.reads_for_dataset("Sample1");
    assert_eq!(sample1_reads.len(), 1);
    assert_eq!(sample1_reads[0].read_id, Some("read1".to_string()));
}