use crate::error::NanogetError;
use std::path::Path;

/// Check if a file exists
pub fn check_file_exists(path: &Path) -> Result<(), NanogetError> {
    if !path.exists() {
        return Err(NanogetError::FileNotFound(
            path.to_string_lossy().to_string(),
        ));
    }
    Ok(())
}

/// Calculate average quality from Phred scores
/// Converts Phred scores to error probabilities, calculates average, then back to Phred
pub fn average_quality(qualities: &[u8]) -> Option<f64> {
    if qualities.is_empty() {
        return None;
    }

    // Skip if all values are 255 (missing quality indicator in BAM)
    if qualities.iter().all(|&q| q == 255) {
        return None;
    }

    // Convert Phred scores to error probabilities
    let error_sum: f64 = qualities
        .iter()
        .map(|&q| 10.0_f64.powf(q as f64 / -10.0))
        .sum();

    let avg_error = error_sum / qualities.len() as f64;

    // Convert back to Phred score
    // Handle edge cases: if avg_error is 0 or very small, log10 returns -inf
    // Clamp to a reasonable maximum quality score (e.g., 60)
    let result = -10.0 * avg_error.log10();
    if result.is_nan() || result.is_infinite() || result > 60.0 {
        Some(60.0) // Cap at Q60 for extremely high quality
    } else if result < 0.0 {
        Some(0.0) // Floor at Q0
    } else {
        Some(result)
    }
}

/// Calculate percent identity from CIGAR operations and reference length
#[allow(dead_code)]
pub fn calculate_percent_identity(matches: u32, total_aligned: u32) -> f64 {
    if total_aligned == 0 {
        0.0
    } else {
        (matches as f64 / total_aligned as f64) * 100.0
    }
}

/// Detect compression type from file extension
#[derive(Debug, Clone, Copy)]
pub enum CompressionType {
    None,
    Gzip,
    Bzip2,
    #[allow(dead_code)]
    Bgzip,
}

impl CompressionType {
    pub fn from_path(path: &Path) -> Self {
        let path_str = path.to_string_lossy().to_lowercase();

        if path_str.ends_with(".gz") {
            // Could be gzip or bgzip, we'll assume gzip for now
            Self::Gzip
        } else if path_str.ends_with(".bz2") {
            Self::Bzip2
        } else {
            Self::None
        }
    }
}

/// Open a file with appropriate decompression
pub fn open_file(path: &Path) -> Result<Box<dyn std::io::Read>, NanogetError> {
    use std::fs::File;
    use std::io::BufReader;

    check_file_exists(path)?;

    let file = File::open(path)?;
    let reader = BufReader::new(file);

    match CompressionType::from_path(path) {
        CompressionType::None => Ok(Box::new(reader)),
        CompressionType::Gzip => {
            use flate2::read::GzDecoder;
            Ok(Box::new(GzDecoder::new(reader)))
        }
        CompressionType::Bzip2 => {
            use bzip2::read::BzDecoder;
            Ok(Box::new(BzDecoder::new(reader)))
        }
        CompressionType::Bgzip => {
            // For now, treat bgzip same as gzip
            use flate2::read::GzDecoder;
            Ok(Box::new(GzDecoder::new(reader)))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_average_quality_basic() {
        // Basic test with typical Nanopore quality scores
        let qualities = vec![30, 35, 40];
        let avg = average_quality(&qualities).unwrap();
        assert!(avg > 30.0 && avg < 40.0);
    }

    #[test]
    fn test_average_quality_empty() {
        let empty: Vec<u8> = vec![];
        assert_eq!(average_quality(&empty), None);
    }

    #[test]
    fn test_average_quality_missing_indicator() {
        // 255 is the missing quality indicator in BAM files
        let missing: Vec<u8> = vec![255, 255, 255];
        assert_eq!(average_quality(&missing), None);
    }

    #[test]
    fn test_average_quality_partial_missing() {
        // Mix of valid and missing quality values - should still calculate
        // Note: 255 values are included in calculation, not skipped individually
        let mixed: Vec<u8> = vec![30, 255, 30];
        let result = average_quality(&mixed);
        assert!(result.is_some());
    }

    #[test]
    fn test_average_quality_very_high_capped() {
        // Very high quality scores should be capped at 60
        let very_high: Vec<u8> = vec![93, 93, 93]; // Illumina max Q93
        let avg_high = average_quality(&very_high).unwrap();
        assert!(avg_high <= 60.0);
        assert!(avg_high >= 59.0); // Should be close to cap
    }

    #[test]
    fn test_average_quality_zero_scores() {
        // Zero quality scores (Q0 = 100% error probability)
        let zeros: Vec<u8> = vec![0, 0, 0];
        let avg_zero = average_quality(&zeros).unwrap();
        assert!((0.0..1.0).contains(&avg_zero));
    }

    #[test]
    fn test_average_quality_single_value() {
        // Single quality value should return approximately that value
        let single: Vec<u8> = vec![20];
        let avg = average_quality(&single).unwrap();
        assert!((avg - 20.0).abs() < 0.01);
    }

    #[test]
    fn test_average_quality_uniform_values() {
        // All same values should return that value
        let uniform: Vec<u8> = vec![25, 25, 25, 25, 25];
        let avg = average_quality(&uniform).unwrap();
        assert!((avg - 25.0).abs() < 0.01);
    }

    #[test]
    fn test_average_quality_typical_nanopore() {
        // Typical Nanopore quality distribution (Q7-Q15 range)
        let nanopore: Vec<u8> = vec![8, 10, 12, 9, 11, 13, 10, 8, 14, 12];
        let avg = average_quality(&nanopore).unwrap();
        assert!(avg > 7.0 && avg < 15.0);
    }

    #[test]
    fn test_average_quality_typical_illumina() {
        // Typical Illumina quality distribution (Q30-Q40 range)
        let illumina: Vec<u8> = vec![35, 38, 40, 37, 36, 39, 38, 40, 35, 37];
        let avg = average_quality(&illumina).unwrap();
        assert!(avg > 30.0 && avg < 40.0);
    }

    #[test]
    fn test_average_quality_low_quality_dominates() {
        // One low quality base should significantly lower the average
        // (Phred averaging gives more weight to low quality)
        let with_low: Vec<u8> = vec![40, 40, 40, 40, 5];
        let without_low: Vec<u8> = vec![40, 40, 40, 40, 40];
        let avg_with = average_quality(&with_low).unwrap();
        let avg_without = average_quality(&without_low).unwrap();
        // The single Q5 base should pull the average down significantly
        assert!(avg_with < avg_without - 5.0);
    }

    #[test]
    fn test_average_quality_result_is_finite() {
        // Test various inputs to ensure result is always finite
        let test_cases: Vec<Vec<u8>> = vec![
            vec![0],
            vec![1],
            vec![10],
            vec![20],
            vec![30],
            vec![40],
            vec![50],
            vec![60],
            vec![0, 60],
            vec![1, 2, 3, 4, 5],
            vec![40, 41, 42, 43, 44],
        ];

        for qualities in test_cases {
            let result = average_quality(&qualities);
            assert!(result.is_some(), "Failed for {:?}", qualities);
            let avg = result.unwrap();
            assert!(
                avg.is_finite(),
                "Non-finite result for {:?}: {}",
                qualities,
                avg
            );
            assert!(avg >= 0.0, "Negative result for {:?}: {}", qualities, avg);
            assert!(
                avg <= 60.0,
                "Result exceeds cap for {:?}: {}",
                qualities,
                avg
            );
        }
    }

    #[test]
    fn test_average_quality_large_input() {
        // Test with a large number of quality scores
        let large: Vec<u8> = vec![20; 10000];
        let avg = average_quality(&large).unwrap();
        assert!((avg - 20.0).abs() < 0.01);
    }

    #[test]
    fn test_average_quality_boundary_values() {
        // Test boundary values
        let min_valid: Vec<u8> = vec![0];
        let max_before_missing: Vec<u8> = vec![254];

        assert!(average_quality(&min_valid).is_some());
        assert!(average_quality(&max_before_missing).is_some());

        // 254 is extremely high quality, should be capped
        let avg_254 = average_quality(&max_before_missing).unwrap();
        assert!(avg_254 <= 60.0);
    }

    #[test]
    fn test_percent_identity() {
        assert_eq!(calculate_percent_identity(95, 100), 95.0);
        assert_eq!(calculate_percent_identity(0, 0), 0.0);
        assert_eq!(calculate_percent_identity(100, 100), 100.0);
    }

    #[test]
    fn test_compression_detection() {
        use std::path::Path;

        assert!(matches!(
            CompressionType::from_path(Path::new("test.fastq")),
            CompressionType::None
        ));
        assert!(matches!(
            CompressionType::from_path(Path::new("test.fastq.gz")),
            CompressionType::Gzip
        ));
        assert!(matches!(
            CompressionType::from_path(Path::new("test.fastq.bz2")),
            CompressionType::Bzip2
        ));
    }
}
