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

    // Convert Phred scores to error probabilities
    let error_sum: f64 = qualities
        .iter()
        .map(|&q| 10.0_f64.powf(q as f64 / -10.0))
        .sum();

    let avg_error = error_sum / qualities.len() as f64;

    // Convert back to Phred score
    Some(-10.0 * avg_error.log10())
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

/// Memory-efficient string interning for read IDs and other repeated strings
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

#[derive(Default)]
#[allow(dead_code)]
pub struct StringInterner {
    strings: HashMap<String, Arc<str>>,
}

impl StringInterner {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            strings: HashMap::new(),
        }
    }

    #[allow(dead_code)]
    pub fn intern(&mut self, s: String) -> Arc<str> {
        if let Some(interned) = self.strings.get(&s) {
            interned.clone()
        } else {
            let arc_str: Arc<str> = s.clone().into();
            self.strings.insert(s, arc_str.clone());
            arc_str
        }
    }
}

/// Thread-safe string interner
#[allow(dead_code)]
pub type ThreadSafeInterner = Arc<Mutex<StringInterner>>;

#[allow(dead_code)]
pub fn create_interner() -> ThreadSafeInterner {
    Arc::new(Mutex::new(StringInterner::new()))
}

/// Progress reporting utilities
use indicatif::{ProgressBar, ProgressStyle};

#[allow(dead_code)]
pub fn create_progress_bar(len: u64, message: &str) -> ProgressBar {
    let pb = ProgressBar::new(len);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{msg} [{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {eta}")
            .unwrap()
            .progress_chars("##-"),
    );
    pb.set_message(message.to_string());
    pb
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_average_quality() {
        let qualities = vec![30, 35, 40]; // High quality scores
        let avg = average_quality(&qualities).unwrap();
        assert!(avg > 30.0 && avg < 40.0);

        let empty: Vec<u8> = vec![];
        assert_eq!(average_quality(&empty), None);
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
