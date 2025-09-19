use clap::ValueEnum;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, ValueEnum, Serialize, Deserialize, PartialEq)]
pub enum FileType {
    /// Standard FASTQ file
    Fastq,
    /// FASTQ file with rich metadata (MinKNOW/Albacore format)
    FastqRich,
    /// Minimal FASTQ processing
    FastqMinimal,
    /// FASTA file
    Fasta,
    /// BAM alignment file
    Bam,
    /// CRAM alignment file  
    Cram,
    /// Unaligned BAM file
    Ubam,
    /// Sequencing summary file
    Summary,
}

impl FileType {
    /// Detect file type from extension
    #[allow(dead_code)]
    pub fn from_extension(path: &std::path::Path) -> Option<Self> {
        let extension = path.extension()?.to_str()?.to_lowercase();
        
        match extension.as_str() {
            "fastq" | "fq" => Some(Self::Fastq),
            "fasta" | "fa" | "fas" => Some(Self::Fasta),
            "bam" => Some(Self::Bam),
            "cram" => Some(Self::Cram),
            "txt" | "tsv" => {
                // Check if filename suggests it's a summary file
                if path.file_name()?.to_str()?.contains("summary") {
                    Some(Self::Summary)
                } else {
                    None
                }
            }
            _ => None,
        }
    }
    
    /// Check if the file type supports parallel processing
    #[allow(dead_code)]
    pub fn supports_parallel(&self) -> bool {
        match self {
            Self::Fastq | Self::FastqRich | Self::FastqMinimal | Self::Fasta => true,
            Self::Bam | Self::Cram | Self::Ubam => true,
            Self::Summary => false, // Summary files are typically processed as a whole
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_file_type_detection() {
        assert_eq!(FileType::from_extension(Path::new("test.fastq")), Some(FileType::Fastq));
        assert_eq!(FileType::from_extension(Path::new("test.bam")), Some(FileType::Bam));
        assert_eq!(FileType::from_extension(Path::new("sequencing_summary.txt")), Some(FileType::Summary));
        assert_eq!(FileType::from_extension(Path::new("test.unknown")), None);
    }
}