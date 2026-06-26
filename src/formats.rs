use crate::error::NanogetError;
use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use std::path::Path;

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
    /// Detect file type by inspecting magic bytes and, when ambiguous, the header content.
    ///
    /// Detection order:
    /// 1. CRAM magic (`CRAM`)
    /// 2. BGZF magic → BAM (or Ubam if the header has no reference sequences)
    /// 3. Plain gzip → fall back to extension
    /// 4. `@` first byte → FASTQ
    /// 5. `>` first byte → FASTA
    /// 6. Tab-separated first line with known summary columns → Summary
    // Public library API (re-exported via `nanoget_rs::FileType`); not yet wired into the binary.
    #[allow(dead_code)]
    pub fn sniff(path: &Path) -> Result<Self, NanogetError> {
        use std::fs::File;
        use std::io::{Read, Seek, SeekFrom};

        let mut f = File::open(path)
            .map_err(|_| NanogetError::FileNotFound(path.to_string_lossy().to_string()))?;

        let mut magic = [0u8; 16];
        let n = f.read(&mut magic).map_err(|e| {
            NanogetError::ParseError(format!("Cannot read {}: {}", path.display(), e))
        })?;
        let magic = &magic[..n];

        if n == 0 {
            return Err(NanogetError::ParseError(format!(
                "Empty file: {}",
                path.display()
            )));
        }

        // CRAM magic: b"CRAM"
        if magic.starts_with(b"CRAM") {
            return Ok(Self::Cram);
        }

        // Gzip / BGZF magic: 0x1f 0x8b
        if magic.starts_with(&[0x1f, 0x8b]) {
            // BGZF adds a BC extra subfield at bytes 12-15: [0x42, 0x43, 0x02, 0x00]
            // (SI1='B', SI2='C', SLEN=2 little-endian) with FEXTRA flag set in byte 3.
            if n >= 16 && (magic[3] & 0x04 != 0) && magic[12..16] == [0x42, 0x43, 0x02, 0x00] {
                return sniff_bam_or_ubam(path);
            }
            // Plain gzip — distinguish FASTQ vs FASTA via extension
            return Self::from_extension(path).ok_or_else(|| {
                NanogetError::ParseError(format!(
                    "Cannot determine format for gzipped file: {} \
                     (use a standard extension like .fastq.gz or .fasta.gz)",
                    path.display()
                ))
            });
        }

        // FASTQ: records start with '@'. Peek the first header line to tell a
        // rich (MinKNOW/albacore metadata) FASTQ apart from a plain one.
        if magic[0] == b'@' {
            f.seek(SeekFrom::Start(0)).map_err(|e| {
                NanogetError::ParseError(format!("Seek error on {}: {}", path.display(), e))
            })?;
            let mut buf = [0u8; 1024];
            let bn = f.read(&mut buf).unwrap_or(0);
            if first_line_looks_rich(&buf[..bn]) {
                return Ok(Self::FastqRich);
            }
            return Ok(Self::Fastq);
        }

        // FASTA: records start with '>'
        if magic[0] == b'>' {
            return Ok(Self::Fasta);
        }

        // Sequencing summary: tab-separated with known column headers
        f.seek(SeekFrom::Start(0)).map_err(|e| {
            NanogetError::ParseError(format!("Seek error on {}: {}", path.display(), e))
        })?;
        let mut header_buf = [0u8; 512];
        let hn = f.read(&mut header_buf).unwrap_or(0);
        if let Ok(text) = std::str::from_utf8(&header_buf[..hn]) {
            let first_line = text.lines().next().unwrap_or("");
            if first_line.contains('\t') {
                let cols: Vec<&str> = first_line.split('\t').collect();
                if cols.contains(&"sequence_length_template")
                    || (cols.contains(&"read_id") && cols.contains(&"channel"))
                {
                    return Ok(Self::Summary);
                }
            }
        }

        Err(NanogetError::ParseError(format!(
            "Cannot determine file format for: {}\n\
             Hint: ensure files have a standard extension (.fastq, .bam, .cram, .fasta) \
             or recognisable content",
            path.display()
        )))
    }

    /// Detect file type from extension, including compressed variants (.gz, .bz2).
    pub fn from_extension(path: &Path) -> Option<Self> {
        let extension = path.extension()?.to_str()?.to_lowercase();

        // Strip one layer of compression to get the inner extension
        if matches!(extension.as_str(), "gz" | "bz2") {
            let stem = path.file_stem()?;
            let inner_ext = Path::new(stem).extension()?.to_str()?.to_lowercase();
            return match inner_ext.as_str() {
                "fastq" | "fq" => Some(Self::Fastq),
                "fasta" | "fa" | "fas" => Some(Self::Fasta),
                "bam" => Some(Self::Bam),
                _ => None,
            };
        }

        match extension.as_str() {
            "fastq" | "fq" => Some(Self::Fastq),
            "fasta" | "fa" | "fas" => Some(Self::Fasta),
            "bam" => Some(Self::Bam),
            "cram" => Some(Self::Cram),
            "txt" | "tsv" => {
                if path.file_name()?.to_str()?.contains("summary") {
                    Some(Self::Summary)
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    /// Returns true for aligned formats (BAM/CRAM).
    // Public library API (re-exported via `nanoget_rs::FileType`); not used by the binary.
    #[allow(dead_code)]
    pub fn is_aligned(&self) -> bool {
        matches!(self, Self::Bam | Self::Cram)
    }

    /// Detect format from the first bytes of a stream (no I/O).
    ///
    /// Used for stdin detection where the caller holds the bytes via `BufReader::fill_buf()`.
    /// Assumes BGZF → aligned BAM (cannot check the SAM header without consuming the stream).
    /// Plain gzip is decompressed (in memory, from the peeked bytes) to reveal the inner format.
    pub fn sniff_stdin_bytes(bytes: &[u8]) -> Result<Self, NanogetError> {
        use std::io::Read;

        if bytes.is_empty() {
            return Err(NanogetError::ParseError("Empty stdin".into()));
        }
        if bytes.starts_with(b"CRAM") {
            return Ok(Self::Cram);
        }
        if bytes.starts_with(&[0x1f, 0x8b]) {
            // BGZF: gzip + FEXTRA + BC subfield → BAM
            if bytes.len() >= 16
                && (bytes[3] & 0x04 != 0)
                && bytes[12..16] == [0x42, 0x43, 0x02, 0x00]
            {
                return Ok(Self::Bam);
            }
            // Plain gzip (FASTQ or FASTA): decompress the first line from the
            // peeked buffer to reveal the inner format (and, for FASTQ, whether
            // the header carries rich metadata).
            let mut gz = flate2::read::GzDecoder::new(bytes);
            let mut head = Vec::new();
            let mut chunk = [0u8; 256];
            loop {
                match gz.read(&mut chunk) {
                    Ok(0) => break,
                    Ok(k) => {
                        head.extend_from_slice(&chunk[..k]);
                        if head.contains(&b'\n') || head.len() >= 1024 {
                            break;
                        }
                    }
                    Err(e) => {
                        return Err(NanogetError::ParseError(format!(
                            "Failed to decompress gzip stdin: {}",
                            e
                        )))
                    }
                }
            }
            return match head.first() {
                Some(b'@') => Ok(if first_line_looks_rich(&head) {
                    Self::FastqRich
                } else {
                    Self::Fastq
                }),
                Some(b'>') => Ok(Self::Fasta),
                _ => Err(NanogetError::ParseError(
                    "Gzip stdin does not appear to be FASTQ or FASTA".into(),
                )),
            };
        }
        if bytes[0] == b'@' {
            if first_line_looks_rich(bytes) {
                return Ok(Self::FastqRich);
            }
            return Ok(Self::Fastq);
        }
        if bytes[0] == b'>' {
            return Ok(Self::Fasta);
        }
        if let Ok(text) = std::str::from_utf8(bytes) {
            let first_line = text.lines().next().unwrap_or("");
            if first_line.contains('\t') {
                let cols: Vec<&str> = first_line.split('\t').collect();
                if cols.contains(&"sequence_length_template")
                    || (cols.contains(&"read_id") && cols.contains(&"channel"))
                {
                    return Ok(Self::Summary);
                }
            }
        }
        Err(NanogetError::ParseError(
            "Cannot determine stdin format from magic bytes — \
             ensure the stream starts with a recognisable header"
                .into(),
        ))
    }
}

/// Extract the first line from a byte buffer and test it for rich-FASTQ metadata.
fn first_line_looks_rich(bytes: &[u8]) -> bool {
    let end = bytes
        .iter()
        .position(|&b| b == b'\n')
        .unwrap_or(bytes.len());
    header_looks_rich(&String::from_utf8_lossy(&bytes[..end]))
}

/// True when a FASTQ header carries MinKNOW/albacore metadata in its description
/// — the same `key=value` (legacy) or `tag:type:value` (MinKNOW >= 26.01) fields
/// the rich-FASTQ reader parses. Used to auto-detect `FastqRich`.
fn header_looks_rich(header: &str) -> bool {
    // The description is everything after the read id (the first whitespace).
    let Some((_id, desc)) = header.split_once(char::is_whitespace) else {
        return false;
    };
    desc.split_whitespace().any(|field| {
        if let Some((key, _)) = field.split_once('=') {
            matches!(key, "ch" | "start_time" | "duration" | "runid")
        } else {
            let mut parts = field.splitn(3, ':');
            matches!(
                (parts.next(), parts.next(), parts.next()),
                (Some("ch"), Some(_), Some(_))
                    | (Some("st"), Some(_), Some(_))
                    | (Some("du"), Some(_), Some(_))
                    | (Some("RG"), Some(_), Some(_))
            )
        }
    })
}

/// Open the BAM header to distinguish aligned BAM from unaligned BAM (no @SQ lines).
fn sniff_bam_or_ubam(path: &Path) -> Result<FileType, NanogetError> {
    use rust_htslib::bam::{self, Read};
    let reader = bam::Reader::from_path(path)
        .map_err(|e| NanogetError::ParseError(format!("Cannot open {}: {}", path.display(), e)))?;
    if reader.header().target_count() == 0 {
        Ok(FileType::Ubam)
    } else {
        Ok(FileType::Bam)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_header_looks_rich() {
        // Legacy MinKNOW/albacore key=value
        assert!(header_looks_rich(
            "@read1 runid=abc ch=42 start_time=2020-01-01T00:00:00Z"
        ));
        assert!(header_looks_rich("@read1 ch=42"));
        // SAM-style tag:type:value (MinKNOW >= 26.01)
        assert!(header_looks_rich(
            "@read1 st:Z:2026-01-01T00:00:00Z ch:i:42"
        ));
        assert!(header_looks_rich("@read1 RG:Z:runid_model@v_barcode"));
        // Plain headers must not be mistaken for rich
        assert!(!header_looks_rich("@read1"));
        assert!(!header_looks_rich("@read1 some free-text description"));
        assert!(!header_looks_rich("@SRR123.1 1 length=1000"));
    }

    #[test]
    fn test_first_line_looks_rich() {
        let rich = b"@read1 ch=42 start_time=2020-01-01T00:00:00Z\nACGT\n+\n!!!!\n";
        assert!(first_line_looks_rich(rich));
        let plain = b"@read1\nACGT\n+\n!!!!\n";
        assert!(!first_line_looks_rich(plain));
    }

    #[test]
    fn test_file_type_detection() {
        assert_eq!(
            FileType::from_extension(Path::new("test.fastq")),
            Some(FileType::Fastq)
        );
        assert_eq!(
            FileType::from_extension(Path::new("test.fastq.gz")),
            Some(FileType::Fastq)
        );
        assert_eq!(
            FileType::from_extension(Path::new("test.bam")),
            Some(FileType::Bam)
        );
        assert_eq!(
            FileType::from_extension(Path::new("test.cram")),
            Some(FileType::Cram)
        );
        assert_eq!(
            FileType::from_extension(Path::new("sequencing_summary.txt")),
            Some(FileType::Summary)
        );
        assert_eq!(FileType::from_extension(Path::new("test.unknown")), None);
    }
}
