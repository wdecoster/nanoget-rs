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

/// Gzip (and BGZF) magic bytes.
const GZIP_MAGIC: &[u8] = &[0x1f, 0x8b];

/// Bzip2 magic bytes.
const BZIP2_MAGIC: &[u8] = b"BZh";

/// How much of a stream to pull in when classifying it by content. Large enough for the
/// header line of a barcoded sequencing summary, which can carry ~38 columns.
const HEAD_BYTES: usize = 8192;

/// Read up to `HEAD_BYTES` from a stream, stopping early once a full line is available.
///
/// Errors are deliberately swallowed: a stream that cannot be read yields a short (or
/// empty) head, which the caller treats as "unrecognisable" and falls back on.
fn read_head<R: std::io::Read>(mut reader: R) -> Vec<u8> {
    let mut head = Vec::with_capacity(256);
    let mut chunk = [0u8; 256];
    while head.len() < HEAD_BYTES {
        match reader.read(&mut chunk) {
            Ok(0) => break,
            Ok(k) => {
                head.extend_from_slice(&chunk[..k]);
                if head.contains(&b'\n') {
                    break;
                }
            }
            Err(_) => break,
        }
    }
    head
}

/// Rewind a file to the start so its content can be re-read through a decoder.
fn rewind(f: &mut std::fs::File, path: &Path) -> Result<(), NanogetError> {
    use std::io::{Seek, SeekFrom};
    f.seek(SeekFrom::Start(0)).map_err(|e| {
        NanogetError::ParseError(format!("Seek error on {}: {}", path.display(), e))
    })?;
    Ok(())
}

/// Classify the (decompressed) head of a text stream: FASTQ — plain or rich — FASTA, or
/// a sequencing summary. `None` when it matches none of them.
///
/// Shared by the file and stdin sniffers so that the same bytes always classify the same
/// way regardless of how they arrived.
fn classify_text_head(head: &[u8]) -> Option<FileType> {
    match head.first()? {
        b'@' => Some(if first_line_looks_rich(head) {
            FileType::FastqRich
        } else {
            FileType::Fastq
        }),
        b'>' => Some(FileType::Fasta),
        _ => {
            // Sequencing summary: tab-separated with known column headers. Compared
            // lossily so a head truncated mid-character still classifies.
            let end = head.iter().position(|&b| b == b'\n').unwrap_or(head.len());
            let first_line = String::from_utf8_lossy(&head[..end]);
            let cols: Vec<&str> = first_line.split('\t').collect();
            if cols.len() > 1
                && (cols.contains(&"sequence_length_template")
                    || (cols.contains(&"read_id") && cols.contains(&"channel")))
            {
                Some(FileType::Summary)
            } else {
                None
            }
        }
    }
}

impl FileType {
    /// Detect file type by inspecting magic bytes and, when ambiguous, the header content.
    ///
    /// Detection order:
    /// 1. CRAM magic (`CRAM`)
    /// 2. BGZF magic → BAM (or Ubam if the header has no reference sequences)
    /// 3. Plain gzip or bzip2 → decompress the head and classify that, falling back
    ///    to the extension if the content is unrecognisable
    /// 4. Otherwise classify the head directly
    // Public library API (re-exported via `nanoget_rs::FileType`); not yet wired into the binary.
    pub fn sniff(path: &Path) -> Result<Self, NanogetError> {
        use std::fs::File;
        use std::io::Read;

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
        if magic.starts_with(GZIP_MAGIC) {
            // BGZF adds a BC extra subfield at bytes 12-15: [0x42, 0x43, 0x02, 0x00]
            // (SI1='B', SI2='C', SLEN=2 little-endian) with FEXTRA flag set in byte 3.
            if n >= 16 && (magic[3] & 0x04 != 0) && magic[12..16] == [0x42, 0x43, 0x02, 0x00] {
                return sniff_bam_or_ubam(path);
            }
            // Plain gzip. Decompress the head so the inner format is read from the
            // content, not guessed from the extension — otherwise a gzipped rich FASTQ
            // is indistinguishable from a plain one, and a gzipped summary file cannot
            // be recognised at all.
            rewind(&mut f, path)?;
            let head = read_head(flate2::read::GzDecoder::new(&mut f));
            return Self::classify_compressed_head(&head, path);
        }

        // Bzip2 magic: b"BZh"
        if magic.starts_with(BZIP2_MAGIC) {
            rewind(&mut f, path)?;
            let head = read_head(bzip2::read::BzDecoder::new(&mut f));
            return Self::classify_compressed_head(&head, path);
        }

        // Uncompressed text: FASTQ (plain or rich), FASTA or sequencing summary.
        rewind(&mut f, path)?;
        let head = read_head(&mut f);
        classify_text_head(&head).ok_or_else(|| {
            NanogetError::ParseError(format!(
                "Cannot determine file format for: {}\n\
                 Hint: ensure files have a standard extension (.fastq, .bam, .cram, .fasta) \
                 or recognisable content",
                path.display()
            ))
        })
    }

    /// Classify the decompressed head of a compressed file, falling back to the
    /// extension when the content is unrecognisable (a truncated or corrupt stream
    /// still yields a usable answer for a conventionally named file).
    fn classify_compressed_head(head: &[u8], path: &Path) -> Result<Self, NanogetError> {
        classify_text_head(head)
            .or_else(|| Self::from_extension(path))
            .ok_or_else(|| {
                NanogetError::ParseError(format!(
                    "Cannot determine format for compressed file: {} \
                     (use a standard extension like .fastq.gz or .fasta.gz)",
                    path.display()
                ))
            })
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
                "txt" | "tsv" => Self::summary_by_name(path),
                _ => None,
            };
        }

        match extension.as_str() {
            "fastq" | "fq" => Some(Self::Fastq),
            "fasta" | "fa" | "fas" => Some(Self::Fasta),
            "bam" => Some(Self::Bam),
            "cram" => Some(Self::Cram),
            "txt" | "tsv" => Self::summary_by_name(path),
            _ => None,
        }
    }

    /// A `.txt`/`.tsv` file is only assumed to be a sequencing summary when its name
    /// says so; the extension alone carries no format information.
    fn summary_by_name(path: &Path) -> Option<Self> {
        if path.file_name()?.to_str()?.contains("summary") {
            Some(Self::Summary)
        } else {
            None
        }
    }

    /// Returns true for aligned formats (BAM/CRAM).
    // Public library API (re-exported via `nanoget_rs::FileType`); not used by the binary.
    pub fn is_aligned(&self) -> bool {
        matches!(self, Self::Bam | Self::Cram)
    }

    /// Detect format from the first bytes of a stream (no I/O).
    ///
    /// Used for stdin detection where the caller holds the bytes via `BufReader::fill_buf()`.
    /// Reports every BGZF stream as `Bam` (the SAM header cannot be checked without
    /// consuming the stream); the caller resolves aligned vs. unaligned from the header
    /// once htslib has the stream — see `extract::process_stdin_alignments`.
    /// Plain gzip is decompressed (in memory, from the peeked bytes) to reveal the inner format.
    pub fn sniff_stdin_bytes(bytes: &[u8]) -> Result<Self, NanogetError> {
        if bytes.is_empty() {
            return Err(NanogetError::ParseError("Empty stdin".into()));
        }
        if bytes.starts_with(b"CRAM") {
            return Ok(Self::Cram);
        }
        if bytes.starts_with(GZIP_MAGIC) {
            // BGZF: gzip + FEXTRA + BC subfield → BAM
            if bytes.len() >= 16
                && (bytes[3] & 0x04 != 0)
                && bytes[12..16] == [0x42, 0x43, 0x02, 0x00]
            {
                return Ok(Self::Bam);
            }
            // Plain gzip: decompress the head from the peeked buffer to reveal the
            // inner format.
            let head = read_head(flate2::read::GzDecoder::new(bytes));
            return classify_text_head(&head).ok_or_else(|| {
                NanogetError::ParseError(
                    "Gzip stdin does not appear to be FASTQ, FASTA or a sequencing summary".into(),
                )
            });
        }
        classify_text_head(bytes).ok_or_else(|| {
            NanogetError::ParseError(
                "Cannot determine stdin format from magic bytes — \
                 ensure the stream starts with a recognisable header"
                    .into(),
            )
        })
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
