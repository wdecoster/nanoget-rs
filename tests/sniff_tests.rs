//! Content-based format detection, including through compression.
//!
//! `FileType::sniff` previously guessed compressed files from their extension alone,
//! which lost rich-FASTQ metadata for `.fastq.gz` and rejected every `.bz2` file and
//! every gzipped summary outright.

use bzip2::write::BzEncoder;
use flate2::write::GzEncoder;
use nanoget_rs::FileType;
use std::io::Write;
use std::path::{Path, PathBuf};

const RICH_FASTQ: &[u8] =
    b"@r1 runid=abc ch=42 start_time=2019-12-23T13:44:31Z\nACGTACGTAC\n+\n5555555555\n";
const PLAIN_FASTQ: &[u8] = b"@r1\nACGTACGTAC\n+\n5555555555\n";
const FASTA: &[u8] = b">s1\nACGTACGTAC\n";
const SUMMARY: &[u8] = b"filename\tread_id\trun_id\tchannel\tstart_time\tduration\t\
sequence_length_template\tmean_qscore_template\nf1\tr1\tabc\t5\t100.0\t1.5\t1000\t12.5\n";

fn write_plain(dir: &Path, name: &str, data: &[u8]) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, data).expect("write");
    path
}

fn write_gz(dir: &Path, name: &str, data: &[u8]) -> PathBuf {
    let path = dir.join(name);
    let mut enc = GzEncoder::new(
        std::fs::File::create(&path).expect("create"),
        flate2::Compression::default(),
    );
    enc.write_all(data).expect("gz write");
    enc.finish().expect("gz finish");
    path
}

fn write_bz2(dir: &Path, name: &str, data: &[u8]) -> PathBuf {
    let path = dir.join(name);
    let mut enc = BzEncoder::new(
        std::fs::File::create(&path).expect("create"),
        bzip2::Compression::default(),
    );
    enc.write_all(data).expect("bz2 write");
    enc.finish().expect("bz2 finish");
    path
}

/// The same bytes must classify the same way plain, gzipped and bzipped — and the same
/// way through the stdin sniffer.
#[test]
fn test_sniff_is_consistent_across_compression() {
    let dir = tempfile::tempdir().expect("tempdir");
    let cases: [(&str, &[u8], FileType); 4] = [
        ("rich.fastq", RICH_FASTQ, FileType::FastqRich),
        ("plain.fastq", PLAIN_FASTQ, FileType::Fastq),
        ("seqs.fasta", FASTA, FileType::Fasta),
        ("sequencing_summary.txt", SUMMARY, FileType::Summary),
    ];

    for (name, data, expected) in cases {
        let plain = write_plain(dir.path(), name, data);
        let gz = write_gz(dir.path(), &format!("{}.gz", name), data);
        let bz2 = write_bz2(dir.path(), &format!("{}.bz2", name), data);

        assert_eq!(FileType::sniff(&plain).unwrap(), expected, "plain {}", name);
        assert_eq!(FileType::sniff(&gz).unwrap(), expected, "gzipped {}", name);
        assert_eq!(FileType::sniff(&bz2).unwrap(), expected, "bzipped {}", name);
        assert_eq!(
            FileType::sniff_stdin_bytes(data).unwrap(),
            expected,
            "stdin {}",
            name
        );
        assert_eq!(
            FileType::sniff_stdin_bytes(&std::fs::read(&gz).unwrap()).unwrap(),
            expected,
            "gzipped stdin {}",
            name
        );
    }
}

/// Rich metadata must survive gzipping. This is the case that mattered in practice:
/// `.fastq.gz` is the common on-disk form, and guessing it from the extension silently
/// dropped channel, start time and run id.
#[test]
fn test_sniff_detects_rich_metadata_through_gzip() {
    let dir = tempfile::tempdir().expect("tempdir");
    let gz = write_gz(dir.path(), "reads.fastq.gz", RICH_FASTQ);
    assert_eq!(FileType::sniff(&gz).unwrap(), FileType::FastqRich);

    // A plain header in the same container must not be upgraded.
    let plain_gz = write_gz(dir.path(), "plain_reads.fastq.gz", PLAIN_FASTQ);
    assert_eq!(FileType::sniff(&plain_gz).unwrap(), FileType::Fastq);
}

/// When the compressed content is unrecognisable, a conventional extension still wins
/// rather than the whole call failing.
#[test]
fn test_sniff_falls_back_to_extension_for_unrecognisable_content() {
    let dir = tempfile::tempdir().expect("tempdir");
    let gz = write_gz(
        dir.path(),
        "mystery.fastq.gz",
        b"not a fastq header at all\n",
    );
    assert_eq!(FileType::sniff(&gz).unwrap(), FileType::Fastq);

    let bz2 = write_bz2(dir.path(), "mystery.fasta.bz2", b"not a fasta header\n");
    assert_eq!(FileType::sniff(&bz2).unwrap(), FileType::Fasta);
}

/// Content that matches nothing and carries no usable extension is an error, not a guess.
#[test]
fn test_sniff_rejects_unclassifiable_input() {
    let dir = tempfile::tempdir().expect("tempdir");

    let empty = write_plain(dir.path(), "empty.fastq", b"");
    assert!(FileType::sniff(&empty).is_err(), "empty file");

    let junk = write_plain(dir.path(), "junk.dat", b"nothing recognisable here\n");
    assert!(FileType::sniff(&junk).is_err(), "unrecognisable plain file");

    let junk_gz = write_gz(dir.path(), "junk.dat.gz", b"nothing recognisable here\n");
    assert!(FileType::sniff(&junk_gz).is_err(), "unrecognisable gzip");

    assert!(FileType::sniff(dir.path()).is_err(), "directory");
    assert!(
        FileType::sniff(&dir.path().join("absent.fastq")).is_err(),
        "missing file"
    );
}

/// A summary header longer than the old 512-byte probe must still be recognised.
#[test]
fn test_sniff_summary_with_many_columns() {
    let dir = tempfile::tempdir().expect("tempdir");
    let mut header: Vec<String> = (0..40)
        .map(|i| format!("padding_column_{:02}", i))
        .collect();
    header.push("sequence_length_template".to_string());
    header.push("mean_qscore_template".to_string());
    let content = format!("{}\n", header.join("\t"));
    assert!(content.len() > 512);

    let path = write_plain(dir.path(), "summary_wide.txt", content.as_bytes());
    assert_eq!(FileType::sniff(&path).unwrap(), FileType::Summary);
}
