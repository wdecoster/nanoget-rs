//! CLI surface: format detection, the supplementary flag, and value-enum validation.

use assert_cmd::Command;
use std::path::Path;

/// A fixture from the `nanotest` submodule.
///
/// Fails loudly rather than skipping: a silent skip made these tests pass on any machine
/// where the fixtures were missing, which is every CI runner that forgets the submodule.
fn nanotest(name: &str) -> String {
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
    path.to_string_lossy().into_owned()
}

fn run(args: &[&str]) -> std::process::Output {
    Command::cargo_bin("nanoget")
        .expect("binary")
        .arg("extract")
        .args(args)
        .output()
        .expect("run nanoget")
}

/// `--file-type` is optional: formats are detected per file, so a mixed set works in one
/// run and each file is classified on its own content.
#[test]
fn test_file_type_is_optional_and_detected_per_file() {
    let (fastq, bam, summary) = (
        nanotest("reads.fastq.gz"),
        nanotest("alignment.bam"),
        nanotest("sequencing_summary.txt"),
    );

    let counts: Vec<usize> = [&fastq, &bam, &summary]
        .iter()
        .map(|f| {
            let out = run(&["-f", "tsv", f]);
            assert!(out.status.success(), "{}: {:?}", f, out.status.code());
            let stdout = String::from_utf8_lossy(&out.stdout);
            stdout
                .lines()
                .find_map(|l| l.strip_prefix("# Total reads: "))
                .and_then(|n| n.parse().ok())
                .unwrap_or_else(|| panic!("no read count for {}", f))
        })
        .collect();

    // All three together must yield exactly the sum of the three individually.
    let out = run(&["-f", "tsv", &fastq, &bam, &summary]);
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    let combined: usize = stdout
        .lines()
        .find_map(|l| l.strip_prefix("# Total reads: "))
        .and_then(|n| n.parse().ok())
        .expect("no combined read count");
    assert_eq!(combined, counts.iter().sum::<usize>());
}

/// Supplementary alignments are kept by default and excluded by `--drop-supplementary`.
#[test]
fn test_drop_supplementary() {
    let bam = nanotest("alignment.bam");

    let count = |args: &[&str]| -> usize {
        let out = run(args);
        assert!(out.status.success());
        String::from_utf8_lossy(&out.stdout)
            .lines()
            .find_map(|l| l.strip_prefix("# Total reads: "))
            .and_then(|n| n.parse().ok())
            .expect("no read count")
    };

    let kept = count(&["-t", "bam", "-f", "tsv", &bam]);
    let dropped = count(&["-t", "bam", "-f", "tsv", "--drop-supplementary", &bam]);
    assert!(
        dropped < kept,
        "--drop-supplementary kept {} of {}",
        dropped,
        kept
    );
}

/// Invalid values for the enum-typed options are rejected by clap with the allowed set,
/// rather than silently falling through to a default or a Debug dump.
#[test]
fn test_invalid_enum_values_are_rejected() {
    let fastq = nanotest("reads.fastq.gz");

    for (args, expected) in [
        (vec!["-f", "xml", &fastq], "json, tsv"),
        (vec!["--combine", "bogus", &fastq], "simple, track"),
        (vec!["--read-type", "3D", &fastq], "1D, 2D, 1D2"),
        (vec!["-t", "nonsense", &fastq], "fastq"),
    ] {
        let out = run(&args);
        assert!(!out.status.success(), "{:?} unexpectedly succeeded", args);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("invalid value") && stderr.contains(expected),
            "{:?} gave: {}",
            args,
            stderr
        );
    }
}
