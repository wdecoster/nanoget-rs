# nanoget-rs

A Rust port of [nanoget](https://github.com/wdecoster/nanoget), a tool for extracting metrics from Oxford Nanopore sequencing data and alignments.

## Features

nanoget-rs can extract metrics from various sequencing file formats:

- **FASTQ files** (standard, rich metadata, minimal processing)
- **FASTA files**
- **BAM/SAM/CRAM files** (aligned reads)
- **uBAM files** (unaligned reads)
- **Sequencing summary files** (from Albacore/Guppy/Dorado)

## Installation

### From source

```bash
git clone https://github.com/wdecoster/nanoget-rs
cd nanoget-rs
cargo build --release
```

The binary will be available at `target/release/nanoget`.

### As a Rust library

Add to your `Cargo.toml`:

```toml
[dependencies]
nanoget-rs = { git = "https://github.com/wdecoster/nanoget-rs" }
# Or from crates.io once published:
# nanoget-rs = "0.1.0"
```

## Usage

### Basic usage

The format of each input file is detected from its content, so you normally just point
nanoget at the files:

```bash
nanoget extract reads.fastq.gz
nanoget extract alignments.bam
nanoget extract sequencing_summary.txt
```

Detection covers FASTQ (plain and rich), FASTA, BAM, CRAM, unaligned BAM and sequencing
summaries, through gzip and bzip2 compression, and works on `-` (stdin) too:

```bash
samtools view -b aln.sam | nanoget extract -
```

Because each file is classified on its own content, a mixed set works in one run:

```bash
nanoget extract reads.fastq.gz alignments.bam sequencing_summary.txt
```

### Overriding detection

`-t` / `--file-type` overrides detection. It is needed only to force a format, or to
select `fastq-minimal`, which is a processing mode (lengths only, no read ids or
qualities) rather than a format and so cannot be detected:

```bash
nanoget extract -t fastq-minimal huge_reads.fastq.gz
```

### Output formats

By default, output is in JSON format. You can also specify TSV:
```bash
nanoget extract reads.fastq -f tsv
```

Save output to a file:
```bash
nanoget extract reads.fastq -o metrics.json
```

### Processing multiple files

Process multiple files and combine results:
```bash
nanoget extract file1.fastq file2.fastq file3.fastq
```

Track datasets separately:
```bash
nanoget extract file1.fastq file2.fastq --combine track --names sample1 sample2
```

### Advanced options

Use multiple threads:
```bash
nanoget extract reads.fastq -j 8
```

For BAM/CRAM, supplementary alignments are included by default. They are hard-clipped
fragments of a read, so counting them inflates read counts and yield; exclude them with:
```bash
nanoget extract alignments.bam --drop-supplementary
```

For summary files, specify read type and barcode analysis:
```bash
nanoget extract sequencing_summary.txt --read-type 1D --barcoded
```

### Strictness

Input that is malformed rather than merely unusual is a hard error, not a warning —
corrupt FASTQ is common and a plausible-looking wrong answer is worse than a stop:

- a FASTQ record whose sequence and quality lines differ in length, or that has no read id
- a record in a rich FASTQ carrying no MinKNOW/albacore metadata in its header
- a sequencing-summary column that is present but holds an unparseable value (a column
  that is simply absent, or a blank cell, is fine)

Zero-length reads are dropped from all formats, matching python nanoget; the number
dropped is reported at `info` log level (`RUST_LOG=info`).

## Library Usage

nanoget-rs can be used as a Rust library for integration into other tools. This is generally **preferred over calling the executable** because it:

- Avoids subprocess overhead
- Provides type-safe access to data structures
- Enables direct manipulation of metrics without JSON parsing
- Allows custom analysis and filtering

### Simple API

```rust
use nanoget_rs::convenience::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Simple one-liners for common use cases
    let metrics = extract_from_fastq("reads.fastq")?;
    
    println!("Found {} reads", metrics.summary.read_count);
    println!("Mean length: {:.0} bp", metrics.summary.length_stats.mean);
    
    // Filter and analyze
    let high_quality = metrics.filter_by_quality(30.0);
    let long_reads = metrics.filter_by_length(1000);
    
    Ok(())
}
```

### Advanced API

```rust
use nanoget_rs::{extract_metrics, ExtractArgs, FileType, MetricsCollection};
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = ExtractArgs {
        files: vec![PathBuf::from("sample1.fastq"), PathBuf::from("sample2.fastq")],
        file_type: FileType::Fastq,
        threads: 8,
        combine: "track".to_string(),
        names: Some(vec!["Control".to_string(), "Treatment".to_string()]),
        // ... other options
    };
    
    let metrics = extract_metrics(&args)?;
    
    // Analyze by dataset
    for dataset in metrics.dataset_names() {
        let reads = metrics.reads_for_dataset(&dataset);
        println!("{}: {} reads", dataset, reads.len());
    }
    
    // Export results
    let json_output = metrics.to_json()?;
    std::fs::write("results.json", json_output)?;
    
    Ok(())
}
```

### When to Use Library vs Executable

**Use the library when:**
- Integrating into existing Rust applications
- Need custom analysis or filtering
- Processing many files programmatically
- Want type-safe access to metrics
- Building pipelines or workflows

**Use the executable when:**
- Simple command-line analysis
- Shell scripting
- One-off data exploration
- Interfacing from non-Rust languages

## Output

nanoget-rs outputs comprehensive metrics including:

- **Read-level metrics**: length, quality scores, alignment statistics
- **Summary statistics**: mean, median, standard deviation, quartiles
- **Distributions**: channel usage, barcode distributions (when applicable)
- **Time-based analysis**: sequencing start times and duration (when available)

Example output structure:
```json
{
  "reads": [
    {
      "read_id": "read_001",
      "length": 1500,
      "quality": 12.5,
      "aligned_length": 1450,
      "mapping_quality": 60,
      "percent_identity": 95.2,
      "channel_id": 100,
      "start_time": "2023-01-01T12:00:00Z",
      "duration": 2.5
    }
  ],
  "summary": {
    "read_count": 10000,
    "length_stats": {
      "mean": 1520.5,
      "median": 1500.0,
      "min": 100.0,
      "max": 50000.0,
      "std_dev": 2500.0
    },
    "quality_stats": { ... },
    "channel_distribution": { ... }
  }
}
```

## Performance

nanoget-rs is designed for high performance with:

- **Parallel processing** for multiple files
- **Memory-efficient** streaming for large files
- **Compressed file support** (gzip, bzip2)
- **Progress reporting** for long-running operations

## Comparison with Python nanoget

nanoget-rs aims to be functionally equivalent to the original Python nanoget while offering:

- **Better performance** through Rust's efficiency
- **Lower memory usage** with streaming and optimized data structures
- **Static typing** for improved reliability
- **Cross-platform** single binary distribution

### Storage

Metrics are held columnar — one array per field rather than one struct per read. A column
that an input format never populates is never allocated, which is most of them for any
given format: a plain FASTQ sets 3 of the 13 fields. Peak memory is 4-9x lower than a row
layout at the same run time.

```rust
let metrics = extract_auto(vec!["reads.fastq.gz"])?;

// Borrow a column outright - no copy, no projection.
let lengths: &[u32] = metrics.reads.lengths();

// Or read row-wise, where that is clearer.
for read in metrics.iter() {
    println!("{:?}\t{}", read.read_id(), read.length());
}

// Filtering selects indices and gathers once.
let long_reads = metrics.filter_by_length(10_000);
```

Float metrics are stored — and serialised — as `f32`, which is ample for a Phred score or
a percent identity. JSON therefore carries about seven significant digits
(`"quality": 7.69542`); against a `f64` pipeline, values can differ in the seventh. Summary
statistics are computed in `f64`, though their `min` and `max` are stored values and so
carry `f32` precision.

### Known differences in output

- **`percent_identity` is gap-compressed, not BLAST-style.** nanoget-rs reports
  `1 - (NM - gap_bases + gap_count) / (matches + gap_count)`, counting each indel once
  regardless of its length — the same definition as the minimap2 `de` tag, which is used
  directly when present. Python nanoget reports BLAST-style identity,
  `1 - NM / (matches + insertions + deletions)`, charging every base of every indel.
  Gap-compressed identity is the more meaningful measure for nanopore reads, whose errors
  are dominated by homopolymer indels, but **the two are not directly comparable**: on the
  nanotest alignment (1115 reads) nanoget-rs reports a mean of 89.96 and a minimum of
  73.84 where python nanoget reports 86.35 and 49.09.
- **`quality` is not reported for aligned BAM/CRAM input.** For aligned reads the percent
  identity above is the more informative measure, so the Phred quality python nanoget also
  extracts is deliberately not computed. FASTQ, uBAM and summary input all report quality
  as usual.
- **`mapping_quality` is `null` when the BAM records 255** (the "unavailable" sentinel);
  python nanoget keeps the literal 255 in the distribution.

## Development

The test suite runs against real ONT files from
[nanotest](https://github.com/wdecoster/nanotest), included as a git submodule:

```bash
git clone --recurse-submodules https://github.com/wdecoster/nanoget-rs.git
# or, in an existing clone:
git submodule update --init --depth 1
cargo test
```

## Contributing

Contributions are welcome! Please feel free to submit issues and enhancement requests.

### Development setup

Install the toolchain components and the git hooks once after cloning:

```bash
make setup          # adds rustfmt & clippy components (and cargo-audit/outdated)
make install-hooks  # installs pre-commit and pre-push hooks into .git/hooks
```

The hooks run `cargo fmt` and `cargo clippy` (and tests on commit) so that
formatting/lint issues are caught locally before they fail CI. Run the same
checks manually with `make pre-commit`, `make pre-push`, or `make ci`.

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This project is a Rust port of the original [nanoget](https://github.com/wdecoster/nanoget) by Wouter De Coster.