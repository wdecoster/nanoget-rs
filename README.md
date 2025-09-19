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

Extract metrics from a FASTQ file:
```bash
nanoget extract -t fastq reads.fastq
```

Extract metrics from a BAM file:
```bash
nanoget extract -t bam alignments.bam
```

Extract metrics from a sequencing summary:
```bash
nanoget extract -t summary sequencing_summary.txt
```

### Output formats

By default, output is in JSON format. You can also specify CSV:
```bash
nanoget extract -t fastq reads.fastq -f csv
```

Save output to a file:
```bash
nanoget extract -t fastq reads.fastq -o metrics.json
```

### Processing multiple files

Process multiple files and combine results:
```bash
nanoget extract -t fastq file1.fastq file2.fastq file3.fastq
```

Track datasets separately:
```bash
nanoget extract -t fastq file1.fastq file2.fastq --combine track --names sample1 sample2
```

### Advanced options

Use multiple threads:
```bash
nanoget extract -t fastq reads.fastq -j 8
```

Process huge files without parallelization:
```bash
nanoget extract -t fastq huge_file.fastq --huge
```

For BAM files, keep supplementary alignments:
```bash
nanoget extract -t bam alignments.bam --keep-supplementary
```

For summary files, specify read type and barcode analysis:
```bash
nanoget extract -t summary sequencing_summary.txt --read-type 1D --barcoded
```

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

## Contributing

Contributions are welcome! Please feel free to submit issues and enhancement requests.

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This project is a Rust port of the original [nanoget](https://github.com/wdecoster/nanoget) by Wouter De Coster.