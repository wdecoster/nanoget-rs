# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-09-09

The first published release. 0.1.0 and 0.1.1 existed in `Cargo.toml` but were never
tagged, so their entries below are recorded for history rather than as shipped versions.

**If you have used this crate before, read the first entry under Fixed.** Every quality
score it reported from FASTQ input was wrong.

### Fixed

- FASTQ quality scores are now Phred+33 decoded. `bio`'s FASTQ reader returns the raw
  ASCII quality line, which was being averaged as if it were already Phred, so every
  FASTQ-derived quality was 33 too high and saturated at the 60 cap. BAM/uBAM/summary
  input was unaffected.
- Unaligned BAM piped on stdin is now detected from the SAM header instead of being
  routed to the aligned extractor, which discarded every record via the unmapped filter
  and reported zero reads with a success exit code.
- Extraction from stdin now errors on zero reads, matching the file path.
- `FileType::sniff` now classifies compressed files by decompressing their head instead
  of guessing from the extension. Gzipped rich FASTQ was reported as plain FASTQ, losing
  channel, start time and run id for the common `.fastq.gz` form; bzip2 files and
  gzipped sequencing summaries were rejected outright. The file and stdin sniffers now
  share one classifier, so the same bytes classify the same way either way.
- `MetricsCollection::reads_above_length_percentile` no longer underflows on an empty
  collection, and clamps out-of-range percentiles instead of indexing past the data.
- Malformed FASTQ is now a hard error naming the record, the read id and both lengths,
  instead of being silently accepted and averaged over the shorter of the sequence and
  quality lines.
- Zero-length reads are dropped from FASTQ, FASTA, BAM, uBAM and summary input, matching
  python nanoget; the number dropped is logged at `info` level.
- A blank `barcode_arrangement` cell in a summary file is treated as absent rather than as
  a literal empty barcode, consistent with the other optional columns.
- `channel_id` of 0 is preserved. It had been used as the absent marker, so a summary file
  carrying a literal 0 silently lost the channel and dropped the read from the channel
  distribution.
- The error for a rich FASTQ record with no metadata now names `--file-type fastq` as the
  way to read a file whose headers are mixed.

### Changed

- **Breaking (library): metrics are stored columnar.** `MetricsCollection::reads` is a
  `ReadColumns` rather than a `Vec<ReadMetrics>`: one array per field, with a column only
  allocated when an input actually populates it. `ReadMetrics` remains the row type, built
  transiently while parsing and never stored; `ReadView` borrows a read back out and
  offers the same fields as accessor methods. Peak memory falls 4-9x — a 2M-read rich
  FASTQ goes from 554 MB to 126 MB, a 4M-row summary from 865 MB to 201 MB — at no
  measurable cost in run time. The JSON and TSV output formats are unchanged.
  - Consumers move from `for r in &c.reads { r.length }` to `for r in c.iter() { r.length() }`,
    or better, borrow the column: `c.reads.lengths()`.
  - Filtering is `c.select(&indices)`; `filter_by_length`, `filter_by_quality` and
    `reads_above_length_percentile` are unchanged in behaviour.
  - `MetricsCollection::from_rows` builds one from a `Vec<ReadMetrics>` for tests and
    callers that already hold rows.
- **Breaking (library):** `MetricsCollection::summary` is now the method `summary()`,
  computed on first call and cached. It is derived from the reads, so exposing it as a
  field invited the two drifting apart; computing it lazily also stops `select` (and so
  every filter and downsample) paying for a full statistics pass a caller may never read.
- Float metrics (quality, aligned quality, percent identity, duration) are stored as
  `f32` and **serialised as `f32`**, so JSON now carries the precision that actually
  exists — `"quality": 7.69542` rather than `7.695419788360596`, which was the shortest
  text round-tripping the widened `f64` and looked exact to 16 digits when only about 7
  were real. TSV output, which is formatted to 2-3 decimals, is unaffected. Against the
  previous `f64` pipeline, 3 of 960 reads in the nanotest alignment shift percent identity
  by 0.001 percentage points and read qualities move by under 1e-6 Phred; summary
  statistics are still computed in `f64` and agree to nine significant figures.
- Output is now streamed rather than built in memory: `MetricsCollection::write_tsv`
  writes a row at a time, and the binary streams TSV and JSON into a `BufWriter`. Peak
  memory is 38-61% lower and no longer depends on the output format; a 2M-read FASTQ went
  from 1311 MB to 554 MB and from 2.94 s to 2.69 s.
- **Breaking (library):** `MetricsCollection::combine` takes the datasets' reads instead
  of whole collections, so summary statistics are computed once over the combined set
  instead of once per input file and then discarded.
- Sequencing summary columns are resolved once from the header rather than through a
  hash map built per row: 35% faster on a 4M-row summary file (4.84 s to 3.16 s). A
  missing required column is now reported from the header instead of the first data row.
- TSV output no longer ends with a trailing blank line.
- **Breaking (library):** `ExtractArgs::file_type` is now `Option<FileType>`, and
  `output_format`, `read_type` and `combine` are the typed enums `OutputFormat`,
  `ReadType` and `CombineMethod` instead of `String`. `MetricsCollection::combine` takes
  a `CombineMethod`. Added `convenience::extract_auto` for the common case.
- `--file-type` is now optional. Each input file's format is detected from its content,
  so a mixed set of files works in one run; pass `-t` only to override detection or to
  select `fastq-minimal`, which is a processing mode rather than a format.
- `--output-format`, `--read-type` and `--combine` are value enums, so an invalid value
  is rejected with the allowed set instead of silently falling back to a default (or, for
  `--output-format`, dumping a Rust `Debug` representation).
- Replaced the non-functional `--keep-supplementary` flag with `--drop-supplementary`.
  Supplementary alignments are kept by default, as before.
- A record in a rich FASTQ with no MinKNOW/albacore metadata in its header is now an
  error rather than a silently metadata-less read.
- Sequencing-summary columns that are present but unparseable are now an error naming the
  column and value; absent columns and blank cells remain silently optional. `start_time`
  now accepts RFC3339 timestamps as well as seconds, matching the rich-FASTQ path.
- Unaligned BAM now uses the BGZF decompression thread pool, like aligned BAM/CRAM; all
  alignment paths share one reader setup.
- Non-finite numbers (`nan`, `inf`) in summary columns are treated as "not measured"
  rather than aborting the run, and are excluded from summary statistics. `StatsSummary`
  fields are now finite by construction, which makes the JSON output round-trip through
  the crate's own `Deserialize` — previously a NaN serialised as `null` and could not be
  read back.
- `channel_distribution` and `barcode_distribution` are `BTreeMap` rather than `HashMap`,
  so JSON output is byte-identical between runs and channels are numerically sorted.
  Previously `HashMap`'s randomised iteration order made output undiffable.
- The release workflow is rebuilt. It had never run successfully: it declared no
  `permissions`, so the token was read-only; it built without `--target` but looked for
  the binary under `target/<target>/`, so four of its five matrix entries could not have
  found their output; it targeted Windows, which htslib does not support; and it used
  `actions/create-release` and `actions/upload-release-asset`, both archived since 2021.
  It now builds Linux gnu, static Linux musl and macOS arm64, and takes its release notes
  from this file.
- The binary now links the library instead of re-declaring every module as a private tree,
  which removed all 14 `#[allow(dead_code)]` attributes and halves what it compiles.
- Test fixtures from [nanotest](https://github.com/wdecoster/nanotest) are included as a
  git submodule, so tests needing real ONT data run everywhere instead of skipping.
- Documented the deliberate output differences from python nanoget in the README:
  gap-compressed rather than BLAST-style `percent_identity`, no `quality` for aligned
  BAM/CRAM, and `null` rather than 255 for an unavailable `mapping_quality`.

## [0.1.1] and [0.1.0] - never released

Recorded for history. Neither version was tagged or published.

### Added

- Initial Rust implementation of nanoget
- Support for FASTQ, FASTA, BAM, CRAM, uBAM, and summary files
- Parallel processing capabilities
- Memory-optimized streaming processing
- Comprehensive test suite
- Both CLI binary and library API
- Feature parity with Python nanoget
- GitHub Actions CI/CD pipeline
- Automated releases on tag push
- Documentation and examples

### Changed

- Complete rewrite from Python to Rust for better performance
- Enhanced error handling and type safety
- Improved memory efficiency

### Fixed

- All compilation warnings resolved
- Proper error propagation throughout codebase
