use crate::cli::CombineMethod;
use crate::columns::{self, ReadColumns, ReadColumnsBuilder, ReadView};
use crate::error::NanogetError;
use chrono::{DateTime, Utc};
use serde::ser::{SerializeSeq, SerializeStruct};
use serde::{Deserialize, Serialize, Serializer};
use std::collections::BTreeMap;
use std::fmt;
use std::io::Write;
use std::sync::OnceLock;

/// Represents the metrics extracted from a single read
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadMetrics {
    /// Read identifier
    pub read_id: Option<String>,

    /// Read length (number of bases)
    pub length: u32,

    /// Average quality score of the read
    pub quality: Option<f64>,

    /// Length of aligned portion (for aligned reads)
    pub aligned_length: Option<u32>,

    /// Average quality of aligned portion
    pub aligned_quality: Option<f64>,

    /// Mapping quality (for aligned reads)
    pub mapping_quality: Option<u8>,

    /// Gap-compressed percent identity to the reference (for aligned reads).
    ///
    /// Each indel counts once regardless of length, matching the minimap2 `de` tag.
    /// Note this is **not** the BLAST-style identity python nanoget reports, so values
    /// are systematically a few points higher and not directly comparable — see
    /// `extract::alignment_stats`.
    pub percent_identity: Option<f64>,

    /// Channel ID (from sequencing summary or rich FASTQ)
    pub channel_id: Option<u16>,

    /// Start time of sequencing
    pub start_time: Option<DateTime<Utc>>,

    /// Duration of sequencing
    pub duration: Option<f64>,

    /// Barcode assignment (for barcoded samples)
    pub barcode: Option<String>,

    /// Run ID
    pub run_id: Option<String>,

    /// Dataset name (when combining multiple files with tracking)
    pub dataset: Option<String>,
}

impl ReadMetrics {
    /// Create a new ReadMetrics with basic information
    pub fn new(read_id: Option<String>, length: u32) -> Self {
        Self {
            read_id,
            length,
            quality: None,
            aligned_length: None,
            aligned_quality: None,
            mapping_quality: None,
            percent_identity: None,
            channel_id: None,
            start_time: None,
            duration: None,
            barcode: None,
            run_id: None,
            dataset: None,
        }
    }

    /// Set quality score
    pub fn with_quality(mut self, quality: f64) -> Self {
        self.quality = Some(quality);
        self
    }

    /// Set alignment information
    pub fn with_alignment(
        mut self,
        aligned_length: u32,
        aligned_quality: Option<f64>,
        mapping_quality: Option<u8>,
        percent_identity: Option<f64>,
    ) -> Self {
        self.aligned_length = Some(aligned_length);
        self.aligned_quality = aligned_quality;
        self.mapping_quality = mapping_quality;
        self.percent_identity = percent_identity;
        self
    }

    /// Set sequencing metadata
    pub fn with_sequencing_metadata(
        mut self,
        channel_id: Option<u16>,
        start_time: Option<DateTime<Utc>>,
        duration: Option<f64>,
    ) -> Self {
        self.channel_id = channel_id;
        self.start_time = start_time;
        self.duration = duration;
        self
    }
}

/// Collection of read metrics with summary statistics.
///
/// Reads are held columnar (see [`ReadColumns`]); `iter()` and `reads.get()` borrow them
/// back as rows. The serialised *shape* is unchanged — JSON still carries a `reads` array
/// of per-read objects — though float fields now render at `f32` precision, which is what
/// the columns hold.
#[derive(Debug)]
pub struct MetricsCollection {
    /// Individual read metrics, one column per field.
    pub reads: ReadColumns,

    /// Summary statistics, computed on first use.
    ///
    /// Not a public field: it is derived from `reads`, so exposing it invited the two
    /// drifting apart. Computing it lazily also matters — filtering and downsampling go
    /// through `select`, and a consumer that only ever reads columns was paying for a
    /// full statistics pass on every intermediate collection.
    summary: OnceLock<MetricsSummary>,
}

impl MetricsCollection {
    /// Create a new collection from columns.
    pub fn new(reads: ReadColumns) -> Self {
        Self {
            reads,
            summary: OnceLock::new(),
        }
    }

    /// Summary statistics over the whole collection, computed on first call and cached.
    pub fn summary(&self) -> &MetricsSummary {
        self.summary
            .get_or_init(|| MetricsSummary::from_columns(&self.reads))
    }

    /// Create a collection from owned rows.
    ///
    /// Convenient for tests and for callers holding a `Vec<ReadMetrics>`; the extraction
    /// path builds columns directly and never materialises the rows.
    pub fn from_rows(rows: Vec<ReadMetrics>) -> Self {
        let mut builder = ReadColumnsBuilder::with_capacity(rows.len());
        for row in rows {
            builder.push(row);
        }
        Self::new(builder.finish())
    }

    /// Number of reads.
    pub fn len(&self) -> usize {
        self.reads.len()
    }

    pub fn is_empty(&self) -> bool {
        self.reads.is_empty()
    }

    /// Iterate over reads as borrowed rows.
    pub fn iter(&self) -> impl Iterator<Item = ReadView<'_>> + '_ {
        self.reads.iter()
    }

    /// Combine the reads of several datasets into one collection.
    ///
    /// Takes the columns rather than whole `MetricsCollection`s so that the summary is
    /// computed exactly once, over the combined set.
    pub fn combine(
        datasets: Vec<ReadColumns>,
        method: CombineMethod,
        names: Option<Vec<String>>,
    ) -> Self {
        let mut datasets = datasets;

        if matches!(method, CombineMethod::Track) {
            for (i, columns) in datasets.iter_mut().enumerate() {
                let name = names
                    .as_ref()
                    .and_then(|n| n.get(i))
                    .cloned()
                    .unwrap_or_else(|| format!("dataset_{}", i));
                columns::set_dataset(columns, &name);
            }
        }

        Self::new(columns::concat(datasets))
    }

    /// Get reads from a specific dataset (when using track mode)
    pub fn reads_for_dataset(&self, dataset_name: &str) -> Vec<ReadView<'_>> {
        self.reads
            .iter()
            .filter(|read| read.dataset() == Some(dataset_name))
            .collect()
    }

    /// Get all unique dataset names
    pub fn dataset_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self
            .reads
            .datasets()
            .map(|d| d.values().to_vec())
            .unwrap_or_default();
        names.sort();
        names.dedup();
        names
    }

    /// Select reads by index, in the order given.
    ///
    /// The building block for filtering, downsampling and reordering: decide on indices,
    /// then gather the columns once.
    pub fn select(&self, indices: &[usize]) -> MetricsCollection {
        MetricsCollection::new(self.reads.select(indices))
    }

    /// Filter reads by minimum length
    pub fn filter_by_length(&self, min_length: u32) -> MetricsCollection {
        let indices: Vec<usize> = self
            .reads
            .lengths()
            .iter()
            .enumerate()
            .filter(|(_, &l)| l >= min_length)
            .map(|(i, _)| i)
            .collect();
        self.select(&indices)
    }

    /// Filter reads by minimum quality
    pub fn filter_by_quality(&self, min_quality: f64) -> MetricsCollection {
        let indices: Vec<usize> = self
            .reads
            .iter()
            .filter(|read| read.quality().is_some_and(|q| q >= min_quality))
            .map(|read| read.index())
            .collect();
        self.select(&indices)
    }

    /// Get reads longer than a percentile threshold
    pub fn reads_above_length_percentile(&self, percentile: f64) -> MetricsCollection {
        if self.reads.is_empty() {
            return MetricsCollection::new(ReadColumns::default());
        }

        let mut lengths: Vec<u32> = self.reads.lengths().to_vec();
        lengths.sort_unstable();

        // `lengths` is non-empty, so `len() - 1` cannot underflow. The percentile is
        // clamped so an out-of-range argument saturates at the ends of the distribution
        // instead of indexing past it (a NaN percentile casts to index 0).
        let last = lengths.len() - 1;
        let index = (percentile.clamp(0.0, 100.0) / 100.0 * last as f64) as usize;
        let threshold = lengths.get(index).copied().unwrap_or(0);

        self.filter_by_length(threshold)
    }

    /// Export to pretty-printed JSON string
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Export to compact JSON string
    pub fn to_json_compact(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string(self)
    }

    /// Export to TSV format.
    ///
    /// Convenience wrapper over [`write_tsv`](Self::write_tsv) that buffers the whole
    /// output in memory. Prefer `write_tsv` for anything large: this string is roughly
    /// 190 bytes per read, which for a full sequencing run exceeds the size of the reads
    /// themselves.
    pub fn to_tsv(&self) -> Result<String, NanogetError> {
        let mut buf = Vec::new();
        self.write_tsv(&mut buf)?;
        String::from_utf8(buf)
            .map_err(|e| NanogetError::ProcessingError(format!("TSV is not valid UTF-8: {}", e)))
    }

    /// Write TSV directly to a sink, one read at a time.
    ///
    /// Keeps memory independent of the read count: nothing larger than a single row is
    /// held at once. Pass a [`BufWriter`](std::io::BufWriter) — this issues many small
    /// writes.
    pub fn write_tsv<W: Write>(&self, mut writer: W) -> Result<(), NanogetError> {
        // Header row for individual reads
        writer.write_all(b"read_id\tlength\tquality\taligned_length\taligned_quality\tmapping_quality\tpercent_identity\tchannel_id\tstart_time\tduration\tbarcode\trun_id\tdataset\n")?;

        // Individual read data
        for read in self.reads.iter() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                read.read_id().unwrap_or(""),
                read.length(),
                OptFixed(read.quality()),
                OptDisplay(read.aligned_length()),
                OptFixed(read.aligned_quality()),
                OptDisplay(read.mapping_quality()),
                OptFixed(read.percent_identity()),
                OptDisplay(read.channel_id()),
                OptRfc3339(read.start_time()),
                OptFixed(read.duration()),
                read.barcode().unwrap_or(""),
                read.run_id().unwrap_or(""),
                read.dataset().unwrap_or("")
            )?;
        }

        // Add summary statistics as a comment section
        writer.write_all(b"\n# Summary Statistics\n")?;
        let summary = self.summary();
        writeln!(writer, "# Total reads: {}", summary.read_count)?;

        write_stats_comment(&mut writer, "Length", &summary.length_stats)?;
        for (label, stats) in [
            ("Quality", &summary.quality_stats),
            ("Mapping quality", &summary.mapping_quality_stats),
            ("Percent identity", &summary.percent_identity_stats),
        ] {
            if let Some(stats) = stats {
                write_stats_comment(&mut writer, label, stats)?;
            }
        }

        Ok(())
    }
}

/// Serialised as `{ "reads": [ {...}, ... ], "summary": {...} }` — the same shape the row
/// layout produced. Rows are materialised one at a time during serialisation and dropped,
/// so the columnar saving is not given back here.
impl Serialize for MetricsCollection {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("MetricsCollection", 2)?;
        state.serialize_field("reads", &RowsSerializer(&self.reads))?;
        state.serialize_field("summary", self.summary())?;
        state.end()
    }
}

struct RowsSerializer<'a>(&'a ReadColumns);

impl Serialize for RowsSerializer<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;
        for view in self.0.iter() {
            seq.serialize_element(&SerializedRead::from(&view))?;
        }
        seq.end()
    }
}

/// A read as it goes onto the wire.
///
/// Float fields are `f32` here because that is what the columns hold. Serialising the
/// widened `f64` instead printed the shortest text that round-trips as `f64` — 16-17
/// digits of which only about 7 carry information, so a quality stored as 7.6954198 came
/// out as 7.695419788360596 and looked exact to a consumer parsing it at full precision.
/// The field names and order match `ReadMetrics`, so the JSON shape is unchanged.
#[derive(Serialize)]
struct SerializedRead<'a> {
    read_id: Option<&'a str>,
    length: u32,
    quality: Option<f32>,
    aligned_length: Option<u32>,
    aligned_quality: Option<f32>,
    mapping_quality: Option<u8>,
    percent_identity: Option<f32>,
    channel_id: Option<u16>,
    start_time: Option<DateTime<Utc>>,
    duration: Option<f32>,
    barcode: Option<&'a str>,
    run_id: Option<&'a str>,
    dataset: Option<&'a str>,
}

impl<'a> From<&ReadView<'a>> for SerializedRead<'a> {
    fn from(view: &ReadView<'a>) -> Self {
        Self {
            read_id: view.read_id(),
            length: view.length(),
            quality: view.quality().map(|v| v as f32),
            aligned_length: view.aligned_length(),
            aligned_quality: view.aligned_quality().map(|v| v as f32),
            mapping_quality: view.mapping_quality(),
            percent_identity: view.percent_identity().map(|v| v as f32),
            channel_id: view.channel_id(),
            start_time: view.start_time(),
            duration: view.duration().map(|v| v as f32),
            barcode: view.barcode(),
            run_id: view.run_id(),
            dataset: view.dataset(),
        }
    }
}

impl<'de> Deserialize<'de> for MetricsCollection {
    /// Only the `reads` array is read; the summary is recomputed from it rather than
    /// trusted, so a document whose `summary` disagrees with its reads (or omits it
    /// entirely) still loads, and loads consistently.
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        // Reading is not the hot path, so the rows are materialised and then packed.
        #[derive(Deserialize)]
        struct Wire {
            reads: Vec<ReadMetrics>,
        }
        let wire = Wire::deserialize(deserializer)?;
        Ok(MetricsCollection::from_rows(wire.reads))
    }
}

/// Write one `# <label> stats - ...` comment line.
fn write_stats_comment<W: Write>(
    mut writer: W,
    label: &str,
    stats: &StatsSummary,
) -> Result<(), NanogetError> {
    writeln!(
        writer,
        "# {} stats - count: {}, mean: {:.2}, median: {:.2}, min: {:.2}, max: {:.2}, \
         std_dev: {:.2}, q25: {:.2}, q75: {:.2}",
        label,
        stats.count,
        stats.mean,
        stats.median,
        stats.min,
        stats.max,
        stats.std_dev,
        stats.q25,
        stats.q75
    )?;
    Ok(())
}

/// `Display` adapters that render an `Option` as either its value or an empty field,
/// writing straight into the output rather than building a `String` per cell.
struct OptFixed(Option<f64>);
struct OptDisplay<T: fmt::Display>(Option<T>);
struct OptRfc3339(Option<DateTime<Utc>>);

impl fmt::Display for OptFixed {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(v) => write!(f, "{:.3}", v),
            None => Ok(()),
        }
    }
}

impl<T: fmt::Display> fmt::Display for OptDisplay<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.0 {
            Some(v) => write!(f, "{}", v),
            None => Ok(()),
        }
    }
}

impl fmt::Display for OptRfc3339 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            // to_rfc3339 allocates; the Display impls chrono offers do not match the
            // previous output, so this keeps the format byte-identical.
            Some(t) => write!(f, "{}", t.to_rfc3339()),
            None => Ok(()),
        }
    }
}

/// Summary statistics for a collection of reads
#[derive(Debug, Serialize, Deserialize)]
pub struct MetricsSummary {
    /// Total number of reads
    pub read_count: usize,

    /// Length statistics
    pub length_stats: StatsSummary,

    /// Quality statistics (if available)
    pub quality_stats: Option<StatsSummary>,

    /// Mapping quality statistics (if available)
    pub mapping_quality_stats: Option<StatsSummary>,

    /// Percent identity statistics (if available)
    pub percent_identity_stats: Option<StatsSummary>,

    /// Channel distribution (if available).
    ///
    /// A `BTreeMap` rather than a `HashMap` so the serialised order is deterministic —
    /// `HashMap` iteration order is randomised per process, which made two runs over the
    /// same input produce byte-different JSON. Channels also come out numerically sorted.
    pub channel_distribution: Option<BTreeMap<u16, usize>>,

    /// Barcode distribution (if available). Ordered, for the same reason.
    pub barcode_distribution: Option<BTreeMap<String, usize>>,
}

impl MetricsSummary {
    /// Calculate summary statistics from columns.
    ///
    /// Each statistic reads its column directly. Where a column is absent, or holds no
    /// finite value, the statistic is `None` rather than a summary of nothing.
    pub fn from_columns(reads: &ReadColumns) -> Self {
        let read_count = reads.len();

        let lengths: Vec<f64> = reads.lengths().iter().map(|&l| l as f64).collect();
        let length_stats = StatsSummary::from_values(&lengths);

        // Non-finite values are skipped alongside absent ones: a NaN is not a
        // measurement, and letting one through would poison every field of the resulting
        // StatsSummary, which serialises it as JSON `null` — unreadable by the type's own
        // Deserialize, since these fields are `f64` and not `Option<f64>`.
        let quality_stats = stats_from_f32(reads.qualities_raw());
        let percent_identity_stats = stats_from_f32(reads.percent_identities_raw());
        let mapping_quality_stats = reads.mapping_qualities_raw().and_then(|column| {
            let values: Vec<f64> = column
                .iter()
                .filter(|&&q| q != 255)
                .map(|&q| q as f64)
                .collect();
            (!values.is_empty()).then(|| StatsSummary::from_values(&values))
        });

        let channel_distribution = reads.channel_ids_raw().and_then(|column| {
            let mut counts: BTreeMap<u16, usize> = BTreeMap::new();
            for &channel in column.iter().filter(|&&c| c != 0) {
                *counts.entry(channel).or_insert(0) += 1;
            }
            (!counts.is_empty()).then_some(counts)
        });

        let barcode_distribution = reads.barcodes().and_then(|_| {
            let mut counts: BTreeMap<String, usize> = BTreeMap::new();
            for read in reads.iter() {
                if let Some(barcode) = read.barcode() {
                    *counts.entry(barcode.to_string()).or_insert(0) += 1;
                }
            }
            (!counts.is_empty()).then_some(counts)
        });

        Self {
            read_count,
            length_stats,
            quality_stats,
            mapping_quality_stats,
            percent_identity_stats,
            channel_distribution,
            barcode_distribution,
        }
    }
}

/// Summarise a float column, skipping absent (`NaN`) entries.
fn stats_from_f32(column: Option<&[f32]>) -> Option<StatsSummary> {
    let column = column?;
    let values: Vec<f64> = column
        .iter()
        .filter(|v| v.is_finite())
        .map(|&v| f64::from(v))
        .collect();
    (!values.is_empty()).then(|| StatsSummary::from_values(&values))
}

/// Basic statistical summary for numerical data.
///
/// Every field is finite: non-finite inputs are filtered out before the summary is
/// built. This is what makes the serialised form round-trip — serde_json writes a
/// non-finite `f64` as `null`, which these non-optional fields cannot read back.
#[derive(Debug, Serialize, Deserialize)]
pub struct StatsSummary {
    pub count: usize,
    pub mean: f64,
    pub median: f64,
    pub min: f64,
    pub max: f64,
    pub std_dev: f64,
    pub q25: f64,
    pub q75: f64,
}

impl StatsSummary {
    /// Calculate statistics from a vector of values
    pub fn from_values(values: &[f64]) -> Self {
        if values.is_empty() {
            return Self {
                count: 0,
                mean: 0.0,
                median: 0.0,
                min: 0.0,
                max: 0.0,
                std_dev: 0.0,
                q25: 0.0,
                q75: 0.0,
            };
        }

        let mut sorted_values = values.to_vec();
        // Use unwrap_or(Equal) to handle NaN values gracefully
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let count = values.len();
        let mean = values.iter().sum::<f64>() / count as f64;
        let median = calculate_percentile(&sorted_values, 50.0);
        let min = sorted_values[0];
        let max = sorted_values[count - 1];
        let q25 = calculate_percentile(&sorted_values, 25.0);
        let q75 = calculate_percentile(&sorted_values, 75.0);

        // Calculate standard deviation
        let variance = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / count as f64;
        let std_dev = variance.sqrt();

        Self {
            count,
            mean,
            median,
            min,
            max,
            std_dev,
            q25,
            q75,
        }
    }
}

/// Calculate percentile from sorted values
fn calculate_percentile(sorted_values: &[f64], percentile: f64) -> f64 {
    if sorted_values.is_empty() {
        return 0.0;
    }

    let index = (percentile / 100.0) * (sorted_values.len() - 1) as f64;
    let lower = index.floor() as usize;
    let upper = index.ceil() as usize;

    if lower == upper {
        sorted_values[lower]
    } else {
        let weight = index - lower as f64;
        sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stats_summary() {
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let stats = StatsSummary::from_values(&values);

        assert_eq!(stats.count, 5);
        assert_eq!(stats.mean, 3.0);
        assert_eq!(stats.median, 3.0);
        assert_eq!(stats.min, 1.0);
        assert_eq!(stats.max, 5.0);
    }

    #[test]
    fn test_read_metrics_builder() {
        let metrics = ReadMetrics::new(Some("read1".to_string()), 1000)
            .with_quality(35.0)
            .with_alignment(950, Some(36.0), Some(60), Some(95.5));

        assert_eq!(metrics.length, 1000);
        assert_eq!(metrics.quality, Some(35.0));
        assert_eq!(metrics.aligned_length, Some(950));
        assert_eq!(metrics.percent_identity, Some(95.5));
    }

    #[test]
    fn test_reads_above_length_percentile_on_empty_collection() {
        // Previously underflowed `lengths.len() - 1` and panicked in debug builds.
        let empty = MetricsCollection::from_rows(Vec::new());
        assert_eq!(empty.reads_above_length_percentile(90.0).reads.len(), 0);
        assert_eq!(empty.reads_above_length_percentile(0.0).reads.len(), 0);
    }

    #[test]
    fn test_reads_above_length_percentile_thresholds() {
        let reads = (1..=5)
            .map(|i| ReadMetrics::new(Some(format!("r{}", i)), i * 100))
            .collect();
        let collection = MetricsCollection::from_rows(reads);

        // Lengths are 100..500; the 50th percentile threshold is 300, keeping 3 reads.
        assert_eq!(
            collection.reads_above_length_percentile(50.0).reads.len(),
            3
        );
        assert_eq!(collection.reads_above_length_percentile(0.0).reads.len(), 5);
        assert_eq!(
            collection.reads_above_length_percentile(100.0).reads.len(),
            1
        );

        // Out-of-range percentiles saturate instead of indexing past the distribution.
        assert_eq!(
            collection.reads_above_length_percentile(-10.0).reads.len(),
            5
        );
        assert_eq!(
            collection.reads_above_length_percentile(150.0).reads.len(),
            1
        );
        assert_eq!(
            collection
                .reads_above_length_percentile(f64::NAN)
                .reads
                .len(),
            5
        );
    }

    #[test]
    fn test_non_finite_values_are_excluded_from_summaries() {
        // ReadMetrics fields are public, so a caller can set a NaN directly. It must not
        // reach the summary, where it would make every field non-finite.
        let reads = vec![
            ReadMetrics::new(Some("ok".into()), 100).with_quality(20.0),
            ReadMetrics::new(Some("nan".into()), 100).with_quality(f64::NAN),
            ReadMetrics::new(Some("inf".into()), 100).with_quality(f64::INFINITY),
        ];
        let collection = MetricsCollection::from_rows(reads);

        let quality = collection
            .summary()
            .quality_stats
            .as_ref()
            .expect("quality stats");
        assert_eq!(
            quality.count, 1,
            "only the finite quality should be counted"
        );
        assert_eq!(quality.mean, 20.0);
        for value in [
            quality.mean,
            quality.median,
            quality.min,
            quality.max,
            quality.std_dev,
            quality.q25,
            quality.q75,
        ] {
            assert!(value.is_finite(), "non-finite field: {}", value);
        }
    }

    #[test]
    fn test_summary_is_absent_when_every_value_is_non_finite() {
        let reads = vec![ReadMetrics::new(Some("nan".into()), 100).with_quality(f64::NAN)];
        let collection = MetricsCollection::from_rows(reads);
        assert!(collection.summary().quality_stats.is_none());
    }

    #[test]
    fn test_json_output_round_trips() {
        // serde_json writes a non-finite f64 as `null`, which StatsSummary's non-optional
        // fields cannot read back — so this only holds because summaries are finite.
        let reads = vec![
            ReadMetrics::new(Some("r1".into()), 1000)
                .with_quality(f64::NAN)
                .with_alignment(950, None, Some(60), Some(95.5)),
            ReadMetrics::new(Some("r2".into()), 2000).with_quality(12.5),
        ];
        let collection = MetricsCollection::from_rows(reads);

        let json = collection.to_json().expect("serialise");
        let parsed: MetricsCollection = serde_json::from_str(&json).expect("round-trip");

        assert_eq!(parsed.summary().read_count, collection.summary().read_count);
        assert_eq!(parsed.reads.len(), 2);
        assert_eq!(parsed.summary().length_stats.mean, 1500.0);
        // The per-read NaN still serialises as null, which Option<f64> reads back as None.
        assert_eq!(parsed.reads.get(0).unwrap().quality(), None);
        assert_eq!(parsed.reads.get(1).unwrap().quality(), Some(12.5));
    }

    #[test]
    fn test_tsv_output() {
        let read1 = ReadMetrics::new(Some("read1".to_string()), 1000).with_quality(35.5);
        let read2 = ReadMetrics::new(Some("read2".to_string()), 2000)
            .with_quality(40.0)
            .with_alignment(1900, Some(41.0), Some(60), Some(95.5));

        let metrics = MetricsCollection::from_rows(vec![read1, read2]);
        let tsv_output = metrics.to_tsv().unwrap();

        // Check that it contains the header
        assert!(tsv_output.contains("read_id\tlength\tquality"));

        // Check that it contains the read data with tabs
        assert!(tsv_output.contains("read1\t1000\t35.500"));
        assert!(tsv_output.contains("read2\t2000\t40.000"));

        // Check that it contains summary statistics
        assert!(tsv_output.contains("# Summary Statistics"));
        assert!(tsv_output.contains("# Total reads: 2"));
        assert!(tsv_output.contains("# Length stats"));
        assert!(tsv_output.contains("# Quality stats"));
    }
}
