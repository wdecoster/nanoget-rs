use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use crate::error::NanogetError;

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

    /// Percent identity to reference (for aligned reads)
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

/// Collection of read metrics with summary statistics
#[derive(Debug, Serialize, Deserialize)]
pub struct MetricsCollection {
    /// Individual read metrics
    pub reads: Vec<ReadMetrics>,

    /// Summary statistics
    pub summary: MetricsSummary,
}

impl MetricsCollection {
    /// Create a new collection from a vector of read metrics
    pub fn new(reads: Vec<ReadMetrics>) -> Self {
        let summary = MetricsSummary::from_reads(&reads);
        Self { reads, summary }
    }

    /// Combine multiple collections
    pub fn combine(collections: Vec<Self>, method: &str, names: Option<Vec<String>>) -> Self {
        let mut all_reads = Vec::new();

        match method {
            "track" => {
                // Add dataset names to reads
                for (i, mut collection) in collections.into_iter().enumerate() {
                    let dataset_name = names
                        .as_ref()
                        .and_then(|n| n.get(i))
                        .cloned()
                        .unwrap_or_else(|| format!("dataset_{}", i));

                    for read in &mut collection.reads {
                        read.dataset = Some(dataset_name.clone());
                    }
                    all_reads.extend(collection.reads);
                }
            }
            _ => {
                // Simple concatenation
                for collection in collections {
                    all_reads.extend(collection.reads);
                }
            }
        }

        Self::new(all_reads)
    }

    /// Get reads from a specific dataset (when using track mode)
    #[allow(dead_code)]
    pub fn reads_for_dataset(&self, dataset_name: &str) -> Vec<&ReadMetrics> {
        self.reads
            .iter()
            .filter(|read| read.dataset.as_deref() == Some(dataset_name))
            .collect()
    }

    /// Get all unique dataset names
    #[allow(dead_code)]
    pub fn dataset_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self
            .reads
            .iter()
            .filter_map(|read| read.dataset.clone())
            .collect();
        names.sort();
        names.dedup();
        names
    }

    /// Filter reads by minimum length
    #[allow(dead_code)]
    pub fn filter_by_length(&self, min_length: u32) -> MetricsCollection {
        let filtered_reads: Vec<ReadMetrics> = self
            .reads
            .iter()
            .filter(|read| read.length >= min_length)
            .cloned()
            .collect();
        MetricsCollection::new(filtered_reads)
    }

    /// Filter reads by minimum quality
    #[allow(dead_code)]
    pub fn filter_by_quality(&self, min_quality: f64) -> MetricsCollection {
        let filtered_reads: Vec<ReadMetrics> = self
            .reads
            .iter()
            .filter(|read| read.quality.map(|q| q >= min_quality).unwrap_or(false))
            .cloned()
            .collect();
        MetricsCollection::new(filtered_reads)
    }

    /// Get reads longer than a percentile threshold
    #[allow(dead_code)]
    pub fn reads_above_length_percentile(&self, percentile: f64) -> MetricsCollection {
        let mut lengths: Vec<u32> = self.reads.iter().map(|r| r.length).collect();
        lengths.sort();

        let index = (percentile / 100.0 * (lengths.len() - 1) as f64) as usize;
        let threshold = lengths.get(index).copied().unwrap_or(0);

        self.filter_by_length(threshold)
    }

    /// Export to JSON string
    /// Export to pretty-printed JSON string
    #[allow(dead_code)]
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Export to compact JSON string
    #[allow(dead_code)]
    pub fn to_json_compact(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string(self)
    }

    /// Export to TSV format
    pub fn to_tsv(&self) -> Result<String, NanogetError> {
        let mut output = String::new();
        
        // Header row for individual reads
        output.push_str("read_id\tlength\tquality\taligned_length\taligned_quality\tmapping_quality\tpercent_identity\tchannel_id\tstart_time\tduration\tbarcode\trun_id\tdataset\n");
        
        // Individual read data
        for read in &self.reads {
            output.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                read.read_id.as_deref().unwrap_or(""),
                read.length,
                read.quality.map(|q| format!("{:.3}", q)).unwrap_or_default(),
                read.aligned_length.map(|l| l.to_string()).unwrap_or_default(),
                read.aligned_quality.map(|q| format!("{:.3}", q)).unwrap_or_default(),
                read.mapping_quality.map(|q| q.to_string()).unwrap_or_default(),
                read.percent_identity.map(|p| format!("{:.3}", p)).unwrap_or_default(),
                read.channel_id.map(|c| c.to_string()).unwrap_or_default(),
                read.start_time.map(|t| t.to_rfc3339()).unwrap_or_default(),
                read.duration.map(|d| format!("{:.3}", d)).unwrap_or_default(),
                read.barcode.as_deref().unwrap_or(""),
                read.run_id.as_deref().unwrap_or(""),
                read.dataset.as_deref().unwrap_or("")
            ));
        }
        
        // Add summary statistics as a comment section
        output.push_str("\n# Summary Statistics\n");
        output.push_str(&format!("# Total reads: {}\n", self.summary.read_count));
        
        // Length statistics
        output.push_str(&format!(
            "# Length stats - count: {}, mean: {:.2}, median: {:.2}, min: {:.2}, max: {:.2}, std_dev: {:.2}, q25: {:.2}, q75: {:.2}\n",
            self.summary.length_stats.count,
            self.summary.length_stats.mean,
            self.summary.length_stats.median,
            self.summary.length_stats.min,
            self.summary.length_stats.max,
            self.summary.length_stats.std_dev,
            self.summary.length_stats.q25,
            self.summary.length_stats.q75
        ));
        
        // Quality statistics if available
        if let Some(quality_stats) = &self.summary.quality_stats {
            output.push_str(&format!(
                "# Quality stats - count: {}, mean: {:.2}, median: {:.2}, min: {:.2}, max: {:.2}, std_dev: {:.2}, q25: {:.2}, q75: {:.2}\n",
                quality_stats.count,
                quality_stats.mean,
                quality_stats.median,
                quality_stats.min,
                quality_stats.max,
                quality_stats.std_dev,
                quality_stats.q25,
                quality_stats.q75
            ));
        }
        
        // Mapping quality statistics if available
        if let Some(mapping_quality_stats) = &self.summary.mapping_quality_stats {
            output.push_str(&format!(
                "# Mapping quality stats - count: {}, mean: {:.2}, median: {:.2}, min: {:.2}, max: {:.2}, std_dev: {:.2}, q25: {:.2}, q75: {:.2}\n",
                mapping_quality_stats.count,
                mapping_quality_stats.mean,
                mapping_quality_stats.median,
                mapping_quality_stats.min,
                mapping_quality_stats.max,
                mapping_quality_stats.std_dev,
                mapping_quality_stats.q25,
                mapping_quality_stats.q75
            ));
        }
        
        // Percent identity statistics if available
        if let Some(percent_identity_stats) = &self.summary.percent_identity_stats {
            output.push_str(&format!(
                "# Percent identity stats - count: {}, mean: {:.2}, median: {:.2}, min: {:.2}, max: {:.2}, std_dev: {:.2}, q25: {:.2}, q75: {:.2}\n",
                percent_identity_stats.count,
                percent_identity_stats.mean,
                percent_identity_stats.median,
                percent_identity_stats.min,
                percent_identity_stats.max,
                percent_identity_stats.std_dev,
                percent_identity_stats.q25,
                percent_identity_stats.q75
            ));
        }
        
        Ok(output)
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

    /// Channel distribution (if available)
    pub channel_distribution: Option<HashMap<u16, usize>>,

    /// Barcode distribution (if available)
    pub barcode_distribution: Option<HashMap<String, usize>>,
}

impl MetricsSummary {
    /// Calculate summary statistics from a collection of reads
    pub fn from_reads(reads: &[ReadMetrics]) -> Self {
        let read_count = reads.len();

        // Length statistics
        let lengths: Vec<f64> = reads.iter().map(|r| r.length as f64).collect();
        let length_stats = StatsSummary::from_values(&lengths);

        // Quality statistics
        let qualities: Vec<f64> = reads.iter().filter_map(|r| r.quality).collect();
        let quality_stats = if !qualities.is_empty() {
            Some(StatsSummary::from_values(&qualities))
        } else {
            None
        };

        // Mapping quality statistics
        let mapping_qualities: Vec<f64> = reads
            .iter()
            .filter_map(|r| r.mapping_quality.map(|q| q as f64))
            .collect();
        let mapping_quality_stats = if !mapping_qualities.is_empty() {
            Some(StatsSummary::from_values(&mapping_qualities))
        } else {
            None
        };

        // Percent identity statistics
        let percent_identities: Vec<f64> =
            reads.iter().filter_map(|r| r.percent_identity).collect();
        let percent_identity_stats = if !percent_identities.is_empty() {
            Some(StatsSummary::from_values(&percent_identities))
        } else {
            None
        };

        // Channel distribution
        let mut channel_counts = HashMap::new();
        for read in reads {
            if let Some(channel) = read.channel_id {
                *channel_counts.entry(channel).or_insert(0) += 1;
            }
        }
        let channel_distribution = if !channel_counts.is_empty() {
            Some(channel_counts)
        } else {
            None
        };

        // Barcode distribution
        let mut barcode_counts = HashMap::new();
        for read in reads {
            if let Some(barcode) = &read.barcode {
                *barcode_counts.entry(barcode.clone()).or_insert(0) += 1;
            }
        }
        let barcode_distribution = if !barcode_counts.is_empty() {
            Some(barcode_counts)
        } else {
            None
        };

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

/// Basic statistical summary for numerical data
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
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

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
    fn test_tsv_output() {
        let read1 = ReadMetrics::new(Some("read1".to_string()), 1000)
            .with_quality(35.5);
        let read2 = ReadMetrics::new(Some("read2".to_string()), 2000)
            .with_quality(40.0)
            .with_alignment(1900, Some(41.0), Some(60), Some(95.5));
        
        let metrics = MetricsCollection::new(vec![read1, read2]);
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
