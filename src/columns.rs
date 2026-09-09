//! Columnar storage for per-read metrics.
//!
//! Metrics are held one column per field rather than one struct per read. Two things
//! follow, and they are the whole point:
//!
//! * **A column that an input format never populates is never allocated.** A plain FASTQ
//!   sets 3 of the 13 fields; a row layout still pays for all 13 on every read.
//! * **Consumers that want a column get a slice, not a copy.** Reading code almost always
//!   wants `&[u32]` of lengths, not a sequence of structs to project.
//!
//! [`ReadMetrics`](crate::ReadMetrics) remains the row type. It is built transiently, one
//! read at a time, and handed to [`ReadColumnsBuilder::push`], which scatters it into the
//! columns; it is never stored. Going the other way, [`ReadView`] borrows one read's
//! worth of columns and reads like the struct it replaces.
//!
//! ## Absence
//!
//! There are two levels of it, and they mean different things:
//!
//! * The whole column is `None` — this input has no such field at all (FASTA has no
//!   qualities).
//! * The column exists but a particular read has no value — encoded in-band with a
//!   sentinel rather than a validity bitmap, because every field already has a natural
//!   one. Floats use `NaN`, which is also how the parsers already treat an unmeasured
//!   value; `mapping_quality` uses 255, which is what SAM itself means by "unavailable".
//!
//! The consequence is that a sentinel value cannot be told apart from absence. For every
//! field this is either exactly the intended reading (`NaN` quality, mapq 255) or a value
//! the parsers cannot produce (`u32::MAX` aligned length, `u16::MAX` channel). It does
//! mean a `ReadMetrics` constructed by hand with one of those values, or JSON carrying
//! one, reads back as `None`.

use crate::metrics::ReadMetrics;
use chrono::{DateTime, TimeZone, Utc};

/// Absent marker for `aligned_length`.
const NO_U32: u32 = u32::MAX;
/// Absent marker for `mapping_quality` — the SAM "unavailable" value.
const NO_MAPQ: u8 = 255;
/// Absent marker for `channel_id`.
///
/// `u16::MAX` rather than 0: ONT channels are numbered from 1, but a summary file can
/// carry a literal 0, and treating that as "absent" silently dropped the read from the
/// channel distribution. No flow cell has 65535 channels.
const NO_CHANNEL: u16 = u16::MAX;
/// Absent marker for `start_time`.
const NO_TIME: i64 = i64::MIN;
/// Absent marker for a dictionary code.
const NO_CODE: u32 = u32::MAX;

/// Read identifiers packed into one buffer, addressed by offset.
///
/// One growing buffer for all ids instead of an allocation per read: a 36-character UUID
/// costs 36 bytes plus an 8-byte offset here, against 24 bytes of `String` plus a 36-byte
/// heap block plus allocator overhead in a row layout.
///
/// Offsets are `u64`. `u32` would be half the width but wraps silently once the packed
/// ids exceed 4 GiB — around 119M UUIDs, which a PromethION run can reach — handing out
/// garbage slices with no error anywhere.
#[derive(Debug, Default, Clone)]
pub struct StringColumn {
    data: String,
    /// `offsets[i]..offsets[i + 1]` bounds read `i`; always one longer than the column.
    offsets: Vec<u64>,
    /// Empty is a legitimate id, so absence needs its own flag.
    present: Vec<bool>,
}

impl StringColumn {
    /// Bytes reserved per read for the id arena. Read ids are typically a 36-character
    /// UUID; this covers most of that up front so the buffer's doubling growth starts from
    /// a useful size rather than from zero.
    const ID_BYTES_HINT: usize = 32;

    fn with_capacity(capacity: usize) -> Self {
        let mut offsets = Vec::with_capacity(capacity + 1);
        offsets.push(0);
        Self {
            data: String::with_capacity(capacity.saturating_mul(Self::ID_BYTES_HINT)),
            offsets,
            present: Vec::with_capacity(capacity),
        }
    }

    fn push(&mut self, value: Option<&str>) {
        match value {
            Some(v) => {
                self.data.push_str(v);
                self.present.push(true);
            }
            None => self.present.push(false),
        }
        self.offsets.push(self.data.len() as u64);
    }

    fn get(&self, index: usize) -> Option<&str> {
        if !*self.present.get(index)? {
            return None;
        }
        let (start, end) = (
            *self.offsets.get(index)? as usize,
            *self.offsets.get(index + 1)? as usize,
        );
        self.data.get(start..end)
    }
}

/// A low-cardinality string column stored as codes into a dictionary.
///
/// Run ids, barcodes and dataset names repeat across every read of a file — one run id
/// per run, one of ~96 barcodes, a handful of dataset names. Storing a code per read and
/// the value once collapses that to 4 bytes per read.
#[derive(Debug, Default, Clone)]
pub struct DictColumn {
    values: Vec<String>,
    codes: Vec<u32>,
}

impl DictColumn {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            values: Vec::new(),
            codes: Vec::with_capacity(capacity),
        }
    }

    fn push(&mut self, value: Option<&str>) {
        let code = match value {
            None => NO_CODE,
            Some(v) => {
                // Linear scan: these columns have a handful of distinct values, so a hash
                // map would cost more than it saves.
                match self.values.iter().position(|existing| existing == v) {
                    Some(i) => i as u32,
                    None => {
                        self.values.push(v.to_string());
                        (self.values.len() - 1) as u32
                    }
                }
            }
        };
        self.codes.push(code);
    }

    fn get(&self, index: usize) -> Option<&str> {
        let code = *self.codes.get(index)?;
        if code == NO_CODE {
            return None;
        }
        self.values.get(code as usize).map(String::as_str)
    }

    /// Overwrite every entry with a single value — how `combine` tags a dataset.
    ///
    /// A zero-length column gets no dictionary entry: otherwise a dataset with no reads
    /// would still be listed by `dataset_names()`.
    fn fill(&mut self, len: usize, value: &str) {
        self.values.clear();
        self.codes.clear();
        if len == 0 {
            return;
        }
        self.values.push(value.to_string());
        self.codes.resize(len, 0);
    }

    /// The distinct values, in first-seen order.
    pub fn values(&self) -> &[String] {
        &self.values
    }
}

/// Per-read metrics, stored one column per field.
#[derive(Debug, Default, Clone)]
pub struct ReadColumns {
    len: usize,
    length: Vec<u32>,
    quality: Option<Vec<f32>>,
    aligned_length: Option<Vec<u32>>,
    aligned_quality: Option<Vec<f32>>,
    mapping_quality: Option<Vec<u8>>,
    percent_identity: Option<Vec<f32>>,
    channel_id: Option<Vec<u16>>,
    start_time: Option<Vec<i64>>,
    duration: Option<Vec<f32>>,
    read_id: Option<StringColumn>,
    barcode: Option<DictColumn>,
    run_id: Option<DictColumn>,
    dataset: Option<DictColumn>,
}

impl ReadColumns {
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Read lengths. Always present, one entry per read.
    pub fn lengths(&self) -> &[u32] {
        &self.length
    }

    /// Raw quality column, `NaN` where a read has no quality.
    pub fn qualities_raw(&self) -> Option<&[f32]> {
        self.quality.as_deref()
    }

    /// Raw aligned-length column, [`u32::MAX`] where absent.
    pub fn aligned_lengths_raw(&self) -> Option<&[u32]> {
        self.aligned_length.as_deref()
    }

    /// Raw percent-identity column, `NaN` where absent.
    pub fn percent_identities_raw(&self) -> Option<&[f32]> {
        self.percent_identity.as_deref()
    }

    /// Raw mapping-quality column, 255 where absent.
    pub fn mapping_qualities_raw(&self) -> Option<&[u8]> {
        self.mapping_quality.as_deref()
    }

    /// Raw start-time column as nanoseconds since the Unix epoch, [`i64::MIN`] where
    /// absent.
    pub fn start_times_raw(&self) -> Option<&[i64]> {
        self.start_time.as_deref()
    }

    /// Raw channel column, 0 where absent.
    pub fn channel_ids_raw(&self) -> Option<&[u16]> {
        self.channel_id.as_deref()
    }

    /// Raw duration column, `NaN` where absent.
    pub fn durations_raw(&self) -> Option<&[f32]> {
        self.duration.as_deref()
    }

    pub fn barcodes(&self) -> Option<&DictColumn> {
        self.barcode.as_ref()
    }

    pub fn datasets(&self) -> Option<&DictColumn> {
        self.dataset.as_ref()
    }

    /// Whether any read carries a value for the field.
    pub fn has_quality(&self) -> bool {
        Self::any_finite(self.quality.as_deref())
    }

    pub fn has_alignment(&self) -> bool {
        self.aligned_length
            .as_deref()
            .is_some_and(|c| c.iter().any(|&v| v != NO_U32))
    }

    pub fn has_time(&self) -> bool {
        self.start_time
            .as_deref()
            .is_some_and(|c| c.iter().any(|&v| v != NO_TIME))
    }

    pub fn has_channel(&self) -> bool {
        self.channel_id
            .as_deref()
            .is_some_and(|c| c.iter().any(|&v| v != NO_CHANNEL))
    }

    fn any_finite(column: Option<&[f32]>) -> bool {
        column.is_some_and(|c| c.iter().any(|v| v.is_finite()))
    }

    /// Borrow one read.
    pub fn get(&self, index: usize) -> Option<ReadView<'_>> {
        (index < self.len).then_some(ReadView {
            columns: self,
            index,
        })
    }

    /// Iterate over reads as borrowed rows.
    pub fn iter(&self) -> impl Iterator<Item = ReadView<'_>> + '_ {
        (0..self.len).map(move |index| ReadView {
            columns: self,
            index,
        })
    }

    /// Build a new collection from the reads at `indices`, in the order given.
    ///
    /// This is how filtering, downsampling and reordering work: choose indices, then
    /// gather once. Only the columns that exist are gathered. An index may repeat, which
    /// duplicates that read.
    ///
    /// Indices outside `0..len()` are skipped, so the result can be shorter than
    /// `indices`. That is checked in debug builds, because a caller producing them is
    /// losing reads it thinks it selected.
    pub fn select(&self, indices: &[usize]) -> Self {
        debug_assert!(
            indices.iter().all(|&i| i < self.len),
            "select() called with an index outside 0..{}",
            self.len
        );
        let mut builder = ReadColumnsBuilder::with_capacity(indices.len());
        for &i in indices {
            if let Some(view) = self.get(i) {
                builder.push_view(&view);
            }
        }
        builder.finish()
    }
}

/// One read, borrowed from the columns.
///
/// Field accessors mirror [`ReadMetrics`](crate::ReadMetrics), so row-wise code — writing
/// a TSV line, evaluating a filter predicate — reads the same as it did against a struct.
#[derive(Debug, Clone, Copy)]
pub struct ReadView<'a> {
    columns: &'a ReadColumns,
    index: usize,
}

impl<'a> ReadView<'a> {
    pub fn index(&self) -> usize {
        self.index
    }

    pub fn read_id(&self) -> Option<&'a str> {
        self.columns.read_id.as_ref()?.get(self.index)
    }

    pub fn length(&self) -> u32 {
        self.columns.length[self.index]
    }

    pub fn quality(&self) -> Option<f64> {
        finite(self.columns.quality.as_ref()?.get(self.index).copied()?)
    }

    pub fn aligned_length(&self) -> Option<u32> {
        let v = *self.columns.aligned_length.as_ref()?.get(self.index)?;
        (v != NO_U32).then_some(v)
    }

    pub fn aligned_quality(&self) -> Option<f64> {
        finite(
            self.columns
                .aligned_quality
                .as_ref()?
                .get(self.index)
                .copied()?,
        )
    }

    pub fn mapping_quality(&self) -> Option<u8> {
        let v = *self.columns.mapping_quality.as_ref()?.get(self.index)?;
        (v != NO_MAPQ).then_some(v)
    }

    pub fn percent_identity(&self) -> Option<f64> {
        finite(
            self.columns
                .percent_identity
                .as_ref()?
                .get(self.index)
                .copied()?,
        )
    }

    pub fn channel_id(&self) -> Option<u16> {
        let v = *self.columns.channel_id.as_ref()?.get(self.index)?;
        (v != NO_CHANNEL).then_some(v)
    }

    pub fn start_time(&self) -> Option<DateTime<Utc>> {
        let v = *self.columns.start_time.as_ref()?.get(self.index)?;
        if v == NO_TIME {
            return None;
        }
        Utc.timestamp_nanos(v).into()
    }

    pub fn duration(&self) -> Option<f64> {
        finite(self.columns.duration.as_ref()?.get(self.index).copied()?)
    }

    pub fn barcode(&self) -> Option<&'a str> {
        self.columns.barcode.as_ref()?.get(self.index)
    }

    pub fn run_id(&self) -> Option<&'a str> {
        self.columns.run_id.as_ref()?.get(self.index)
    }

    pub fn dataset(&self) -> Option<&'a str> {
        self.columns.dataset.as_ref()?.get(self.index)
    }

    /// Materialise this read as an owned [`ReadMetrics`].
    pub fn to_owned_metrics(&self) -> ReadMetrics {
        ReadMetrics {
            read_id: self.read_id().map(str::to_string),
            length: self.length(),
            quality: self.quality(),
            aligned_length: self.aligned_length(),
            aligned_quality: self.aligned_quality(),
            mapping_quality: self.mapping_quality(),
            percent_identity: self.percent_identity(),
            channel_id: self.channel_id(),
            start_time: self.start_time(),
            duration: self.duration(),
            barcode: self.barcode().map(str::to_string),
            run_id: self.run_id().map(str::to_string),
            dataset: self.dataset().map(str::to_string),
        }
    }
}

fn finite(value: f32) -> Option<f64> {
    value.is_finite().then(|| f64::from(value))
}

/// Appends reads into columns, allocating each column only when a value first appears.
///
/// A column that first sees a value at read *n* is back-filled with *n* absent markers,
/// so a field that only some reads carry still costs nothing until it does.
#[derive(Debug, Default)]
pub struct ReadColumnsBuilder {
    capacity: usize,
    columns: ReadColumns,
}

impl ReadColumnsBuilder {
    pub fn new() -> Self {
        Self::with_capacity(0)
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            capacity,
            columns: ReadColumns {
                length: Vec::with_capacity(capacity),
                ..ReadColumns::default()
            },
        }
    }

    pub fn len(&self) -> usize {
        self.columns.len
    }

    pub fn is_empty(&self) -> bool {
        self.columns.len == 0
    }

    /// Append one read, consuming the transient row struct.
    pub fn push(&mut self, read: ReadMetrics) {
        let len = self.columns.len;
        let capacity = self.capacity;

        self.columns.length.push(read.length);

        push_opt_f32(&mut self.columns.quality, len, capacity, read.quality);
        push_opt_u32(
            &mut self.columns.aligned_length,
            len,
            capacity,
            read.aligned_length,
        );
        push_opt_f32(
            &mut self.columns.aligned_quality,
            len,
            capacity,
            read.aligned_quality,
        );
        push_opt_mapq(
            &mut self.columns.mapping_quality,
            len,
            capacity,
            read.mapping_quality,
        );
        push_opt_f32(
            &mut self.columns.percent_identity,
            len,
            capacity,
            read.percent_identity,
        );
        push_opt_channel(&mut self.columns.channel_id, len, capacity, read.channel_id);
        push_opt_time(&mut self.columns.start_time, len, capacity, read.start_time);
        push_opt_f32(&mut self.columns.duration, len, capacity, read.duration);

        push_string(
            &mut self.columns.read_id,
            len,
            capacity,
            read.read_id.as_deref(),
        );
        push_dict(
            &mut self.columns.barcode,
            len,
            capacity,
            read.barcode.as_deref(),
        );
        push_dict(
            &mut self.columns.run_id,
            len,
            capacity,
            read.run_id.as_deref(),
        );
        push_dict(
            &mut self.columns.dataset,
            len,
            capacity,
            read.dataset.as_deref(),
        );

        self.columns.len += 1;
    }

    /// Append a read borrowed from another collection, without materialising a row.
    pub fn push_view(&mut self, view: &ReadView<'_>) {
        let len = self.columns.len;
        let capacity = self.capacity;

        self.columns.length.push(view.length());

        push_opt_f32(&mut self.columns.quality, len, capacity, view.quality());
        push_opt_u32(
            &mut self.columns.aligned_length,
            len,
            capacity,
            view.aligned_length(),
        );
        push_opt_f32(
            &mut self.columns.aligned_quality,
            len,
            capacity,
            view.aligned_quality(),
        );
        push_opt_mapq(
            &mut self.columns.mapping_quality,
            len,
            capacity,
            view.mapping_quality(),
        );
        push_opt_f32(
            &mut self.columns.percent_identity,
            len,
            capacity,
            view.percent_identity(),
        );
        push_opt_channel(
            &mut self.columns.channel_id,
            len,
            capacity,
            view.channel_id(),
        );
        push_opt_time(
            &mut self.columns.start_time,
            len,
            capacity,
            view.start_time(),
        );
        push_opt_f32(&mut self.columns.duration, len, capacity, view.duration());

        push_string(&mut self.columns.read_id, len, capacity, view.read_id());
        push_dict(&mut self.columns.barcode, len, capacity, view.barcode());
        push_dict(&mut self.columns.run_id, len, capacity, view.run_id());
        push_dict(&mut self.columns.dataset, len, capacity, view.dataset());

        self.columns.len += 1;
    }

    pub fn finish(self) -> ReadColumns {
        self.columns
    }
}

// Column appenders. Each lazily creates the column on the first present value,
// back-filling the reads that came before with the absent marker.

fn push_opt_f32(column: &mut Option<Vec<f32>>, len: usize, capacity: usize, value: Option<f64>) {
    match (column.as_mut(), value) {
        (Some(c), v) => c.push(v.map_or(f32::NAN, |v| v as f32)),
        (None, None) => {}
        (None, Some(v)) => {
            let mut c = Vec::with_capacity(capacity.max(len + 1));
            c.resize(len, f32::NAN);
            c.push(v as f32);
            *column = Some(c);
        }
    }
}

fn push_opt_u32(column: &mut Option<Vec<u32>>, len: usize, capacity: usize, value: Option<u32>) {
    match (column.as_mut(), value) {
        (Some(c), v) => c.push(v.unwrap_or(NO_U32)),
        (None, None) => {}
        (None, Some(v)) => {
            let mut c = Vec::with_capacity(capacity.max(len + 1));
            c.resize(len, NO_U32);
            c.push(v);
            *column = Some(c);
        }
    }
}

fn push_opt_mapq(column: &mut Option<Vec<u8>>, len: usize, capacity: usize, value: Option<u8>) {
    match (column.as_mut(), value) {
        (Some(c), v) => c.push(v.unwrap_or(NO_MAPQ)),
        (None, None) => {}
        (None, Some(v)) => {
            let mut c = Vec::with_capacity(capacity.max(len + 1));
            c.resize(len, NO_MAPQ);
            c.push(v);
            *column = Some(c);
        }
    }
}

fn push_opt_channel(
    column: &mut Option<Vec<u16>>,
    len: usize,
    capacity: usize,
    value: Option<u16>,
) {
    match (column.as_mut(), value) {
        (Some(c), v) => c.push(v.unwrap_or(NO_CHANNEL)),
        (None, None) => {}
        (None, Some(v)) => {
            let mut c = Vec::with_capacity(capacity.max(len + 1));
            c.resize(len, NO_CHANNEL);
            c.push(v);
            *column = Some(c);
        }
    }
}

fn push_opt_time(
    column: &mut Option<Vec<i64>>,
    len: usize,
    capacity: usize,
    value: Option<DateTime<Utc>>,
) {
    // Nanoseconds since the epoch spans 1677-2262, which covers any sequencing run. A
    // timestamp outside that range is stored as absent rather than wrapping.
    let nanos = value.and_then(|t| t.timestamp_nanos_opt());
    match (column.as_mut(), nanos) {
        (Some(c), v) => c.push(v.unwrap_or(NO_TIME)),
        (None, None) => {}
        (None, Some(v)) => {
            let mut c = Vec::with_capacity(capacity.max(len + 1));
            c.resize(len, NO_TIME);
            c.push(v);
            *column = Some(c);
        }
    }
}

fn push_string(
    column: &mut Option<StringColumn>,
    len: usize,
    capacity: usize,
    value: Option<&str>,
) {
    match (column.as_mut(), value) {
        (Some(c), v) => c.push(v),
        (None, None) => {}
        (None, Some(v)) => {
            let mut c = StringColumn::with_capacity(capacity.max(len + 1));
            for _ in 0..len {
                c.push(None);
            }
            c.push(Some(v));
            *column = Some(c);
        }
    }
}

fn push_dict(column: &mut Option<DictColumn>, len: usize, capacity: usize, value: Option<&str>) {
    match (column.as_mut(), value) {
        (Some(c), v) => c.push(v),
        (None, None) => {}
        (None, Some(v)) => {
            let mut c = DictColumn::with_capacity(capacity.max(len + 1));
            for _ in 0..len {
                c.push(None);
            }
            c.push(Some(v));
            *column = Some(c);
        }
    }
}

/// Tag every read with a dataset name, replacing whatever was there.
pub(crate) fn set_dataset(columns: &mut ReadColumns, name: &str) {
    let len = columns.len;
    columns
        .dataset
        .get_or_insert_with(DictColumn::default)
        .fill(len, name);
}

/// Concatenate columns, preserving order.
pub(crate) fn concat(parts: Vec<ReadColumns>) -> ReadColumns {
    if parts.len() == 1 {
        return parts.into_iter().next().unwrap_or_default();
    }
    let total: usize = parts.iter().map(ReadColumns::len).sum();
    let mut builder = ReadColumnsBuilder::with_capacity(total);
    for part in &parts {
        for view in part.iter() {
            builder.push_view(&view);
        }
    }
    builder.finish()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;

    fn row(length: u32) -> ReadMetrics {
        ReadMetrics::new(None, length)
    }

    fn build(rows: Vec<ReadMetrics>) -> ReadColumns {
        let mut b = ReadColumnsBuilder::with_capacity(rows.len());
        for r in rows {
            b.push(r);
        }
        b.finish()
    }

    /// Every allocated column must be exactly as long as the collection, whenever the
    /// first present value appears. This is the invariant the whole layout rests on.
    #[test]
    fn columns_stay_aligned_when_a_field_appears_late() {
        for first_present in [0usize, 1, 7, 63] {
            let rows: Vec<ReadMetrics> = (0..64)
                .map(|i| {
                    let mut r = row(100 + i as u32);
                    if i == first_present {
                        r.quality = Some(20.0);
                        r.channel_id = Some(5);
                        r.read_id = Some(format!("read{}", i));
                        r.run_id = Some("run".into());
                    }
                    r
                })
                .collect();
            let c = build(rows);

            assert_eq!(c.len(), 64);
            assert_eq!(c.lengths().len(), 64, "first_present={}", first_present);
            assert_eq!(c.qualities_raw().unwrap().len(), 64);
            assert_eq!(c.channel_ids_raw().unwrap().len(), 64);

            for i in 0..64 {
                let v = c.get(i).unwrap();
                assert_eq!(v.length(), 100 + i as u32);
                let expected = i == first_present;
                assert_eq!(v.quality().is_some(), expected, "i={}", i);
                assert_eq!(v.channel_id().is_some(), expected, "i={}", i);
                assert_eq!(v.read_id().is_some(), expected, "i={}", i);
                assert_eq!(v.run_id().is_some(), expected, "i={}", i);
            }
        }
    }

    /// A field no read carries must not allocate a column at all.
    #[test]
    fn absent_fields_allocate_nothing() {
        let c = build(vec![row(100), row(200)]);
        assert!(c.qualities_raw().is_none());
        assert!(c.aligned_lengths_raw().is_none());
        assert!(c.mapping_qualities_raw().is_none());
        assert!(c.start_times_raw().is_none());
        assert!(c.barcodes().is_none());
        assert!(!c.has_quality() && !c.has_alignment() && !c.has_time() && !c.has_channel());
    }

    /// An empty id is a legitimate value and must not read back as absent.
    #[test]
    fn string_column_distinguishes_empty_from_absent() {
        let mut a = row(1);
        a.read_id = Some(String::new());
        let mut b = row(2);
        b.read_id = Some("xyz".into());
        let c = build(vec![a, b, row(3)]);

        assert_eq!(c.get(0).unwrap().read_id(), Some(""));
        assert_eq!(c.get(1).unwrap().read_id(), Some("xyz"));
        assert_eq!(c.get(2).unwrap().read_id(), None);
    }

    /// Repeated values share one dictionary entry; distinct ones do not collide.
    #[test]
    fn dict_column_deduplicates() {
        let rows: Vec<ReadMetrics> = ["a", "b", "a", "a", "b"]
            .iter()
            .map(|v| {
                let mut r = row(1);
                r.run_id = Some((*v).into());
                r
            })
            .collect();
        let c = build(rows);

        let got: Vec<Option<&str>> = c.iter().map(|v| v.run_id()).collect();
        assert_eq!(
            got,
            vec![Some("a"), Some("b"), Some("a"), Some("a"), Some("b")]
        );
    }

    /// Concatenating parts with different column sets must widen, not lose data.
    #[test]
    fn concat_handles_differing_column_sets() {
        let mut with_qual = row(10);
        with_qual.quality = Some(12.5);
        let mut with_aln = row(30);
        with_aln.aligned_length = Some(25);

        let c = concat(vec![
            build(vec![with_qual]),
            build(vec![row(20)]),
            build(vec![with_aln]),
        ]);

        assert_eq!(c.len(), 3);
        assert_eq!(c.qualities_raw().unwrap().len(), 3);
        assert_eq!(c.aligned_lengths_raw().unwrap().len(), 3);
        assert_eq!(c.get(0).unwrap().quality(), Some(12.5));
        assert_eq!(c.get(1).unwrap().quality(), None);
        assert_eq!(c.get(2).unwrap().aligned_length(), Some(25));
        assert_eq!(c.get(0).unwrap().aligned_length(), None);
    }

    #[test]
    fn select_duplicates_reorders_and_skips_out_of_range() {
        let c = build((0..5).map(|i| row(i * 10)).collect());

        let picked = c.select(&[3, 0, 3]);
        assert_eq!(
            picked.lengths(),
            &[30, 0, 30],
            "duplicates and order must be honoured"
        );

        assert_eq!(c.select(&[]).len(), 0);

        // Out of range is skipped; the debug_assert fires only in debug builds, so this
        // is checked for the release behaviour.
        #[cfg(not(debug_assertions))]
        assert_eq!(c.select(&[1, 99]).lengths(), &[10]);
    }

    /// Tagging a dataset must not invent a name for a part with no reads.
    #[test]
    fn empty_tracked_dataset_has_no_name() {
        let mut empty = ReadColumns::default();
        set_dataset(&mut empty, "only");
        assert_eq!(empty.len(), 0);
        assert!(empty.datasets().unwrap().values().is_empty());
    }

    #[test]
    fn set_dataset_replaces_existing_values() {
        let rows: Vec<ReadMetrics> = ["x", "y"]
            .iter()
            .map(|v| {
                let mut r = row(1);
                r.dataset = Some((*v).into());
                r
            })
            .collect();
        let mut c = build(rows);
        set_dataset(&mut c, "z");

        assert_eq!(c.datasets().unwrap().values(), &["z".to_string()]);
        assert!(c.iter().all(|v| v.dataset() == Some("z")));
    }

    /// Channel 0 occurs in real summary files and must survive; the sentinel is u16::MAX.
    #[test]
    fn channel_zero_is_a_real_value() {
        let mut r = row(1);
        r.channel_id = Some(0);
        let c = build(vec![r, row(2)]);

        assert_eq!(c.get(0).unwrap().channel_id(), Some(0));
        assert_eq!(c.get(1).unwrap().channel_id(), None);
    }

    /// Values that collide with a sentinel read back as absent. Documented behaviour,
    /// pinned here so a change to it is deliberate.
    #[test]
    fn sentinel_values_read_back_as_absent() {
        let mut r = row(1);
        r.mapping_quality = Some(255); // SAM's own "unavailable"
        r.aligned_length = Some(u32::MAX);
        r.quality = Some(f64::NAN);
        let c = build(vec![r]);

        let v = c.get(0).unwrap();
        assert_eq!(v.mapping_quality(), None);
        assert_eq!(v.aligned_length(), None);
        assert_eq!(v.quality(), None);
    }

    #[test]
    fn timestamps_survive_a_round_trip_at_nanosecond_precision() {
        let t = Utc.timestamp_opt(3733, 25_749_999).single().unwrap();
        let mut r = row(1);
        r.start_time = Some(t);
        let c = build(vec![r, row(2)]);

        assert_eq!(c.get(0).unwrap().start_time(), Some(t));
        assert_eq!(c.get(1).unwrap().start_time(), None);
        assert!(c.has_time());
    }

    #[test]
    fn get_and_iter_agree_and_stop_at_the_end() {
        let c = build((0..3).map(row).collect());
        assert!(c.get(3).is_none());
        assert_eq!(c.iter().count(), 3);
        let via_get: Vec<u32> = (0..c.len()).map(|i| c.get(i).unwrap().length()).collect();
        let via_iter: Vec<u32> = c.iter().map(|v| v.length()).collect();
        assert_eq!(via_get, via_iter);
    }

    #[test]
    fn view_round_trips_through_owned_metrics() {
        let mut r = row(1234);
        r.read_id = Some("id".into());
        r.quality = Some(12.5);
        r.channel_id = Some(7);
        r.barcode = Some("bc01".into());
        let c = build(vec![r]);

        let owned = c.get(0).unwrap().to_owned_metrics();
        assert_eq!(owned.read_id.as_deref(), Some("id"));
        assert_eq!(owned.length, 1234);
        assert_eq!(owned.quality, Some(12.5));
        assert_eq!(owned.channel_id, Some(7));
        assert_eq!(owned.barcode.as_deref(), Some("bc01"));

        // And back again, unchanged.
        let again = build(vec![owned]);
        assert_eq!(again.get(0).unwrap().read_id(), Some("id"));
        assert_eq!(again.get(0).unwrap().quality(), Some(12.5));
    }

    #[test]
    fn concat_of_one_part_is_the_part() {
        let c = concat(vec![build(vec![row(1), row(2)])]);
        assert_eq!(c.lengths(), &[1, 2]);
        assert_eq!(concat(vec![]).len(), 0);
    }
}
