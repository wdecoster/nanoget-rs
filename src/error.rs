use thiserror::Error;

#[derive(Error, Debug)]
pub enum NanogetError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("BAM/SAM parsing error: {0}")]
    Htslib(#[from] rust_htslib::errors::Error),

    #[error("CSV parsing error: {0}")]
    Csv(#[from] csv::Error),

    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),

    #[error("File not found: {0}")]
    FileNotFound(String),

    #[error("Unsupported file format: {0}")]
    UnsupportedFormat(String),

    // No longer constructed: the read-type validation that raised it moved into the
    // `ReadType` value enum, which clap rejects before extraction starts. Kept because
    // removing a variant from a public error enum is a breaking change.
    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Processing error: {0}")]
    ProcessingError(String),
}
