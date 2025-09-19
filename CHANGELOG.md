# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial Rust implementation of nanoget
- Support for FASTQ, FASTA, BAM, CRAM, uBAM, and summary files
- Parallel processing capabilities
- Memory-optimized streaming processing
- Comprehensive test suite
- Both CLI binary and library API
- GitHub Actions CI/CD pipeline
- Automated releases on tag push
- Cross-platform binary builds (Linux, macOS, Windows)
- Documentation and examples

### Changed
- Complete rewrite from Python to Rust for better performance
- Enhanced error handling and type safety
- Improved memory efficiency

### Fixed
- All compilation warnings resolved
- Proper error propagation throughout codebase

## [0.1.0] - TBD

### Added
- Initial release of nanoget-rs
- Feature parity with Python nanoget
- Performance improvements through Rust implementation
- Library API for integration with other Rust tools