# Contributing to nanoget-rs

We welcome contributions to nanoget-rs! This document provides guidelines for contributing to the project.

## Development Setup

1. **Prerequisites:**
   - Rust 1.70+ (stable)
   - Git
   - System dependencies: `libclang-dev`, `libbz2-dev`, `liblzma-dev`, `zlib1g-dev`, `libcurl4-openssl-dev`, `libssl-dev`, `pkg-config`

2. **Setup:**
   ```bash
   git clone https://github.com/wdecoster/nanoget-rs.git
   cd nanoget-rs
   cargo build
   cargo test
   ```

## Development Workflow

1. **Fork the repository** on GitHub
2. **Create a feature branch** from `master`:
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. **Make your changes** following the coding standards below
4. **Test your changes:**
   ```bash
   cargo test
   cargo fmt --check
   cargo clippy -- -D warnings
   ```
5. **Commit your changes** with clear, descriptive messages
6. **Push to your fork** and create a Pull Request

## Coding Standards

- **Formatting:** Use `cargo fmt` to format code
- **Linting:** Fix all `cargo clippy` warnings
- **Testing:** Add tests for new functionality
- **Documentation:** Update documentation for public APIs
- **Error Handling:** Use appropriate error types and handle errors gracefully

## Testing

- **Unit tests:** Test individual functions and modules
- **Integration tests:** Test complete workflows with real data
- **Performance tests:** Ensure changes don't regress performance

Run tests with:
```bash
# All tests
cargo test

# Specific test
cargo test test_name

# With output
cargo test -- --nocapture
```

## Pull Request Guidelines

- **Clear description:** Explain what your PR does and why
- **Link issues:** Reference any related issues
- **Small PRs:** Keep changes focused and reviewable
- **Update docs:** Include documentation updates
- **Add tests:** Ensure new code is tested

## Performance Considerations

nanoget-rs is designed for high-performance processing of large genomic files:

- **Memory efficiency:** Avoid loading entire files into memory
- **Parallel processing:** Use Rayon for CPU-intensive operations
- **Streaming:** Process data incrementally when possible
- **Benchmarking:** Profile performance-critical changes

## Code Review Process

1. **Automated checks** must pass (CI, tests, formatting)
2. **Manual review** by maintainers
3. **Performance review** for changes affecting core algorithms
4. **Documentation review** for public API changes

## Release Process

1. **Version bump** in `Cargo.toml`
2. **Update CHANGELOG.md**
3. **Create release tag:** `git tag v0.x.y`
4. **Push tag:** Automated release will be created

## Getting Help

- **Issues:** Open an issue for bugs or feature requests
- **Discussions:** Use GitHub Discussions for questions
- **Documentation:** Check the README and code comments

## License

By contributing, you agree that your contributions will be licensed under the MIT License.