use clap::Parser;
use std::fs::File;
use std::io::{BufWriter, Write};

use nanoget_rs::{extract_metrics, Cli, Commands, NanogetError, OutputFormat};

fn main() -> Result<(), NanogetError> {
    env_logger::init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Extract(args) => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(args.threads)
                .build()
                .map_err(|e| NanogetError::ProcessingError(e.to_string()))?;

            let metrics = pool.install(|| extract_metrics(&args))?;

            // Stream straight to the destination. Materialising the whole output as a
            // String first costs more memory than the reads themselves for a full run.
            let mut writer: BufWriter<Box<dyn Write>> = BufWriter::new(match &args.output {
                Some(path) => Box::new(File::create(path)?),
                None => Box::new(std::io::stdout().lock()),
            });

            match args.output_format {
                OutputFormat::Json => {
                    serde_json::to_writer_pretty(&mut writer, &metrics)?;
                    writeln!(writer)?;
                }
                OutputFormat::Tsv => metrics.write_tsv(&mut writer)?,
            }

            writer.flush()?;
        }
    }

    Ok(())
}
