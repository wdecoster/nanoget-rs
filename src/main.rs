use clap::Parser;

mod cli;
mod error;
mod extract;
mod formats;
mod metrics;
mod utils;

use crate::cli::{Cli, Commands};
use crate::error::NanogetError;

fn main() -> Result<(), NanogetError> {
    env_logger::init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Extract(args) => {
            let metrics = extract::extract_metrics(&args)?;

            // Generate output based on format
            let output = match args.output_format.as_str() {
                "json" => serde_json::to_string_pretty(&metrics)?,
                "tsv" => metrics.to_tsv()?,
                _ => format!("{:#?}", metrics),
            };

            // Write to file or stdout
            if let Some(output_path) = &args.output {
                std::fs::write(output_path, output)?;
            } else {
                println!("{}", output);
            }
        }
    }

    Ok(())
}
