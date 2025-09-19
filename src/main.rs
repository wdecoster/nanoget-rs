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
            
            // Output results based on format
            match args.output_format.as_str() {
                "json" => println!("{}", serde_json::to_string_pretty(&metrics)?),
                "csv" => {
                    // TODO: Implement CSV output
                    println!("CSV output not yet implemented");
                }
                _ => println!("{:#?}", metrics),
            }
        }
    }
    
    Ok(())
}