// Simple example using the convenience API

use nanoget_rs::convenience::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Simple one-liner for common use cases
    let fastq_metrics = extract_from_fastq("reads.fastq")?;
    let bam_metrics = extract_from_bam("alignments.bam")?;
    let fasta_metrics = extract_from_fasta("sequences.fasta")?;
    
    // Print some basic stats
    println!("FASTQ: {} reads, mean length: {:.0}", 
             fastq_metrics.summary.read_count,
             fastq_metrics.summary.length_stats.mean);
    
    println!("BAM: {} reads, mean length: {:.0}", 
             bam_metrics.summary.read_count,
             bam_metrics.summary.length_stats.mean);
    
    println!("FASTA: {} reads, mean length: {:.0}", 
             fasta_metrics.summary.read_count,
             fasta_metrics.summary.length_stats.mean);
    
    Ok(())
}