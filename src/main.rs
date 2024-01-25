use  std::io::BufWriter;

use clap::{Args, Parser, Subcommand, ValueEnum};
use needletail::{parse_fastx_file, Sequence, FastxReader};
use simdutf8::basic::from_utf8;

#[derive(Parser)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Adds files to myapp
    Sanitize { filename: String, output_base: String },
}

fn main() {
    let cli = Cli::parse();

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Commands::Sanitize { filename , output_base} => {
            sanitze(filename, output_base);
        }
    }
}

fn sanitze(filename: &str, output_base: &str) {
    let mut reader = parse_fastx_file(&filename).expect("invalid path/file");

    let mut id_translation: Vec<(String, String)> = Vec::new();

    let output_fasta = std::fs::File::create(format!("{}.fasta", output_base)).expect("Unable to create file");
    let mut output_fasta = BufWriter::new(output_fasta);

    let output_translation_table = std::fs::File::create(format!("{}.translation_table.tsv", output_base)).expect("Unable to create file");
    let mut output_translation_table = BufWriter::new(output_translation_table);

    let mut record_number = 0;

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.seq();
        let id = from_utf8(record.id()).unwrap();
        let new_id = format!("{}", record_number);
        
        id_translation.push((id.to_string(), new_id.to_string()));

        
        
    }


}