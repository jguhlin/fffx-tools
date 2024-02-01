use  std::io::{BufWriter, Write};

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
    Chunk { filename: String, chunk_size: u64 }
}

fn main() {
    let cli = Cli::parse();

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Commands::Sanitize { filename , output_base} => {
            sanitze(filename, output_base);
        },
        Commands::Chunk { filename, chunk_size } => {
            chunk(filename, chunk_size);
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
        let seq = record.normalize(false);
        let id = from_utf8(record.id()).unwrap();
        let new_id = format!("c{}", record_number);

        // Split id at first space or "|"
        let id = id.split(|c| c == ' ' || c == '|').next().unwrap();
        
        id_translation.push((new_id.to_string(), id.to_string()));

        output_fasta.write(format!(">seq{}\n", new_id).as_bytes()).expect("Unable to write data");
        output_fasta.write_all(&seq).expect("Unable to write data");
        output_fasta.write(b"\n").expect("Unable to write data");

        record_number += 1;
    }

    // Check for duplicates of the normal ids
    let mut id_translation_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    for (id, new_id) in id_translation.clone() {
        if id_translation_map.contains_key(&id) {
            panic!("Duplicate id: {}", id);           
        }
        id_translation_map.insert(id, new_id);
    }

    for (id, new_id) in id_translation {
        output_translation_table.write(format!("{}\t{}\n", id, new_id).as_bytes()).expect("Unable to write data");
    }
}

fn chunk(filename: &str, chunk_size: &u64) {
    let mut reader = parse_fastx_file(&filename).expect("invalid path/file");

    let mut record_number = 0;
    let mut chunk_number = 0;

    let mut output_fasta = std::fs::File::create(format!("{}.chunk{}.fasta", filename, chunk_number)).expect("Unable to create file");
    let mut output_fasta = BufWriter::new(output_fasta);

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.normalize(false);
        let id = from_utf8(record.id()).unwrap();

        if record_number % chunk_size == 0 {
            output_fasta = BufWriter::new(std::fs::File::create(format!("{}.chunk{}.fasta", filename, chunk_number)).expect("Unable to create file"));
            chunk_number += 1;
        }

        output_fasta.write(format!(">{}_{}\n", id, record_number).as_bytes()).expect("Unable to write data");
        output_fasta.write_all(&seq).expect("Unable to write data");
        output_fasta.write(b"\n").expect("Unable to write data");

        record_number += 1;
    }
}
