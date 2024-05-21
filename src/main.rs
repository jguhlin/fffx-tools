use std::io::{BufRead, BufWriter, Write};

use clap::{Args, Parser, Subcommand, ValueEnum};
use needletail::{parse_fastx_file, FastxReader, Sequence};
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
    #[command(
        about = "Sanitizes a fasta file. Removes all non-standard characters from the sequence and renumbers the sequences."
    )]
    Sanitize {
        filename: String,
        output_base: String,
    },

    #[command(about = "Chunks a fasta file into smaller files. Useful for parallel processing.")]
    Chunk { filename: String, chunk_size: u64 },

    #[command(
        about = "Filters sequences by length. Can be used to filter out short or long sequences."
    )]
    LengthFilter {
        filename: String,
        output: String,
        min_length: Option<u64>,
        max_length: Option<u64>,
    },

    #[command(
        about = "Removes sequences by matching IDs found in a keyfile. Useful for fast filtering of files."
    )]
    RemoveByIdKeyfile {
        filename: String,
        output: String,
        ids: String,
    },

    #[command(
        about = "Renumbers the sequences in a fasta file - No translation table is provided. Useful for dealing with duplicates."
    )]
    Renumber { filename: String, output: String },

    #[command(
        about = "Removes sequences by matching landmarks (scaffold, contigs, etc) in the sequence. Useful for removing multiple sequences at once."
    )]
    Remove {
        filename: String,
        landmarks: Vec<String>,
    },
}

fn main() {
    let cli = Cli::parse();

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Commands::Sanitize {
            filename,
            output_base,
        } => {
            sanitze(filename, output_base);
        }
        Commands::Chunk {
            filename,
            chunk_size,
        } => {
            chunk(filename, chunk_size);
        }
        Commands::LengthFilter {
            filename,
            output,
            min_length,
            max_length,
        } => {
            length_filter(filename, output, min_length, max_length);
        }
        Commands::Renumber { filename, output } => {
            renumber(filename, output);
        }
        Commands::RemoveByIdKeyfile {
            filename,
            output,
            ids,
        } => {
            remove_by_id_keyfile(filename, output, ids);
        }
        Commands::Remove {
            filename,
            landmarks,
        } => {
            remove_by_landmarks(filename, landmarks);
        }
    }
}

// Output to stdout
fn remove_by_landmarks(filename: &str, landmarks: &Vec<String>) {
    let mut reader = parse_fastx_file(&filename).expect("invalid path/file");

    let mut output_fasta = std::io::stdout();

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.normalize(false);
        let id = from_utf8(record.id()).unwrap();

        let mut skip = false;
        for landmark in landmarks {
            if id == *landmark {
                skip = true;
                break;
            }
        }

        if skip {
            continue;
        }

        output_fasta
            .write(format!(">{}\n", id).as_bytes())
            .expect("Unable to write data");
        output_fasta.write_all(&seq).expect("Unable to write data");
        output_fasta.write(b"\n").expect("Unable to write data");
    }
}

fn remove_by_id_keyfile(filename: &str, output: &str, keyfile: &str) {
    let mut reader = std::fs::File::open(keyfile).expect("Unable to open keyfile");
    let mut reader = std::io::BufReader::new(reader);
    // Keyfiles should be one ID per line
    let mut ids: std::collections::HashSet<String> = std::collections::HashSet::new();
    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        // Clear any whitespace
        let line = line.trim();
        ids.insert(line.to_string());
    }

    let mut reader = parse_fastx_file(&filename).expect("invalid path/file");

    let mut output_fasta =
        std::fs::File::create(format!("{}", output)).expect("Unable to create file");
    let mut output_fasta = BufWriter::new(output_fasta);

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.seq();
        let id = from_utf8(record.id()).unwrap();

        if ids.contains(id) {
            continue;
        }

        output_fasta
            .write(format!(">{}\n", id).as_bytes())
            .expect("Unable to write data");
        output_fasta.write_all(&seq).expect("Unable to write data");
        output_fasta.write(b"\n").expect("Unable to write data");
    }
}

fn sanitze(filename: &str, output_base: &str) {
    let mut reader = parse_fastx_file(&filename).expect("invalid path/file");

    let mut id_translation: Vec<(String, String)> = Vec::new();

    let output_fasta =
        std::fs::File::create(format!("{}.fasta", output_base)).expect("Unable to create file");
    let mut output_fasta = BufWriter::new(output_fasta);

    let output_translation_table =
        std::fs::File::create(format!("{}.translation_table.tsv", output_base))
            .expect("Unable to create file");
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

        output_fasta
            .write(format!(">seq{}\n", new_id).as_bytes())
            .expect("Unable to write data");
        output_fasta.write_all(&seq).expect("Unable to write data");
        output_fasta.write(b"\n").expect("Unable to write data");

        record_number += 1;
    }

    // Check for duplicates of the normal ids
    let mut id_translation_map: std::collections::HashMap<String, String> =
        std::collections::HashMap::new();
    for (id, new_id) in id_translation.clone() {
        if id_translation_map.contains_key(&id) {
            panic!("Duplicate id: {}", id);
        }
        id_translation_map.insert(id, new_id);
    }

    for (id, new_id) in id_translation {
        output_translation_table
            .write(format!("{}\t{}\n", id, new_id).as_bytes())
            .expect("Unable to write data");
    }
}

fn renumber(filename: &str, output: &str) {
    let mut reader = parse_fastx_file(&filename).expect("invalid path/file");

    let output_fasta =
        std::fs::File::create(format!("{}", &output)).expect("Unable to create file");
    let mut output_fasta = BufWriter::new(output_fasta);

    let mut record_number = 0;

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.normalize(false);
        let id = from_utf8(record.id()).unwrap();
        let new_id = format!("c{}", record_number);

        output_fasta
            .write(format!(">seq{}\n", new_id).as_bytes())
            .expect("Unable to write data");
        output_fasta.write_all(&seq).expect("Unable to write data");
        output_fasta.write(b"\n").expect("Unable to write data");

        record_number += 1;
    }
}

fn length_filter(filename: &str, output: &str, min_length: &Option<u64>, max_length: &Option<u64>) {
    let mut reader = parse_fastx_file(&filename).expect("invalid path/file");

    let mut output_fasta =
        std::fs::File::create(format!("{}", output)).expect("Unable to create file");
    let mut output_fasta = BufWriter::new(output_fasta);

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.normalize(false);
        let id = from_utf8(record.id()).unwrap();

        if let Some(min_length) = min_length {
            if seq.len() < *min_length as usize {
                continue;
            }
        }

        if let Some(max_length) = max_length {
            if seq.len() > *max_length as usize {
                continue;
            }
        }

        output_fasta
            .write(format!(">{}\n", id).as_bytes())
            .expect("Unable to write data");
        output_fasta.write_all(&seq).expect("Unable to write data");
        output_fasta.write(b"\n").expect("Unable to write data");
    }
}

fn chunk(filename: &str, chunk_size: &u64) {
    let mut reader = parse_fastx_file(&filename).expect("invalid path/file");

    let mut record_number = 0;
    let mut chunk_number = 0;

    let mut output_fasta =
        std::fs::File::create(format!("{}.chunk{}.fasta", filename, chunk_number))
            .expect("Unable to create file");
    let mut output_fasta = BufWriter::new(output_fasta);

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.normalize(false);
        let id = from_utf8(record.id()).unwrap();

        if record_number % chunk_size == 0 {
            output_fasta = BufWriter::new(
                std::fs::File::create(format!("{}.chunk{}.fasta", filename, chunk_number))
                    .expect("Unable to create file"),
            );
            chunk_number += 1;
        }

        output_fasta
            .write(format!(">{}_{}\n", id, record_number).as_bytes())
            .expect("Unable to write data");
        output_fasta.write_all(&seq).expect("Unable to write data");
        output_fasta.write(b"\n").expect("Unable to write data");

        record_number += 1;
    }
}
