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
    #[command(
        about = "Compute FASTQ stats: number of reads, total length, GC content, N content, length distribution, mean, medium, s.d., of quality score."
    )]
    FastqStats {
        filename: String,
    },
    #[command(
        about = "Compute FASTA stats: number of reads, total length, GC content, N content, length distribution."
    )]
    FastaStats {
        filename: String,
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
        Commands::FastqStats { filename } => {
            fastq_stats(filename);
        }
        Commands::FastaStats { filename } => {
            fasta_stats(filename);
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

// about = "Compute FASTQ stats: number of reads, total length, GC content, N content, length distribution, mean, medium, s.d., of quality score."
fn fastq_stats(filename: &str) {
    let mut reader = parse_fastx_file(&filename).expect("invalid path/file");

    let mut total_length = 0;
    let mut total_gc = 0;
    let mut total_n = 0;
    let mut total_reads = 0;
    let mut length_distribution: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();
    let mut quality_distribution: std::collections::HashMap<u8, usize> =
        std::collections::HashMap::new();

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.seq();
        let qual = record.qual();

        total_reads += 1;
        total_length += seq.len();

        let mut gc = 0;
        let mut n = 0;
        for base in seq.iter() {
            match base {
                b'G' | b'C' => gc += 1,
                b'N' => n += 1,
                _ => (),
            }
        }

        total_gc += gc;
        total_n += n;

        *length_distribution.entry(seq.len()).or_insert(0) += 1;

        if let Some(qual) = qual {
            for q in qual.iter() {
                *quality_distribution.entry(*q - 33).or_insert(0) += 1;
            }
        }
    }

    let mean = total_length as f64 / total_reads as f64;
    let mut median = 0;
    let mut median_count = 0;
    let mut median_found = false;
    for (length, count) in length_distribution.iter() {
        median_count += count;
        if median_count >= total_reads / 2 {
            median = *length;
            median_found = true;
            break;
        }
    }

    if !median_found {
        panic!("Unable to find median");
    }

    let mut quality_sum = 0;
    let mut quality_count = 0;
    for (quality, count) in quality_distribution.iter() {
        quality_sum += (*quality as usize) * count;
        quality_count += count;
    }

    let quality_mean = quality_sum as f64 / quality_count as f64;

    let mut quality_variance = 0;
    for (quality, count) in quality_distribution.iter() {
        quality_variance += (*quality as usize - quality_mean as usize).pow(2) * count;
    }

    let quality_sd = (quality_variance as f64 / quality_count as f64).sqrt();

    println!("Total reads: {}", total_reads);
    println!("Total length: {}", total_length);
    println!("GC content: {}", total_gc as f64 / total_length as f64);
    println!("N content: {}", total_n as f64 / total_length as f64);
    println!("Length distribution: {:?}", length_distribution);
    println!("Mean: {}", mean);
    println!("Median: {}", median);
    println!("Quality mean: {}", quality_mean);
    println!("Quality median: {}", quality_distribution.len() / 2);
    println!("Quality s.d.: {}", quality_sd);
}

// about = "Compute FASTA stats: number of reads, total length, GC content, N content, length distribution."
fn fasta_stats(filename: &str) {
    let mut reader = parse_fastx_file(&filename).expect("invalid path/file");

    let mut total_length = 0;
    let mut total_gc = 0;
    let mut total_n = 0;
    let mut total_reads = 0;
    let mut length_distribution: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();
    
    // Calculate N50, N90, etc.
    let mut lengths: Vec<usize> = Vec::new();

    // Calculate Per contig and Per Scaffold stats
    let mut contig_lengths: Vec<usize> = Vec::new();
    let mut scaffold_lengths: Vec<usize> = Vec::new();
    
    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.seq();

        total_reads += 1;
        total_length += seq.len();

        let mut gc = 0;
        let mut n = 0;
        for base in seq.iter() {
            match base {
                b'G' | b'C' => gc += 1,
                b'N' => n += 1,
                _ => (),
            }
        }

        total_gc += gc;
        total_n += n;

        *length_distribution.entry(seq.len()).or_insert(0) += 1;

        lengths.push(seq.len());

        // Scaffolds have at least one gap that is represented by a N
        // and that gap is 10bp long
        // todo fix this
        if seq.contains(&b'N') {
            scaffold_lengths.push(seq.len());
        } else {
            contig_lengths.push(seq.len());
        }
    }

    println!("Total reads: {}", total_reads);
    println!("Total length: {}", total_length);
    println!("GC content: {}", total_gc as f64 / total_length as f64);
    println!("N content: {}", total_n as f64 / total_length as f64);
    println!("Length distribution: {:?}", length_distribution);

    // Calculate N50, N90, etc.
    lengths.sort();
    let mut total = 0;
    let mut n50 = 0;
    let mut n90 = 0;

    for length in lengths.iter().rev() {
        total += length;
        if total >= total_length / 2 {
            n50 = *length;
            break;
        }
    }

    total = 0;
    for length in lengths.iter().rev() {
        total += length;
        if total >= (total_length as f64 * 0.9) as usize {
            n90 = *length;
            break;
        }
    }

    println!("N50: {}", n50);
    println!("N90: {}", n90);

    // Calculate Per contig and Per Scaffold stats
    let mut total_contig_length = 0;
    let mut total_scaffold_length = 0;
    for length in contig_lengths.iter() {
        total_contig_length += length;
    }

    for length in scaffold_lengths.iter() {
        total_scaffold_length += length;
    }

    let mean_contig_length = total_contig_length as f64 / contig_lengths.len() as f64;

    let mean_scaffold_length = total_scaffold_length as f64 / scaffold_lengths.len() as f64;

    println!("Mean contig length: {}", mean_contig_length);
    println!("Mean scaffold length: {}", mean_scaffold_length);



}