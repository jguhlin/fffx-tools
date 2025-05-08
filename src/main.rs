use std::io::{BufRead, BufWriter, Write};

use clap::{Args, Parser, Subcommand, ValueEnum};
use needletail::{parse_fastx_file, FastxReader, Sequence};
use simdutf8::basic::from_utf8;
use rayon::prelude::*;
use bytecount::count;
use memchr::memmem::Finder;

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

    #[command(about = "Reverses the sanitization process.")]
    Desanitize { file_base: String, output: String },

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
    FastqStats { filename: String },
    #[command(
        about = "Compute FASTA stats: number of reads, total length, GC content, N content, length distribution."
    )]
    FastaStats { 
        filename: Vec<String>, 

        #[arg(short, long, default_value_t = false)]
        header: bool,
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
        Commands::Desanitize { file_base, output } => {
            desanitize(file_base, output);
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
        Commands::FastaStats { filename, header } => {
            fasta_stats(*header, filename);
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
    if &format!("{}.fasta", output_base) == filename {
        panic!(
            "Output file cannot be the same as the input file - Set base name to something else"
        );
    }

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

// todo untested
fn desanitize(file_base: &str, output: &str) {
    if &format!("{}.fasta", output) == file_base {
        panic!(
            "Output file cannot be the same as the input file - Set base name to something else"
        );
    }

    let mut reader = parse_fastx_file(&format!("{}.fasta", file_base)).expect("invalid path/file");

    let mut translation_table: std::collections::HashMap<String, String> =
        std::collections::HashMap::new();

    let mut translation_table_reader =
        std::fs::File::open(format!("{}.translation_table.tsv", file_base))
            .expect("Unable to open translation table");
    let mut translation_table_reader = std::io::BufReader::new(translation_table_reader);

    for line in translation_table_reader.lines() {
        let line = line.expect("Unable to read line");
        let mut line = line.split('\t');
        let id = line.next().unwrap();
        let new_id = line.next().unwrap();
        translation_table.insert(new_id.to_string(), id.to_string());
    }

    let mut output_fasta =
        std::fs::File::create(format!("{}", output)).expect("Unable to create file");
    let mut output_fasta = BufWriter::new(output_fasta);

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.normalize(false);
        let id = from_utf8(record.id()).unwrap();

        let new_id = translation_table.get(id).unwrap();

        output_fasta
            .write(format!(">{}\n", new_id).as_bytes())
            .expect("Unable to write data");
        output_fasta.write_all(&seq).expect("Unable to write data");
        output_fasta.write(b"\n").expect("Unable to write data");
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

    let mut total_gc = 0;
    let mut total_n = 0;
    let mut lengths: Vec<usize> = Vec::new();
    let mut average_read_qualities: Vec<f64> = Vec::new();

    while let Some(record) = reader.next() {
        let record = record.expect("Invalid record");
        let seq = record.seq();
        let qual = record.qual();

        lengths.push(seq.len());

        let gc = count(seq.as_ref(), b'G') + count(seq.as_ref(), b'C');
        let n = count(seq.as_ref(), b'N');

        total_gc += gc;
        total_n += n;

        if let Some(qual) = qual {
            // Get average
            let total_qual = qual.iter().map(|&x| x as f64).sum::<f64>();
            let average_qual = (total_qual / qual.len() as f64) - 33.0;
            average_read_qualities.push(average_qual as f64);
        }
    }

    let total_length = lengths.iter().sum::<usize>();
    let total_reads = lengths.len();

    let read_min = *lengths.iter().min().unwrap();
    let read_max = *lengths.iter().max().unwrap();

    let mean = total_length as f64 / total_reads as f64;
    let median = lengths[lengths.len() / 2];

    // Stored as f64
    let read_qual_min = *average_read_qualities
        .iter()
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let read_qual_max = *average_read_qualities
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    let quality_mean =
        average_read_qualities.iter().sum::<f64>() / average_read_qualities.len() as f64;

    // Calculate N50
    let mut total = 0;
    let mut n50 = 0;

    // sort
    lengths.sort();

    for length in lengths.iter().rev() {
        total += length;
        if total >= total_length / 2 {
            n50 = *length;
            break;
        }
    }

    println!("Total reads: {}", total_reads);
    println!("Total length: {}", total_length);
    println!("GC content: {:.2}", total_gc as f64 / total_length as f64);
    println!("N content: {}", total_n as f64 / total_length as f64);
    println!("Length N50 min/max: {} -  {} / {}", n50, read_min, read_max);
    println!("Mean: {:.2}", mean);
    println!("Median: {}", median);
    println!(
        "Quality mean (min/max): {:.2} ({:.2}/{:.2})",
        quality_mean, read_qual_min, read_qual_max
    );
}

struct FastaStats {
    /// The input file
    filename: String,
    entries: usize,
    length: usize,
    gc: f32,
    n: usize,
    n50: usize,
    n90: usize,
    mean_contig: f32,
    mean_scaffold: f32,
}

// about = "Compute FASTA stats: number of reads, total length, GC content, N content, length distribution."
fn fasta_stats(header: bool, filenames: &Vec<String>) {
    // Let's print it out as a TSV entry
    if header {
        println!("Filename\tEntries\tLength\tGC\tN\tN50\tN90\tMean contig\tMean scaffold");
    }

    let scaffold_finder = Finder::new(&[b'N'; 10]);  // build once

    let stats: Vec<FastaStats> = filenames.par_iter().map(|file| {
        let mut reader = parse_fastx_file(&file).expect("invalid path/file");

        let mut total_length = 0;
        let mut total_gc = 0;
        let mut total_n = 0;
        let mut total_landmarks = 0;

        // Calculate N50, N90, etc.
        let mut lengths: Vec<usize> = Vec::new();

        // Calculate Per contig and Per Scaffold stats
        let mut contig_lengths: Vec<usize> = Vec::new();
        let mut scaffold_lengths: Vec<usize> = Vec::new();

        while let Some(record) = reader.next() {
            let record = record.expect("Invalid record");
            let seq = record.seq();

            total_landmarks += 1;
            total_length += seq.len();

            let gc = count(seq.as_ref(), b'G') + count(seq.as_ref(), b'C');
            let n = count(seq.as_ref(), b'N');
    
            total_gc += gc;
            total_n += n;

            lengths.push(seq.len());

            // Scaffolds have at least one gap that is represented by a N
            // and that gap is 10bp long
            // todo fix this
            
            if scaffold_finder.find(seq.as_ref()).is_some() {
                scaffold_lengths.push(seq.len()); // Changed from len to seq.len()
            } else {
                contig_lengths.push(seq.len()); // Changed from len to seq.len()
            }
        }

        // Calculate N50, N90, etc.
        lengths.sort_unstable();
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

        FastaStats {
            filename: file.to_string(),
            entries: total_landmarks,
            length: total_length,
            gc: total_gc as f32 / total_length as f32,
            n: total_n,
            n50,
            n90,
            mean_contig: mean_contig_length as f32,
            mean_scaffold: mean_scaffold_length as f32,
        }
    }).collect();

    for stat in stats {
        println!(
            "{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{:.2}\t{:.2}",
            stat.filename,
            stat.entries,
            stat.length,
            stat.gc,
            stat.n,
            stat.n50,
            stat.n90,
            stat.mean_contig,
            stat.mean_scaffold
        );
    }

}

#[cfg(test)]
mod test {

    // Test file
    // test_data/test_input.fasta

    use std::fs;

    #[test]
    fn test_sanitize() {

        fs::create_dir("test_data/test").unwrap();
        let filename = "test_data/test_input.fasta";
        let output_base = "test_data/test/test_output";
        super::sanitze(filename, output_base);

        let output_fasta = std::fs::read_to_string(format!("{}.fasta", output_base)).unwrap();
        let output_translation_table =
            std::fs::read_to_string(format!("{}.translation_table.tsv", output_base)).unwrap();

        let expected_fasta = std::fs::read_to_string("test_data/test_output.fasta").unwrap();
        let expected_translation_table =
            std::fs::read_to_string("test_data/test/test_output.translation_table.tsv").unwrap();

        assert_eq!(output_fasta, expected_fasta);
        assert_eq!(output_translation_table, expected_translation_table);

        // Delete files
        std::fs::remove_file(format!("{}.fasta", output_base)).unwrap();
        std::fs::remove_file(format!("{}.translation_table.tsv", output_base)).unwrap();
    }
}
