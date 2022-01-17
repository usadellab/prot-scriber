#[macro_use]
extern crate lazy_static;

use annotation_process::AnnotationProcess;
use clap::{App, Arg};
use seq_family_reader::parse_seq_families_file;

mod annotation_process;
mod default;
mod generate_hrd_associated_funcs;
mod model_funcs;
mod output_writer;
mod playground;
mod query;
mod seq_family;
mod seq_family_reader;
mod seq_sim_table_reader;

fn main() {
    let matches = App::new("Prot-Scriber")
        .version("0.1.0")
        .about("Assigns short human readable descriptions (HRD) to query biological sequences using reference canditate descriptions. In this, prot-scriber consumes sequence similarity search (Blast or Diamond or similar) results in tabular format. A customized lexical analysis is carried out on the descriptions of these Blast Hits and a resulting HRD is assigned to the query sequences.\nprot-scriber can also apply the same methodology to produce HRDs for sets of biological sequences, i.e. gene families.\nSee below on how to run your favorite sequence similarity search tool, so that it produces tabular results in the format prot-scriber needs them.")
        .arg(
            Arg::new("output")
            .required(true)
            .takes_value(true)
            .short('o')
            .long("output")
            .help("Filename in which the tabular output will be stored."),
        )
        .arg(
            Arg::new("seq-sim-table")
            .required(true)
            .short('s')
            .takes_value(true)
            .long("seq-sim-table")
            .multiple_occurrences(true)
            .help("File in which to find sequence similarity search results in tabular format. Use e.g. Blast or Diamond to produce them. Required columns are: 'qacc sacc qlen qstart qend slen sstart send bitscore stitle'. If given in different order than shown here use the --header arg."),
        )
        .arg(
            Arg::new("header columns of sequence similarity search result table")
            .short('h')
            .takes_value(true)
            .long("header")
            .multiple_occurrences(true)
            .help("Header of the --seq-sim-table (-s) arg. Separated by space (' ') the names of the columns as they appear in the respective table. Required and default columns are 'qacc sacc qlen qstart qend slen sstart send bitscore stitle'. You can have additional columns that will be ignored. If multiple --seq-sim-table (-s) args are provided make sure the --header (-h) args appear in the correct order, e.g. the first -h arg will be used for the first -s arg, the second -h will be used for the second -s and so on."),
        )
        .arg(
            Arg::new("biological sequence families")
            .short('f')
            .takes_value(true)
            .long("seq-families")
            .help("A file in which families of biological sequences are stored, one family per line. Each line must have format 'fam-name TAB gene1,gene2,gene3'. Make sure no gene appears in more than one family."),
        )
        .arg(
            Arg::new("debug")
            .short('d')
            .long("debug")
            .help("Turn debugging information on."),
        )
        .get_matches();

    // Create a new AnnotationProcess instance and provide it with the necessary input data:
    let mut annotation_process = AnnotationProcess::new();

    // Add biological sequence families information, if provided as input by the user:
    if let Some(seq_families) = matches.value_of("seq-families") {
        parse_seq_families_file(seq_families, &mut annotation_process);
    }

    // Set the input sequence similarity search result (SSSR) tables (Blast or Diamond):
    annotation_process.seq_sim_search_tables = matches
        .values_of("seq-sim-table")
        .unwrap()
        .map(|x| (*x).to_string())
        .collect();

    // Execute the Annotation-Process:
    annotation_process.run();

    // Save output:
    if let Some(o) = matches.value_of("output") {
        match output_writer::write_output_table(
            o.to_string(),
            annotation_process.human_readable_descriptions,
        ) {
            Ok(()) => println!("output written to file {}!", o),
            Err(e) => eprintln!(
                "We are sorry, an error occurred when attempting to write output to file {} \n{}",
                o, e
            ),
        };
    }
}
