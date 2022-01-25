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
    let matches = App::new("prot-scriber")
        .version("version 0.1.0")
        .about("\n1. Summary\n----------\n`prot-scriber` uses reference descriptions (`stitle` in Blast terminology) from sequence similarity search results (Blast Hits) to assign short human readable descriptions (HRD) to query biological sequences or sets of them (a.k.a gene, or sequence, families). In this, prot-scriber consumes sequence similarity search (Blast, Diamond, or similar) results in tabular format. A customized lexical analysis is carried out on the descriptions (`stitle` in Blast terminology) of these Blast Hits and a resulting HRD is assigned to the query sequences or query families, respectively.\n\n2. `prot-scriber` input preparation\n-----------------------------------\nThis sections explains how to run your favorite sequence similarity search tool, so that it produces tabular results in the format prot-scriber needs them. You can run sequence similarity searches with Blast [McGinnis, S. & Madden, T. L. BLAST: at the core of a powerful and diverse set of sequence analysis tools. Nucleic Acids Res 32, W20–W25 (2004).] or Diamond [Buchfink, B., Xie, C. & Huson, D. H. Fast and sensitive protein alignment using DIAMOND. Nat Meth 12, 59–60 (2015).]. Note that there are other tools to carry out sequence similarity searches which can be used to generate the input for prot-scriber. As long as you have a tabular text file with the three required columns holding the query identifier, the subject (`Hit`) identifier, and the subject (`Hit`) description (`stitle` in Blast terminology) prot-scriber will accept this as input.\nDepending on the type of your query sequences the search method and searched reference databases vary. For amino acid queries search protein reference databases, for nucleotide query sequences search nucleotide reference databases. If you have protein coding nucleotide query sequences you can choose to either search protein reference databases using translated nucleotide ueries with `blastx` or `diamond blastx` or search reference nucleotide databases with `blastn` or `diamond blastn`. Note, that before carrying out any sequence similarity searches you need to format your reference databases. This is achieved by either the `makeblastdb` (Blast) or `makedb` (Diamond) commands, respectively. Please see the respective tool's (Blast or Diamond) manual for details on how to format your reference sequence database.\n\n2.1 Which reference databases to search\n---------------------------------------\nFor amino acid (protein) or protein coding nucleotide query sequences we recommend searching UniProt's Swissprot and trEMBL. For nucleotide sequences UniRef100 and, or UniParc might be good choices. Note that you can search _any_ database you deem to hold valuable reference sequences.\n\n2.2 Example Blast or Diamond commands\n-------------------------------------\nNote that the following instructions on how to execute your sequence similarity searches with Blast or Diamond only include the information - in terms of selected output table columns - absolutely required by `prot-scriber`. You are welcome, of course, to have more columns in your tabular output, e.g. `bitscore` or `evalue` etc.\n\n2.2.1 Blast\n-----------\nGenerate prot-scriber input with Blast as follows. Note that you need to search each of your reference databases (`-db`) with a separate Blast command. The following example uses `blastp`, replace it, if your query sequence type makes that necessary with `blastn` or `blastx`.\n\nblastp -db <reference_database.fasta> -query <your_query_sequences.fasta> -num_threads <how-many-do-you-want-to-use> -out <queries_vs_reference_db_name_blastout.txt> -outfmt \"6 delim=<TAB> qacc sacc stitle\"\n\nIt is important to note, that in the above `outfmt` argument the `delim` set to `<TAB>` means you need to actually type in a TAB character. (We write `<TAB>` here, so you see something, not only whitespace.) Typically you can type it by hitting Ctrl+Tab in the terminal.\n\n2.2.2 Diamond\n-------------\n")
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
            Arg::new("header")
            .short('h')
            .takes_value(true)
            .long("header")
            .multiple_occurrences(true)
            .help("Header of the --seq-sim-table (-s) arg. Separated by space (' ') the names of the columns as they appear in the respective table. Required and default columns are 'qacc sacc qlen qstart qend slen sstart send bitscore stitle'. You can have additional columns that will be ignored. If multiple --seq-sim-table (-s) args are provided make sure the --header (-h) args appear in the correct order, e.g. the first -h arg will be used for the first -s arg, the second -h will be used for the second -s and so on."),
        )
        .arg(
            Arg::new("seq-families")
            .short('f')
            .takes_value(true)
            .long("seq-families")
            .help("A file in which families of biological sequences are stored, one family per line. Each line must have format 'fam-name TAB gene1,gene2,gene3'. Make sure no gene appears in more than one family."),
        )
        .arg(
            Arg::new("annotate-non-family-queries")
            .short('a')
            .takes_value(false)
            .long("annotate-non-family-queries")
            .help("If this flag is given, queries for which there are sequence similarity search (Blast) results but that are NOT member of a sequence family will receive an annotation (human readable description) in the output file, too. Default value of this setting is 'OFF' (false)."),
        )
        .arg(
            Arg::new("debug")
            .short('d')
            .long("debug")
            .help("Turn debugging information on."),
        )
        .arg(
            Arg::new("n-threads")
            .short('n')
            .takes_value(true)
            .long("n-threads")
            .help("The maximum number of parallel threads to use. Default is the number of logical cores. Note that at most one thread is used per input sequence similarity search result (Blast table) file. After parsing these annotation may use up to this number of threads to generate human readable descriptions."),
        ).get_matches();

    // Create a new AnnotationProcess instance and provide it with the necessary input data:
    let mut annotation_process = AnnotationProcess::new();

    // Set number of parallel processes to use:
    if let Some(n_threads) = matches.value_of("n-threads") {
        annotation_process.n_threads = n_threads
            .parse()
            .expect("Could not parse argument 'n-threads' ('n') into a positive integer");
    }
    // Set the number of parallel processes to be used by `rayon` (see
    // `AnnotationProcess::process_rest_data`):
    rayon::ThreadPoolBuilder::new()
        .num_threads(annotation_process.n_threads)
        .build_global()
        .expect("Could not set the number of parallel processes to be used to generate human readable descriptions (AnnotationProcess::process_rest_data).");

    // Add biological sequence families information, if provided as input by the user:
    if let Some(seq_families) = matches.value_of("seq-families") {
        parse_seq_families_file(seq_families, &mut annotation_process);
        println!(
            "Loaded {:?} sequence families from {:?}",
            &annotation_process.seq_families.len(),
            &seq_families
        );
    }
    // Shall non family queries also be annotated?
    if matches.is_present("annotate-non-family-queries") {
        annotation_process.annotate_lonely_queries = true;
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
