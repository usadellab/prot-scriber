#[macro_use]
extern crate lazy_static;

use annotation_process::{run, AnnotationProcess};
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
        .about("\n1. Summary\n----------\n'prot-scriber' uses reference descriptions ('stitle' in Blast terminology) from sequence similarity search results (Blast Hits) to assign short human readable descriptions (HRD) to query biological sequences or sets of them (a.k.a gene, or sequence, families). In this, prot-scriber consumes sequence similarity search (Blast, Diamond, or similar) results in tabular format. A customized lexical analysis is carried out on the descriptions ('stitle' in Blast terminology) of these Blast Hits and a resulting HRD is assigned to the query sequences or query families, respectively.\n\n2. 'prot-scriber' input preparation\n-----------------------------------\nThis sections explains how to run your favorite sequence similarity search tool, so that it produces tabular results in the format prot-scriber needs them. You can run sequence similarity searches with Blast [McGinnis, S. & Madden, T. L. BLAST: at the core of a powerful and diverse set of sequence analysis tools. Nucleic Acids Res 32, W20–W25 (2004).] or Diamond [Buchfink, B., Xie, C. & Huson, D. H. Fast and sensitive protein alignment using DIAMOND. Nat Meth 12, 59–60 (2015).]. Note that there are other tools to carry out sequence similarity searches which can be used to generate the input for prot-scriber. As long as you have a tabular text file with the three required columns holding the query identifier, the subject ('Hit') identifier, and the subject ('Hit') description ('stitle' in Blast terminology) prot-scriber will accept this as input.\nDepending on the type of your query sequences the search method and searched reference databases vary. For amino acid queries search protein reference databases, for nucleotide query sequences search nucleotide reference databases. If you have protein coding nucleotide query sequences you can choose to either search protein reference databases using translated nucleotide queries with 'blastx' or 'diamond blastx' or search reference nucleotide databases with 'blastn' or 'diamond blastn'. Note, that before carrying out any sequence similarity searches you need to format your reference databases. This is achieved by either the 'makeblastdb' (Blast) or 'makedb' (Diamond) commands, respectively. Please see the respective tool's (Blast or Diamond) manual for details on how to format your reference sequence database.\n\n2.1 Which reference databases to search\n---------------------------------------\nFor amino acid (protein) or protein coding nucleotide query sequences we recommend searching UniProt's Swissprot and trEMBL. For nucleotide sequences UniRef100 and, or UniParc might be good choices. Note that you can search _any_ database you deem to hold valuable reference sequences.\n\n2.2 Example Blast or Diamond commands\n-------------------------------------\nNote that the following instructions on how to execute your sequence similarity searches with Blast or Diamond only include the information - in terms of selected output table columns - absolutely required by 'prot-scriber'. You are welcome, of course, to have more columns in your tabular output, e.g. 'bitscore' or 'evalue' etc. Note that you need to search each of your reference databases with a separate Blast or Diamond command, respectively.\n\n2.2.1 Blast\n-----------\nGenerate prot-scriber input with Blast as follows. The following example uses 'blastp', replace it, if your query sequence type makes that necessary with 'blastn' or 'blastx'.\n\nblastp -db <reference_database.fasta> -query <your_query_sequences.fasta> -num_threads <how-many-do-you-want-to-use> -out <queries_vs_reference_db_name_blastout.txt> -outfmt \"6 delim=<TAB> qacc sacc stitle\"\n\nIt is important to note, that in the above 'outfmt' argument the 'delim' set to '<TAB>' means you need to actually type in a TAB character. (We write '<TAB>' here, so you see something, not only whitespace.) Typically you can type it by hitting Ctrl+Tab in the terminal.\n\n2.2.2 Diamond\n-------------\nGenerate prot-scriber input with Diamond as follows. The following example uses 'blastp', replace it, if your query sequence type makes that necessary with 'blastn' or 'blastx'.\n\ndiamond blastp -p <how-many-threads-do-you-want-to-use> --quiet -d <reference-database.dmnd> -q <your_query_sequences.fasta> -o <queries_vs_reference_db_name_diamondout.txt> -f 6 qseqid sseqid stitle\n\nNote that diamond by default uses the '<TAB>' character as a field-separator for its output tables.\n\n2.3 Gene Family preparation and analysis\n----------------------------------------\nAssume you have the proteomes of eight crucifer plant species and want to cluster the respective amino acid sequences into gene families. Note that the following example provides code to be executed in a BASH Shell (also available on Windows). We provide a very basic procedure to perform the clustering:\n\n(i) \"All versus all\" Blast or Diamond\n\nAssume all amino acid sequences of the eight example proteomes stored in a single file 'all_proteins.fasta'\nRun:\n\ndiamond makedb --in all_proteins.fasta -d all_proteins.fasta\n\ndiamond blastp --quiet -p <how-many-threads-do-you-want-to-use?> -d all_proteins.fasta.dmnd -q all_proteins.fasta -o all_proteins_vs_all.txt -f 6 qseqid sseqid pident\n\n(ii) Run markov clustering\n\nNote that 'mcl' is a command line tool implementing the original Markov Clustering algorithm [Stijn van Dongen, A cluster algorithm for graphs. Technical Report INS-R0010, National Research Institute for Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000]. On most systems you can install the 'mcl' binary using the respective package manager, e.g. 'sudo apt-get update && sudo apt-get install -y mcl' (Debian / Ubuntu).\n\nmcl all_proteins_vs_all.txt -o all_proteins_gene_clusters.txt --abc -I 2.0\n\n(iii) Add gene family names to mcl output and filter out singleton clusters\n\nNote that we use the GNU tools 'sed' and 'awk' to do some basic post-processing of the 'mcl' output.\n\nsed -e 's/\\t/,/g' all_proteins_gene_clusters.txt | awk -F \",\" 'BEGIN{i=1}{if (NF > 1){print \"Seq-Fam_\" i \"\\t\" $0; i=i+1}}' > all_proteins_gene_families.txt\n\nCongratulations! You now have clustered your eight plant crucifer proteomes into gene families (file 'all_proteins_gene_families.txt').\n\n(iv) Run prot-scriber\n\nWe assume that you ran either 'blastp' or 'diamond blastp' (see section 2.2 for details) to search your selected reference databases with the 'all_proteins.fasta' queries. Here, we assume you have searched UniProt's Swissprot and trEMBL databases.\n\nprot-scriber -f all_proteins_gene_families.txt -s all_proteins_vs_Swissprot_blastout.txt -s all_proteins_vs_trEMBL_blastout.txt -o all_proteins_gene_families_HRDs.txt\n\n\n3. Technical manual\n-------------------\n")
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
            .help("File in which to find sequence similarity search results in tabular format (SSST). Use e.g. Blast or Diamond to produce them. Required columns are: 'qacc sacc stitle' (Blast) or 'qseqid sseqid stitle' (Diamond). (See section '2. prot-scriber input preparation' for more details.) If the required columns, or more, appear in different order than shown here you must use the --header (-e) argument. If any of the input SSSTs uses a different field-separator than the '<TAB>' character, you must provide the --field-separator (-p) argument. You can provide multiple SSSTs, simply by repeating the -s argument, e.g. '-s queries_vs_swissprot_diamond_out.txt -s queries_vs_trembl_diamond_out.txt'. Providing multiple --seq-sim-table (-s) arguments might imply the order in which you give other arguments like --header (-e) and --field-separator (-p). See there for more details."),
        )
        .arg(
            Arg::new("header")
            .short('e')
            .takes_value(true)
            .long("header")
            .multiple_occurrences(true)
            .help("Header of the --seq-sim-table (-s) arg. Separated by space (' ') the names of the columns in order of appearance in the respective table. Required and default columns are 'qacc sacc stitle'. Note that this option only understands Blast terminology, i.e. even if you ran Diamond, please provide 'qacc' instead of 'qseqid' and 'sacc' instead of 'sseqid'. Luckily 'stitle' is 'stitle' in Diamond, too. You can have additional columns that will be ignored, as long as the required columns appear in the correct order. Consider this example: 'qacc sacc evalue bitscore stitle'. If multiple --seq-sim-table (-s) args are provided make sure the --header (-e) args appear in the correct order, e.g. the first -e arg will be used for the first -s arg, the second -e will be used for the second -s and so on. Set to 'default' to use the hard coded default."),
        )
        .arg(
            Arg::new("field-separator")
            .short('p')
            .takes_value(true)
            .long("field-separator")
            .multiple_occurrences(true)
            .help("Field-Separator of the --seq-sim-table (-s) arg. The default value is the '<TAB>' character. Consider this example: '-p @'. If multiple --seq-sim-table (-s) args are provided make sure the --field-separator (-p) args appear in the correct order, e.g. the first -p arg will be used for the first -s arg, the second -p will be used for the second -s and so on. You can provide '-p default' to use the hard coded default (TAB)."),
        )
        .arg(
            Arg::new("seq-families")
            .short('f')
            .takes_value(true)
            .long("seq-families")
            .help("A file in which families of biological sequences are stored, one family per line. Each line must have format 'fam-name TAB gene1,gene2,gene3'. Make sure no gene appears in more than one family."),
        )
        .arg(
            Arg::new("seq-family-id-genes-separator")
            .short('i')
            .takes_value(true)
            .long("seq-family-id-genes-separator")
            .help("A string used as separator in the argument --seq-families (-f) gene families file. This string separates the gene-family-identifier (name) from the gene-identifier list that family comprises. Default is '<TAB>' (\"\\t\")."),
        )
        .arg(
            Arg::new("seq-family-gene-ids-separator")
            .short('g')
            .takes_value(true)
            .long("seq-family-gene-ids-separator")
            .help("A regular expression (Rust syntax) used to split the list of gene-identifiers in the argument --seq-families (-f) gene families file. Default is '(\\s*,\\s*|\\s+)'."),
        )
        .arg(
            Arg::new("annotate-non-family-queries")
            .short('a')
            .takes_value(false)
            .long("annotate-non-family-queries")
            .help("Use this option only in combination with --seq-families (-f), i.e. when prot-scriber is used to generate human readable descriptions for gene families. If in that context this flag is given, queries for which there are sequence similarity search (Blast) results but that are NOT member of a sequence family will receive an annotation (human readable description) in the output file, too. Default value of this setting is 'OFF' (false)."),
        )
        .arg(
            Arg::new("verbose")
            .short('v')
            .long("verbose")
            .help("Print informative messages about the annotation process."),
        )
        .arg(
            Arg::new("n-threads")
            .short('n')
            .takes_value(true)
            .long("n-threads")
            .help("The maximum number of parallel threads to use. Default is the number of logical cores. Required minimum is two (2). Note that at most one thread is used per input sequence similarity search result (Blast table) file. After parsing these annotation may use up to this number of threads to generate human readable descriptions."),
        ).get_matches();

    // Create a new AnnotationProcess instance and provide it with the necessary input data:
    let mut annotation_process = AnnotationProcess::new();

    // Does the user want informative messages printed out?
    if matches.is_present("verbose") {
        annotation_process.verbose = true;
    }

    // Set number of parallel processes to use:
    if let Some(n_threads) = matches.value_of("n-threads") {
        annotation_process.n_threads = n_threads
            .trim()
            .parse()
            .expect("Could not parse argument '--n-threads' ('-n') into a positive integer");
    }
    // Set the number of parallel processes to be used by `rayon` (see
    // `AnnotationProcess::process_rest_data`):
    rayon::ThreadPoolBuilder::new()
        .num_threads(annotation_process.n_threads)
        .build_global()
        .expect("Could not set the number of parallel processes to be used to generate human readable descriptions (AnnotationProcess::process_rest_data).");

    // Add biological sequence families information, if provided as input by the user:
    if let Some(seq_families) = matches.value_of("seq-families") {
        // What is the character that separates a gene-family-identifier from its list of
        // gene-identifiers?
        if matches.is_present("seq-family-id-genes-separator") {
            annotation_process.seq_family_id_genes_separator = matches
                .value_of("seq-family-id-genes-separator")
                .unwrap()
                .trim()
                .to_string();
        }
        // What is the regular expression (string representation) that shall be used to split the list
        // of gene-identifiers a gene-family comprises?
        if matches.is_present("seq-family-gene-ids-separator") {
            annotation_process.seq_family_gene_ids_separator = matches
                .value_of("seq-family-gene-ids-separator")
                .unwrap()
                .trim()
                .to_string();
        }
        parse_seq_families_file(seq_families, &mut annotation_process);
        if annotation_process.verbose {
            println!(
                "Loaded {:?} sequence families from {:?}",
                &annotation_process.seq_families.len(),
                &seq_families
            );
        }
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
    // For each of the above to be parsed SSSR tables set their column mappings, if given by the
    // user:
    if matches.is_present("header") {
        for header_arg in matches.values_of("header").unwrap() {
            annotation_process.add_ssst_columns(header_arg);
        }
    }
    // For each of the above to be parsed SSSR tables set their their respective field-separator,
    // if given by the user:
    if matches.is_present("field-separator") {
        for field_separator in matches.values_of("field-separator").unwrap() {
            annotation_process.add_ssst_field_separator(field_separator);
        }
    }

    // Execute the Annotation-Process:
    annotation_process = run(annotation_process);

    // Save output:
    if let Some(o) = matches.value_of("output") {
        match output_writer::write_output_table(
            o.to_string(),
            annotation_process.human_readable_descriptions,
        ) {
            Ok(()) => {
                if annotation_process.verbose {
                    println!("output written to file {:?}.", o);
                }
            }
            Err(e) => eprintln!(
                "We are sorry, an error occurred when attempting to write output to file {:?} \n{:?}",
                o, e
            ),
        };
    }
}
