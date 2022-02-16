use super::default::{
    BLACKLIST_STITLE_REGEXS, CAPTURE_REPLACE_DESCRIPTION_PAIRS, FILTER_REGEXS,
    NON_INFORMATIVE_WORDS_REGEXS, SEQ_SIM_TABLE_COLUMNS, SPLIT_DESCRIPTION_REGEX,
    SPLIT_GENE_FAMILY_GENES_REGEX, SPLIT_GENE_FAMILY_ID_FROM_GENE_SET, SSSR_TABLE_FIELD_SEPARATOR,
};
use super::model_funcs::{parse_regex_file, parse_regex_replace_tuple_file};
use super::query::Query;
use super::seq_family::SeqFamily;
use super::seq_sim_table_reader::parse_table;
use num_cpus;
use rayon::prelude::*;
use regex::Regex;
use std::collections::HashMap;
use std::sync::{mpsc, Arc, Mutex};
use std::thread;

/// An instance of AnnotationProcess represents exactly what its name suggest, the assignment of
/// human readable descriptions, i.e. the annotation of queries or sets of these (biological
/// sequence families) with short and concise textual descriptions.
#[derive(Debug, Clone)]
pub struct AnnotationProcess {
    /// The valid file paths to tabular sequence similarity search results, the input.
    pub seq_sim_search_tables: Vec<String>,
    /// The order of columns (`qacc`, `sacc`, and `stitle`) in the above `seq_sim_search_tables`:
    pub ssst_columns: Vec<HashMap<String, usize>>,
    /// The field-separators used in the above `seq_sim_search_tables`:
    pub ssst_field_separators: Vec<char>,
    /// For each sequence similarity search result table (ssst) the list of regular expressions
    /// used to identify to be discarded descriptions (`stitle`) - note that the vector-index is
    /// used to pair input ssst with its blacklist regexs.
    pub ssst_blacklist_regexs: Vec<Vec<Regex>>,
    /// For each sequence similarity search result table (ssst) the list of regular expressions
    /// used to delete, i.e. filter out, matching sub-strings from (`stitle`) - note that the
    /// vector-index is used to pair input ssst with its filter regexs.
    pub ssst_filter_regexs: Vec<Vec<Regex>>,
    /// For each sequence similarity search result table (ssst) the list of
    /// "capture-replace-pairs", tuples of regular expressions and replace strings, is held here.
    /// These pairs are used to transform matching sub-strings from Blast Hit descriptions
    /// (`stitle`) - note that the vector-index is used to pair input ssst with its filter regexs.
    pub ssst_capture_replace_pairs: Vec<Vec<(Regex, String)>>,
    /// The in memory database of parsed sequence similarity search results in terms of Queries
    /// with their respective Hits.
    pub queries: HashMap<String, Query>,
    /// The in memory database of biological sequence families, i.e. sets of query identifiers, to
    /// be annotated with human readable descriptions. Keys are the families identifier and values
    /// are the SeqFamily instances.
    pub seq_families: HashMap<String, SeqFamily>,
    /// This string separates the gene-family-identifier (name) from the gene-identifier list
    /// that family comprises.
    pub seq_family_id_genes_separator: String,
    /// A regular expression (Rust syntax) represented as String to satisfy the `Default` trait.
    /// This regex is used to split the list of gene-identifiers in the gene families file.
    pub seq_family_gene_ids_separator: String,
    /// An in memory index from Query identifier to SeqFamily identifier:
    pub query_id_to_seq_family_id_index: HashMap<String, String>,
    /// A regular expression used to split descriptions (`stitle` in Blast terminology) into words.
    pub description_split_regex: Regex,
    /// The path to the optional argument file holding regular expression, one per line, used to
    /// recognize non informative words. If not given, the
    /// `default::BLACKLIST_DESCRIPTION_WORDS_REGEXS` is used.
    pub non_informative_words_regexs: Vec<Regex>,
    /// The human readable descriptions (HRDs) generated for the queries, i.e. either single query
    /// sequences or families (sets of query sequences). Stored here using the query identifier as
    /// key and the generated HRD as values.
    pub human_readable_descriptions: HashMap<String, String>,
    /// The number of parallel threads to use.
    pub n_threads: usize,
    /// In mode FamilyAnnotation also annotate lonely queries, i.e. queries not comprised in a
    /// sequence family?
    pub annotate_lonely_queries: bool,
    /// Does the user want informative messages about the annotation process printed out?
    pub verbose: bool,
}

/// Representation of the mode an instance of AnnotationProcess runs in. Can be either (i)
/// annotation of single biological query sequences `SequenceAnnotation`, or (ii) annotation of
/// sets of such query sequences `FamilyAnnotation`. Annotation means the generation of human
/// readable descriptions for either (i) single queries, or (ii) whole sets of biological
/// sequences.
#[derive(Debug, Clone)]
pub enum AnnotationProcessMode {
    SequenceAnnotation,
    FamilyAnnotation,
}

/// The central function that runs an annotation process. Note that because this function uses an
/// `Arc` `Mutex` reference to the argument `annotation_process`, this could not be implemented as
/// an instance function, i.e. using the argument `&mut self`, because then the Rust compiler
/// required the reference to have the `static` lifetime. So, this functions takes ownership of the
/// argument `annotation_process` and returns it in a modified form.
///
/// # Arguments
///
/// * `annotation_process` - The instance of `AnnotationProcess` to run.
pub fn run(mut annotation_process: AnnotationProcess) -> AnnotationProcess {
    // Are we printing information verbosely? (Note that by copying this boolean, we avoid
    // running into problems with the borrow-checker in the threads' println! statement:
    let verbose = annotation_process.verbose;

    // Validate input; if invalid panic! with a comprehensive error message:
    annotation_process.validate_fields();

    // If there are more input tables than the annotation_process.n_threads, only use n_threads
    // parallel processes.
    let n = if annotation_process.seq_sim_search_tables.len() <= annotation_process.n_threads {
        annotation_process.seq_sim_search_tables.len()
    } else {
        annotation_process.n_threads
    };

    // Setup communication between threads:
    let (tx, rx) = mpsc::channel();

    // Enable the threads to access the input sequence similarity search result tables, including
    // their index (position) in the original command line argument call:
    let sssts_mutex = Arc::new(Mutex::new(
        annotation_process
            .seq_sim_search_tables
            .iter()
            .map(|x| x.clone())
            .enumerate()
            .collect::<Vec<(usize, String)>>(),
    ));

    // Enable the threads to access the input column order in each input sequence similarity search
    // result table:
    let ssst_cols_mutex = Arc::new(Mutex::new(annotation_process.ssst_columns.clone()));

    // Enable the threads to access the input blacklist-regex-lists used to identify to be
    // discarded Blast Hit descriptions (`stitle`) in each input sequence similarity search result
    // table:
    let ssst_blacklist_regexs_mutex =
        Arc::new(Mutex::new(annotation_process.ssst_blacklist_regexs.clone()));

    // Enable the threads to access the input filter-regex-lists used to identify to be deleted
    // sub-strings in the parsed Blast Hit descriptions (`stitle`) in each input sequence
    // similarity search result table:
    let ssst_filter_regexs_mutex =
        Arc::new(Mutex::new(annotation_process.ssst_filter_regexs.clone()));

    // Enable the threads to access the input filter-regex-lists used to identify to be deleted
    // sub-strings in the parsed Blast Hit descriptions (`stitle`) in each input sequence
    // similarity search result table:
    let ssst_capture_replace_pairs_mutex = Arc::new(Mutex::new(
        annotation_process.ssst_capture_replace_pairs.clone(),
    ));

    // Enable the threads to access the input field separator used in the respective input
    // sequence similarity search result tables:
    let ssst_field_seps_mutex =
        Arc::new(Mutex::new(annotation_process.ssst_field_separators.clone()));

    // Prepare `n` threads for sequence similarity search parsing, each thread will parse a table
    // not yet processed until no tables are left to be processed:
    for _ in 0..n {
        let tx_i = tx.clone();

        // Start this sss_tbl's dedicated threat -
        // ... prepare thread local variables:
        let sssts_mutex_clone = sssts_mutex.clone();
        let ssst_cols_mutex_clone = ssst_cols_mutex.clone();
        let ssst_blacklist_regexs_mutex_clone = ssst_blacklist_regexs_mutex.clone();
        let ssst_filter_regexs_mutex_clone = ssst_filter_regexs_mutex.clone();
        let ssst_capture_replace_pairs_mutex_clone = ssst_capture_replace_pairs_mutex.clone();
        let ssst_field_seps_mutex_clone = ssst_field_seps_mutex.clone();

        // ... start the thread:
        thread::spawn(move || {
            // Field-Separator in Sequence Similarity Search (Blast) Result rows (lines):
            let mut field_separator = *SSSR_TABLE_FIELD_SEPARATOR;
            // Sequence Similarity Search (Blast) Result column indices:
            let mut qacc_col: usize = (*SEQ_SIM_TABLE_COLUMNS).get("qacc").unwrap().clone();
            let mut sacc_col: usize = (*SEQ_SIM_TABLE_COLUMNS).get("sacc").unwrap().clone();
            let mut stitle_col: usize = (*SEQ_SIM_TABLE_COLUMNS).get("stitle").unwrap().clone();

            loop {
                let mut ssst = sssts_mutex_clone.lock().unwrap();

                // Stop, if all input sequence similarity search tables have been parsed
                // already:
                if ssst.is_empty() {
                    break;
                }

                // Get the current input sequence similarity search table and its index:
                let (i, sss_tbl) = ssst.pop().unwrap();
                // Free the lock, so other threads may access `ssst_arc_mutex`:
                drop(ssst);

                // Did the user provide values for the column positions in the argument `sss_tbl`?
                let ssst_columns = ssst_cols_mutex_clone.lock().unwrap();
                if !ssst_columns.is_empty() {
                    let ssst_cols_i = &ssst_columns[i];
                    qacc_col = ssst_cols_i.get("qacc").unwrap().clone();
                    sacc_col = ssst_cols_i.get("sacc").unwrap().clone();
                    stitle_col = ssst_cols_i.get("stitle").unwrap().clone();
                }
                // Enable other threads to access `annotation_process.ssst_columns`:
                drop(ssst_columns);

                // Did the user provide values for the blacklist regex list to be used to identify
                // to be discarded Blast Hit descriptions (`stitle`) in the argument `sss_tbl`?
                let ssst_blacklist_regexs = ssst_blacklist_regexs_mutex_clone.lock().unwrap();
                let mut blacklist_regexs_i = (*BLACKLIST_STITLE_REGEXS).clone();
                if !ssst_blacklist_regexs.is_empty() {
                    blacklist_regexs_i = ssst_blacklist_regexs[i].clone();
                }
                // Enable other threads to access `annotation_process.ssst_blacklist_regexs`:
                drop(ssst_blacklist_regexs);

                // Did the user provide values for the filter regex list to be used to identify to
                // be deleted sub-strings in the Blast Hit descriptions (`stitle`) in the argument
                // `sss_tbl`?
                let ssst_filter_regexs = ssst_filter_regexs_mutex_clone.lock().unwrap();
                let mut filter_regexs_i = (*FILTER_REGEXS).clone();
                if !ssst_filter_regexs.is_empty() {
                    filter_regexs_i = ssst_filter_regexs[i].clone();
                }
                // Enable other threads to access `annotation_process.ssst_filter_regexs`:
                drop(ssst_filter_regexs);

                // Did the user provide values for the capture-replace-pairs to be used to identify
                // to be transformed sub-strings in the Blast Hit descriptions (`stitle`) in the
                // argument `sss_tbl`?
                let ssst_capture_replace_pairs =
                    ssst_capture_replace_pairs_mutex_clone.lock().unwrap();
                let mut capture_replace_pairs_i = (*CAPTURE_REPLACE_DESCRIPTION_PAIRS).clone();
                if !ssst_capture_replace_pairs.is_empty() {
                    capture_replace_pairs_i = ssst_capture_replace_pairs[i].clone();
                }
                // Enable other threads to access `annotation_process.ssst_capture_replace_pairs`:
                drop(ssst_capture_replace_pairs);

                // Did the user provide a custom field-separator for the argument `sss_tbl`?
                let ssst_field_separators = ssst_field_seps_mutex_clone.lock().unwrap();
                if !ssst_field_separators.is_empty() {
                    field_separator = ssst_field_separators[i];
                }
                // Enable other threads to access `annotation_process.ssst_field_separators`:
                drop(ssst_field_separators);

                parse_table(
                    &sss_tbl,
                    &field_separator,
                    &qacc_col,
                    &sacc_col,
                    &stitle_col,
                    &blacklist_regexs_i,
                    &filter_regexs_i,
                    Some(&capture_replace_pairs_i),
                    // Because we are in a `loop` we need to clone the cloned sender:
                    tx_i.clone(),
                );

                // Inform user, if requested:
                if verbose {
                    println!("Finished parsing {:?}", &sss_tbl);
                }
            }
        });
    }
    // Because of the above for loop tx needs to be cloned into tx_i's. tx needs to be dropped,
    // otherwise the below receiver loop will wait forever for tx to send some messages.
    drop(tx);

    // Process messages sent by the above threads. Note that this might trigger the annotation of
    // some queries or sequence families, if their data has been parsed completely:
    for (qacc, query) in rx {
        annotation_process.insert_query(qacc, query);
    }

    // Make sure all queries or sequence families are annotated:
    annotation_process.process_rest_data();

    // Return the modified `annotation_process`:
    annotation_process
}

impl AnnotationProcess {
    /// Creates a default instance of struct AnnotationProcess and returns it.
    pub fn new() -> AnnotationProcess {
        let nt = if num_cpus::get() < 2 {
            2
        } else {
            num_cpus::get()
        };
        AnnotationProcess {
            seq_sim_search_tables: vec![],
            ssst_columns: vec![],
            ssst_blacklist_regexs: vec![],
            ssst_filter_regexs: vec![],
            ssst_capture_replace_pairs: vec![],
            ssst_field_separators: vec![],
            queries: HashMap::new(),
            seq_families: HashMap::new(),
            seq_family_id_genes_separator: (*SPLIT_GENE_FAMILY_ID_FROM_GENE_SET).to_string(),
            seq_family_gene_ids_separator: (*SPLIT_GENE_FAMILY_GENES_REGEX).to_string(),
            description_split_regex: (*SPLIT_DESCRIPTION_REGEX).clone(),
            non_informative_words_regexs: (*NON_INFORMATIVE_WORDS_REGEXS).clone(),
            query_id_to_seq_family_id_index: HashMap::new(),
            human_readable_descriptions: HashMap::new(),
            n_threads: nt,
            annotate_lonely_queries: false,
            verbose: false,
        }
    }

    /// Processes the sequence similarity search result (SSSR) data parsed for the argument
    /// `query`. Adds the data to existing one, of data already has been parsed for the argument
    /// `query` from a different SSSR file or inserts the new data into the in-memory database
    /// `self.queries`. If for the argument query all input SSSR tables
    /// (`self.seq_sim_search_tables`) have produced data, this data will be processed by invoking
    /// `self.process_query_data_complete`.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A mutable reference to the current instance of AnnotationProcess, which
    ///                 serves as an in memory database into which to insert the parsed query.
    /// * `qacc: String` - The identifier of the argument query, i.e. the to be key in
    ///                    self.queries.
    /// * `query: Query` - A reference to the query to be inserted into the in memory database.
    pub fn insert_query(&mut self, qacc: String, query: Query) {
        // panic! if query.id already in results, this means the input SSSR files were not sorted
        // by query identifiers (`qacc` in Blast terminology):
        if self.human_readable_descriptions.contains_key(&qacc) {
            panic!( "Found an unexpected occurrence of query {:?} while parsing input files. Make sure your sequence similarity search result tables are sorted by query identifiers, i.e. `qacc` in Blast terminology. Use GNU sort, e.g. `sort -k <qacc-col-no> <your-blast-out-table>`.", &qacc);
        }
        if !self.queries.contains_key(&qacc) {
            self.queries.insert(qacc.clone(), query);
        }
        let mut stored_query = self.queries.get_mut(&qacc).unwrap();
        stored_query.n_parsed_from_sssr_tables += 1;
        // Have all input SSSR files provided data for the argument `query`?
        if stored_query.n_parsed_from_sssr_tables == self.seq_sim_search_tables.len() as u16 {
            drop(stored_query);
            // If yes, then process the parsed data:
            self.process_query_data_complete(qacc);
        }
    }

    /// Inserts the argument `seq_family: SeqFamily` into this AnnotationProcess instance's
    /// `seq_families`, while also updating the in memory index of biological query sequence
    /// identifiers pointing to their respective sequence family (see
    /// `query_id_to_seq_family_id_index`).
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A mutable reference to the current instance of AnnotationProcess, which
    ///                 serves as an in memory database into which to insert the argument
    ///                 biological sequence family.
    pub fn insert_seq_family(&mut self, seq_family_id: String, seq_family: SeqFamily) {
        for query_id in &seq_family.query_ids {
            if self.query_id_to_seq_family_id_index.contains_key(query_id) {
                let other_family_id = self.query_id_to_seq_family_id_index.get(query_id).unwrap();
                if *other_family_id != seq_family_id {
                    panic!("Biological sequence {:?} already set as member of family {:?}. But found {:?} again declared as member of another family {:?}.\nMake sure each biological sequence appears in one and only one family to avoid this problem.", query_id, other_family_id, query_id, seq_family_id);
                }
            }
            self.query_id_to_seq_family_id_index
                .insert((*query_id).clone(), seq_family_id.clone());
        }
        self.seq_families.insert(seq_family_id, seq_family);
    }

    /// Informs and returns the mode an AnnotationProcess (argument `&self`) is running in, can be
    /// either (i) AnnotationProcessMode::SequenceAnnotation or (ii)
    /// AnnotationProcessMode::FamilyAnnotation.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A mutable reference to the current instance of AnnotationProcess, which
    ///                 serves as an in memory database into which to insert the parsed query.
    pub fn mode(&self) -> AnnotationProcessMode {
        if self.seq_families.len() > 0 {
            AnnotationProcessMode::FamilyAnnotation
        } else {
            AnnotationProcessMode::SequenceAnnotation
        }
    }

    /// Function generates a human readable description (HRD) for the argument `query_id`. The
    /// resulting HRD is stored in `self.human_readable_descriptions` and thus the query is marked
    /// as processed. In order to optimize memory footprint the query and all of its contained
    /// sequence similarity search result (Hits in Blast terminology) data is removed.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A mutable reference to the current instance of AnnotationProcess, which
    ///                 serves as an in memory database into which to insert the parsed query.
    /// * `query_id: String` - An instance of `String` representing the query identifier
    pub fn annotate_query(&mut self, query_id: String) {
        // Generate the desired result, i.e. a human readable description for the Query:
        let hrd = self.queries.get(&query_id).unwrap().annotate(
            &self.description_split_regex,
            &self.non_informative_words_regexs,
        );
        // Add the new result to the in memory database, i.e.
        // `self.human_readable_descriptions`:
        self.human_readable_descriptions
            .insert(query_id.clone(), hrd);
        // Free memory by removing the parsed input data, no longer required:
        self.queries.remove(&query_id);
    }

    /// Function generates a human readable description (HRD) for the argument `seq_family_id`. The
    /// resulting HRD is stored in `self.human_readable_descriptions` and thus the family is marked
    /// as processed. In order to optimize memory footprint the family and all of its contained
    /// query data is removed.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A mutable reference to the current instance of AnnotationProcess, which
    ///                 serves as an in memory database into which to insert the parsed query.
    /// * `seq_family_id: &String` - A reference to a `String` representing the biological sequence
    ///                              family's (`SeqFamily`) identifier.
    pub fn annotate_seq_family(&mut self, seq_family_id: &String) {
        // Generate the desired result, i.e. a human readable description for the SeqFamily:
        let seq_family = self.seq_families.get(seq_family_id).unwrap();
        let hrd = seq_family.annotate(
            &self.queries,
            &self.description_split_regex,
            &self.non_informative_words_regexs,
        );
        // Add the new result to the in memory database, i.e.
        // `self.human_readable_descriptions`:
        self.human_readable_descriptions
            .insert((*seq_family_id).clone(), hrd);
        // need to clone, otherwise had problems with the compiler (E0599):
        let query_ids = seq_family.query_ids.clone();
        // Free memory by removing the parsed input data, no longer required:
        for query_id in query_ids.iter() {
            self.queries.remove(query_id);
            self.query_id_to_seq_family_id_index.remove(query_id);
            self.seq_families.remove(seq_family_id);
        }
    }

    /// Invoked whenever a query instance has been supplied with results from _all_ sequence
    /// similarity search result (SSSR) files, implying that for that particular instance of
    /// `Query` no more SSSR results (Hits in Blast terminology) can be parsed. Thus, that query
    /// can be processed and a human readable description can be generated for it. If this instance
    /// of `AnnotationProcess` (`self`) is run in `AnnotationProcessMode::FamilyAnnotation` a
    /// similar approach is triggered for the SeqFamily that contains the argument `query_id`. If
    /// that family has SSSR result data for _all_ of its contained queries, the family will be
    /// processed and annotated.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A mutable reference to the current instance of AnnotationProcess, which
    ///                 serves as an in memory database into which to insert the parsed query.
    pub fn process_query_data_complete(&mut self, query_id: String) {
        let mode = self.mode();
        match mode {
            // Handle annotation of single biological sequences:
            AnnotationProcessMode::SequenceAnnotation => {
                self.annotate_query(query_id);
            }
            // Handle annotation of sets of biological sequences, so called "Gene Families":
            AnnotationProcessMode::FamilyAnnotation => {
                // Get SeqFamily for query sequence identifier (`query_id`):
                if self.query_id_to_seq_family_id_index.contains_key(&query_id) {
                    let seq_fam_id = self
                        .query_id_to_seq_family_id_index
                        .get(&query_id)
                        .unwrap()
                        .clone();
                    if self.seq_families.contains_key(&seq_fam_id) {
                        // Tell Family that parsing of Blast results for argument `query_id` has
                        // been completed:
                        let seq_fam = self.seq_families.get_mut(&seq_fam_id).unwrap();
                        seq_fam.mark_query_id_with_complete_data(&query_id);
                        // Ask SeqFamily if all queries have complete data:
                        if seq_fam.all_query_data_complete() {
                            self.annotate_seq_family(&seq_fam_id);
                        }
                    }
                } else if self.annotate_lonely_queries {
                    // If no family for qs_id can be found, annotate query as in
                    // SequenceAnnotation:
                    self.annotate_query(query_id);
                }
            }
        }
    }

    /// Invoked after all parsing of input sequence similarity search result (SSSR) files has
    /// finished to annotate those queries or biological sequence families that have not yet been
    /// annotated. These are those that do not have results in each separate SSSR file.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A mutable reference to the current instance of AnnotationProcess, which
    ///                 serves as an in memory database into which to insert the parsed query.
    pub fn process_rest_data(&mut self) {
        // Note that below `par_iter` is used to process the data _in parallel_. To make this work
        // the parallel processes need to be independent and cannot write results of annotations
        // (HRDs) into the current instance of AnnotationProcess without using something like an
        // Mutex. Thus results are collected in terms of tuples containing the annotee identifier
        // and the generated human readable description.
        let mode = self.mode();
        let hrd_tuples: Vec<(String, String)>;
        match mode {
            // Handle annotation of single biological sequences:
            AnnotationProcessMode::SequenceAnnotation => {
                // Process queries that might have gotten parsed results only from a subset of the input
                // sequence similarity search result (SSSR) files:
                hrd_tuples = self
                    .queries
                    .keys()
                    .cloned()
                    .collect::<Vec<String>>()
                    .par_iter()
                    .map(|query_id| {
                        let query = self.queries.get(query_id).unwrap();
                        let hrd = query.annotate(
                            &self.description_split_regex,
                            &self.non_informative_words_regexs,
                        );
                        ((*query_id).to_string(), hrd)
                    })
                    .collect();
            }
            // Handle annotation of sets of biological sequences, so called "Gene Families":
            AnnotationProcessMode::FamilyAnnotation => {
                // Process seq families that might have queries that got no blast hits in some
                // input blast tables:
                hrd_tuples = self
                    .seq_families
                    .keys()
                    .cloned()
                    .collect::<Vec<String>>()
                    .par_iter()
                    .map(|seq_fam_id| {
                        let seq_fam = self.seq_families.get(seq_fam_id).unwrap();
                        let hrd = seq_fam.annotate(
                            &self.queries,
                            &self.description_split_regex,
                            &self.non_informative_words_regexs,
                        );
                        ((*seq_fam_id).to_string(), hrd)
                    })
                    .collect();
            }
        }

        // Free memory:
        self.queries = Default::default();
        self.seq_families = Default::default();
        self.query_id_to_seq_family_id_index = Default::default();

        // Set the human readable descriptions generated in parallel:
        for i_tpl in hrd_tuples {
            self.human_readable_descriptions.insert(i_tpl.0, i_tpl.1);
        }
    }

    /// Parses the command line argument `header` into a HashMap<String, usize> in which the
    /// sequence similarity search result (Blast) table (SSST) column names are mapped to their
    /// respective position in the to be parsed SSST. Inserts the parsed HashMap into
    /// `self.ssst_columns`, and uses the `default::SEQ_SIM_TABLE_COLUMNS` HashMap if the argument
    /// `header_arg` equals `"default"` (case insensitive).
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A reference to a mutable instance of AnnotationProcess.
    /// * `header_arg: &str` - The passed header argument
    pub fn add_ssst_columns(&mut self, header_arg: &str) {
        let mut seq_sim_table_cols: HashMap<String, usize>;
        if header_arg.trim().to_lowercase() == "default" {
            seq_sim_table_cols = (*SEQ_SIM_TABLE_COLUMNS).clone();
        } else {
            seq_sim_table_cols = HashMap::new();
            for (i, col_name) in header_arg
                .trim()
                .split(" ")
                .filter(|x| !x.is_empty())
                .enumerate()
            {
                seq_sim_table_cols.insert(col_name.to_string(), i as usize);
            }
        }
        self.ssst_columns.push(seq_sim_table_cols);
    }

    /// Parses one (of potentially many) command line argument `blacklist-regexs` into a
    /// `Vec<Regex>` in which the regular expressions are stored used to identify to be discarded
    /// `stitle`. These `stitle` strings are parsed from the input sequence similarity search
    /// result (Blast) table (SSST). Inserts the parsed `Vec<Regex>` into
    /// `self.ssst_blacklist_regexs`, and uses the `default::BLACKLIST_STITLE_REGEXS` `Vec<Regex>`
    /// if the argument `blacklist_regexs_arg` equals `"default"` (case insensitive).
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A reference to a mutable instance of AnnotationProcess.
    /// * `blacklist_regexs_arg: &str` - The passed command line argument
    pub fn add_ssst_blacklist_regexs(&mut self, blacklist_regexs_arg: &str) {
        let blacklist_regexs: Vec<Regex> =
            if blacklist_regexs_arg.trim().to_lowercase() == "default" {
                (*BLACKLIST_STITLE_REGEXS).clone()
            } else {
                parse_regex_file(blacklist_regexs_arg)
            };
        self.ssst_blacklist_regexs.push(blacklist_regexs);
    }

    /// Parses one (of potentially many) command line argument `filter-regexs` into a `Vec<Regex>`
    /// in which the regular expressions are stored used to identify to be deleted, filtered out
    /// matching sub-strings in the Blast Hit descriptions (`stitle`). These `stitle` strings are
    /// parsed from the input sequence similarity search result (Blast) table (SSST). Inserts the
    /// parsed `Vec<Regex>` into `self.ssst_filter_regexs`, and uses the `default::FILTER_REGEXS`
    /// `Vec<Regex>` if the argument `filter_regexs_arg` equals `"default"` (case insensitive).
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A reference to a mutable instance of AnnotationProcess.
    /// * `filter_regexs_arg: &str` - The passed command line argument
    pub fn add_ssst_filter_regexs(&mut self, filter_regexs_arg: &str) {
        let filter_regexs: Vec<Regex> = if filter_regexs_arg.trim().to_lowercase() == "default" {
            (*FILTER_REGEXS).clone()
        } else {
            parse_regex_file(filter_regexs_arg)
        };
        self.ssst_filter_regexs.push(filter_regexs);
    }

    /// Parses one (of potentially many) command line argument `capture-replace-pairs` into a
    /// `Vec<(Regex,String)>` in which the regular expressions and the replacements are stored used
    /// to identify to be transformed matching sub-strings in the Blast Hit descriptions
    /// (`stitle`). This transformation is done using the standard regular expression `replace`
    /// method. Note that the `stitle` strings (descriptions) are parsed from the input sequence
    /// similarity search result (Blast) table (SSST). This function inserts the parsed
    /// `Vec<(Regex,String)>` into `self.ssst_capture_replace_pairs`, and uses the
    /// `default::REPLACE_REGEXS_DESCRIPTION` `Vec<(Regex, String)>` if the argument
    /// `capture_replace_pairs_arg` equals `"default"` (case insensitive).
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A reference to a mutable instance of AnnotationProcess.
    /// * `capture_replace_pairs_arg: &str` - The passed command line argument
    pub fn add_ssst_capture_replace_pairs(&mut self, capture_replace_pairs_arg: &str) {
        let capture_replace_pairs: Vec<(Regex, String)> =
            if capture_replace_pairs_arg.trim().to_lowercase() == "default" {
                (*CAPTURE_REPLACE_DESCRIPTION_PAIRS).clone()
            } else {
                parse_regex_replace_tuple_file(capture_replace_pairs_arg)
            };
        self.ssst_capture_replace_pairs.push(capture_replace_pairs);
    }

    /// Parses the command line argument `field-separator` into a `char` used to split a line (row)
    /// in a sequence similarity search result table into fields, i.e. a Blast Hit record. If the
    /// argument `field_separator_arg` equals `"default"` (case insensitive) the value of
    /// `default::SSSR_TABLE_FIELD_SEPARATOR` is used.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A reference to a mutable instance of AnnotationProcess.
    /// * `field_separator_arg: &str` - The passed field_separator argument
    pub fn add_ssst_field_separator(&mut self, field_separator_arg: &str) {
        let mut seq_sim_table_field_separator = *SSSR_TABLE_FIELD_SEPARATOR;
        if field_separator_arg.trim().to_lowercase() != "default" {
            seq_sim_table_field_separator = field_separator_arg.chars().next().unwrap();
        }
        self.ssst_field_separators
            .push(seq_sim_table_field_separator);
    }

    /// Function validates the AnnotationProcess's fields and checks whether they are valid and
    /// complete to start `run`. If invalid the function panics! with a comprehensive error
    /// message.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A reference to a mutable instance of AnnotationProcess.
    pub fn validate_fields(&mut self) {
        let n_ssst = self.seq_sim_search_tables.len();

        // --header
        if !self.ssst_columns.is_empty() && self.ssst_columns.len() != n_ssst {
            let n_ssst_cols = self.ssst_columns.len();
            panic!("Cannot run Annotation-Process, because got {} sequence similarity search result tables (SSSTs), but only {} column mappings. Please provide either no column mappings, causing the default to be used for all SSSTs, or provide one --field_separator (-e) argument for each of your input SSSTs. See --help for more details.", n_ssst_cols, n_ssst);
        }
        if !self.ssst_columns.is_empty() {
            let required_cols = vec!["qacc", "sacc", "stitle"];
            for (i, ssst_cols_i) in self.ssst_columns.iter().enumerate() {
                for col_i in &required_cols {
                    if !ssst_cols_i.contains_key(&col_i.to_string()) {
                        panic!(
                            "Cannot run Annotation-Process, because --header (-e) argument number {} does not contain required column {:?}!",
                            i + 1, col_i
                        );
                    }
                }
            }
        }

        // --blacklist-regexs
        if !self.ssst_blacklist_regexs.is_empty() && self.ssst_blacklist_regexs.len() != n_ssst {
            let n_blacklist_regexs = self.ssst_blacklist_regexs.len();
            panic!("Cannot run Annotation-Process, because got {} sequence similarity search result tables (SSSTs), but {} --blacklist-regexs (-b). Please provide either no --blacklist-regexs, causing the default to be used for all SSSTs, or provide one --blacklist-regexs (-b) argument for each of your input SSSTs. See --help for more details.", n_ssst, n_blacklist_regexs);
        }

        // --filter-regexs
        if !self.ssst_filter_regexs.is_empty() && self.ssst_filter_regexs.len() != n_ssst {
            let n_ssst_filter_regexs = self.ssst_filter_regexs.len();
            panic!("Cannot run Annotation-Process, because got {} sequence similarity search result tables (SSSTs), but {} --filter-regexs (-l). Please provide either no --filter-regexs, causing the default to be used for all SSSTs, or provide one --filter-regexs (-l) argument for each of your input SSSTs. See --help for more details.", n_ssst, n_ssst_filter_regexs);
        }

        // --capture-replace-pairs
        if !self.ssst_capture_replace_pairs.is_empty()
            && self.ssst_capture_replace_pairs.len() != n_ssst
        {
            let n_ssst_capture_replace_pairs = self.ssst_capture_replace_pairs.len();
            panic!("Cannot run Annotation-Process, because got {} sequence similarity search result tables (SSSTs), but {} --capture-replace-pairs (-c). Please provide either no --capture-replace-pairs, causing the default to be used for all SSSTs, or provide one --capture-replace-pairs (-c) argument for each of your input SSSTs. See --help for more details.", n_ssst, n_ssst_capture_replace_pairs);
        }

        // --field-separator
        if !self.ssst_field_separators.is_empty() && self.ssst_field_separators.len() != n_ssst {
            let n_ssst_field_seps = self.ssst_field_separators.len();
            panic!("Cannot run Annotation-Process, because got {} sequence similarity search result tables (SSSTs), but {} field-separators. Please provide either no field-separators, causing the default to be used for all SSSTs, or provide one --field-separator (-p) argument for each of your input SSSTs. See --help for more details.", n_ssst, n_ssst_field_seps);
        }

        // --n-threads
        if self.n_threads < 2 {
            panic!("Cannot run Annotation-Process, because option '--n-threads' ('-n') must at least be minimum of two (2)!");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn new_annotation_process_initializes_fields() {
        let ap = AnnotationProcess::new();
        assert_eq!(ap.seq_sim_search_tables.len(), 0);
        assert_eq!(ap.seq_families.len(), 0);
    }

    #[test]
    fn annotation_process_mode_detected_correctly() {
        let mut ap = AnnotationProcess::new();
        assert!(matches!(
            ap.mode(),
            AnnotationProcessMode::SequenceAnnotation
        ));
        // meaningless empty SeqFamily, but none the less...
        ap.seq_families
            .insert("Family1".to_string(), SeqFamily::new());
        assert!(matches!(ap.mode(), AnnotationProcessMode::FamilyAnnotation));
    }

    #[test]
    fn insert_query_works() {
        let mut ap = AnnotationProcess::new();
        // let mut nq1 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());
        let mut nq1 = Query::new();
        let h1 = ("hit_One","sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1");
        let h2 = ("hit_Two","sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1");
        nq1.hits.insert(h1.0.to_string(), h1.1.to_string());
        nq1.hits.insert(h2.0.to_string(), h2.1.to_string());
        // Test insert_query
        let qacc = "Soltu.DM.02G015700.1".to_string();
        ap.insert_query(qacc.clone(), nq1);

        // check if query got inserted correctly
        assert!(ap.queries.contains_key(&qacc));
        assert_eq!(ap.queries.get(&qacc).unwrap().hits.len(), 2);

        // check if the hits for the query are present
        let sq = ap.queries.get(&qacc).unwrap();
        assert!(sq.hits.contains_key(h1.0));
        assert!(sq.hits.contains_key(h2.0));
    }

    #[test]
    #[should_panic]
    fn insert_query_panics_in_case_of_unsorted_blast_table() {
        let mut ap = AnnotationProcess::new();
        let nq1 = Query::new();
        let qacc = "Soltu.DM.02G015700.1".to_string();

        // Mark nq1 as already processed:
        ap.human_readable_descriptions
            .insert(qacc.clone(), "Unknown protein".to_string());
        // Should panic:
        ap.insert_query(qacc, nq1);
    }

    #[test]
    fn insert_family_works() {
        let mut ap = AnnotationProcess::new();
        let mut sf1 = SeqFamily::new();
        sf1.query_ids = vec![
            "Query1".to_string(),
            "Query2".to_string(),
            "Query3".to_string(),
        ];
        let sf_id1 = "SeqFamily1".to_string();
        ap.insert_seq_family(sf_id1.clone(), sf1);
        assert!(ap.seq_families.contains_key("SeqFamily1"));
        assert_eq!(
            *ap.query_id_to_seq_family_id_index.get("Query1").unwrap(),
            sf_id1
        );
        assert_eq!(
            *ap.query_id_to_seq_family_id_index.get("Query2").unwrap(),
            sf_id1
        );
        assert_eq!(
            *ap.query_id_to_seq_family_id_index.get("Query3").unwrap(),
            sf_id1
        );
        let mut sf2 = SeqFamily::new();
        sf2.query_ids = vec![
            "Query4".to_string(),
            "Query5".to_string(),
            "Query6".to_string(),
        ];
        let sf_id2 = "SeqFamily2".to_string();
        ap.insert_seq_family(sf_id2.clone(), sf2);
        assert!(ap.seq_families.contains_key("SeqFamily2"));
        assert_eq!(
            *ap.query_id_to_seq_family_id_index.get("Query4").unwrap(),
            sf_id2
        );
        assert_eq!(
            *ap.query_id_to_seq_family_id_index.get("Query5").unwrap(),
            sf_id2
        );
        assert_eq!(
            *ap.query_id_to_seq_family_id_index.get("Query6").unwrap(),
            sf_id2
        );
    }

    #[test]
    #[should_panic]
    fn double_assignment_of_seq_id_to_different_families_causes_panic() {
        let mut ap = AnnotationProcess::new();
        let mut sf1 = SeqFamily::new();
        sf1.query_ids = vec![
            "Query1".to_string(),
            "Query2".to_string(),
            "Query3".to_string(),
        ];
        let sf_id1 = "SeqFamily1".to_string();
        ap.insert_seq_family(sf_id1.clone(), sf1);
        let mut sf2 = SeqFamily::new();
        sf2.query_ids = vec!["Query1".to_string(), "Query4".to_string()];
        let sf_id2 = "SeqFamily2".to_string();
        ap.insert_seq_family(sf_id2.clone(), sf2);
    }

    // This test also tests the functions
    // * annotate_query
    // * annotate_seq_family
    // implicitly
    #[test]
    fn process_query_data_complete_works() {
        // Test queries:
        let mut ap = AnnotationProcess::new();
        ap.seq_sim_search_tables = vec!["blast_out_table.txt".to_string()];
        // let mut nq1 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());
        let mut nq1 = Query::new();
        let qacc = "Soltu.DM.02G015700.1".to_string();
        ap.insert_query(qacc.clone(), nq1);
        // Query should have been annotated:
        assert!(ap.human_readable_descriptions.contains_key(&qacc));
        assert!(!ap.queries.contains_key(&qacc));
        // Test families:
        ap = AnnotationProcess::new();
        ap.seq_sim_search_tables = vec!["blast_out_table.txt".to_string()];
        let mut sf1 = SeqFamily::new();
        sf1.query_ids = vec!["Soltu.DM.02G015700.1".to_string()];
        let sf_id1 = "SeqFamily1".to_string();
        ap.insert_seq_family(sf_id1.clone(), sf1);
        nq1 = Query::new();
        ap.insert_query(qacc.clone(), nq1);
        // Family should have been annotated:
        assert!(ap.human_readable_descriptions.contains_key(&sf_id1));
        assert!(!ap.queries.contains_key(&qacc));
        assert!(!ap.seq_families.contains_key(&sf_id1));
        assert!(!ap.query_id_to_seq_family_id_index.contains_key(&qacc));
    }

    #[test]
    fn process_rest_data_works() {
        // Test Queries:
        let mut ap = AnnotationProcess::new();
        // let mut nq1 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());

        let mut nq1 = Query::new();
        let qacc = "Soltu.DM.02G015700.1".to_string();
        let h1 = ("hit_One","sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1");
        let h2 = ("hit_Two","sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1");
        nq1.hits.insert(h1.0.to_string(), h1.1.to_string());
        nq1.hits.insert(h2.0.to_string(), h2.1.to_string());

        ap.insert_query(qacc.clone(), nq1);
        ap.process_rest_data();
        // Query should have been annotated:
        assert!(ap.human_readable_descriptions.contains_key(&qacc));
        assert!(!ap.queries.contains_key(&qacc));
        // Test Families:
        ap = AnnotationProcess::new();
        let mut sf1 = SeqFamily::new();
        sf1.query_ids = vec![qacc.clone()];
        let sf_id1 = "SeqFamily1".to_string();
        ap.insert_seq_family(sf_id1.clone(), sf1);
        nq1 = Query::new();
        ap.insert_query(qacc.clone(), nq1);
        let mut sf2 = SeqFamily::new();
        sf2.query_ids = vec!["The protein without known relatives".to_string()];
        let sf_id2 = "SeqFamily2".to_string();
        ap.insert_seq_family(sf_id2.clone(), sf2);
        ap.process_rest_data();
        // Families should have been annotated:
        assert!(ap.human_readable_descriptions.contains_key(&sf_id1));
        assert!(!ap.queries.contains_key(&qacc));
        assert!(!ap.seq_families.contains_key(&sf_id1));
        assert!(!ap.query_id_to_seq_family_id_index.contains_key(&qacc));
        assert!(ap.human_readable_descriptions.contains_key(&sf_id2));
        assert!(!ap.seq_families.contains_key(&sf_id2));
    }

    #[test]
    fn run_annotates_queries() {
        let mut ap = AnnotationProcess::new();
        ap.seq_sim_search_tables.push(
            Path::new("misc")
                .join("Twelve_Proteins_vs_Swissprot_blastp.txt")
                .to_str()
                .unwrap()
                .to_string(),
        );
        ap.seq_sim_search_tables.push(
            Path::new("misc")
                .join("Twelve_Proteins_vs_trembl_blastp.txt")
                .to_str()
                .unwrap()
                .to_string(),
        );
        ap = run(ap);
        let hrds = ap.human_readable_descriptions;
        assert!(hrds.len() > 0);
        let queries_with_expected_result = vec![
            "Soltu.DM.01G022510.1".to_string(),
            "Soltu.DM.01G045390.1".to_string(),
            "Soltu.DM.02G015700.1".to_string(),
            "Soltu.DM.02G020600.1".to_string(),
            "Soltu.DM.03G011280.1".to_string(),
            "Soltu.DM.03G026010.1".to_string(),
            "Soltu.DM.04G035790.1".to_string(),
            "Soltu.DM.07G016620.1".to_string(),
            "Soltu.DM.09G022410.3".to_string(),
            "Soltu.DM.10G003150.1".to_string(),
            "Soltu.DM.S001650.1".to_string(),
        ];
        for qid in queries_with_expected_result {
            assert!(hrds.contains_key(&qid))
        }
        for (_, v) in hrds {
            assert!(!v.is_empty());
        }
    }

    #[test]
    fn run_annotates_families() {
        let mut ap = AnnotationProcess::new();
        ap.seq_sim_search_tables.push(
            Path::new("misc")
                .join("Twelve_Proteins_vs_Swissprot_blastp.txt")
                .to_str()
                .unwrap()
                .to_string(),
        );
        ap.seq_sim_search_tables.push(
            Path::new("misc")
                .join("Twelve_Proteins_vs_trembl_blastp.txt")
                .to_str()
                .unwrap()
                .to_string(),
        );
        let mut sf1 = SeqFamily::new();
        let sf1_id = "SeqFamily1".to_string();
        sf1.query_ids = vec![
            "Soltu.DM.01G022510.1".to_string(),
            "Soltu.DM.01G045390.1".to_string(),
            "Soltu.DM.02G015700.1".to_string(),
            "Soltu.DM.02G020600.1".to_string(),
            "Soltu.DM.03G011280.1".to_string(),
            "Soltu.DM.03G026010.1".to_string(),
            "Soltu.DM.04G035790.1".to_string(),
        ];
        let mut sf2 = SeqFamily::new();
        let sf2_id = "SeqFamily2".to_string();
        sf2.query_ids = vec![
            "Soltu.DM.07G016620.1".to_string(),
            "Soltu.DM.09G022410.3".to_string(),
            "Soltu.DM.10G003150.1".to_string(),
            "Soltu.DM.S001650.1".to_string(),
            "The_Protein_Without_Blast_hits".to_string(),
        ];
        ap.insert_seq_family(sf1_id.clone(), sf1);
        ap.insert_seq_family(sf2_id.clone(), sf2);
        ap = run(ap);
        let hrds = ap.human_readable_descriptions;
        assert_eq!(hrds.len(), 2);
        let queries_with_expected_result = vec![sf1_id, sf2_id];
        for qid in queries_with_expected_result {
            assert!(hrds.contains_key(&qid))
        }
        for (_, v) in hrds {
            assert!(!v.is_empty());
        }
    }
}
