use super::default::{SEQ_SIM_TABLE_COLUMNS, SSSR_TABLE_FIELD_SEPARATOR};
use super::query::Query;
use super::seq_family::SeqFamily;
use super::seq_sim_table_reader::parse_table;
use std::collections::HashMap;
use std::sync::mpsc::Sender;
use std::sync::{mpsc, Arc, Mutex};
use std::{mem, thread};

/// An instance of AnnotationProcess represents exactly what its name suggest, the assignment of
/// human readable descriptions, i.e. the annotation of queries or sets of these (biological
/// sequence families) with short and concise textual descriptions.
#[derive(Debug, Clone, Default)]
pub struct AnnotationProcess {
    /// The valid file paths to tabular sequence similarity search results, the input.
    pub seq_sim_search_tables: Vec<String>,
    /// The in memory database of parsed sequence similarity search results in terms of Queries
    /// with their respective Hits.
    pub queries: HashMap<String, Query>,
    /// The in memory database of biological sequence families, i.e. sets of query identifiers, to
    /// be annotated with human readable descriptions. Keys are the families identifier and values
    /// are the SeqFamily instances.
    pub seq_families: HashMap<String, SeqFamily>,
    /// An in memory index from Query identifier to SeqFamily identifier:
    pub query_id_to_seq_family_id_index: HashMap<String, String>,
    /// The human readable descriptions (HRDs) generated for the queries, i.e. either single query
    /// sequences or families (sets of query sequences). Stored here using the query identifier as
    /// key and the generated HRD as values.
    pub human_readable_descriptions: HashMap<String, String>,
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

/// The central function that runs an annotation process. Note that this is implemented as a
/// "static" function (in Java terminology - apologies offered), because using a reference to an
/// instance of `AnnotationProcess` would cause lifetime issues. This is, because a
/// Arc<Mutex<AnnotationProcess>> is constructed and shared between threads. Thus this function
/// that actually executes an instance of `AnnotationProcess` takes ownership of an argument
/// `annotation_process`. Returns the results of running the Annotation-Process in the form of a
/// `HashMap<String, String>`.
///
/// # Arguments
///
/// * `annotation_process: AnnotationProcess` - The instance of `AnnotationProcess` to run.
pub fn run(annotation_process: AnnotationProcess) -> AnnotationProcess {
    // Prepare:
    let sssr_tables: Vec<String> = annotation_process
        .seq_sim_search_tables
        .iter()
        .map(|x| x.clone())
        .collect();
    let (tx, rx) = mpsc::channel();
    let ap_arc_mutex: Arc<Mutex<AnnotationProcess>> = Arc::new(Mutex::new(annotation_process));

    // Parse and process each sequence similarity search result table in a dedicated
    // thread:
    for sss_tbl in sssr_tables {
        let tx_i = tx.clone();
        let ap_i = ap_arc_mutex.clone();

        // Start this sss_tbl's dedicated threat:
        thread::spawn(move || {
            parse_table(
                sss_tbl,
                *SSSR_TABLE_FIELD_SEPARATOR,
                &(*SEQ_SIM_TABLE_COLUMNS),
                ap_i,
                tx_i,
            );
        });
    }
    // In the unlikely case not a single messages is send through any tx_i, send one, so this does
    // not become an endless loop:
    let final_message = "PROT-SCRIBER THREADS HAVE FINISHED".to_string();
    tx.send(final_message.clone()).unwrap();
    drop(tx);

    // Process queries whose Blast result data has been parsed from all input tables:
    for received in rx {
        if received != final_message {
            let ap_i = ap_arc_mutex.clone();
            ap_i.lock().unwrap().process_query_data_complete(received);
        }
    }

    // How to take ownership of argument `annotation_process` back from having it moved into an
    // Arc<Mutex<AnnotationProcess>>; see:
    // https://stackoverflow.com/questions/29177449/how-to-take-ownership-of-t-from-arcmutext
    let mut ap = Arc::try_unwrap(ap_arc_mutex).unwrap().into_inner().unwrap();

    // Make sure all queries or sequence families are annotated:
    ap.process_rest_data();

    // Return modified version of input argument `annotation_process`:
    ap
}

impl AnnotationProcess {
    /// Creates an empty (`Default`) instance of struct AnnotationProcess.
    pub fn new() -> AnnotationProcess {
        Default::default()
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
    /// * `query: &mut Query` - A reference to the query to be inserted into the in memory
    ///                         database.
    /// * `transmitter: std::sync::mpsc::Sender<String>` - An instance of transmitter to tell the
    /// parent thread when a query has been completely parsed.
    pub fn insert_query(&mut self, query: &mut Query, transmitter: Sender<String>) {
        // panic! if query.id already in results, this means the input SSSR files were not sorted
        // by query identifiers (`qacc` in Blast terminology):
        if self.human_readable_descriptions.contains_key(&query.id) {
            panic!( "Found an unexpected occurrence of query {:?} while parsing input files. Make sure your sequence similarity search result tables are sorted by query identifiers, i.e. `qacc` in Blast terminology. Use GNU sort, e.g. `sort -k <qacc-col-no> <your-blast-out-table>`.", &query.id);
        }
        let stored_query: &mut Query;
        if self.queries.contains_key(&query.id) {
            // Insert sequence similarity search results (SSSR) in the format of a Query, where this
            // Query has been parsed from another input SSSR file before:
            stored_query = self.queries.get_mut(&query.id).unwrap();
            // Add the parsed results to the already existing data:
            stored_query.add_hits(query);
        } else {
            // Insert new Query that has NOT been parsed from another input SSSR file before:
            self.queries.insert(query.id.clone(), query.clone());
            stored_query = self.queries.get_mut(&query.id).unwrap();
        }
        stored_query.n_parsed_from_sssr_tables += 1;
        // Have all input SSSR files provided data for the argument `query`?
        if stored_query.n_parsed_from_sssr_tables == self.seq_sim_search_tables.len() as u16 {
            let sqi = stored_query.id.clone();
            // If yes, then tell the parent thread that all Blast result data for query `sqi` has
            // been completely parsed:
            transmitter.send(sqi).unwrap();
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
        let hrd = self.queries.get(&query_id).unwrap().annotate();
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
        let hrd = seq_family.annotate(&self.queries);
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
                } else {
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
        // Cannot get away without cloning the keys. Rust compiler complains about an immutable
        // borrow and at the same time a mutable borrow of `self` in the call of
        // process_rest_data or annotate_seq_family, respectively.

        // Process seq families that might have queries that got no blast hits in any input blast
        // table:
        for seq_family_id in self
            .seq_families
            .iter()
            .map(|(k, _v)| k.clone())
            .collect::<Vec<String>>()
            .iter()
        {
            self.annotate_seq_family(seq_family_id);
        }
        // Process queries that might have gotten parsed results only from a subset of the input
        // sequence similarity search result (SSSR) files:
        for query_id in self
            .queries
            .iter()
            .map(|(k, _v)| k.clone())
            .collect::<Vec<String>>()
            .iter()
        {
            self.process_query_data_complete(query_id.to_string());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hit::Hit;
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
        let (tx, rx) = mpsc::channel();
        let mut ap = AnnotationProcess::new();
        let mut nq1 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
        );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
        );
        nq1.add_hit(&h1);
        nq1.add_hit(&h2);
        // Test insert_query
        ap.insert_query(&mut nq1, tx.clone());
        assert!(ap.queries.contains_key(&nq1.id));
        assert_eq!(ap.queries.get(&nq1.id).unwrap().hits.len(), 2);
        // New query, but for the same `qacc`, supposedly parsed from another Blast result table:
        let mut nq2 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());
        let h3 = Hit::new(
            "Hit_Three", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
        );
        let h4 = Hit::new(
            "Hit_Four", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
        );
        nq2.add_hit(&h3);
        nq2.add_hit(&h4);
        // Test inserting the same `qacc` parsed from another input Blast result:
        ap.insert_query(&mut nq2, tx);
        for received in rx {
            ap.process_query_data_complete(received);
        }
        assert!(ap.queries.contains_key(&nq2.id));
        assert_eq!(ap.queries.len(), 1);
        // Assert that the hits from nq1 and nq2 were in fact joined in the stored query:
        assert_eq!(ap.queries.get(&nq2.id).unwrap().hits.len(), 4);
        let sq = ap.queries.get(&nq2.id).unwrap();
        assert!(sq.hits.contains_key(&h1.id));
        assert!(sq.hits.contains_key(&h2.id));
        assert!(sq.hits.contains_key(&h3.id));
        assert!(sq.hits.contains_key(&h4.id));
    }

    #[test]
    #[should_panic]
    fn insert_query_panics_in_case_of_unsorted_blast_table() {
        let (tx, _) = mpsc::channel();
        let mut ap = AnnotationProcess::new();
        let mut nq1 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());
        // Mark nq1 as already processed:
        ap.human_readable_descriptions
            .insert(nq1.id.clone(), "Unknown protein".to_string());
        // Should panic:
        ap.insert_query(&mut nq1, tx);
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
        let (tx, rx) = mpsc::channel();
        let mut ap = AnnotationProcess::new();
        ap.seq_sim_search_tables = vec!["blast_out_table.txt".to_string()];
        let mut nq1 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());
        ap.insert_query(&mut nq1, tx);
        for received in rx {
            ap.process_query_data_complete(received);
        }
        // Query should have been annotated:
        assert!(ap.human_readable_descriptions.contains_key(&nq1.id));
        assert!(!ap.queries.contains_key(&nq1.id));
        // Test families:
        let (tx1, rx1) = mpsc::channel();
        ap = AnnotationProcess::new();
        ap.seq_sim_search_tables = vec!["blast_out_table.txt".to_string()];
        let mut sf1 = SeqFamily::new();
        sf1.query_ids = vec!["Soltu.DM.02G015700.1".to_string()];
        let sf_id1 = "SeqFamily1".to_string();
        ap.insert_seq_family(sf_id1.clone(), sf1);
        nq1 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());
        ap.insert_query(&mut nq1, tx1);
        for received in rx1 {
            ap.process_query_data_complete(received);
        }
        // Family should have been annotated:
        assert!(ap.human_readable_descriptions.contains_key(&sf_id1));
        assert!(!ap.queries.contains_key(&nq1.id));
        assert!(!ap.seq_families.contains_key(&sf_id1));
        assert!(!ap.query_id_to_seq_family_id_index.contains_key(&nq1.id));
    }

    #[test]
    fn process_rest_data_works() {
        // Test Queries:
        let (tx, _) = mpsc::channel();
        let mut ap = AnnotationProcess::new();
        let mut nq1 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
        );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
        );
        nq1.add_hit(&h1);
        nq1.add_hit(&h2);
        ap.insert_query(&mut nq1, tx.clone());
        ap.process_rest_data();
        // Query should have been annotated:
        assert!(ap.human_readable_descriptions.contains_key(&nq1.id));
        assert!(!ap.queries.contains_key(&nq1.id));
        // Test Families:
        ap = AnnotationProcess::new();
        let mut sf1 = SeqFamily::new();
        sf1.query_ids = vec!["Soltu.DM.02G015700.1".to_string()];
        let sf_id1 = "SeqFamily1".to_string();
        ap.insert_seq_family(sf_id1.clone(), sf1);
        nq1 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());
        ap.insert_query(&mut nq1, tx);
        let mut sf2 = SeqFamily::new();
        sf2.query_ids = vec!["The protein without known relatives".to_string()];
        let sf_id2 = "SeqFamily2".to_string();
        ap.insert_seq_family(sf_id2.clone(), sf2);
        ap.process_rest_data();
        // Families should have been annotated:
        assert!(ap.human_readable_descriptions.contains_key(&sf_id1));
        assert!(!ap.queries.contains_key(&nq1.id));
        assert!(!ap.seq_families.contains_key(&sf_id1));
        assert!(!ap.query_id_to_seq_family_id_index.contains_key(&nq1.id));
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
            "The_Protein_Without_Blast_Hits".to_string(),
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
