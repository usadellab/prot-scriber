use super::query::Query;
use super::seq_family::SeqFamily;
use std::collections::HashMap;
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::thread;

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

impl AnnotationProcess {
    /// Creates and empty (`Default`) instance of struct AnnotationProcess.
    pub fn new() -> AnnotationProcess {
        Default::default()
    }

    /// The central function that runs an annotation process.
    ///
    /// # Arguments
    ///
    /// * `&self` - A reference to an instance of AnnotationProcess
    pub fn run(&mut self) {
        // let (tx, rx) = mpsc::channel();
        // let qhd_arc: Arc<Mutex<HashMap<String, Vec<String>>>> =
        //     Arc::new(Mutex::new((*self).queries));

        // // Parse and process each sequence similarity search result table in a dedicated
        // // thread:
        // for sss_tbl in (*self).seq_sim_search_tables {
        //     let tx_i = tx.clone();
        //     let qhd_arc_i = qhd_arc.clone();

        //     // Start this sss_tbl's dedicated threat:
        //     thread::spawn(move || {
        //         let vals = vec![
        //             format!("({}) hello", i),
        //             env!("CARGO_MANIFEST_DIR").to_string(),
        //             String::from("from"),
        //             String::from("the"),
        //             String::from("thread"),
        //         ];

        //         for val in vals {
        //             let mut vec_i = qhd_arc_i.lock().unwrap();
        //             vec_i.push(val);
        //         }
        //         tx_i.send(format!("thread {}'s message", i)).unwrap();
        //     });
        // }

        // // Receiver thread will wait eternally, as long as tx has not been dropped:
        // drop(tx);

        // // Listen to messages coming from any thread:
        // for received in rx {
        //     println!("Got: {}", received);
        // }
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
    pub fn insert_query(&mut self, query: &mut Query) {
        // To Do: panic! if query.id already in results, i.e. self.human_readable_descriptions
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
            // If yes, then process the parsed data:
            self.process_query_data_complete(sqi);
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
                    panic!("Biological sequence '{:?}' already set as member of family '{:?}'. But found '{:?}' again declared as member of another family '{:?}'.\nMake sure each biological sequence appears in one and only one family to avoid this problem.", query_id, other_family_id, query_id, seq_family_id);
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

    pub fn process_query_data_complete(&mut self, query_id: String) {
        let mode = self.mode();
        match mode {
            // Handle annotation of single biological sequences:
            AnnotationProcessMode::SequenceAnnotation => {
                // Generate the desired result, i.e. a human readable description for the Query:
                let hrd = self.queries.get(&query_id).unwrap().annotate();
                // Add the new result to the in memory database, i.e.
                // `self.human_readable_descriptions`:
                self.human_readable_descriptions
                    .insert(query_id.clone(), hrd);
                // Free memory by removing the parsed input data, no longer required:
                self.queries.remove(&query_id);
            }
            // Handle annotation of sets of biological sequences, so called "Gene Families":
            AnnotationProcessMode::FamilyAnnotation => {
                // Get SeqFamily for query sequence identifier (qs_id)
                // - if no family for qs_id can be found, annotate query as in SequenceAnnotation
                // mode
                // Tell SeqFamily about the completed data
                // Ask SeqFamily if all queries have complete data ...
                // ... and if so, annotate,
                // ... add resulting HRD to in memory database of results,
                // ... and remove data for memory efficiency.
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        // To Do
        assert!(false);
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
}
