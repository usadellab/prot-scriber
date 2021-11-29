use std::collections::HashMap;
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::thread;

#[derive(Debug, Clone, Default)]
pub struct AnnotationProcess {
    pub seq_sim_search_tables: Vec<String>,
    pub queriesHitDescriptions: HashMap<String, Vec<String>>,
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
        let (tx, rx) = mpsc::channel();
        let qhd_arc: Arc<Mutex<HashMap<String, Vec<String>>>> =
            Arc::new(Mutex::new((*self).queriesHitDescriptions));

        // Parse and process each sequence similarity search result table in a dedicated
        // thread:
        for sss_tbl in (*self).seq_sim_search_tables {
            let tx_i = tx.clone();
            let qhd_arc_i = qhd_arc.clone();

            // Start this sss_tbl's dedicated threat:
            thread::spawn(move || {
                let vals = vec![
                    format!("({}) hello", i),
                    env!("CARGO_MANIFEST_DIR").to_string(),
                    String::from("from"),
                    String::from("the"),
                    String::from("thread"),
                ];

                for val in vals {
                    let mut vec_i = qhd_arc_i.lock().unwrap();
                    vec_i.push(val);
                }
                tx_i.send(format!("thread {}'s message", i)).unwrap();
            });
        }

        // Receiver thread will wait eternally, as long as tx has not been dropped:
        drop(tx);

        // Listen to messages coming from any thread:
        for received in rx {
            println!("Got: {}", received);
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
    }
}
