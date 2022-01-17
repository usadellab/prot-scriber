use super::default::{SPLIT_DESCRIPTION_REGEX, UNKNOWN_PROTEIN_DESCRIPTION};
use crate::generate_hrd_associated_funcs::generate_human_readable_description;
use std::collections::HashMap;

/// A sequence similarity search is executed for a query sequence, which is represented by `Query`.
#[derive(Debug, Clone, Default)]
pub struct Query {
    /// The sequence similarity search results (Blast Hits)
    pub hits: HashMap<String, String>,
    /// A counter of how many times this query was parsed in sequence similarity search results
    pub n_parsed_from_sssr_tables: u16,
}

/// Representation of a query in a sequence similarity search (SSS), e.g. Blast or Diamond.
impl Query {
    /// Returns a new and initialized instance of struct `Query`.
    pub fn new() -> Query {
        Query {
            hits: HashMap::<String, String>::new(),
            n_parsed_from_sssr_tables: 0,
        }
    }

    /// Generates and returns a human readable description (`String`) for this biological query
    /// sequence.
    ///
    /// # Arguments
    ///
    /// * `&self` - A mutable reference to self, this instance of Query
    pub fn annotate(&self) -> String {
        let mut hrd: String = (*UNKNOWN_PROTEIN_DESCRIPTION).to_string();
        if self.hits.len() > 0 {
            let hit_descriptions = self
                .hits
                .values()
                .map(|hit_desc| (*hit_desc).clone())
                .collect();
            let hrd_option =
                generate_human_readable_description(&hit_descriptions, &(*SPLIT_DESCRIPTION_REGEX));
            match hrd_option {
                Some(hum_read_desc) => {
                    hrd = hum_read_desc;
                }
                None => {}
            }
        }
        hrd
    }
}
