use crate::generate_hrd_associated_funcs::generate_human_readable_description;
use regex::Regex;
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
    /// * `split_regex` - A reference to a regular expression used to split descriptions (`stitle`
    /// in Blast terminology) into words.
    /// * `non_informative_words_regexs` - A reference to a vector holding regular expressions used
    /// to identify non informative words, that receive only a minimum score.
    /// * `center_at_quantile` - A real value between zero and one used to center the inverse
    /// information content scores.
    pub fn annotate(
        &self,
        split_regex: &Regex,
        non_informative_words_regexs: &Vec<Regex>,
        center_at_quantile: &f64,
    ) -> Option<String> {
        if self.hits.len() > 0 {
            let hit_descriptions = self
                .hits
                .values()
                .map(|hit_desc| (*hit_desc).clone())
                .collect();
            generate_human_readable_description(
                &hit_descriptions,
                split_regex,
                non_informative_words_regexs,
                center_at_quantile,
            )
        } else {
            None
        }
    }
}
