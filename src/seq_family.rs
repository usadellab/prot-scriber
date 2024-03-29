use super::query::Query;
use crate::generate_hrd_associated_funcs::generate_human_readable_description;
use regex::Regex;
use std::collections::HashMap;

/// Represenation of a set of biological sequences, e.g. a gene family generated by sequence
/// similarity based clustering.
#[derive(Debug, Clone, Default)]
pub struct SeqFamily {
    /// The biological sequence identifiers this SeqFamily comprises:
    pub query_ids: Vec<String>,
    /// Those Query-Identifiers for which all input sequence similarity search result files have
    /// produced data. So in other words those Query-IDs that are ready to be used as input for the
    /// generation of a human readable description:
    pub query_ids_with_complete_data: Vec<usize>,
}

impl SeqFamily {
    /// Creates and empty (`Default::default()`) instance of struct SeqFamily.
    pub fn new() -> SeqFamily {
        Default::default()
    }

    /// Returns `true` if and only if all query identifiers in argument `self.query_ids` are
    /// contained in `self.query_ids_with_complete_data`, false otherwise.
    ///
    /// # Arguments
    ///
    /// * `&self` - a reference to an instance of SeqFamily
    pub fn all_query_data_complete(&self) -> bool {
        (0..(self.query_ids.len()))
            .into_iter()
            .all(|indx| self.query_ids_with_complete_data.contains(&indx))
    }

    /// Stores the information that the Query of argument `query_id` has been parsed completely, so
    /// all of its associated sequence similarity search result data has been successfully read
    /// from the provided input. To implement this the argument `query_id`'s index in the field
    /// `query_ids` is stored in field `query_ids_with_complete_data`.
    ///
    /// # Arguments
    ///
    /// * `&self` - a reference to an instance of SeqFamily
    /// * `query_id: &String` - a reference to the Query Identifier to be marked as successfully
    ///                         and completely parsed.
    pub fn mark_query_id_with_complete_data(&mut self, query_id: &String) {
        let query_indx = self
            .query_ids
            .iter()
            .position(|qid| *qid == *query_id)
            .unwrap();
        self.query_ids_with_complete_data.push(query_indx);
    }

    /// Generates and returns a human readable description (`String`) for this set (family) of
    /// biological query sequences.
    ///
    /// # Arguments
    ///
    /// * `&self` - A mutable reference to self, this instance of SeqFamily
    /// * `queries: &HashMap<String, Query>` - A constant reference to the in memory database of
    /// `Query` instances. This is used to extract the `Hit.description`s from.
    /// * `split_regex` - A reference to a regular expression used to split descriptions (`stitle`
    /// in Blast terminology) into words.
    /// * `non_informative_words_regexs` - A reference to a vector holding regular expressions used
    /// to identify non informative words, that receive only a minimum score.
    /// * `center_at_quantile` - A real value between zero and one used to center the inverse
    /// information content scores.
    pub fn annotate(
        &self,
        queries: &HashMap<String, Query>,
        split_regex: &Regex,
        non_informative_words_regexs: &Vec<Regex>,
        center_at_quantile: &f64,
    ) -> Option<String> {
        let mut hit_descriptions: Vec<String> = vec![];
        // Gather all Hit descriptions of all queries belonging to this sequence family. This
        // means collecting all queries' hit-descriptions:
        for qid in self.query_ids.iter() {
            // If the searches found hits of significant similarity for the query sequence:
            if queries.contains_key(qid) {
                for (_, hit_desc) in &queries.get(qid).unwrap().hits {
                    hit_descriptions.push(hit_desc.clone());
                }
            }
        }
        if hit_descriptions.len() > 0 {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_query_data_complete_works() {
        let mut sf1 = SeqFamily::new();
        assert!(sf1.all_query_data_complete());
        sf1.query_ids.push("Query1".to_string());
        sf1.query_ids.push("Query2".to_string());
        sf1.query_ids.push("Query3".to_string());
        assert!(!sf1.all_query_data_complete());
        sf1.query_ids_with_complete_data.push(0);
        sf1.query_ids_with_complete_data.push(1);
        assert!(!sf1.all_query_data_complete());
        sf1.query_ids_with_complete_data.push(2);
        assert!(sf1.all_query_data_complete());
    }

    #[test]
    fn mark_query_id_with_complete_data_works() {
        let mut sf1 = SeqFamily::new();
        assert!(sf1.all_query_data_complete());
        sf1.query_ids.push("Query1".to_string());
        sf1.query_ids.push("Query2".to_string());
        sf1.query_ids.push("Query3".to_string());
        assert!(!sf1.all_query_data_complete());
        sf1.mark_query_id_with_complete_data(&"Query1".to_string());
        sf1.mark_query_id_with_complete_data(&"Query2".to_string());
        assert!(!sf1.all_query_data_complete());
        sf1.mark_query_id_with_complete_data(&"Query3".to_string());
        assert!(sf1.all_query_data_complete());
    }
}
