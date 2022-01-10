use super::default::UNKNOWN_PROTEIN_DESCRIPTION;
use super::hit::*;
use crate::generate_hrd_associated_funcs::generate_human_readable_description;
use std::collections::HashMap;

// NOTE that unit tests for struct Query are located in `./src/query_tests.rs`.

/// A sequence similarity search is executed for a query sequence, which is represented by `Query`.
#[derive(Debug, Clone, Default)]
pub struct Query {
    /// The sequence identifier
    pub id: String,
    /// The sequence similarity search results (Blast Hits)
    pub hits: HashMap<String, Hit>,
    /// A counter of how many times this query was parsed in sequence similarity search results
    pub n_parsed_from_sssr_tables: u16,
}

/// Representation of a query in a sequence similarity search (SSS), e.g. Blast or Diamond.
impl Query {
    /// Returns a new and initialized instance of struct `Query`.
    pub fn new() -> Query {
        Query {
            id: Default::default(),
            hits: HashMap::<String, Hit>::new(),
            n_parsed_from_sssr_tables: 0,
        }
    }

    /// Creates an empty Query with an initialized empty HashMap (`hits`) and initialized `Default`
    /// `bitscore_scaling_factor`. Sets the `id` field to the provided argument.
    ///
    /// # Arguments
    ///
    /// * id - The `qacc` value, i.e. the sequence accession (identifier) of the query sequence
    /// as provided to the used sequence similarity search tool (e.g. Blast or Diamond).
    pub fn from_qacc(id: String) -> Query {
        Query {
            id,
            hits: HashMap::<String, Hit>::new(),
            n_parsed_from_sssr_tables: 0,
        }
    }

    /// Adds a hit instance to the query's (`&self`) hits HashMap field. The hit is added only, if
    /// it either does not yet exist among `hits` or if its biscore is higher than the recorded hit
    /// of identical `id`. In the latter case the argument hit replaces the previously stored one.
    /// Function returns the number of hits in the query's `hits` HashMap.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - mutable reference to self (the query)
    /// * `hit` - a reference to the hit instance to be added
    pub fn add_hit(&mut self, hit: &Hit) -> usize {
        if !self.hits.contains_key(&hit.id)
            || self.hits.get(&hit.id).unwrap().bitscore < hit.bitscore
        {
            self.hits.insert(hit.id.clone(), hit.clone());
        }
        self.hits.len()
    }

    /// Simple helper method that iteratively adds the Hit instances of argument `query` to
    /// self.hits.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - mutable reference to self (the query)
    /// * `&query` - reference to a Query whose hits to add to self
    pub fn add_hits(&mut self, query: &Query) {
        for h in query.hits.values() {
            self.add_hit(&h);
        }
    }

    /// Generates and returns a human readable description (`String`) for this biological query
    /// sequence.
    ///
    /// # Arguments
    ///
    /// * `&self` - A mutable reference to self, this instance of Query
    pub fn annotate(&self) -> String {
        if self.hits.len() > 0 {
            let hit_descriptions = self.hits.values().map(|h| h.description.clone()).collect();
            generate_human_readable_description(&hit_descriptions)
        } else {
            (*UNKNOWN_PROTEIN_DESCRIPTION).to_string()
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::query::*;

    #[test]
    fn new_query_is_correctly_initialized() {
        let nq1 = Query::from_qacc("Soltu.DM.02G015700.1".to_string());
        assert_eq!(nq1.n_parsed_from_sssr_tables, 0);
        assert_eq!(nq1.id, "Soltu.DM.02G015700.1".to_string());
        let nq2 = Query::new();
        assert_eq!(nq2.n_parsed_from_sssr_tables, 0);
    }

    #[test]
    fn query_add_hit_only_uses_higest_bitscore() {
        let high = Hit::new(
            "hit_One", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let low = Hit::new(
            "hit_One", "1.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let highest = Hit::new(
            "hit_One", "666.6",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let other = Hit::new(
            "hit_Two", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.add_hit(&high);
        assert_eq!(query.hits.len(), 1);
        query.add_hit(&low);
        assert_eq!(query.hits.len(), 1);
        assert_eq!(*query.hits.get("hit_One").unwrap(), high);
        query.add_hit(&other);
        assert_eq!(query.hits.len(), 2);
        assert_eq!(*query.hits.get("hit_One").unwrap(), high);
        assert_eq!(*query.hits.get("hit_Two").unwrap(), other);
        query.add_hit(&highest);
        assert_eq!(query.hits.len(), 2);
        assert_eq!(*query.hits.get("hit_One").unwrap(), highest);
        assert_eq!(*query.hits.get("hit_Two").unwrap(), other);
    }
}
