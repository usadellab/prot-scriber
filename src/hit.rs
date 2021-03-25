use super::default::FILTER_REGEXS;
use super::default::SPLIT_DESCRIPTION_REGEX;
use super::model_funcs::{filter_stitle, overlap_with_query};
use eq_float::F64;
use std::cmp::{max, min, Ordering};
use std::collections::HashSet;

/// A sequence similarity search result is called a "Hit" and is represented by its namesake.
#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Hit {
    pub id: String,
    pub qstart: u32,
    pub qend: u32,
    pub bitscore: F64,
    pub overlap_with_query: F64,
    pub description: String,
    // the "similarity score" between the Query and this Hit:
    pub query_similarity_score: F64,
}

impl Ord for Hit {
    /// Sorting of hits should rely on their respective query_similarity_scores. If scores are
    /// equal the one with more hit sequences is greater than the other.
    ///
    /// # Arguments
    ///
    /// * `&self` - A reference to the cluster to be compared with another
    /// * `other: &Hit` - A reference to another cluster `&self` should be compared to
    fn cmp(&self, other: &Hit) -> Ordering {
        self.query_similarity_score
            .cmp(&other.query_similarity_score)
    }
}

impl PartialOrd for Hit {
    /// Sorting of hits should rely on their respective query_similarity_scores. If scores are
    /// equal the one with more hit sequences is greater than the other. _Note_, that the result is
    /// an `Option<Ordering>` which makes it an implementation of partial ordering.
    ///
    /// # Arguments
    ///
    /// * `&self` - A reference to the cluster to be compared with another
    /// * `other: &Hit` - A reference to another cluster `&self` should be compared to
    fn partial_cmp(&self, other: &Hit) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Hit {
    /// Constructor to generate a new instance of structure `Hit`. The arguments are mostly named
    /// after their namesakes in the output table as produced e.g. by Blast or Diamond (see e.g.
    /// their respective manuals for more details).
    ///
    /// # Arguments
    ///
    /// * `id: &str` - The hit sequence's accession (identifier) as stored in the `sacc` column
    /// * `qlen: &str` - The query sequence's length
    /// * `qstart: &str` - The start of the local alignment in the query sequence
    /// * `qend: &str` - The end of the local alignment in the query sequence
    /// * `slen: &str` - The length of the hit sequence
    /// * `sstart: &str` - The start of the local alignment in the hit sequence
    /// * `send: &str` - The end of the local alignment in the hit sequence
    /// * `bitscore: &str` - The bitscore of the local alignment
    /// * `stitle: &str` - The full title line of the hit in the reference database (stored in
    ///                    fasta format)
    pub fn new(
        id: &str,
        qlen: &str,
        qstart: &str,
        qend: &str,
        slen: &str,
        sstart: &str,
        send: &str,
        bitscore: &str,
        stitle: &str,
    ) -> Hit {
        let slen_prsd = slen.parse().unwrap();
        let sstart_prsd = sstart.parse().unwrap();
        let send_prsd = send.parse().unwrap();
        let qlen_prsd = qlen.parse().unwrap();
        let qstart_prsd = qstart.parse().unwrap();
        let qend_prsd = qend.parse().unwrap();
        Hit {
            id: String::from(id),
            qstart: qstart_prsd,
            qend: qend_prsd,
            bitscore: F64(bitscore.parse().unwrap()),
            overlap_with_query: overlap_with_query(
                sstart_prsd,
                send_prsd,
                slen_prsd,
                qstart_prsd,
                qend_prsd,
                qlen_prsd,
            ),
            description: filter_stitle(&stitle, &(*FILTER_REGEXS)),
            query_similarity_score: Default::default(),
        }
    }

    /// Overlap between hits `o.hits( h1, h2 )` is defined as the overlap between the query's
    /// region the two respective hits aligned to. _Note_: Make sure both hits belong to the same
    /// query! Returns the overlap (`f64`) calculated as follows:
    ///
    /// o.hits( h1, h2 ) =
    ///   ( min(qend.h1, qend.h2) - max(qstart.h1, qstart.h2) + 1.0 ) / qlen
    ///
    /// _Note_ that the above `qlen` is of course identical between the two hits (`self` and
    /// `with`).
    ///
    /// # Arguments
    ///
    /// * `&self` - reference to an instance of Hit
    /// * `with: &Hit` - the hit which which to calculate the overlap
    /// * `qlenL &f64` - the query sequence's length
    pub fn overlap_with_hit(&self, with: &Hit, qlen: &f64) -> f64 {
        if *qlen <= 0.0 {
            panic!(
                "Function Hit.overlap_with_hit was given invalid query length (qlen) argument <= 0: {}",
                qlen
            );
        }
        (min(self.qend, with.qend) as f64 - max(self.qstart, with.qstart) as f64 + 1.0) / qlen
    }

    /// Splits the Hit's description into words using the argument regular expression.
    ///
    /// # Arguments
    ///
    /// * `&self` - The Hit whose description to split into words
    pub fn description_words(&self) -> HashSet<&str> {
        (*SPLIT_DESCRIPTION_REGEX)
            .split(&self.description)
            .into_iter()
            .collect::<HashSet<&str>>()
    }

    /// Computes a numerical similarity ranging from zero to one between two respective descriptions.
    /// The formula is number of shared words divided by the shorter description's length (as in
    /// number of words). _Note_ that double appearances of words are not counted.
    ///
    /// # Arguments
    ///
    /// * `&self` - The Hit whose description to compare to the other
    /// * `to: &Hit` - The other Hit whose description to compute the similarity to
    pub fn description_similarity(&self, to: &Hit) -> f64 {
        let self_words = self.description_words();
        let to_words = to.description_words();
        let intersection: HashSet<_> = self_words.intersection(&to_words).collect();
        intersection.len() as f64 / min(self_words.len(), to_words.len()) as f64
    }

    /// Computes the similarity between two hits (`&self` and `to: &Hit`) based on a linear
    /// combination of the overlap between the two hits respective local alignment with the
    /// original query and the similarity between the hits' descriptions (see `description_similarity`
    /// for details).
    ///
    /// # Arguments
    ///
    /// * `&self` - The instance of Hit itself
    /// * `to: &Hit` - The other instance of Hit taken from sequence similarity searches for the
    ///                same query.
    /// * `qlen: &f64` - The length of the query the two argument hits were found for.
    pub fn similarity(&self, to: &Hit, qlen: &f64) -> f64 {
        let o = self.overlap_with_hit(to, qlen);
        let dd = self.description_similarity(to);
        (o + dd) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn parse_hit_from_strs() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        assert_eq!(h1.id, String::from("Hit_One"));
        assert_eq!(h1.bitscore, F64(123.4));
        assert_eq!(h1.qstart, 1);
        assert_eq!(h1.qend, 50);
        assert_approx_eq!(h1.overlap_with_query.0, 0.37, 0.01);
        assert_eq!(
            h1.description,
            "Probable LRR receptor-like serine/threonine-protein kinase At3g47570"
        );
    }

    #[test]
    fn calculate_overlap_between_hits() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "26", "50", "200", "76", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let o = h1.overlap_with_hit(&h2, &100.0);
        assert_eq!(o, 0.25);
    }

    #[test]
    fn split_description_into_words() {
        let h = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let expected_words: HashSet<&str> = vec![
            "Probable",
            "LRR",
            "receptor-like",
            "serine/threonine-protein",
            "kinase",
            "At3g47570",
        ]
        .into_iter()
        .collect();
        assert_eq!(h.description_words(), expected_words);
    }

    #[test]
    fn similarity_between_hits() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "26", "50", "200", "76", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let d1 = h1.similarity(&h2, &100.0);
        assert_eq!(d1, 0.625);
        let h3 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h4 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let d2 = h3.similarity(&h4, &100.0);
        // println!(
        //     "h4.description_similarity(h3) = {}",
        //     h4.description_similarity(&h3)
        // );
        assert_eq!(d2, 0.75);
    }

    #[test]
    fn test_similarity_between_hits_is_never_infinite() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h3 = Hit::new(
            "Hit_Three", "100", "51", "100", "300", "201", "300", "50.0",
            "sp|P15538|C11B1_HUMAN Cytochrome P450 11B1, mitochondrial OS=Homo sapiens OX=9606 GN=CYP11B1 PE=1 SV=5"
            );
        // let mut query = Query::from_qacc("Query_Curious".to_string());
        // query.qlen = F64(100.0);
        // query.add_hit(&h1);
        // query.add_hit(&h2);
        let sim = h1.similarity(&h2, &100.0);
        assert!(!sim.is_infinite());
        assert_eq!(sim, 0.75);
        let sim_zero = h2.similarity(&h3, &100.0);
        assert!(!sim_zero.is_infinite());
        assert_eq!(sim_zero, 0.0);
        // Test similarity between hits after setting the bit-score scaling factor:
        let h3 = Hit::new(
            "Hit_Three", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h4 = Hit::new(
            "Hit_Four", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h5 = Hit::new(
            "Hit_Five", "100", "51", "100", "300", "201", "300", "10.0",
            "sp|P15538|C11B1_HUMAN Cytochrome P450 11B1, mitochondrial OS=Homo sapiens OX=9606 GN=CYP11B1 PE=1 SV=5"
            );
        let mut query_frivolous = Query::from_qacc("Query_Frivolous".to_string());
        query_frivolous.add_hit(&h3);
        query_frivolous.add_hit(&h4);
        query_frivolous.add_hit(&h5);
        query_frivolous.find_bitscore_scaling_factor();
        let h3_h5 = h3.similarity(&h5, &100.0);
        // println!(
        //     "Similarity between Hits that align to DISJOINT regions of the query and have no shared words in their respective descriptions is: {:?}",
        //     h3_h5
        // );
        assert!(!h3_h5.is_nan());
        assert!(!h3_h5.is_infinite());
        assert_eq!(h3_h5, 0.0);
    }
}
