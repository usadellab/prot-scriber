use super::default::FILTER_REGEXS;
use super::default::SPLIT_DESCRIPTION_REGEX;
use super::seq_sim_clustering::cluster;
use super::seq_sim_clustering::get_clusters;
use eq_float::F64;
use ndarray::{Array2, ArrayBase};
use regex::Regex;
use std::cmp::Ordering;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::collections::HashSet;

/// Fasta entries have a long title in which the sequence identifier and often taxonomic
/// information is given along with a short human readable protein description. We are only
/// interested in the latter. This function extracts the short description using regular
/// expressions.
///
/// # Arguments
///
/// * stitle - The sequence title line as found in the original Fasta file.
/// * regexs - A vector of regular expressions to be applied in series to the argument stitle to
///            extract the desired short description.
pub fn filter_stitle(stitle: &str, regexs: &Vec<Regex>) -> String {
    String::from(
        regexs
            .iter()
            .fold(stitle.to_string(), |accumulated, current| {
                current.replace_all(&accumulated, "").to_string()
            })
            .trim(),
    )
}

/// Calculates the overlap between a query and one of its hits as produced by sequence similarity
/// searches (e.g. Blast or Diamond). These search algorithms produce local alignments and the
/// arguments to this function. Overlap is calculated as
///
/// `((qend - qstart + 1) + (send - sstart + 1))`
///
/// # Arguments
///
/// * sstart - The start position of the local alignment in the hit
/// * send - The end position of the local alignment in the hit
/// * slen - The overall sequence length of the hit
/// * qstart - The start position of the local alignment in the query
/// * qend - The end position of the local alignment in the query
/// * qlen - The overall sequence length of the query
pub fn overlap_with_query(
    sstart: i64,
    send: i64,
    slen: i64,
    qstart: i64,
    qend: i64,
    qlen: i64,
) -> F64 {
    F64((((qend - qstart + 1) + (send - sstart + 1)) as f64) / ((qlen + slen) as f64))
}

/// A sequence similarity search result is called a "Hit" and is represented by its namesake.
#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Hit {
    pub id: String,
    pub qstart: i64,
    pub qend: i64,
    pub bitscore: F64,
    pub overlap_with_query: F64,
    pub description: String,
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
    fn overlap_with_hit(&self, with: &Hit, qlen: &f64) -> f64 {
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
    fn description_words(&self) -> HashSet<&str> {
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
    fn description_similarity(&self, to: &Hit) -> f64 {
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
    fn similarity(&self, to: &Hit, qlen: &f64) -> f64 {
        let o = self.overlap_with_hit(to, qlen);
        let dd = self.description_similarity(to);
        (o + dd) / 2.0
    }
}

/// A sequence similarity search is executed for a query sequence, which is represented by `Query`.
#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Query {
    pub id: String,
    pub qlen: F64,
    pub hits: HashMap<String, Hit>,
    // The calculation of similarity requires scaling the bit-scores by division by the maximum
    // bit-score found in all Hits of a respective Query instance:
    pub bitscore_scaling_factor: F64,
    // The dimension names of the resulting markov clustered similarity matrix:
    pub seq_sim_mtrx_node_names: Vec<String>,
    // The resulting clusters, represented by the respective sequence identifiers (`String`):
    pub clusters: Vec<Cluster>,
}

impl Query {
    /// Function translates the query's sequence similarity search results into a similarity matrix
    /// and markov clusters it. The results are stored in the query (`&self.clusters`) itself.
    /// Additionally, the cluster containing the query is marked in `self.querys_cluster_ind`.
    /// _Note_ that the query's clusters will be *sorted* by their respective score _descending_,
    /// i.e. the best scoring cluster will be the first.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A mutable reference to the query itself.
    pub fn cluster_hits(&mut self) {
        if self.hits.len() > 0 {
            let mcl_matrix = cluster(self);
            let mcl_clusters = get_clusters(&mcl_matrix);
            self.clusters = mcl_clusters
                .iter()
                .map(|clstr| {
                    let hits = clstr
                        .iter()
                        .map(|seq_node| {
                            self.seq_sim_mtrx_node_names.get(*seq_node).unwrap().clone()
                        })
                        .collect();
                    let score = self.calculate_cluster_score(&hits);
                    Cluster { hits, score }
                })
                .collect();
            // Sort by score descending:
            self.clusters.sort_by(|a, b| b.cmp(&a));
        }
    }

    /// Calculates the score of a cluster of `Hit`s generated for a Query (`&self`). The score is
    /// the arithmetic mean of two scoring measures. The first measure is the maximum similarity -
    /// see function `Query::similarity(&self, to: &Hit)` - between the Query (`&self`) and any hit
    /// contained in the cluster (`custer_index`). The second scoring measure is the fraction of
    /// all the Query's (`&self`) hits contained in the cluster. Thus, the cluster score is
    /// calculated as:
    /// score( cluster ) =
    ///   mean(
    ///     max( Query.similarity( hit_i ) ),
    ///     cluster.len() / Query.hits.len()
    ///   )
    /// Note that the score is returned as an instance of `F64`.
    ///
    /// # Arguments
    ///
    /// * `&self` - A reference to the query instance itself
    /// * `hit_ids: &Vec<String>` - The Hit sequences identifiers.
    pub fn calculate_cluster_score(&self, hit_ids: &Vec<String>) -> F64 {
        let max_similarity_between_query_and_hit = hit_ids.iter().fold(0.0, |acc, curr_hit_id| {
            let curr_query_hit_sim = self.similarity(&self.hits.get(curr_hit_id).unwrap());
            if acc > curr_query_hit_sim {
                acc
            } else {
                curr_query_hit_sim
            }
        });
        let cluster_hit_coverage = hit_ids.len() as f64 / self.hits.len() as f64;
        F64((max_similarity_between_query_and_hit + cluster_hit_coverage) / 2.0)
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
            qlen: Default::default(),
            bitscore_scaling_factor: Default::default(),
            hits: HashMap::<String, Hit>::new(),
            seq_sim_mtrx_node_names: Default::default(),
            clusters: Default::default(),
        }
    }

    /// Creates and empty (`Default`) instance of struct Query.
    pub fn new() -> Query {
        Default::default()
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

    /// Identifies the bitscore-scaling-factor as the maximum of all bitscores of the query
    /// (`&self`) hits.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - mutable reference to self (the query)
    pub fn find_bitscore_scaling_factor(&mut self) {
        // find maximum bitscore
        let bs_scl_fct = F64(self.hits.iter().fold(0.0, |accumulated, (_, hit)| {
            if hit.bitscore.0 > accumulated {
                hit.bitscore.0
            } else {
                accumulated
            }
        }));
        if bs_scl_fct == F64(0.0) {
            panic!("bitscore scaling factor for Query '{:?}' is zero. This would cause division by zero.", self.id);
        } else {
            self.bitscore_scaling_factor = bs_scl_fct;
        }
    }

    /// Calculates the similarity between the query (`self`) and the argument hit (`to`) as the mean
    /// of overlap_with_query and bitscore, where the bitscore is scaled by dividing it with the max. bitscore
    /// found for all the query's hits. Returns the similarity ranging from zero to one, inclusive.
    ///
    /// # Arguments
    ///
    /// * `&self` - the query
    /// * `to: &Hit` - the hit (one of the query's) to which the calculate the similarity to
    pub fn similarity(&self, to: &Hit) -> f64 {
        if !self.hits.contains_key(&to.id) {
            panic!(
                "The query '{}' does not contain the hit '{}' to which calculate the similarity to.",
                self.id, to.id
            );
        }
        (to.overlap_with_query.0 + to.bitscore.0 / self.bitscore_scaling_factor.0) / 2.0
    }

    /// Calculates similarities between the query (`&self`) its `hits` and in between the hits.
    /// Stores those in a quadratic two dimensional similarity matrix (`Array2<f64>`). Note that
    /// similarities of a node (matrix-cell) to itself is not yet initialized. Use function
    /// seq_sim_clustering::add_self_loops to do that. Returns a tupel, `(Vec<String>,
    /// Array2<f64>)`, with first entry the sequence identifiers, i.e. the row and column names of
    /// the second entry, which is the similarity matrix. The sequence IDs consist of the
    /// alphabetically sorted hit IDs followed by the terminal Query ID. This order is guaranteed.
    ///
    /// # Arguments
    ///
    /// * `&self` - the query
    pub fn to_similarity_matrix(&mut self) -> (Vec<String>, Array2<f64>) {
        // find maximum bit-score and use it as scaling factor for all Hits' bit-scores:
        self.find_bitscore_scaling_factor();
        // the names of the rows and columns in the similarity matrix:
        let mut seq_ids: Vec<String> = self.hits.keys().map(|k: &String| (*k).clone()).collect();
        seq_ids.sort();
        // remember the sequence identifiers for future reference:
        self.seq_sim_mtrx_node_names = seq_ids.clone();
        let mut mtrx: Array2<f64> = ArrayBase::zeros((seq_ids.len(), seq_ids.len()));
        // Save time and iterate over the similarities of the upper triangle including the
        // diagonal, but set each calculated similarity also in the lower one. Note that the inner
        // loop starts at the outer loop's index and not at zero.
        for row_ind in 0..seq_ids.len() {
            // Similarities between hit "row_ind" and ...
            for col_ind in row_ind..self.hits.len() {
                // ... and hit "col_i"
                // Note, that similarity to self will be set later using
                // seq_sim_clustering::add_self_loops
                if row_ind != col_ind {
                    let hit_row_i = self.hits.get(&seq_ids[row_ind]).unwrap();
                    let hit_col_i = self.hits.get(&seq_ids[col_ind]).unwrap();
                    let d = hit_row_i.similarity(hit_col_i, &self.qlen.0);
                    mtrx[[row_ind, col_ind]] = d;
                    mtrx[[col_ind, row_ind]] = d;
                }
            }
        }
        (seq_ids, mtrx)
    }
}

/// A cluster of sequence similarity search hits is represented by this struct.
#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Cluster {
    pub hits: Vec<String>,
    pub score: F64,
}

impl Cluster {
    /// Checks whether the cluster (`&self`) hit identifier vector contains the argument `hit_id`.
    /// Returns true if so, false otherwise.
    ///
    /// # Arguments
    ///
    /// * `&self` - The reference to the cluster to be queried
    /// * `hit_id: &String` - A reference to the hit identifier to check whether it is contained in
    ///                       this (`&self`) cluster.
    pub fn contains(&self, hit_id: &String) -> bool {
        self.hits.contains(hit_id)
    }
}

impl Ord for Cluster {
    /// Sorting of clusters should rely on their respective scores. If scores are equal the one
    /// with more hit sequences is greater than the other.
    ///
    /// # Arguments
    ///
    /// * `&self` - A reference to the cluster to be compared with another
    /// * `other: &Cluster` - A reference to another cluster `&self` should be compared to
    fn cmp(&self, other: &Cluster) -> Ordering {
        let o = self.score.cmp(&other.score);
        if o == Ordering::Equal {
            self.hits.len().cmp(&other.hits.len())
        } else {
            o
        }
    }
}

impl PartialOrd for Cluster {
    /// Sorting of clusters should rely on their respective scores. If scores are equal the one
    /// with more hit sequences is greater than the other. _Note_, that the result is an
    /// `Option<Ordering>` which makes it an implementation of partial ordering.
    ///
    /// # Arguments
    ///
    /// * `&self` - A reference to the cluster to be compared with another
    /// * `other: &Cluster` - A reference to another cluster `&self` should be compared to
    fn partial_cmp(&self, other: &Cluster) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use ndarray::arr2;

    #[test]
    fn default_filter_regexs_extract_uni_prot_descriptions() {
        let t1 = "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1";
        assert_eq!(
            filter_stitle(t1, &(*FILTER_REGEXS)),
            "Probable LRR receptor-like serine/threonine-protein kinase At3g47570"
        );
    }

    #[test]
    fn query_add_hit_only_uses_higest_bitscore() {
        let high = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let low = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "1.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let highest = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "666.6",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let other = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.add_hit(&high);
        assert_eq!(query.hits.len(), 1);
        query.add_hit(&low);
        assert_eq!(query.hits.len(), 1);
        assert_eq!(*query.hits.get("Hit_One").unwrap(), high);
        query.add_hit(&other);
        assert_eq!(query.hits.len(), 2);
        assert_eq!(*query.hits.get("Hit_One").unwrap(), high);
        assert_eq!(*query.hits.get("Hit_Two").unwrap(), other);
        query.add_hit(&highest);
        assert_eq!(query.hits.len(), 2);
        assert_eq!(*query.hits.get("Hit_One").unwrap(), highest);
        assert_eq!(*query.hits.get("Hit_Two").unwrap(), other);
    }

    #[test]
    fn test_calculate_bitscore() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "666.6",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.add_hit(&h1);
        query.add_hit(&h2);
        query.find_bitscore_scaling_factor();
        assert_eq!(query.bitscore_scaling_factor.0, 666.6);
    }

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

    #[test]
    fn similarity_between_query_and_hit() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.add_hit(&h1);
        query.add_hit(&h2);
        query.find_bitscore_scaling_factor();
        assert_eq!(h2.overlap_with_query.0, 0.5);
        assert_eq!(query.bitscore_scaling_factor.0, 500.0);
        let d = query.similarity(&h2);
        // similarity( query, h2 ) =
        // mean( overlap[ .5 ] + scaled_bitscore[ 1/5 ] )
        assert_eq!(d, (0.5 + 1.0 / 5.0) / 2.0);
    }

    #[test]
    pub fn test_to_similarity_matrix() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.qlen = F64(100.0);
        query.add_hit(&h1);
        query.add_hit(&h2);
        query.find_bitscore_scaling_factor();
        let sim_mtrx = query.to_similarity_matrix();
        // println!("Query.to_similarity_matrix() yields:\n{:?}\n", sim_mtrx);
        let expected = arr2(&[
            [0.0, h1.similarity(&h2, &query.qlen.0)],
            [h1.similarity(&h2, &query.qlen.0), 0.0],
        ]);
        assert_eq!(sim_mtrx.1, expected);
        assert!(query.seq_sim_mtrx_node_names.len() == 2);
        assert_eq!(query.seq_sim_mtrx_node_names, vec!["Hit_One", "Hit_Two"]);
        // Test query that has hits aligning to disjoint regions of the query sequence and do not
        // share words in their respective descriptions:
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
        query_frivolous.qlen = F64(100.0);
        query_frivolous.add_hit(&h3);
        query_frivolous.add_hit(&h4);
        query_frivolous.add_hit(&h5);
        let sim_mtrx = query_frivolous.to_similarity_matrix();
        // println!(
        //     "Query with hits that align to disjoint regions and do not share words in their respective descriptions yield similarity matrix:\n{:?}",
        //     sim_mtrx
        // );
        let expected_sim_mtrx = arr2(&[[0.0, 0.0, 0.0], [0.0, 0.0, 0.75], [0.0, 0.75, 0.0]]);
        assert_eq!(sim_mtrx.1, expected_sim_mtrx);
    }

    #[test]
    pub fn test_cluster_hits() {
        // Test query that should result in a single cluster:
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.qlen = F64(100.0);
        query.add_hit(&h1);
        query.add_hit(&h2);
        query.cluster_hits();
        assert_eq!(query.clusters.len(), 1);
        let querys_cluster = &query.clusters.get(0).unwrap();
        // println!("query.cluster_hits() yields:\n{:?}", querys_cluster);
        assert!(querys_cluster.contains(&"Hit_One".to_string()));
        assert!(querys_cluster.contains(&"Hit_Two".to_string()));
        // Test query that has NO hits:
        let mut query_no_hits = Query::from_qacc("Query_So_Lonely".to_string());
        query_no_hits.qlen = F64(123.0);
        query_no_hits.cluster_hits();
        assert_eq!(query_no_hits.clusters.len(), 0);
        // Test query that should produce TWO clusters:
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
        query_frivolous.qlen = F64(100.0);
        query_frivolous.add_hit(&h3);
        query_frivolous.add_hit(&h4);
        query_frivolous.add_hit(&h5);
        query_frivolous.cluster_hits();
        // println!(
        //     "Query with three Hits yields these clusters:\n{:?}",
        //     query_frivolous.clusters
        // );
        // Test includes assertion of correct sort order of clusters:
        assert_eq!(query_frivolous.clusters.len(), 2);
        assert!(query_frivolous
            .clusters
            .get(0)
            .unwrap()
            .contains(&"Hit_Three".to_string()));
        assert!(query_frivolous
            .clusters
            .get(0)
            .unwrap()
            .contains(&"Hit_Four".to_string()));
        assert!(query_frivolous
            .clusters
            .get(1)
            .unwrap()
            .contains(&"Hit_Five".to_string()));
    }

    #[test]
    fn test_calculate_cluster_score() {
        let h3 = Hit::new(
            "Hit_Three", "100", "1", "50", "100", "51", "100", "500.0",
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
        query_frivolous.qlen = F64(100.0);
        query_frivolous.add_hit(&h3);
        query_frivolous.add_hit(&h4);
        query_frivolous.add_hit(&h5);
        query_frivolous.cluster_hits();
        // println!(
        //     "Query with three Hits yields these clusters:\n{:?}",
        //     query_frivolous.clusters
        // );
        let h3_h4_indx = query_frivolous
            .clusters
            .iter()
            .position(|clstr| clstr.contains(&"Hit_Three".to_string()))
            .unwrap();
        let cluster_h3_h4_score = query_frivolous.calculate_cluster_score(
            &query_frivolous
                .clusters
                .get(h3_h4_indx)
                .unwrap()
                .hits
                .clone(),
        );
        let expected_h3_h4_clstr_score = F64((query_frivolous.similarity(&h3) + 2.0 / 3.0) / 2.0);
        // println!(
        //     "Cluster of Hits 'Hit_Three' and 'Hit_Four' of Query '{}' has index {} and cluster-score {}",
        //     query_frivolous.id, h3_h4_indx, cluster_h3_h4_score
        // );
        assert_eq!(cluster_h3_h4_score, expected_h3_h4_clstr_score);
        let h5_indx = query_frivolous
            .clusters
            .iter()
            .position(|clstr| clstr.contains(&"Hit_Five".to_string()))
            .unwrap();
        let cluster_h5_score = query_frivolous
            .calculate_cluster_score(&query_frivolous.clusters.get(h5_indx).unwrap().hits.clone());
        let expected_h5_clstr_score = F64((query_frivolous.similarity(&h5) + 1.0 / 3.0) / 2.0);
        // println!(
        //     "Cluster of Hit 'Hit_Five' of Query '{}' has index {} and cluster-score {}",
        //     query_frivolous.id, h5_indx, cluster_h5_score
        // );
        assert_eq!(cluster_h5_score, expected_h5_clstr_score);
    }
}
