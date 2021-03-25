use super::cluster::*;
use super::hit::*;
use super::seq_sim_clustering::*;
use eq_float::F64;
use ndarray::{Array2, ArrayBase};
use regex::Regex;
use std::collections::{HashMap, HashSet};

// NOTE that unit tests for struct Query are located in `./src/query_tests.rs`.

/// A sequence similarity search is executed for a query sequence, which is represented by `Query`.
#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Query {
    pub id: String,
    pub qlen: u32,
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
                    let mut hits = clstr
                        .iter()
                        .map(|seq_node| {
                            self.seq_sim_mtrx_node_names.get(*seq_node).unwrap().clone()
                        })
                        .collect();
                    let score = self.calculate_cluster_score(&hits);
                    let aligned_query_region = self.cluster_aligned_query_region(&hits);
                    // Sort hit identifier (in `hits`) by their respective
                    // Hit.query_similarity_score descending:
                    hits.sort_by(|a, b| {
                        let a_hit = self.hits.get(a).unwrap();
                        let b_hit = self.hits.get(b).unwrap();
                        // See Ord implementation for struct Query:
                        b_hit.cmp(a_hit)
                    });
                    Cluster {
                        hits,
                        score,
                        aligned_query_region,
                    }
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
    /// Note that the score is returned as an instance of `F64`. Also _note_ that this function
    /// initializes the field `query_similarity_score` of each instance of struct Hit referenced by
    /// `hit_ids`.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - A mutable reference to the query instance itself
    /// * `hit_ids: &Vec<String>` - The Hit sequences identifiers.
    pub fn calculate_cluster_score(&mut self, hit_ids: &Vec<String>) -> F64 {
        let mut max_similarity_between_query_and_hit = 0.0;
        for curr_hit_id in hit_ids.iter() {
            let ch_score = self.similarity(self.hits.get(curr_hit_id).unwrap());
            if ch_score > max_similarity_between_query_and_hit {
                max_similarity_between_query_and_hit = ch_score;
            }
            // Remember the similarity score between the Query (`self`) this Hit (`curr_hit_id`)
            // for future references:
            self.hits
                .get_mut(curr_hit_id)
                .unwrap()
                .query_similarity_score = F64(ch_score);
        }
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
                    let d = hit_row_i.similarity(hit_col_i, &self.qlen);
                    mtrx[[row_ind, col_ind]] = d;
                    mtrx[[col_ind, row_ind]] = d;
                }
            }
        }
        (seq_ids, mtrx)
    }

    /// Generates a consensus description for a cluster (argument `cluster_indx`) of this query
    /// (`&self`). Hits' descriptions, contained in the cluster, are split using argument regular
    /// expression. The returned consensus description is composed of words in the intersection of
    /// all Hits' descriptions in the order they appear in the best scoring (see
    /// `Query::similarity(&self, to: &Hit)`) hit's description.
    ///
    /// # Arguments
    ///
    /// * `&self` - A reference to the Query containing the cluster `cluster_indx`
    /// * `cluster_indx` - A usize indicating for which of the Query's cluster to generate the
    ///                    consensus description.
    /// * `splitting_axe` - A reference to the regular expression used to split Hit descriptions
    ///                     into sets or vectors of words.
    pub fn cluster_consensus_description(
        &self,
        cluster_indx: usize,
        splitting_axe: &Regex,
    ) -> String {
        let cluster = self.clusters.get(cluster_indx).unwrap();
        let best_scoring_hit_description: Vec<String> = splitting_axe
            .split(
                &self
                    .hits
                    .get(cluster.hits.get(0).unwrap())
                    .unwrap()
                    .description,
            )
            .into_iter()
            .map(|x| String::from(x))
            .collect();
        let hits_description_words: Vec<HashSet<String>> = cluster
            .hits
            .iter()
            .map(|hit_id| {
                splitting_axe
                    .split(&self.hits.get(hit_id).unwrap().description)
                    .into_iter()
                    .map(|x| String::from(x))
                    .collect()
            })
            .collect();
        let mut shared_words = hits_description_words.get(0).unwrap().clone();
        for i in 1..hits_description_words.len() {
            shared_words = shared_words
                .intersection(hits_description_words.get(i).unwrap())
                .into_iter()
                .map(|x| String::from(x))
                .collect();
        }
        best_scoring_hit_description
            .into_iter()
            .filter(|x| shared_words.contains(x))
            .collect::<Vec<String>>()
            .join(" ")
    }

    /// Function identifies the region in the pairwise local alignments between the hits and the
    /// query. The region is returned as that stretch of the query sequence that is covered by all
    /// of the pairwise local alignments between the hit_i and the query. Note that an Option<(u32,
    /// u32)>` is returned. In case there is _no_ query sequence region a cluster aligns to, None
    /// is returned.
    ///
    /// # Arguments
    ///
    /// * `&self` - A reference to the instance of Query.
    /// * `hit_ids: &Vec<String>` - A reference to a vector of Hit identifiers that should be
    ///                             present among this query's hits. The region returned by this
    ///                             function will be calculated for these argument hit identifiers.
    pub fn cluster_aligned_query_region(&self, hit_ids: &Vec<String>) -> Option<(u32, u32)> {
        let mut max_qstart = 0u32;
        let mut min_qend = self.qlen;
        for hit_id in hit_ids {
            let hit_i = self.hits.get(hit_id).unwrap();
            if hit_i.qstart > max_qstart {
                max_qstart = hit_i.qstart;
            }
            if hit_i.qend < min_qend {
                min_qend = hit_i.qend;
            }
        }
        if max_qstart < min_qend {
            Some((max_qstart, min_qend))
        } else {
            None
        }
    }
}

// For unit tests of struct Query see separate module `./src/query_tests.rs`. This file grew too
// large, so the tests were moved.
