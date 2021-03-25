use eq_float::F64;
use std::cmp::Ordering;

/// A cluster of sequence similarity search hits is represented by this struct.
#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Cluster {
    pub hits: Vec<String>,
    pub score: F64,
    pub aligned_query_region: Option<(u32, u32)>,
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
