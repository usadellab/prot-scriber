use eq_float::F64;
use std::cmp::Ordering;

/// A cluster of sequence similarity search hits is represented by this struct.
#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Cluster {
    pub hits: Vec<String>,
    pub score: F64,
    pub aligned_query_region: AlignedQueryRegion,
}

/// An instance of AlignedQueryRegion informs about the region, or sub-sequence, of the query to
/// which all a cluster's hits' pairwise local sequence align. A region can be "strict" in the
/// sense, that _all_ hits' align to this region, or "relaxed", i.e. non-strict, if just a subset
/// of the hits align to parts of the region. See the following graphical examples:
///
/// * strict (`all_hits_overlap == true`)
///
/// query  =============
/// hit_1     ----
/// hit_2      ----
/// hit_3    ----
/// region     **
///
/// * relaxed (`all_hits_overlap == false`)
///
/// query  =============
/// hit_1    ----
/// hit_2        ---
/// hit_3      ---
/// region   *******
#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct AlignedQueryRegion {
    /// The position in the query sequence the region starts
    pub qstart: u32,
    /// The position in the query sequence the region ends
    pub qend: u32,
    /// Is the region a strict overlap, i.e. do all hits have pairwise alignment in this region?
    pub all_hits_overlap: bool,
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
