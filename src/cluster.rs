use eq_float::F64;
use std::cmp::{max, min, Ordering};

/// A cluster of sequence similarity search hits is represented by this struct.
#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Cluster {
    pub hits: Vec<String>,
    pub score: F64,
    pub aligned_query_region: AlignedQueryRegion,
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

impl AlignedQueryRegion {
    /// Checks whether there is an overlap between two instances of `AlignedQueryRegion`; if so
    /// returns `false`. In case of _no_ overlap the function returns `true`.
    ///
    /// # Arguments
    ///
    /// * `&self` - A reference to an instance of AlignedQueryRegion
    /// * `other: &AlignedQueryRegion` - A reference to the "other" instance of AlignedQueryRegion
    ///                                  for which to check whether it aligns at least partially to
    ///                                  the same query region, i.e. sub-sequence.
    pub fn is_disjoint(&self, other: &AlignedQueryRegion) -> bool {
        let max_qstart = max(self.qstart, other.qstart);
        let min_qend = min(self.qend, other.qend);
        // Note that we need to convert to i64 to allow for negative subtraction results, otherwise
        // we can cause an overflow panic:
        (min_qend as i64 - max_qstart as i64) < 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aligned_query_region_is_disjoint() {
        let aqr_1 = AlignedQueryRegion {
            qstart: 1,
            qend: 10,
            all_hits_overlap: true,
        };
        let aqr_2 = AlignedQueryRegion {
            qstart: 11,
            qend: 20,
            all_hits_overlap: true,
        };
        let aqr_1_2_disjoint = aqr_1.is_disjoint(&aqr_2);
        let aqr_2_1_disjoint = aqr_2.is_disjoint(&aqr_1);
        assert_eq!(aqr_1_2_disjoint, aqr_2_1_disjoint);
        assert!(aqr_1_2_disjoint);
        let aqr_3 = AlignedQueryRegion {
            qstart: 5,
            qend: 8,
            all_hits_overlap: true,
        };
        let aqr_1_3_disjoint = aqr_1.is_disjoint(&aqr_3);
        let aqr_3_1_disjoint = aqr_3.is_disjoint(&aqr_1);
        assert_eq!(aqr_1_3_disjoint, aqr_3_1_disjoint);
        assert!(!aqr_1_3_disjoint);
        let aqr_4 = AlignedQueryRegion {
            qstart: 10,
            qend: 18,
            all_hits_overlap: true,
        };
        let aqr_1_4_disjoint = aqr_1.is_disjoint(&aqr_4);
        let aqr_4_1_disjoint = aqr_4.is_disjoint(&aqr_1);
        assert_eq!(aqr_1_4_disjoint, aqr_4_1_disjoint);
        assert!(!aqr_1_4_disjoint);
        let aqr_5 = AlignedQueryRegion {
            qstart: 1,
            qend: 1,
            all_hits_overlap: true,
        };
        let aqr_1_5_disjoint = aqr_1.is_disjoint(&aqr_5);
        let aqr_5_1_disjoint = aqr_5.is_disjoint(&aqr_1);
        assert_eq!(aqr_1_5_disjoint, aqr_5_1_disjoint);
        assert!(!aqr_1_5_disjoint);
    }
}
