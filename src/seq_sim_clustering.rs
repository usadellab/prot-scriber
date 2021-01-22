use super::models::Query;
use ndarray::{arr2, Array2};

pub fn cluster_hits(query: &Query) -> Array2<f64> {
    // dummy
    arr2(&[[6f64, 5f64, 4f64], [3f64, 2f64, 1f64]])
}
