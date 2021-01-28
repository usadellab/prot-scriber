//! `seq_sim_clustering` implements the popular Markov Clustering algorithm (`mcl`) to cluster a
//! Query and its Hits so that from the cluster the query is member of a short human readable
//! protein description can be generated. For details on `mcl` see https://micans.org/mcl/ or
//!
//! * Stijn van Dongen, Graph Clustering by Flow Simulation, PhD thesis, University of Utrecht, May
//! 2000.  ( http://www.library.uu.nl/digiarchief/dip/diss/1895620/inhoud.htm )
//!
//! * Stijn van Dongen, A cluster algorithm for graphs, Technical Report INS-R0010, National
//! Research Institute for Mathematics and Computer Science in the Netherlands, Amsterdam, May
//! 2000.  ( http://www.cwi.nl/ftp/CWIreports/INS/INS-R0010.ps.Z )
//!
//! * Stijn van Dongen, Graph clustering via a discrete uncoupling process, Siam Journal on Matrix
//! Analysis and Applications 30-1, p121-141, 2008.
//!
//! * Stijn van Dongen and Cei Abreu-Goodger, Using MCL to extract clusters from networks, in
//! Bacterial Molecular Networks: Methods and Protocols, Methods in Molecular Biology, Vol 804,
//! pages 281—295 (2012). PMID 22144159.

use super::matrix::Triang2dMatrix;
use super::models::Query;
use eq_float::F64;

/// Function clusters a query and its hits to find the cluster of which the query is member of and
/// use that as a basis to generate a short human readable protein function description.
///
/// # Arguments
///
/// * query - The query including its hits to be subjected to clustering
pub fn cluster_hits(query: &Query) -> Triang2dMatrix {
    // To Do !!!
    Triang2dMatrix::new(
        vec![F64(1.0), F64(2.0), F64(3.0), F64(4.0), F64(5.0), F64(6.0)],
        vec!["A".to_string(), "B".to_string(), "C".to_string()],
    )
}

/// Clusters a two dimensional square stochastic adjacency matrix. "stochastic" means that the
/// cells contain probabilities of moving from one node to another in a random walk. Also each row
/// sums up to one. Clustering is done by simulation of random walks. In this, use matrix
/// multiplication, then amplify the signal by taking hadamard power with the inflation parameter,
/// followed by normalization, so that each row sums up to one, and thus cells represent true
/// probabilities.  See Markov Clustering for more details.
///
/// # Arguments
///
/// * distance_matrix - The two dimensional stochastic matrix to cluster
/// * inflation - The inflation parameter (I), should range from 0.0 < I <= 5.0
/// * delta - The maximum numeric cell-wise difference between the last iteration and the current
///           that is considered still worth continueing the clustering.
/// * n - The number of the current iteration
/// * max_iter - The maximum number of interation steps to carry out
pub fn mcl(
    distance_matrix: &Triang2dMatrix,
    inflation: &f64,
    delta: &f64,
    n: i8,
    max_iter: &i8,
) -> Triang2dMatrix {
    // expansion: matrix multipliction (squaring):
    let mut m = distance_matrix.square();
    // inflation: hadamard power of the matrix
    m = m.hadamard_raise_to_the_power(&inflation);
    // normalize: each row should sum up to one so that the cells represent probabilities:
    m = m.stochastic_normalize();
    // estimate max difference between the input and processed matrices:
    let max_abs_diff = distance_matrix.max_abs_cellwise_difference(&m);
    // Have we reached the max number of iterations or hasn't this iteration yielded significant
    // differences in the probabilities?
    if n == *max_iter || max_abs_diff <= *delta {
        m
    } else {
        // Do another iteration recursively:
        mcl(&m, inflation, delta, n + 1, max_iter)
    }
}

/// Wrapper function around the recursive implementation of the Markov Clustering algorithm (see
/// function `mcl` for more details).
///
/// # Arguments
///
/// * distance_matrix - The two dimensional stochastic matrix to cluster
/// * inflation - The inflation parameter (I), should range from 0.0 < I <= 5.0
/// * delta - The maximum numeric cell-wise difference between the last iteration and the current
///           that is considered still worth continueing the clustering.
/// * max_iter - The maximum number of interation steps to carry out
/// * round_digits - The number of digits to round the resulting clustered two dimensional
///                  stochastic adjacency matrix to. Set to negative value, if no rounding is
///                  wanted. Typical value should range between one and four.
pub fn markov_cluster(
    distance_matrix: &Triang2dMatrix,
    inflation: &f64,
    delta: &f64,
    max_iter: &i8,
    round_digits: &i32,
) -> Triang2dMatrix {
    let clustered_mtrx = mcl(distance_matrix, inflation, delta, 0, max_iter);
    if *round_digits > -1 {
        clustered_mtrx.round(round_digits)
    } else {
        clustered_mtrx
    }
}

/// Unit tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mcl() {
        let dm = Triang2dMatrix::new(
            vec![
                F64(0.0),
                F64(0.65),
                F64(0.25),
                F64(0.05),
                F64(0.05),
                F64(0.0),
                F64(0.25),
                F64(0.05),
                F64(0.05),
                F64(0.0),
                F64(0.0),
                F64(0.0),
                F64(0.0),
                F64(0.9),
                F64(0.0),
            ],
            vec![
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "D".to_string(),
                "E".to_string(),
            ],
        );
        let cm = mcl(&dm, &5.0, &0.001, 0, &10).round(&4i32);
        println!("Testing markov clustering. Input matrix for nodes 1-5 is:\n{:?}\n, and clustered matrix is:\n{:?}",
           dm, cm);
        let expected = Triang2dMatrix::new(
            vec![
                F64(1.0),
                F64(0.0),
                F64(0.0),
                F64(0.0),
                F64(0.0),
                F64(1.0),
                F64(0.0),
                F64(0.0),
                F64(0.0),
                F64(1.0),
                F64(0.0),
                F64(0.0),
                F64(1.0),
                F64(0.0),
                F64(1.0),
            ],
            vec![
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "D".to_string(),
                "E".to_string(),
            ],
        );
        assert_eq!(cm, expected);
    }
}
