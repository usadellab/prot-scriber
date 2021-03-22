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

use super::default::*;
use super::models::*;
use ndarray::Array2;
use std::collections::HashSet;

/// Function clusters a query and its hits to find the cluster of which the query is member of and
/// use that as a basis to generate a short human readable protein function description.
///
/// # Arguments
///
/// * query - The query including its hits to be subjected to clustering
pub fn cluster(query: &mut Query) -> Array2<f64> {
    let m = query.to_similarity_matrix();
    println!("Query '{}' yields similarity matrix:\n{:?}", query.id, m);
    // To Do: Use user arguments and the defaults only if no user args have been supplied:
    markov_cluster(
        &m.1,
        &(*MCL_INFLATION),
        &(*MCL_DELTA),
        &(*MCL_MAX_ITERATIONS),
        &(*ROUND_DECIMAL_DIGITS),
    )
}

/// Function extracts the clusters from the result matrix produced by markov clusterin (MCL).
/// Returns a vector of HashSets of node indices (`usize`) starting from zero. In short, the result
/// is a vector of mathematic sets in which each set represents a cluster. A cluster indicates
/// contained nodes by their indices.
///
/// # Arguments
///
/// * mcl_matrix - The result matrix produced by markov clustering
pub fn get_clusters(mcl_matrix: &Array2<f64>) -> Vec<HashSet<usize>> {
    let mut clusters: Vec<HashSet<usize>> = Vec::new();
    // Iterate over rows. Create a vector of clusters (one per row), where each cluster is a
    // HashSet of node-indices, indicating nodes that have edges between them.
    for row_i in 0..mcl_matrix.shape()[0] {
        let mut node_i_cluster: HashSet<usize> = HashSet::new();
        node_i_cluster.insert(row_i);
        for (cell_i, &cell) in mcl_matrix.row(row_i).indexed_iter() {
            if cell > 0.0 {
                node_i_cluster.insert(cell_i);
            }
        }
        clusters.push(node_i_cluster);
    }
    // Join clusters by iterating over node-indices and joining (set-union) all clusters that
    // contain the current node-index:
    for node_i in 0..mcl_matrix.shape()[0] {
        let clusters_clone = clusters.clone();
        let mut cluster_inds = clusters_clone
            .iter()
            .enumerate()
            .filter(|(_, v)| v.contains(&node_i)) // cluster has current node i?
            .map(|(i, _)| i)
            .collect::<Vec<usize>>();
        // Sort descending by comparing not a with b, but inverted argument order b.cmp(a).
        // Note, that this sort order is required to enable iterative removal of vector elements
        // without affecting lower indices.
        cluster_inds.sort_by(|a, b| b.cmp(a));
        let mut joined_cluster: HashSet<usize> = HashSet::new();
        for i in cluster_inds {
            let cluster_i = clusters_clone.get(i).unwrap();
            // join clusters that contain node i:
            joined_cluster.extend(cluster_i);
            // and remove the just joined cluster:
            clusters.remove(i);
        }
        // Replace the clusters that contained node i with the set-union of them:
        clusters.push(joined_cluster);
    }
    clusters // return
}

/// Normalizes a square two dimensional matrix so that all rows sum up to one, and thus a row's
/// cells can be interpreted as probabilities. See Markov Clustering for more details on this step.
///
/// # Arguments
///
/// * mtrx - The square two-dimensional matrix to normalize.
pub fn normalize(mtrx: &Array2<f64>) -> Array2<f64> {
    let mut norm_mtrx = mtrx.clone();
    for i in 0..mtrx.shape()[0] {
        let row_sum = mtrx.row(i).scalar_sum();
        for k in 0..mtrx.shape()[1] {
            norm_mtrx[[i, k]] = mtrx[[i, k]] / row_sum;
        }
    }
    norm_mtrx
}

/// Rounds the cells of a two dimensional matrix of `f64` values. Each value is rounded to the
/// argument number of digits.
///
/// # Arguments
///
/// * mtrx - The two dimensional matrix of floats to be rounded
/// * n_digits - The number of digits to round the values to
pub fn round_matrix(mtrx: &Array2<f64>, n_digits: &i32) -> Array2<f64> {
    let fctr = 10f64.powi(*n_digits);
    mtrx.map(|x: &f64| (x * fctr).round() / fctr)
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
/// * stochastic_matrix - The two dimensional stochastic matrix to cluster
/// * inflation - The inflation parameter (I), should range from 0.0 < I <= 5.0
/// * delta - The maximum numeric cell-wise difference between the last iteration and the current
///           that is considered still worth continueing the clustering.
/// * n - The number of the current iteration
/// * max_iter - The maximum number of interation steps to carry out
pub fn mcl(
    stochastic_matrix: &Array2<f64>,
    inflation: &f64,
    delta: &f64,
    n: i8,
    max_iter: &i8,
) -> Array2<f64> {
    // expansion: matrix multipliction (squaring):
    let mut m = stochastic_matrix.dot(stochastic_matrix);
    // inflation: hadamard power of the matrix
    m = m.map(|x: &f64| x.powf(*inflation));
    // normalize: each row should sum up to one so that the cells represent probabilities:
    m = normalize(&m);
    // estimate max difference between the input and processed matrices:
    let max_abs_diff = (&m - stochastic_matrix).iter().fold(0.0, |acc, cur| {
        let abs_cur_diff = cur.abs();
        if acc < abs_cur_diff {
            abs_cur_diff
        } else {
            acc
        }
    });
    // Have we reached the max number of iterations or hasn't this iteration yielded significant
    // differences in the probabilities?
    if n == *max_iter || max_abs_diff <= *delta {
        m
    } else {
        // Do another iteration recursively:
        mcl(&m, inflation, delta, n + 1, max_iter)
    }
}

/// Within each row _i_ of the stochastic quadratic matrix `stochastic_matrix` the numeric maximum
/// is identified and set as the probability to move from the current node _i_ back to itself, thus
/// setting self-loops to max likelihood. Why this is done is explained in the publications
/// concerning the Markov Clustering algorithm. Returns a new matrix with added self-loops which
/// still requires normalization before being a stochastic matrix again.
///
/// # Arguments
///
/// * `stochastic_matrix: &Array2<f64>` - A reference to the stochastic matrix to add self loops
pub fn add_self_loops(stochastic_matrix: &Array2<f64>) -> Array2<f64> {
    let mut loops_added = stochastic_matrix.clone();
    for row_i in 0..loops_added.shape()[0] {
        let row_max = stochastic_matrix
            .row(row_i)
            .fold(0.0, |acc, curr| if *curr > acc { *curr } else { acc });
        loops_added[[row_i, row_i]] = row_max;
    }
    loops_added
}

/// Wrapper function around the recursive implementation of the Markov Clustering algorithm (see
/// function `mcl` for more details).
///
/// # Arguments
///
/// * stochastic_matrix - The two dimensional stochastic matrix to cluster
/// * inflation - The inflation parameter (I), should range from 1.0 < I <= 5.0
/// * delta - The maximum numeric cell-wise difference between the last iteration and the current
///           that is considered still worth continueing the clustering.
/// * max_iter - The maximum number of interation steps to carry out
/// * round_digits - The number of digits to round the resulting clustered two dimensional
///                  stochastic adjacency matrix to. Set to negative value, if no rounding is
///                  wanted. Typical value should range between one and four.
pub fn markov_cluster(
    stochastic_matrix: &Array2<f64>,
    inflation: &f64,
    delta: &f64,
    max_iter: &i8,
    round_digits: &i32,
) -> Array2<f64> {
    let looped_normalized_mtrx = normalize(&add_self_loops(stochastic_matrix));
    let clustered_mtrx = mcl(&looped_normalized_mtrx, inflation, delta, 0, max_iter);
    if *round_digits > -1 {
        round_matrix(&clustered_mtrx, round_digits)
    } else {
        clustered_mtrx
    }
}

/// Unit tests
#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr2;

    #[test]
    fn test_add_self_loops() {
        let dm = arr2(&[
            [0.0, 0.65, 0.25, 0.05, 0.05],
            [0.65, 0.0, 0.25, 0.05, 0.05],
            [0.5, 0.5, 0.0, 0.0, 0.0],
            [0.05, 0.05, 0.0, 0.0, 0.9],
            [0.05, 0.05, 0.0, 0.9, 0.0],
        ]);
        let looped_m = add_self_loops(&dm);
        //println!(
        //    "Added self-loops to matrix. Input is:\n{:?}\n.Result is:\n{:?}",
        //    dm, looped_m
        //);
        let expected = arr2(&[
            [0.65, 0.65, 0.25, 0.05, 0.05],
            [0.65, 0.65, 0.25, 0.05, 0.05],
            [0.5, 0.5, 0.5, 0.0, 0.0],
            [0.05, 0.05, 0.0, 0.9, 0.9],
            [0.05, 0.05, 0.0, 0.9, 0.9],
        ]);
        // Test diagonal
        assert_eq!(looped_m[[0, 0]], 0.65);
        assert_eq!(looped_m[[1, 1]], 0.65);
        assert_eq!(looped_m[[2, 2]], 0.50);
        assert_eq!(looped_m[[3, 3]], 0.90);
        assert_eq!(looped_m[[4, 4]], 0.90);
        // And test the whole thing - I know that in this case we do not need the above tests,
        // because they are included in the following. We just like to waste electrical energy and
        // heat the planet. Feeling cold! Sorry!
        assert_eq!(looped_m, expected);
    }

    #[test]
    fn test_mcl() {
        let dm = arr2(&[
            [0.0, 0.65, 0.25, 0.05, 0.05],
            [0.65, 0.0, 0.25, 0.05, 0.05],
            [0.5, 0.5, 0.0, 0.0, 0.0],
            [0.05, 0.05, 0.0, 0.0, 0.9],
            [0.05, 0.05, 0.0, 0.9, 0.0],
        ]);
        let cm = round_matrix(&mcl(&dm, &5.0, &0.0001, 0, &10), &4i32);
        // println!(
        // "Testing mcl. - Without adding self-loops! - Input matrix for nodes 1-5 is:\n{:?}\n, and clustered matrix is:\n{:?}",
        // dm, cm
        // );
        let expected = arr2(&[
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);
        assert_eq!(cm, expected);
    }

    #[test]
    fn test_markov_cluster() {
        let dm = arr2(&[
            [0.0, 0.65, 0.25, 0.05, 0.05],
            [0.65, 0.0, 0.25, 0.05, 0.05],
            [0.5, 0.5, 0.0, 0.0, 0.0],
            [0.05, 0.05, 0.0, 0.0, 0.9],
            [0.05, 0.05, 0.0, 0.9, 0.0],
        ]);
        let cm = markov_cluster(&dm, &5.0, &0.0001, &10, &4);
        // println!("Testing markov clustering. Input matrix for nodes 1-5 is:\n{:?}\n, and clustered matrix is:\n{:?}",
        //    dm, cm);
        let expected = arr2(&[
            [0.5, 0.5, 0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5, 0.5],
            [0.0, 0.0, 0.0, 0.5, 0.5],
        ]);
        assert_eq!(cm, expected);
    }

    #[test]
    fn test_get_clusters() {
        let mcl_matrix = arr2(&[
            [0.0, 0.5, 0.5, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
        ]);
        let clusters = get_clusters(&mcl_matrix);
        // println!(
        //     "Testing get_clusters from the result of Markov Clustering.\nInput is:\n{:?}\nClusters are:\n{:?}\n",
        //     &mcl_matrix, &clusters
        // );
        let cluster_0: HashSet<usize> = vec![0, 1, 2].iter().cloned().collect();
        let cluster_1: HashSet<usize> = vec![3, 4].iter().cloned().collect();
        assert!(clusters.contains(&cluster_0));
        assert!(clusters.contains(&cluster_1));
    }
}
