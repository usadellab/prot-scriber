//! Implementation of a two dimensional triangular square matrix
use eq_float::F64;
use std::fmt;

/// Representation of a triangular square two dimensional matrix. Triangular means that the
/// diagonal mirrors cells, i.e. matrix cells [i,k] are identical to [k,i]. _Note_ that cells are
/// stored in the "C"-style manner, i.e. by row, and not in the "R" or "Fortran" style by column.
#[derive(PartialEq, Eq, Clone, Default)]
pub struct Triang2dMatrix {
    pub cells: Vec<F64>,
    pub row_n_col_names: Vec<String>,
}

impl Triang2dMatrix {
    /// Constructor to generate a new instance of Triang2dMatrix.
    ///
    /// # Arguments
    ///
    /// * `cells: Vec<F64>` - The content of the cells. Only provide the upper triangle, row by
    ///                       row.
    /// * `row_n_col_names: Vec<String>` - The row and column names. Because this is a square
    ///                                    matrix and it is supposed to be used as an adjacency
    ///                                    matrix the row and column names are _identical_.
    pub fn new(cells: Vec<F64>, row_n_col_names: Vec<String>) -> Triang2dMatrix {
        let axis_len = row_n_col_names.len() as i32;
        let expected_n_cells = axis_len.pow(2) - (1..axis_len as i32).sum::<i32>();
        if cells.len() as i32 != expected_n_cells {
            panic!("Dimensions do not match! Provided {} cell values and {} row and column names for the construction of a triangular square matrix.", cells.len(), axis_len);
        }
        Triang2dMatrix {
            cells,
            row_n_col_names,
        }
    }

    /// Returns the number of cells per axis, i.e. number of cells per row and column. Remember,
    /// this is a square matrix.
    ///
    /// # Arguments
    ///
    /// * `&self` - The instance of Triang2dMatrix
    pub fn axis_len(&self) -> usize {
        self.row_n_col_names.len()
    }

    /// Returns a reference to the value indicated by the argument coordinates. Note, that
    /// coordinates start at zero. Internally two dimensional coordinates get transformed to one
    /// dimensional indices of field `self.cells`. The cell index (`i`) is calculated by adding
    /// argument `col` to an `offset` depending on the argument `row`, in which the offset is
    /// computes as the sum over index `k` starting from 1 until and including `row` of the terms
    /// `self.axis_len() - k`.
    ///
    /// # Arguments
    ///
    /// * `&self` - the instance of Triang2dMatrix
    /// * `row` - The row coordinate starting with zero
    /// * `col` - The col coordinate starting with zero
    pub fn get(&self, row: usize, col: usize) -> &F64 {
        if row >= self.axis_len() || col >= self.axis_len() {
            panic!("Cell index (row {}, col {}) is out of bounds. Matrix has {} rows and cols. Note that row and col indices start at 0.",
                row, col, self.axis_len());
        }
        if row > col {
            // This implementation of a triangular matrix only stores cells of the upper triangle.
            self.get(col, row)
        } else {
            let offset: i32 = (1..=(row as i32)).map(|i| self.axis_len() as i32 - i).sum();
            let i = offset as usize + col;
            // println!("get({}, {}) use offset {} and i {}", row, col, offset, i);
            &self.cells[i]
        }
    }

    /// Performs matrix multipliction ("mm") with itself, i.e. `&self` "mm" `&self` and returns a
    /// new instance of Triang2dMatrix as the result.
    ///
    /// # Arguments
    ///
    /// * `&self` - the instance of Triang2dMatrix
    pub fn square(&self) -> Triang2dMatrix {
        let mut sqrd_cells: Vec<F64> = Vec::with_capacity(self.axis_len());
        for i in 0..self.axis_len() {
            for k in i..self.axis_len() {
                let cell: F64 = F64((0..self.axis_len())
                    .map(|x| self.get(i, x).0 * self.get(x, k).0)
                    .sum());
                sqrd_cells.push(cell);
            }
        }
        Triang2dMatrix::new(sqrd_cells, self.row_n_col_names.clone())
    }

    /// Raises each cell (`F64`) to the power of argument `p` and returns a new instance of
    /// Triang2dMatrix as result.
    ///
    /// # Arguments
    ///
    /// * `&self` - the instance of Triang2dMatrix
    /// * `p: &f64` - a reference to the argument power to which to raise each element of `&self`
    ///               to.
    pub fn hadamard_raise_to_the_power(&self, p: &f64) -> Triang2dMatrix {
        let raised_cells = self.cells.iter().map(|x| F64(x.0.powf(*p))).collect();
        Triang2dMatrix::new(raised_cells, self.row_n_col_names.clone())
    }

    /// Normalizes the square triangular matrix so that each row sums up to unit value (`1.0`).
    /// Thus the cells (i,k) can be seen as probabilities to move from node i to node k in a random
    /// walk, if the matrix is interpreted as an adjacency matrix. A new instance of Triang2dMatrix
    /// is returned as result.
    ///
    /// # Arguments
    ///
    /// * `&self` - the instance of Triang2dMatrix
    pub fn stochastic_normalize(&self) -> Triang2dMatrix {
        let mut stoch_norm_cells: Vec<F64> = Vec::new();
        for i in 0..self.axis_len() {
            // The Isley Brothers - Summer Breeze
            // is a great song! Just and simply great!
            let row_summer =
                (0..self.axis_len()).fold(0.0, |acc, iter_col| acc + self.get(i, iter_col).0);
            for k in i..self.axis_len() {
                stoch_norm_cells.push(F64(self.get(i, k).0 / row_summer));
            }
        }
        Triang2dMatrix::new(stoch_norm_cells, self.row_n_col_names.clone())
    }

    /// Comutes the maximum absolute difference between cells of argument `&self` and `to`
    /// comparing cells of the same index only. Returns the difference as `f64`.
    ///
    /// # Arguments
    ///
    /// * `&self` - the instance of Triang2dMatrix
    /// * `to: &Triang2dMatrix` - the other instance of Triang2dMatrix to compute the max absolute
    ///                           cellwise difference to.
    pub fn max_abs_cellwise_difference(&self, to: &Triang2dMatrix) -> f64 {
        if self.axis_len() != to.axis_len() {
            panic!("Cannot compute cellwise absolute differences on matrices of different sizes.");
        }
        (0..self.cells.len())
            .map(|i| F64((self.cells[i].0 - to.cells[i].0).abs()))
            .max()
            .unwrap()
            .0
    }

    /// Rounds the cells of a two dimensional matrix of `f64` values. Each value is rounded to the
    /// argument number of digits.
    ///
    /// # Arguments
    ///
    /// * `&self` - the instance of Triang2dMatrix
    /// * `n_digits: &i32` - The number of digits to round the `F64` cells to
    pub fn round(&self, n_digits: &i32) -> Triang2dMatrix {
        let fctr = 10f64.powi(*n_digits);
        Triang2dMatrix::new(
            self.cells
                .iter()
                .map(|x| F64((x.0 * fctr).round() / fctr))
                .collect(),
            self.row_n_col_names.clone(),
        )
    }
}

/// Triang2dMatrix derives from Debug
impl fmt::Debug for Triang2dMatrix {
    /// Manual implementation of Debug formatting displays the triangular square stoachastic matrix
    /// in a grid to ease human readibility.
    ///
    /// # Arguments
    ///
    /// * `&self` - the instance of Triang2dMatrix
    /// * `f: &mut fmt::Formatter<'_>` - the formatter instance used to construct the string
    ///                                  representation of `&self`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut cells_str = "\n".to_string();
        for row_i in 0..self.axis_len() {
            for col_k in 0..self.axis_len() {
                cells_str.push_str(
                    &format!(" {cell:^6.*} ", 4, cell = self.get(row_i, col_k)).to_string(),
                );
            }
            cells_str.push_str("\n");
        }
        f.debug_struct("Triang2dMatrix")
            .field("row_n_col_names", &self.row_n_col_names)
            .field("cells", &format_args!("{}", cells_str))
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_validates_dimensions() {
        let c = vec![F64(1.0), F64(2.0), F64(3.0), F64(4.0), F64(5.0), F64(6.0)];
        let n = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let m1 = Triang2dMatrix::new(c, n);
        assert_eq!(
            m1.cells,
            vec![F64(1.0), F64(2.0), F64(3.0), F64(4.0), F64(5.0), F64(6.0)]
        );
        assert_eq!(
            m1.row_n_col_names,
            vec!["a".to_string(), "b".to_string(), "c".to_string()]
        );
        assert_eq!(m1.axis_len(), 3);
    }

    #[test]
    #[should_panic]
    fn new_panics_with_wrong_dimensions() {
        let c = vec![F64(1.0), F64(2.0), F64(3.0)];
        let n = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        Triang2dMatrix::new(c, n);
    }

    #[test]
    #[should_panic]
    fn get_panics_with_invalid_indices() {
        let c = vec![F64(1.0), F64(2.0), F64(3.0), F64(4.0), F64(5.0), F64(6.0)];
        let n = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let m = Triang2dMatrix::new(c, n);
        m.get(3, 3);
    }

    #[test]
    #[should_panic]
    fn get_panics_with_invalid_row_index() {
        let c = vec![F64(1.0), F64(2.0), F64(3.0), F64(4.0), F64(5.0), F64(6.0)];
        let n = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let m = Triang2dMatrix::new(c, n);
        m.get(3, 0);
    }

    #[test]
    #[should_panic]
    fn get_panics_with_invalid_col_index() {
        let c = vec![F64(1.0), F64(2.0), F64(3.0), F64(4.0), F64(5.0), F64(6.0)];
        let n = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let m = Triang2dMatrix::new(c, n);
        println!("{:?}", m);
        m.get(0, 3);
    }

    #[test]
    fn get_correct_cells() {
        let c = vec![
            F64(0.0),
            F64(1.0),
            F64(2.0),
            F64(3.0),
            F64(4.0),
            F64(5.0),
            F64(6.0),
            F64(7.0),
            F64(8.0),
            F64(9.0),
        ];
        let n = vec![
            "a".to_string(),
            "b".to_string(),
            "c".to_string(),
            "d".to_string(),
        ];
        let m = Triang2dMatrix::new(c, n);
        //println!("{:?}", m);
        assert_eq!(m.get(0, 0).0, 0.0);
        assert_eq!(m.get(0, 1).0, 1.0);
        assert_eq!(m.get(0, 2).0, 2.0);
        assert_eq!(m.get(0, 3).0, 3.0);
        assert_eq!(m.get(1, 0).0, 1.0);
        assert_eq!(m.get(1, 1).0, 4.0);
        assert_eq!(m.get(1, 2).0, 5.0);
        assert_eq!(m.get(1, 3).0, 6.0);
        assert_eq!(m.get(2, 0).0, 2.0);
        assert_eq!(m.get(2, 1).0, 5.0);
        assert_eq!(m.get(2, 2).0, 7.0);
        assert_eq!(m.get(2, 3).0, 8.0);
        assert_eq!(m.get(3, 0).0, 3.0);
        assert_eq!(m.get(3, 1).0, 6.0);
        assert_eq!(m.get(3, 2).0, 8.0);
        assert_eq!(m.get(3, 3).0, 9.0);
    }

    #[test]
    fn square_is_matrix_multiplication_with_self() {
        let c = vec![F64(1.0), F64(2.0), F64(3.0), F64(4.0), F64(5.0), F64(6.0)];
        let n = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let m = Triang2dMatrix::new(c, n);
        let sqrd_m = m.square();
        assert_eq!(sqrd_m.axis_len(), 3);
        assert_eq!(
            sqrd_m.cells,
            vec![
                F64(14.0),
                F64(25.0),
                F64(31.0),
                F64(45.0),
                F64(56.0),
                F64(70.0)
            ]
        );
    }

    #[test]
    fn test_hadamard_raise_to_the_power() {
        let c = vec![F64(1.0), F64(2.0), F64(3.0), F64(4.0), F64(5.0), F64(6.0)];
        let n = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let m = Triang2dMatrix::new(c, n);
        let raised_m = m.hadamard_raise_to_the_power(&2.0);
        assert_eq!(raised_m.axis_len(), 3);
        assert_eq!(
            raised_m.cells,
            vec![
                F64(1.0),
                F64(4.0),
                F64(9.0),
                F64(16.0),
                F64(25.0),
                F64(36.0)
            ]
        )
    }

    #[test]
    fn test_stochastic_normalize() {
        let c = vec![F64(1.0), F64(2.0), F64(3.0), F64(4.0), F64(5.0), F64(6.0)];
        let n = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let m = Triang2dMatrix::new(c, n);
        let stoch_norm = m.stochastic_normalize();
        assert_eq!(
            stoch_norm.cells,
            vec![
                F64(1.0 / 6.0),
                F64(2.0 / 6.0),
                F64(3.0 / 6.0),
                F64(4.0 / 11.0),
                F64(5.0 / 11.0),
                F64(6.0 / 14.0)
            ]
        );
    }

    #[test]
    fn test_max_absolute_difference() {
        let c1 = vec![F64(1.0), F64(2.0), F64(3.0), F64(4.0), F64(5.0), F64(6.0)];
        let n1 = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let m1 = Triang2dMatrix::new(c1, n1);
        let c2 = vec![F64(2.0), F64(1.0), F64(4.0), F64(3.0), F64(6.0), F64(8.0)];
        let n2 = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let m2 = Triang2dMatrix::new(c2, n2);
        assert_eq!(m1.max_abs_cellwise_difference(&m2), 2.0);
    }
}
