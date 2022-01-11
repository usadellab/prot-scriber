use super::default::FILTER_REGEXS;
use super::model_funcs::filter_stitle;
use eq_float::F64;

/// A sequence similarity search result is called a "Hit" and is represented by its namesake.
#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Hit {
    pub id: String,
    pub bitscore: F64,
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
    /// * `bitscore: &str` - The bitscore of the local alignment
    /// * `stitle: &str` - The full title line of the hit in the reference database (stored in
    ///                    fasta format)
    pub fn new(id: &str, bitscore: &str, stitle: &str) -> Hit {
        Hit {
            id: String::from(id),
            bitscore: F64(bitscore.parse().unwrap()),
            description: filter_stitle(&stitle, &(*FILTER_REGEXS)),
        }
    }

    /// Parses a line in the respective sequence similarity search result table. The line is already
    /// split into cells (columns). Returns a new instance of `Hit`.
    ///
    /// # Arguments
    ///
    /// * `record_cols: &Vec<String>` - The record parsed from splitting the respective line
    /// * `sacc_col_ind: &usize` - The index of argument `record_cols` where the `sacc` is stored.
    /// * `bitscore_col_ind: &usize` - The index of argument `record_cols` where the `bitscore` is
    ///                               stored.
    /// * `stitle_col_ind: &usize` - The index of argument `record_cols` where the `stitle` is stored.
    pub fn parse_hit(
        record_cols: &Vec<String>,
        sacc_col_ind: &usize,
        bitscore_col_ind: &usize,
        stitle_col_ind: &usize,
    ) -> Hit {
        Hit::new(
            &record_cols[*sacc_col_ind],
            &record_cols[*bitscore_col_ind],
            &record_cols[*stitle_col_ind],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_hit_from_strs() {
        let h1 = Hit::new(
            "hit_One", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        assert_eq!(h1.id, String::from("hit_One"));
        assert_eq!(h1.bitscore, F64(123.4));
        assert_eq!(
            h1.description,
            "lrr receptor serine/threonine-protein kinase at3g47570"
        );
    }

    #[test]
    fn parses_hit_from_record_line() {
        let record_cols = vec![
            "Soltu.DM.02G015700.1".to_string(),
            "sp|C0LGP4|Y3475_ARATH".to_string(),
            "2209".to_string(),
            "1284".to_string(),
            "2199".to_string(),
            "1010".to_string(),
            "64".to_string(),
            "998".to_string(),
            "580".to_string(),
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1".to_string()
        ];
        let hit = Hit::parse_hit(&record_cols, &1, &8, &9);
        let expected = Hit::new( "sp|C0LGP4|Y3475_ARATH", "580",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1");
        assert_eq!(hit, expected);
    }
}
