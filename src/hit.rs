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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_hit_from_strs() {
        let h1 = Hit::new(
            "Hit_One", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        assert_eq!(h1.id, String::from("Hit_One"));
        assert_eq!(h1.bitscore, F64(123.4));
        assert_eq!(
            h1.description,
            "Probable LRR receptor-like serine/threonine-protein kinase At3g47570"
        );
    }
}
