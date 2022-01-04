use eq_float::F64;
use regex::Regex;

/// Applies a list of regular expressions to a string. If any of them match returns true,
/// otherwise returns false.
/// This function is used to exclude non-informative descriptions from a human readable
/// sequence title and to check for non-informative words to be excluded from scoring
/// 
/// # Arguments
/// 
/// * stitle - The sequence title string
/// * regexs - A vector of regular expression to be applied to the stitle argument.
pub fn matches_blacklist(stitle: &str, regexs: &Vec<Regex>) -> bool {
    regexs
        .iter()
        .any(|x| x.is_match(&stitle.to_string()))
}

/// Fasta entries have a long title in which the sequence identifier and often taxonomic
/// information is given along with a short human readable protein description. We are only
/// interested in the latter. This function extracts the short description using regular
/// expressions.
///
/// # Arguments
///
/// * stitle - The sequence title line as found in the original Fasta file.
/// * regexs - A vector of regular expressions to be applied in series to the argument stitle to
///            extract the desired short description.
pub fn filter_stitle(stitle: &str, regexs: &Vec<Regex>) -> String {
    String::from(
        regexs
            .iter()
            .fold(stitle.to_string(), |accumulated, current| {
                current.replace_all(&accumulated, "").to_string()
            })
            .trim(),
    )
}

/// Calculates the overlap between a query and one of its hits as produced by sequence similarity
/// searches (e.g. Blast or Diamond). These search algorithms produce local alignments and the
/// arguments to this function. Overlap is calculated as
///
/// `((qend - qstart + 1) + (send - sstart + 1))`
///
/// # Arguments
///
/// * sstart - The start position of the local alignment in the hit
/// * send - The end position of the local alignment in the hit
/// * slen - The overall sequence length of the hit
/// * qstart - The start position of the local alignment in the query
/// * qend - The end position of the local alignment in the query
/// * qlen - The overall sequence length of the query
pub fn overlap_with_query(
    sstart: u32,
    send: u32,
    slen: u32,
    qstart: u32,
    qend: u32,
    qlen: u32,
) -> F64 {
    F64((((qend - qstart + 1) + (send - sstart + 1)) as f64) / ((qlen + slen) as f64))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::default::*;

    #[test]
    fn default_filter_regexs_extract_uni_prot_descriptions() {
        let t1 = "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1";
        assert_eq!(
            filter_stitle(t1, &(*FILTER_REGEXS)),
            "Probable LRR receptor-like serine/threonine-protein kinase At3g47570"
        );
    }

    #[test]
    fn default_matches_blacklist_regexs(){
        let t1 = "LRR receptor-like serine/threonine-protein kinase EFR";
        assert_eq!(false,matches_blacklist(t1, &(*BLACKLIST_STITLE_REGEXS)));

        let t2 = "Probable LRR receptor-like serine/threonine-protein kinase At3g47570";
        assert_eq!(true,matches_blacklist(t2, &(*BLACKLIST_STITLE_REGEXS)));
        
        let t3 = "Putative receptor-like protein kinase At3g47110";
        assert_eq!(true,matches_blacklist(t3, &(*BLACKLIST_STITLE_REGEXS)));

        let t4 = "hypothetical receptor-like protein kinase At3g47110";
        assert_eq!(true,matches_blacklist(t4, &(*BLACKLIST_STITLE_REGEXS)));

        let t5 = "whole Genome shotgun Sequence";
        assert_eq!(true,matches_blacklist(t5, &(*BLACKLIST_STITLE_REGEXS)));

        let t6 = "predicted Receptor-like protein kinase";
        assert_eq!(true,matches_blacklist(t6, &(*BLACKLIST_STITLE_REGEXS)));

    }
}
