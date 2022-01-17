use regex::Regex;

/// Applies a list of regular expressions to a string. If any of them match returns true,
/// otherwise returns false.
/// This function is used to exclude non-informative descriptions from a human readable
/// sequence title and to check for non-informative words to be excluded from scoring
///
/// # Arguments
///
/// * testee - The text to be tested for any matching argument regular expression (`regexs`)
/// * regexs - A vector of regular expression to be applied to the testee argument.
pub fn matches_blacklist(testee: &str, regexs: &Vec<Regex>) -> bool {
    regexs.iter().any(|x| x.is_match(&testee.to_string()))
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
            .trim()
            .to_lowercase(),
    )
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
            "lrr receptor serine/threonine-protein kinase at3g47570"
        );
    }

    #[test]
    fn default_matches_blacklist_regexs() {
        let t1 = "LRR receptor-like serine/threonine-protein kinase EFR";
        assert_eq!(false, matches_blacklist(t1, &(*BLACKLIST_STITLE_REGEXS)));

        let t2 = "Probable LRR receptor-like serine/threonine-protein kinase At3g47570";
        assert_eq!(true, matches_blacklist(t2, &(*BLACKLIST_STITLE_REGEXS)));

        let t3 = "Putative receptor-like protein kinase At3g47110";
        assert_eq!(true, matches_blacklist(t3, &(*BLACKLIST_STITLE_REGEXS)));

        let t4 = "hypothetical receptor-like protein kinase At3g47110";
        assert_eq!(true, matches_blacklist(t4, &(*BLACKLIST_STITLE_REGEXS)));

        let t5 = "whole Genome shotgun Sequence";
        assert_eq!(true, matches_blacklist(t5, &(*BLACKLIST_STITLE_REGEXS)));

        let t6 = "predicted Receptor-like protein kinase";
        assert_eq!(true, matches_blacklist(t6, &(*BLACKLIST_STITLE_REGEXS)));
    }
}
