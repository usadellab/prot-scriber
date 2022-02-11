use regex::Regex;
use std::fs::File;
use std::io::{BufRead, BufReader};

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

/// Reads in and parses a file specified by argument `path` and converts each line into an instance
/// of `Regex`. Returns a vector of the so instantiated regular expressions.
///
/// # Arguments
///
/// * `path` - A `&str` representing the path to the file containing on regular expression per
/// line.
pub fn parse_regex_file(path: &str) -> Vec<Regex> {
    // Open stream to the gene families input file
    let file_path = path.to_string();
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);

    // read file line by line
    let mut parsed_regexs = vec![];
    for (i, line) in reader.lines().enumerate() {
        let regex_line = line.unwrap();

        match Regex::new(&regex_line) {
            Ok(regex) => {
                parsed_regexs.push(regex);
            }
            Err(e) => panic!("{:?} in file {:?} line <{:?}>. Could not parse the line into a Rust regular expression", e, file_path, i),
        }
    }
    parsed_regexs
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
            "lrr receptor serine/threonine- kinase"
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
