use super::default::MAX_MATCH_REPLACE_ITERATIONS;
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

/// Fasta entries (`stitle` in sequence similarity search, e.g. Blast output) have a long title in
/// which the sequence identifier and often taxonomic information is given along with a short human
/// readable protein description. We are only interested in the latter. This function extracts the
/// short description using regular expressions.
///
/// # Arguments
///
/// * stitle - The sequence title line as found in the original Fasta file.
/// * regexs - A vector of regular expressions to be applied in series to the argument stitle to
///            extract the desired short description.
/// * `capture_replace_pairs` - An `Option` of a vector of tuples, pairing a regular expression
///                             (see crate fancy-regex for details on the syntax) and the
///                             capture-group replacement string. These are iteratively applied and
///                             the argument descriptions to prepare it for final splitting into
///                             words (see `split_descriptions` for details).
pub fn filter_stitle(
    stitle: &str,
    regexs: &Vec<Regex>,
    capture_replace_pairs: Option<&Vec<(fancy_regex::Regex, String)>>,
) -> String {
    let mut desc = regexs
        .iter()
        .fold(stitle.to_string(), |accumulated, current| {
            current.replace_all(&accumulated, "").to_string()
        })
        .to_lowercase();
    apply_capture_replace_pairs(&mut desc, capture_replace_pairs);
    // Remove preceding and trailing whitespaces, and return:
    desc.trim().to_string()
}

/// Iteratively applies argument pairs of regular expressions (fancy-regex) and replace
/// instructions (strings) to change the string referenced by argument `s`.
///
/// # Arguments
///
/// * s - A reference to a String to be modified by iterative application of the argument
/// capture-replace-pairs.
/// * capture_replace_pairs - An `Option` containing a vector of tuples, within each the first
/// entry is a regular expression (fancy-regex) and a replace instruction (string).
pub fn apply_capture_replace_pairs(
    s: &mut String,
    capture_replace_pairs: Option<&Vec<(fancy_regex::Regex, String)>>,
) {
    // Use regular expressions and replace with capture groups, if argument is given:
    if let Some(rr_tuples) = capture_replace_pairs {
        for rr_tpl in rr_tuples {
            for _ in 0..(*MAX_MATCH_REPLACE_ITERATIONS) {
                if rr_tpl.0.is_match(s).unwrap() {
                    // `s` is a mutable reference. Set the value it points to to the string
                    // produced by `replace`:
                    *s = rr_tpl.0.replace(s, &rr_tpl.1).to_string();
                } else {
                    break;
                }
            }
        }
    }
}

/// Reads in and parses a file specified by argument `path` and converts each line into an instance
/// of `Regex`. Returns a vector of the so instantiated regular expressions.
///
/// # Arguments
///
/// * `path` - A `&str` representing the path to the file containing one regular expression per
/// line.
pub fn parse_regex_file(path: &str) -> Vec<Regex> {
    // Open stream to the file
    let file_path = path.to_string();
    let file = File::open(path).expect(format!("No such file {:?}", path).as_str());
    let reader = BufReader::new(file);

    // read file line by line
    let mut parsed_regexs = vec![];
    for (i, line) in reader.lines().enumerate() {
        let regex_line = line.unwrap();

        match Regex::new(&regex_line) {
            Ok(regex) => {
                parsed_regexs.push(regex);
            }
            Err(e) => panic!("\n\n{:?} in file {:?} line <{:?}>. Could not parse the line into a Rust regular expression\n\n", e, file_path, i),
        }
    }
    parsed_regexs
}

/// Reads in and parses a file specified by argument `path` and converts each pair of lines into a
/// tupel `(fancy_regex::Regex, String)`. Returns a vector of the so instantiated tupels.
///
/// # Arguments
///
/// * `path` - A `&str` representing the path to the file containing pairs of lines. The first
/// always going to be parsed into a regular expression, the second returned as instance of
/// `String`.
pub fn parse_regex_replace_tuple_file(path: &str) -> Vec<(fancy_regex::Regex, String)> {
    // Open stream to the file
    let file_path = path.to_string();
    let file = File::open(file_path).expect(format!("No such file {:?}", path).as_str());
    let reader = BufReader::new(file);

    // Parse tuples, i.e. pairs of lines:
    let mut regex_replace_tuples: Vec<(fancy_regex::Regex, String)> = vec![];
    let mut is_regex_line = true;
    let mut regex_i: fancy_regex::Regex = fancy_regex::Regex::new("").unwrap();
    let mut n_lines = 0;
    for line in reader.lines() {
        let line_str = line.unwrap();
        if is_regex_line {
            regex_i = fancy_regex::Regex::new(&line_str).expect(
                format!(
                    "Could not parse line {:?} as a regular expression (Rust syntax).",
                    line_str
                )
                .as_str(),
            );
        } else {
            regex_replace_tuples.push((regex_i.clone(), line_str));
        }
        is_regex_line = !is_regex_line;
        n_lines += 1;
    }

    // If we have an odd number of lines, the file is malformed:
    if n_lines & 1 == 1 {
        panic!(
            "\n\n--capture-replace-pairs (-c) argument file {:?} has {:?} lines. But to construct pairs we need an even number of lines. See --help (-h) for more details.\n\n",
            path, n_lines
        );
    }

    regex_replace_tuples
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::default::*;

    #[test]
    fn default_filter_regexs_extract_uni_prot_descriptions() {
        // Test 1:
        let t1 = "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1";
        assert_eq!(
            filter_stitle(t1, &(*FILTER_REGEXS), None),
            "lrr receptor serine/threonine-protein kinase"
        );

        // Test 2 - using `default::REPLACE_REGEXS_DESCRIPTION`:
        let mut hit_words = "sp|C0LGP4|Y3475_ARATH receptor-like protein eix2 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1".to_string();
        let mut expected = "receptor protein eix";
        assert_eq!(
            expected,
            filter_stitle(
                &hit_words,
                &(*FILTER_REGEXS),
                Some(&(*CAPTURE_REPLACE_DESCRIPTION_PAIRS))
            )
        );

        // Test 3 - using `default::REPLACE_REGEXS_DESCRIPTION`:
        hit_words = "sp|C0LGP4|Y3475_ARATH subtilisin-like protease sbt4.15 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1".to_string();
        expected = "subtilisin protease sbt";
        assert_eq!(
            expected,
            filter_stitle(
                &hit_words,
                &(*FILTER_REGEXS),
                Some(&(*CAPTURE_REPLACE_DESCRIPTION_PAIRS))
            )
        );

        // Test 4 - using `default::REPLACE_REGEXS_DESCRIPTION`:
        hit_words = "sp|C0LGP4|Y3475_ARATH duf4228 domain protein OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1".to_string();
        expected = "duf~4228 domain protein";
        assert_eq!(
            expected,
            filter_stitle(
                &hit_words,
                &(*FILTER_REGEXS),
                Some(&(*CAPTURE_REPLACE_DESCRIPTION_PAIRS))
            )
        );

        // Test 5 - removes "8-7" from input stitle "sp|Q6YZZ2|GL87_ORYSJ Germin-like protein 8-7 OS=Oryza sativa subsp. japonica OX=39947 GN=GER6 PE=2 SV=1":
        hit_words = "sp|Q6YZZ2|GL87_ORYSJ Germin-like protein 8-7 OS=Oryza sativa subsp. japonica OX=39947 GN=GER6 PE=2 SV=1".to_string();
        expected = "germin protein";
        assert_eq!(
            expected,
            filter_stitle(
                &hit_words,
                &(*FILTER_REGEXS),
                Some(&(*CAPTURE_REPLACE_DESCRIPTION_PAIRS))
            )
        );

        // Test 6 remove duplicated words using default capture replace pairs:
        hit_words = "sp|Q6YZZ2|GL87_ORYSJ WRKY-like wrky-domain protein OS=Oryza sativa subsp. japonica OX=39947 GN=GER6 PE=2 SV=1".to_string();
        expected = "wrky domain protein";
        assert_eq!(
            expected,
            filter_stitle(
                &hit_words,
                &(*FILTER_REGEXS),
                Some(&(*CAPTURE_REPLACE_DESCRIPTION_PAIRS))
            )
        );

        // Test 7 checks that the identifier is filtered out
        hit_words = "sp|Q9C8M9|SRF6_ARATH Protein STRUBBELIG-RECEPTOR FAMILY 6 OS=Arabidopsis thaliana OX=3702 GN=SRF6 PE=1 SV=1".to_string();
        expected = "protein strubbelig receptor family";
        assert_eq!(
            expected,
            filter_stitle(
                &hit_words,
                &(*FILTER_REGEXS),
                Some(&(*CAPTURE_REPLACE_DESCRIPTION_PAIRS))
            )
        );

        // Test 8 also checks that no additional letters are deleted
        hit_words = "sp|Q6R2K2|SRF4_ARATH n Transferase Domain Containing Protein OS=Arabidopsis thaliana OX=3702 GN=SRF4 PE=2 SV=1".to_string();
        expected = "n transferase domain containing protein";
        assert_eq!(
            expected,
            filter_stitle(
                &hit_words,
                &(*FILTER_REGEXS),
                Some(&(*CAPTURE_REPLACE_DESCRIPTION_PAIRS))
            )
        );

        // Test 9 also checks that no additional letters are deleted
        hit_words = "sp|Q6R2K2|SRF4_ARATH P Transferase Domain Containing Protein OS=Arabidopsis thaliana OX=3702 GN=SRF4 PE=2 SV=1".to_string();
        expected = "p transferase domain containing protein";
        assert_eq!(
            expected,
            filter_stitle(
                &hit_words,
                &(*FILTER_REGEXS),
                Some(&(*CAPTURE_REPLACE_DESCRIPTION_PAIRS))
            )
        );

        // Test 10 checks that the Drosophila specific HRD description prefix 'LOW QUALITY
        // PROTEIN:' is removed:
        hit_words = "tr|A0A6P4E2J9|A0A6P4E2J9_DRORH LOW QUALITY PROTEIN: muscarinic acetylcholine receptor DM1 OS=Drosophila rhopaloa OX=1041015 GN=LOC108039593 PE=3 SV=1".to_string();
        expected = "muscarinic acetylcholine receptor dm";
        assert_eq!(
            expected,
            filter_stitle(
                &hit_words,
                &(*FILTER_REGEXS),
                Some(&(*CAPTURE_REPLACE_DESCRIPTION_PAIRS))
            )
        );

        // Test 11 checks that the Drosophila specific HRD description prefix 'Blast:' is removed:
        hit_words = "tr|A0A3B0K592|A0A3B0K592_DROGU Blast:Homeobox protein abdominal-A OS=Drosophila guanche OX=7266 GN=DGUA_6G017991 PE=3 SV=1".to_string();
        expected = "homeobox protein abdominal a";
        assert_eq!(
            expected,
            filter_stitle(
                &hit_words,
                &(*FILTER_REGEXS),
                Some(&(*CAPTURE_REPLACE_DESCRIPTION_PAIRS))
            )
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
