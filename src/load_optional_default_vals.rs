use regex::Regex;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug)]
pub struct OptionalsUserRegexLists {
    pub family_gene_separator_regex: Regex,
}

impl OptionalsUserRegexLists {
    pub fn new() -> OptionalsUserRegexLists {
        OptionalsUserRegexLists {
            family_gene_separator_regex: Regex::new(r"").unwrap(),
        }
    }

    pub fn parse_family_gene_separator_regex(path: &str) -> Regex {
        // Open stream to the gene separator regex file input file
        let file_path = path.to_string();
        let file = File::open(path).unwrap();
        let reader = BufReader::new(file);
        // read file line by line
        for (i, line) in reader.lines().enumerate() {
            let regex_line = line.unwrap();
            if regex_line.len() > 0 {
                let family_gene_separator_regex = Regex::new(&regex_line).unwrap();
            }
        }
        family_gene_separator_regex
    }
}

// Helper function to get the first index of the vector
// pub fn first<T>(v: &Vec<T>) -> Option<&T> {
//     v.first()
// }

// Load it to our struct that can be made available globally.
// This will the be checked if it has value or not later and decide whether to implement the defaults or not.
