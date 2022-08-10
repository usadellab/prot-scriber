//! Code used to parse sequence similarity search result tables is implemented in this module.
use super::model_funcs::{filter_stitle, matches_blacklist};
use super::query::*;
use regex::Regex;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::sync::mpsc::Sender;

/// Finds a tabular file (`path`) and parses it in a stream approach, i.e. line by line. Every time
/// an instance of Query is successfully and completely parsed it is send using the argument
/// `transmitter` to the respective registered receiver.
///
/// # Arguments
///
/// * `path: String` - The path to the tabular sequence similarity search result file to parse
/// * `field_separator: char` - The separator to use to split a line into an array of columns
/// * `qacc_col: &usize` - The column index in which to find the `qacc`
/// * `sacc_col: &usize` - The column index in which to find the `sacc`
/// * `stitle_col: &usize` - The column index in which to find the `stitle`
/// * `blacklist_regexs: &Vec<Regex>` - The list of regular expressions used to identify to be
/// discarded descriptions (`stitle`) parsed from the argument `path` sequence similarity search
/// result table.
/// * `filter_regexs: &Vec<Regex>` - The list of regular expressions used to identify to be deleted
/// matching sub-strings in the descriptions (`stitle`) parsed from the argument `path` sequence
/// similarity search result table.
/// * `capture_replace_pairs` - An `Option` of a vector of tuples, pairing a regular expression
/// and the capture-group replacement string. These are iteratively applied and the argument
/// descriptions to prepare it for final splitting into words (see `split_descriptions` for
/// details).
/// * `transmitter: Sender<Query>` - Used to send instances of `Query` to any receiver.
pub fn parse_table(
    path: &String,
    field_separator: &char,
    qacc_col: &usize,
    sacc_col: &usize,
    stitle_col: &usize,
    blacklist_regexs: &Vec<Regex>,
    filter_regexs: &Vec<Regex>,
    capture_replace_pairs: Option<&Vec<(fancy_regex::Regex, String)>>,
    transmitter: Sender<(String, Query)>,
) {
    let lines =
        read_lines(&path).expect(format!("An error occurred reading file {:?}", &path).as_str());
    let mut last_qacc = String::new();
    let mut curr_query = Query::new();
    for line_rslt in lines {
        match line_rslt {
            Ok(line) => {
                let cols: Vec<&str> = line.trim().split(*field_separator).collect();
                let qacc = cols[*qacc_col];
                let sacc = cols[*sacc_col];
                let stitle = cols[*stitle_col];

                if qacc != last_qacc && !last_qacc.is_empty() {
                    transmitter.send((last_qacc, curr_query)).unwrap();
                    curr_query = Query::new();
                }

                if !curr_query.hits.contains_key(&sacc.to_string())
                    && !matches_blacklist(stitle, blacklist_regexs)
                {
                    let desc = filter_stitle(stitle, filter_regexs, capture_replace_pairs)
                        .trim()
                        .to_lowercase();
                    if !desc.is_empty() {
                        curr_query.hits.insert(sacc.to_string(), desc);
                    }
                }

                last_qacc = qacc.to_string();
            }
            Err(e) => {
                eprintln!(
                    "\nAn error occurred while parsing {:?}:\n{:?}\nContinuing anyway!\n",
                    path, e
                );
            }
        }
    }

    // Send last parsed query:
    if curr_query.hits.len() > 0 && !last_qacc.is_empty() {
        transmitter.send((last_qacc, curr_query)).unwrap();
    }
}

/// The output is wrapped in a Result to allow matching on errors Returns an Iterator to the Reader
/// of the lines of the file.
///
/// # Arguments
///
/// * `filename` The path to the file to open a `BufReader` for.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
