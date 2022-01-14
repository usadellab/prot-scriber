//! Code used to parse sequence similarity search result tables is implemented in this module.
use super::default::{BLACKLIST_STITLE_REGEXS, FILTER_REGEXS};
use super::model_funcs::{filter_stitle, matches_blacklist};
use super::query::*;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::sync::mpsc::Sender;

/// Finds a tabular file (`path`) and parses it in a stream approach, i.e. line by line. Every time
/// an instance of Query is successfully and completely pares it is send using the argument
/// `transmitter` to the respective registered receiver.
///
/// # Arguments
///
/// * `path: String` - The path to the tabular sequence similarity search result file to parse
/// * `separator: char` - The separator to use to split a line into an array of columns
/// * `table_cols: &HashMap<String, usize>` - The header information, i.e. the column names and
/// their respective position in the table (`path`).
/// * `transmitter: Sender<Query>` - Used to send instances of `Query` to any receiver.
pub fn parse_table(path: String, transmitter: Sender<(String, Query)>) {
    let lines =
        read_lines(&path).expect(format!("An error occurred reading file '{}'", &path).as_str());
    let mut last_qacc = String::new();
    let mut curr_query = Query::new();
    for line_rslt in lines {
        match line_rslt {
            Ok(line) => {
                let cols: Vec<&str> = line.trim().split('\t').collect();
                let qacc = cols[0];
                let sacc = cols[1];
                let stitle = cols[9];

                if qacc != last_qacc && !last_qacc.is_empty() {
                    transmitter.send((last_qacc, curr_query)).unwrap();
                    curr_query = Query::new();
                }

                if !curr_query.hits.contains_key(&sacc.to_string())
                    && !matches_blacklist(stitle, &(*BLACKLIST_STITLE_REGEXS))
                {
                    let desc = filter_stitle(stitle, &(*FILTER_REGEXS))
                        .trim()
                        .to_lowercase();
                    if !desc.is_empty() {
                        curr_query.hits.insert(sacc.to_string(), desc);
                    }
                }

                last_qacc = qacc.to_string();
            }
            Err(e) => {
                eprintln!("An error occurred while parsing '{}':\n{:?}", path, e);
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
