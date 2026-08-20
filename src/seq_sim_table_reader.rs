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

    // Send last parsed query. Note this must NOT require `curr_query.hits` to be non-empty: the
    // qacc-change branch above always sends `curr_query` regardless of whether any hits survived
    // blacklist/filtering, so a query whose hits are all blacklisted (e.g. every Hit is
    // "hypothetical protein") must be sent here too, or it is silently dropped from the
    // annotation process entirely whenever it happens to be the last query in the file -- instead
    // of being registered with zero hits and annotated as "unknown protein", like an otherwise
    // identical query positioned anywhere else in the (sorted) input file:
    if !last_qacc.is_empty() {
        transmitter.send((last_qacc, curr_query)).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::default::{BLACKLIST_STITLE_REGEXS, FILTER_REGEXS};
    use std::collections::HashMap;
    use std::sync::mpsc;

    // Regression test for a bug where a query is silently dropped -- never sent to the receiver
    // at all -- if (a) every one of its Hit descriptions matches a blacklist regex (e.g.
    // "hypothetical protein") AND (b) it happens to be the last `qacc` block in the input file.
    // Every other position in the file is unaffected by (a) alone, because the qacc-change branch
    // inside the parsing loop always sends the accumulated `curr_query` regardless of whether any
    // of its hits survived blacklist/filtering; only the trailing send after the loop, for the
    // very last query block, additionally requires `curr_query.hits` to be non-empty -- which an
    // all-blacklisted last query never is. The practical impact: such a query does not even show
    // up as "unknown protein" in the output. It just vanishes, and whether it does so depends
    // purely on its position in the (correctly sorted) input file, not on its data.
    #[test]
    fn parse_table_sends_last_query_even_if_all_its_hits_are_blacklisted() {
        let path = Path::new("misc")
            .join("tmp_test_parse_table_last_query_all_blacklisted.txt")
            .to_str()
            .unwrap()
            .to_string();
        std::fs::write(
            &path,
            concat!(
                "Query1\tHit1\tsome informative kinase domain\n",
                "Query2\tHit2\thypothetical protein\n"
            ),
        )
        .unwrap();

        let (tx, rx) = mpsc::channel();
        parse_table(
            &path,
            &'\t',
            &0,
            &1,
            &2,
            &BLACKLIST_STITLE_REGEXS,
            &FILTER_REGEXS,
            None,
            tx,
        );
        std::fs::remove_file(&path).unwrap();

        let received: HashMap<String, Query> = rx.into_iter().collect();
        assert!(received.contains_key("Query1"));
        // Query2 is the last block in the file and all its hits are blacklisted ("hypothetical
        // protein"); it must still be reported, with zero surviving hits, not silently dropped:
        assert!(received.contains_key("Query2"));
        assert!(received.get("Query2").unwrap().hits.is_empty());
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
