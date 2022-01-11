//! Code used to parse sequence similarity search result tables is implemented in this module.
use super::default::BLACKLIST_STITLE_REGEXS;
use super::hit::*;
use super::model_funcs::matches_blacklist;
use super::query::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::sync::mpsc::SyncSender;

pub fn row_to_cells(row: &str, butt_splitter_of_doom: &char) -> Vec<String> {
    row.split(*butt_splitter_of_doom)
        .into_iter()
        .map(|x| (*x).to_string())
        .collect()
}

pub fn parse_file(
    path: &String,
    butt_splitter_of_doom: &char,
    qacc_col: &usize,
    tx: SyncSender<(Option<Vec<String>>, Option<String>)>,
) {
    let file = File::open(path).unwrap();
    let mut reader = BufReader::new(file);
    let mut buf = String::new();
    let mut curr_query_id = String::new();
    loop {
        match reader.read_line(&mut buf) {
            Err(_) | Ok(0) => {
                // The file has been parsed, or an error occurred, mark the last query as parsed
                // completely:
                tx.send((None, Some(curr_query_id))).unwrap();
                break;
            }
            Ok(_) => {
                let row_cells = row_to_cells(&buf, &butt_splitter_of_doom);
                // Check, whether the currently processed query has been parsed completely:
                let parsed_query_option =
                    if !curr_query_id.is_empty() && curr_query_id != row_cells[*qacc_col] {
                        Some(curr_query_id)
                    } else {
                        None
                    };
                curr_query_id = row_cells[*qacc_col].clone();
                tx.send((Some(row_cells), parsed_query_option)).unwrap();
                buf.clear();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::default::SEQ_SIM_TABLE_COLUMNS;
    use std::path::Path;
    use std::sync::mpsc;

    #[test]
    fn parses_seq_sim_result_table() {
        let p = Path::new("misc")
            .join("Two_Potato_Proteins_vs_trEMBL_blastpout.txt")
            .to_str()
            .unwrap()
            .to_string();
        let (tx, rx) = mpsc::channel();
        parse_table(p, '\t', &(*SEQ_SIM_TABLE_COLUMNS), tx);
        let mut h: HashMap<String, Query> = HashMap::new();
        for received in rx {
            h.insert(received.id.clone(), received);
        }
        assert_eq!(h.len(), 2);
        assert!(h.contains_key("Soltu.DM.10G003150.1"));
        assert_eq!(h.get("Soltu.DM.10G003150.1").unwrap().hits.len(), 4);
        assert_eq!(
            h.get("Soltu.DM.10G003150.1")
                .unwrap()
                .hits
                .get("sp|P15538|C11B1_HUMAN")
                .unwrap()
                .bitscore
                .0,
            32.3
        );
        assert!(h.contains_key("Soltu.DM.02G015700.1"));
        assert_eq!(h.get("Soltu.DM.02G015700.1").unwrap().hits.len(), 500);
        assert_eq!(
            h.get("Soltu.DM.02G015700.1")
                .unwrap()
                .hits
                .get("sp|Q9FIZ3|GSO2_ARATH")
                .unwrap()
                .bitscore
                .0,
            560.0
        );
    }
}
