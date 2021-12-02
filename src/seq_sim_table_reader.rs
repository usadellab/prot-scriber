//! Code used to parse sequence similarity search result tables is implemented in this module.
use super::annotation_process::AnnotationProcess;
use super::hit::*;
use super::query::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::sync::{Arc, Mutex};

/// Finds a tabuar file (`path`) and parses it in a stream approach, i.e. line by line. Returns a
/// in memory database of the respective parsed queries and their hits. Note that this function is
/// making use of a thread safe reference to the argument `annotation_process`.
///
/// # Arguments
///
/// * `path: &str` - The path to the tabular sequence similarity search result file to parse
/// * `separator: char` - The separator to use to split a line into an array of columns
/// * `table_cols: &HashMap<String, usize>` - The header information, i.e. the column names and
///                                           their respective position in the table (`path`).
/// * `annotation_process: AnnotationProcess` - A mutable reference to an instance of
///                                                  AnnotationProcess in which to gather the
///                                                  parsed sequence similarity search results,
///                                                  i.e. the Queries and the Hits.
pub fn parse_table(
    path: &str,
    separator: char,
    table_cols: &HashMap<String, usize>,
    annotation_process: Arc<Mutex<&mut AnnotationProcess>>,
) {
    // Open stream to the sequence similarity search result table:
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);

    // The current query for which to read in hits:
    let mut curr_query = Query::new();
    // Process line by line in streamed read from table:
    for line in reader.lines() {
        let record_line = line.unwrap();
        let record_cols: Vec<&str> = record_line.trim().split(separator).collect();
        // Obtain the current query accession (ID) knowing its column number:
        let qacc_i = record_cols[*table_cols.get("qacc").unwrap()];
        // Did we hit a new block, i.e. a new query?:
        if qacc_i != curr_query.id {
            // Is it the very first line?:
            if curr_query.id != String::new() {
                // Insert the results collected for the last query:
                let mut ap = annotation_process.lock().unwrap();
                ap.insert_query(&mut curr_query);
                // release lock:
                drop(ap);
            }
            // Prepare gathering of results for the next query:
            curr_query = Query::from_qacc(qacc_i.to_string());
        }
        // Process the current hit:
        let hit_i = parse_hit(&record_cols, table_cols);
        curr_query.add_hit(&hit_i);
    }
    // Insert the last query, because we reached the end of the file:
    let mut ap = annotation_process.lock().unwrap();
    ap.insert_query(&mut curr_query);
    // release lock:
    drop(ap);
}

/// Parses a line in the respective sequence similarity search result table. The line is already
/// split into cells (columns). Returns a new instance of `Hit`.
///
/// # Arguments
///
/// * `record_cols: &Vec<&str>` - The record parsed from splitting the respective line
/// * `table_cols: &HashMap<String, usize>` - The table header, i.e. the column names and their
///                                           position in the argument `record_cols`.
pub fn parse_hit(record_cols: &Vec<&str>, table_cols: &HashMap<String, usize>) -> Hit {
    Hit::new(
        record_cols[*table_cols.get("sacc").unwrap()],
        record_cols[*table_cols.get("qlen").unwrap()],
        record_cols[*table_cols.get("qstart").unwrap()],
        record_cols[*table_cols.get("qend").unwrap()],
        record_cols[*table_cols.get("slen").unwrap()],
        record_cols[*table_cols.get("sstart").unwrap()],
        record_cols[*table_cols.get("send").unwrap()],
        record_cols[*table_cols.get("bitscore").unwrap()],
        record_cols[*table_cols.get("stitle").unwrap()],
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::default::SEQ_SIM_TABLE_COLUMNS;
    use std::path::Path;

    #[test]
    fn parses_hit_from_record_line() {
        let record_cols = vec![
            "Soltu.DM.02G015700.1", "sp|C0LGP4|Y3475_ARATH", "2209", "1284", "2199", "1010", "64", "998", "580",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
        ];
        let hit = parse_hit(&record_cols, &(*SEQ_SIM_TABLE_COLUMNS));
        let expected = Hit::new( "sp|C0LGP4|Y3475_ARATH",
            "2209", "1284", "2199", "1010", "64", "998", "580",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        assert_eq!(hit, expected);
    }

    #[test]
    fn parses_seq_sim_result_table() {
        let p = Path::new("misc").join("Two_Potato_Proteins_vs_trEMBL_blastpout.txt");
        let mut ap = AnnotationProcess::new();
        parse_table(
            p.to_str().unwrap(),
            '\t',
            &(*SEQ_SIM_TABLE_COLUMNS),
            &mut ap,
        );
        let h = ap.queries;
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