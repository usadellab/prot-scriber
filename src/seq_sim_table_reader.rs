use super::default::FILTER_REGEXS;
use super::default::SEQ_SIM_TABLE_COLUMNS;
use super::models::Hit;
use super::models::Query;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

pub fn parse_table(
    path: &str,
    separator: char,
    table_cols: &HashMap<String, usize>,
) -> HashMap<String, Query> {
    println!("path is {}", path);
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);

    let mut h = HashMap::<String, Query>::new();
    let mut curr_query = Query::new();
    for line in reader.lines() {
        let record_line = line.unwrap();
        let record_cols: Vec<&str> = record_line.trim().split(separator).collect();
        println!("record is {:?}", record_cols);
        let qacc_i = record_cols[*table_cols.get("qacc").unwrap()];
        if qacc_i != curr_query.id {
            if curr_query.id != String::new() {
                h.insert(curr_query.id.clone(), curr_query);
            }
            curr_query = Query::from_qacc(qacc_i.to_string());
        }
        let hit_i = parse_hit(&record_cols, separator, table_cols);
        curr_query.hits.insert(hit_i.id.clone(), hit_i);
    }
    h
}

pub fn parse_hit(
    record_cols: &Vec<&str>,
    separator: char,
    table_cols: &HashMap<String, usize>,
) -> Hit {
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
    use std::path::Path;

    #[test]
    fn parses_hit_from_record_line() {
        let record_cols = vec![
            "Soltu.DM.02G015700.1", "sp|C0LGP4|Y3475_ARATH", "2209", "1284", "2199", "1010", "64", "998", "580",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
        ];
        let hit = parse_hit(&record_cols, '\t', &(*SEQ_SIM_TABLE_COLUMNS));
        let expected = Hit::new( "sp|C0LGP4|Y3475_ARATH",
            "2209", "1284", "2199", "1010", "64", "998", "580",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        assert_eq!(hit, expected);
    }

    #[test]
    fn parses_seq_sim_result_table() {
        let p = Path::new("misc").join("Soltu.DM.02G015700.1_AA_vs_trEMBL_blastpout.txt");
        let h = parse_table(p.to_str().unwrap(), '\t', &(*SEQ_SIM_TABLE_COLUMNS));
    }
}
