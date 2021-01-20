use regex::Regex;
use std::collections::HashMap;

lazy_static! {
    pub static ref FILTER_REGEXS: Vec<Regex> = vec![
        Regex::new(r"\sOS=.*$").unwrap(),
        Regex::new(r"(?i)OS.*[.].*protein").unwrap(),
        Regex::new(r"(?i)^H0.*protein").unwrap(),
        Regex::new(r"(?i)contains.*").unwrap(),
        Regex::new(r"IPR.*").unwrap(),
        Regex::new(r"\w{2,}\d{1,2}(g|G)\d+(\.\d)*\s+").unwrap(),
        Regex::new(r"\b\[.*").unwrap(),
        Regex::new(r"\b\S+\|\S+\|\S+").unwrap(),
        Regex::new(r"\(\s*Fragment\s*\)").unwrap(),
        Regex::new(r"^(\s|/|\(|\)|-|\+|\*|,|;|\.|:|\||\d)+$").unwrap(),
    ];
    pub static ref SEQ_SIM_TABLE_COLUMNS: HashMap<String, usize> = {
        let mut h = HashMap::new();
        h.insert("sacc".to_string(), 0);
        h.insert("qacc".to_string(), 1);
        h.insert("bitscore".to_string(), 2);
        h.insert("qstart".to_string(), 3);
        h.insert("qend".to_string(), 4);
        h.insert("qlen".to_string(), 5);
        h.insert("sstart".to_string(), 6);
        h.insert("send".to_string(), 7);
        h.insert("slen".to_string(), 8);
        h.insert("stitle".to_string(), 9);
        h
    };
}
