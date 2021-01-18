use regex::Regex;

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
}
