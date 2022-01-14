//! Default values and global constants are kept in this module.
use regex::Regex;
use std::collections::HashMap;

lazy_static! {

    /// The score assigned to non informative words:
    pub static ref NON_INFORMATIVE_WORD_SCORE : f32 = 0.000000001;

    /// The default Blacklist of regular expressions used to check for non-informative words
    /// in the description to be excluded from scoring. If ANY of these expression matches
    /// the word is considered as non-informative
    pub static ref BLACKLIST_DESCRIPTION_WORDS_REGEXS: Vec<Regex> = vec![
        Regex::new(r"(?i)\bmember\b").unwrap(),
        Regex::new(r"(?i)\bprotein\b").unwrap(),
        Regex::new(r"(?i)\bgene\b").unwrap(),
        Regex::new(r"(?i)\btair\b").unwrap(),
        Regex::new(r"(?i)\bfragment\b").unwrap(),
        Regex::new(r"(?i)\bhomolog\b").unwrap(),
        Regex::new(r"(?i)\bcontig\b").unwrap(),
    ];

    /// The default Blacklist of regular expressions used to filter out Hit title (`stitle`) fields
    /// if they match ANY of these expressions.
    pub static ref BLACKLIST_STITLE_REGEXS: Vec<Regex> = vec![
        Regex::new(r"(?i)^similar\s+to").unwrap(),
        Regex::new(r"(?i)^probable\b").unwrap(),
        Regex::new(r"(?i)^putative\b").unwrap(),
        Regex::new(r"(?i)^predicted\b").unwrap(),
        Regex::new(r"(?i)^uncharacterized\b").unwrap(),
        Regex::new(r"(?i)^unknown\b").unwrap(),
        Regex::new(r"(?i)^hypothetical\b").unwrap(),
        Regex::new(r"(?i)^unnamed\b").unwrap(),
        Regex::new(r"(?i)^whole\s+genome\s+shotgun\s+sequence\b").unwrap(),
        Regex::new(r"(?i)^clone\b").unwrap(),
    ];

    /// The default regular expressions used to filter a Hit title (`stitle`) and retain the short
    /// human readable description.
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
        Regex::new(r"(?i)\bunknown\b").unwrap(),
        Regex::new(r"(?i)-?\blike\b").unwrap(),
        Regex::new(r"(?i)\bsimilar\b").unwrap(),
        Regex::new(r"(?i)\bpredicted\b").unwrap(),
        Regex::new(r"(?i)\bputative\b").unwrap(),
        Regex::new(r"(?i)\buncharacterized\b").unwrap(),
        Regex::new(r"(?i)\bprobable\b").unwrap(),
        Regex::new(r"(?i)\bcontig\b").unwrap(),
        Regex::new(r"(?i)\brelated\b").unwrap(),
        Regex::new(r"(?i)\bremark\b").unwrap(),
        Regex::new(r"(?i)\b\w?orf(\w?|\d+)\b").unwrap(),
    ];

    /// The default header definition of sequence similarity search result tables, i.e. mapping
    /// column names to their factual position in the to be parsed table.
    pub static ref SEQ_SIM_TABLE_COLUMNS: HashMap<String, usize> = {
        let mut h = HashMap::new();
        // Default header is 'qacc sacc bitscore stitle'
        h.insert("qacc".to_string(), 0);
        h.insert("sacc".to_string(), 1);
        h.insert("stitle".to_string(), 9);
        h
    };

    /// A Hit's description is split into words using this default regular expression.
    //pub static ref SPLIT_DESCRIPTION_REGEX: Regex = Regex::new(r"([-/|/\\;,':().\s+]+)").unwrap();
    pub static ref SPLIT_DESCRIPTION_REGEX: Regex = Regex::new(r" ").unwrap();

    /// Default sequence similarity search result table field separator:
    pub static ref SSSR_TABLE_FIELD_SEPARATOR: char = '\t';

    /// The default short description to be used for queries for which no reasonable description
    /// can be generated
    pub static ref UNKNOWN_PROTEIN_DESCRIPTION: &'static str = "unknown protein";

    /// The default short description to be used for sequence families for which no reasonable
    /// description can be generated
    pub static ref UNKNOWN_FAMILY_DESCRIPTION: &'static str = "unknown sequence family";

    /// The default regular expression to split gene family genes
    pub static ref SPLIT_GENE_FAMILY_GENES_REGEX: Regex = Regex::new(r"(\s*,\s*|\s+)").unwrap();
}
