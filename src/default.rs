//! Default values and global constants are kept in this module.
use regex::Regex;
use std::collections::HashMap;

lazy_static! {

    /// The default Blacklist of regular expressions used to check for non-informative words
    /// in the description to be excluded from scoring. If ANY of these expression matches
    /// the word is considered as non-informative
    pub static ref BLACKLIST_DESCRIPTION_WORDS_REGEXS: Vec<Regex> = vec![
        Regex::new(r"(?i)\bunknown\b").unwrap(),
        Regex::new(r"(?i)\bmember\b").unwrap(),
        Regex::new(r"(?i)\blike\b").unwrap(),
        Regex::new(r"(?i)\bassociated\b").unwrap(),
        Regex::new(r"(?i)\bcontaining\b").unwrap(),
        Regex::new(r"(?i)\bactivated\b").unwrap(),
        Regex::new(r"(?i)\bfamily\b").unwrap(),
        Regex::new(r"(?i)\bsubfamily\b").unwrap(),
        Regex::new(r"(?i)\binteracting\b").unwrap(),
        Regex::new(r"(?i)\bactivity\b").unwrap(),
        Regex::new(r"(?i)\bsimilar\b").unwrap(),
        Regex::new(r"(?i)\bproduct\b").unwrap(),
        Regex::new(r"(?i)\bexpressed\b").unwrap(),
        Regex::new(r"(?i)\bpredicted\b").unwrap(),
        Regex::new(r"(?i)\bputative\b").unwrap(),
        Regex::new(r"(?i)\buncharacterized\b").unwrap(),
        Regex::new(r"(?i)\bprobable\b").unwrap(),
        Regex::new(r"(?i)\bprotein\b").unwrap(),
        Regex::new(r"(?i)\bgene\b").unwrap(),
        Regex::new(r"(?i)\btair\b").unwrap(),
        Regex::new(r"(?i)\bfragment\b").unwrap(),
        Regex::new(r"(?i)\bhomolog\b").unwrap(),
        Regex::new(r"(?i)\bcontig\b").unwrap(),
        Regex::new(r"(?i)\brelated\b").unwrap(),
        Regex::new(r"(?i)\bremark\b").unwrap(),
        Regex::new(r"(?i)\b\w?orf(\w?|\d+)\b").unwrap(),
    ];

    /// The default Blacklist of regular expressions used to filter out Hit title (`stitle`) fields
    /// if they match ANY of these expressions.
    pub static ref BLACKLIST_STITLE_REGEXS: Vec<Regex> = vec![
        Regex::new(r"(?i)^similar\s+to").unwrap(),
        Regex::new(r"(?i)^probable ").unwrap(),
        Regex::new(r"(?i)^putative ").unwrap(),
        Regex::new(r"(?i)^predicted ").unwrap(),
        Regex::new(r"(?i)^uncharacterized").unwrap(),
        Regex::new(r"(?i)^unknown").unwrap(),
        Regex::new(r"(?i)^hypothetical").unwrap(),
        Regex::new(r"(?i)^unnamed").unwrap(),
        Regex::new(r"(?i)^whole\s+genome\s+shotgun\s+sequence").unwrap(),
        Regex::new(r"(?i)^clone").unwrap(),
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
    ];

    /// The default header definition of sequence similarity search result tables, i.e. mapping
    /// column names to their factual position in the to be parsed table.
    pub static ref SEQ_SIM_TABLE_COLUMNS: HashMap<String, usize> = {
        let mut h = HashMap::new();
        // Default header is 'qacc sacc qlen qstart qend slen sstart send bitscore stitle'
        h.insert("qacc".to_string(), 0);
        h.insert("sacc".to_string(), 1);
        h.insert("qlen".to_string(), 2);
        h.insert("qstart".to_string(), 3);
        h.insert("qend".to_string(), 4);
        h.insert("slen".to_string(), 5);
        h.insert("sstart".to_string(), 6);
        h.insert("send".to_string(), 7);
        h.insert("bitscore".to_string(), 8);
        h.insert("stitle".to_string(), 9);
        h
    };

    /// A Hit's description is split into words using this default regular expression.
    pub static ref SPLIT_DESCRIPTION_REGEX: Regex = Regex::new(r"([-/|/\\;,':().\s+]+)").unwrap();

    /// Default sequence similarity search result table field separator:
    pub static ref SSSR_TABLE_FIELD_SEPARATOR: char = '\t';

    /// The default inflation parameter (I) for Markov Clustering
    pub static ref MCL_INFLATION: f64 = 5.0;

    /// The default maximum number of MCL iterations to be done
    pub static ref MCL_MAX_ITERATIONS: i8 = 10;

    /// The minimum delta to accept as no change between Markov Clustering iterations
    pub static ref MCL_DELTA: f64 = 0.0001;

    /// The default number of digits to round stochastic matrices' decimal values to
    pub static ref ROUND_DECIMAL_DIGITS: i32 = 4;

    /// The default short description to be used for queries for which no reasonable description
    /// can be generated
    pub static ref UNKNOWN_PROTEIN_DESCRIPTION: &'static str = "unknown protein";

    /// The default short description to be used for sequence families for which no reasonable
    /// description can be generated
    pub static ref UNKNOWN_FAMILY_DESCRIPTION: &'static str = "unknown sequence family";

    /// The default string is used to collapse (join) consensus descriptions of disjoint
    /// hit-clusters:
    pub static ref CLUSTER_CONSENSUS_DESCRIPTIONS_JOIN: &'static str = "; ";

    /// The default regular expressions used to label / filter uninformative words.
    pub static ref UNINFORMATIVE_REGEXS: Vec<Regex> = vec![
        Regex::new(r"(?i)\bprotein\b").unwrap(),
        Regex::new(r"(?i)\bterminal\b").unwrap(),
        Regex::new(r"(?i)\bc\b").unwrap(),
    ];

    /// The default regular expression to split gene family genes
    pub static ref SPLIT_GENE_FAMILY_GENES_REGEX: Regex = Regex::new(r"(\s*,\s*|\s+)").unwrap();
}
