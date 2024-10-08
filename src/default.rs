//! Default values and global constants are kept in this module.
use regex::Regex;
use std::collections::HashMap;

lazy_static! {

    /// The score assigned to non informative words:
    pub static ref NON_INFORMATIVE_WORD_SCORE : f64 = 0.000001;

    /// The default argument `AnnotationProcess.center_iic_at_quantile` to be used for centering
    /// the inverse information content values of words. The literal 50.0 indicates centering at
    /// the mean and not actually a quantile:
    pub static ref CENTER_INVERSE_INFORMATION_CONTENT_AT_QUANTILE : f64 = 50.0;

    /// The default Blacklist of regular expressions used to check for non-informative words
    /// in the description to be excluded from scoring. If ANY of these expression matches
    /// the word is considered as non-informative
    pub static ref NON_INFORMATIVE_WORDS_REGEXS: Vec<Regex> = vec![
        Regex::new(r"(?i)\bthe\b").unwrap(),
        Regex::new(r"(?i)\band\b").unwrap(),
        Regex::new(r"(?i)\bor\b").unwrap(),
        Regex::new(r"(?i)\bfrom\b").unwrap(),
        Regex::new(r"(?i)\bto\b").unwrap(),
        Regex::new(r"(?i)\bmember\b").unwrap(),
        Regex::new(r"(?i)\bprotein\b").unwrap(),
        Regex::new(r"(?i)\bgene\b").unwrap(),
        Regex::new(r"(?i)\btair\b").unwrap(),
        Regex::new(r"\b\d+\b").unwrap(),
    ];

    /// The default Blacklist of regular expressions used to filter out Hit title (`stitle`) fields
    /// if they match ANY of these expressions.
    pub static ref BLACKLIST_STITLE_REGEXS: Vec<Regex> = vec![
        Regex::new(r"(?i)\bsimilar\s+to").unwrap(),
        Regex::new(r"(?i)\bprobable\b").unwrap(),
        Regex::new(r"(?i)\bputative\b").unwrap(),
        Regex::new(r"(?i)\bpredicted\b").unwrap(),
        Regex::new(r"(?i)\buncharacterized\b").unwrap(),
        Regex::new(r"(?i)\bunknown\b").unwrap(),
        Regex::new(r"(?i)\bhypothetical\b").unwrap(),
        Regex::new(r"(?i)\bunnamed\b").unwrap(),
        Regex::new(r"(?i)\bfragment\b").unwrap(),
        Regex::new(r"(?i)\bwhole\s+genome\s+shotgun\s+sequence\b").unwrap(),
        Regex::new(r"(?i)\bclone\b").unwrap(),
    ];

    /// The default regular expressions used to filter a Hit title (`stitle`) and retain the short
    /// human readable description.
    pub static ref FILTER_REGEXS: Vec<Regex> = vec![
        Regex::new(r"(?i)low\s+quality\s+protein:").unwrap(),
        Regex::new(r"\sOS=.*$").unwrap(),
        Regex::new(r"(?i)OS.*[.].*protein").unwrap(),
        Regex::new(r"(?i)^H0.*protein").unwrap(),
        Regex::new(r"(?i)contains.*").unwrap(),
        Regex::new(r"IPR.*").unwrap(),
        Regex::new(r"\b\S+\|\S+\|\S+").unwrap(),
        Regex::new(r"\w{2,}\d{1,2}[gGmMcC]\d+(\.\d+)*").unwrap(),
        Regex::new(r"\b\[.*").unwrap(),
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
        Regex::new(r"(?i)\bisoform\b").unwrap(),
        Regex::new(r"(?i)\bhomolog\b").unwrap(),
        Regex::new(r"(?i)\bfragment\b").unwrap(),
        Regex::new(r"(?i)\bcontig\b").unwrap(),
        Regex::new(r"(?i)\bblast:\b").unwrap(),
        ];

    /// The default header definition of sequence similarity search result tables, i.e. mapping
    /// column names to their factual position in the to be parsed table.
    pub static ref SEQ_SIM_TABLE_COLUMNS: HashMap<String, usize> = {
        let mut h = HashMap::new();
        // Default header is 'qacc sacc bitscore stitle'
        h.insert("qacc".to_string(), 0);
        h.insert("sacc".to_string(), 1);
        h.insert("stitle".to_string(), 2);
        h
    };

    /// A Hit's description is split into words using this default regular expression.
    pub static ref SPLIT_DESCRIPTION_REGEX: Regex = Regex::new(r"([()~_\-/|\\;,':.\s]+)").unwrap();

    /// The default vector of regular expressions _with_ match-groups to be used to split
    /// descriptions (parsed `stitle`) into separate words by replacing the matched region with
    /// the first and second captures:
    pub static ref CAPTURE_REPLACE_DESCRIPTION_PAIRS: Vec<(fancy_regex::Regex, String)> = {
        let mut rrd : Vec<(fancy_regex::Regex, String)> = vec![];
        rrd.push(
            (
                // Protects InterPro, PANTHER, Pfam annotations from being mangled by subsequent
                // tuples:
                fancy_regex::Regex::new(r"(?i)\b(?P<first>duf|pf|ipr|pthr|go|kegg|ec)(?P<second>[0-9:]+)\b").unwrap(),
                r"$first~$second".to_string()
            )
        );
        rrd.push(
            (
                // Transforms e.g. "eix2" into "eix" or "SBT4.15" into "SBT" - case insensitive:
                fancy_regex::Regex::new(r"(?i)\b(?P<first>[a-z]{2,})[-.,\d]+\b").unwrap(),
                "$first ".to_string()
            )
        );
        rrd.push(
            (
                // Deletes numbers like "4", "4.12", or "12-4":
                fancy_regex::Regex::new(r"(^|\s+)[-.\d]+(\s+|$)").unwrap(),
                " ".to_string()
            )
        );
        rrd.push(
            (
                // Deletes multiple occurrences of words using the extended fancy-regex crate
                // syntax. Only the first mention of a word occurring multiple times is retained:
                fancy_regex::Regex::new(r"(?i)\b(?P<first>\b\w+\b)(?P<spacer>.*)\b\k<first>\b").unwrap(),
                r"$first$spacer".to_string()
            )
        );
        rrd.push(
            (
                // Replace multiple adjacent whitespae characters with a single one:
                fancy_regex::Regex::new(r"\s{2,}").unwrap(),
                r" ".to_string()
            )
        );
        rrd
    };

    /// The default vector of regular expressions _with_ match-groups to be used to post-process
    /// ("polish") assigned human readable descriptions before using them as final output:
    pub static ref POLISH_CAPTURE_REPLACE_PAIRS: Vec<(fancy_regex::Regex, String)> = {
        let mut rrd : Vec<(fancy_regex::Regex, String)> = vec![];
        rrd.push(
            (
                // Deletes trailing non-informative words like 'and', 'or', 'the' etc:
                fancy_regex::Regex::new(r"(?i)\s*\b(and|or|the|from|to)\b\s*$").unwrap(),
                r"".to_string()
            )
        );
        rrd
    };

    /// The maximum number of applying one tuple of regular expression and match-group replacing in
    /// generate_hrd_associated_funcs::split_descriptions( ..., `replace_regexs`) (see above
    /// `REPLACE_REGEXS_DESCRIPTION`):
    pub static ref MAX_MATCH_REPLACE_ITERATIONS: u8 = u8::MAX;

    /// Default sequence similarity search result table field separator:
    pub static ref SSSR_TABLE_FIELD_SEPARATOR: char = '\t';

    /// The default short description to be used for queries for which no reasonable description
    /// can be generated
    pub static ref UNKNOWN_PROTEIN_DESCRIPTION: &'static str = "unknown protein";

    /// The default short description to be used for sequence families for which no reasonable
    /// description can be generated
    pub static ref UNKNOWN_FAMILY_DESCRIPTION: &'static str = "unknown sequence family";

    /// The default regular expression to split gene family genes
    pub static ref SPLIT_GENE_FAMILY_GENES_REGEX: &'static str = r"(\s*,\s*|\s+)";

    /// The default character used to split gene-family-identifiers from the set of genes the
    /// respective family is comprised of:
    pub static ref SPLIT_GENE_FAMILY_ID_FROM_GENE_SET: &'static str = "\t";
}
