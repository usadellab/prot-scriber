// extern crate linked_hash_set;
// use linked_hash_set::LinkedHashSet;
use regex::Regex;
use std::collections::HashMap;


// Given filtered candidate descriptions it splits each word and returns a vector. 
pub fn split_candidates_descripts(candidate_words:&str)-> Vec<&str> {    
    let mut word_vec : Vec<&str> = candidate_words.split(|c : char| c.is_whitespace()).collect();
    word_vec.push(candidate_words); // add the original candidate word
    word_vec
}

// creates a map/dictionary containing individual candidate descriptions(keys) and if its informative (value : bool).
// The matches are converted by the if condition to match the logic. bool (True) is assigned to an informative word.
pub fn word_classification_map(candidate_word_vec:Vec<&str>,re:Vec<Regex>) -> HashMap<&str,bool> {
    let mut word_classif = HashMap::new();
    for w in candidate_word_vec{
        // let mut tag = false;
        if matches_uninformative_list(&w, &re) == true {
            let tag = false;
            word_classif.insert(w,  tag);
        }else{
            let tag = true;
            word_classif.insert(w,  tag);
        }
    };
    word_classif
}
// When parsed a set of regex strings in a vector it matches and returns true if match.
pub fn matches_uninformative_list(word: &str, regex: &Vec<Regex>) -> bool {
    regex.iter().any(|x| x.is_match(&word.to_string()))
}



pub struct PhraseInfo {
    phrases : Vec<String>,
    word_classification: HashMap<String,bool>,
}

impl PhraseInfo {
    pub fn new(hdr : Vec<String>, input_word_classification : HashMap<String,bool>) -> PhraseInfo {
        PhraseInfo {
            phrases : hdr,
            word_classification : input_word_classification,
        }
    }
    pub fn get_word_universe(&self) {
        let universe = (*self).word_classification.keys();
        
    }
}

pub fn phrasesFromCandidateDescription( word_vec : Vec<String>, word_classif : HashMap<String,bool>) -> PhraseInfo{
    let mut Phrase_info = PhraseInfo::new(word_vec, word_classif);
    let universe = Phrase_info.get_word_universe();
    Phrase_info   
}

// phrasesFromCandidateDescription(word_vec , world_classif);


#[cfg(test)]
mod tests {
    use super::*;
    // use crate::default::*;
    #[test]
    fn test_match_uninformative_regex_list(){
        let re = vec![Regex::new(r"(?i)\bterminal\b").unwrap()];
        let word = "alcohol dehydrogenase terminal";
        assert_eq!(true, matches_uninformative_list(&word, &re));

    }
    #[test]
    fn test_split_candidates_words(){
        let candidate_words = "alcohol dehydrogenase c terminal";
        let result = vec!["alcohol", "dehydrogenase", "c", "terminal", "alcohol dehydrogenase c terminal"];
        assert_eq!(result, split_candidates_descripts(candidate_words));
         
    }
    #[test]
    fn test_word_classification_map(){
        let candidate_words = "alcohol dehydrogenase c terminal";
        let re = vec![Regex::new(r"(?i)\bterminal\b").unwrap(),Regex::new(r"(?i)\bc\b").unwrap() ];
        let result = HashMap::from([("alcohol", true),("c", false), ("dehydrogenase", true), ("terminal", false), ("alcohol dehydrogenase c terminal", false)]);
        assert_eq!(result, word_classification_map(split_candidates_descripts(candidate_words),re))

    }
    // #[test]

}
