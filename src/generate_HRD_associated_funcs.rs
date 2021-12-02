use regex::Regex;
use std::collections::HashMap;
use std::collections::HashSet;
use std::iter::FromIterator;

/// Given filtered candidate descriptions it splits each word and returns a vector. 
pub fn split_candidates_descripts(candidate_words:&str)-> Vec<&str> {    
    let mut word_vec : Vec<&str> = candidate_words.split(|c : char| c.is_whitespace()).collect();
    word_vec.push(candidate_words); // add the original candidate word
    word_vec
}

/// creates a map/dictionary containing individual candidate descriptions(keys) and if its informative (value : bool).
/// The matches are converted by the if condition to match the logic. bool (True) is assigned to an informative word.
pub fn word_classification_map(candidate_word_vec:Vec<&str>,re:Vec<Regex>) -> HashMap<&str,bool> {
    let mut word_classif = HashMap::new();
    for w in candidate_word_vec{
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
/// When parsed a set of regex strings in a vector it matches and returns true if match.
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


/// Powersets any slice
pub fn powerset<T: Clone>(slice: &[T]) -> Vec<Vec<T>> {
    let mut v: Vec<Vec<T>> = Vec::new();

    for mask in 0..(1 << slice.len()) {
        let mut ss: Vec<T> = vec![];
        let mut bitset = mask;
        while bitset > 0 {
            // isolate the rightmost bit to select one item
            let rightmost: u64 = bitset & !(bitset - 1);
            // turn the isolated bit into an array index
            let idx = rightmost.trailing_zeros();
            let item = (*slice.get(idx as usize).unwrap()).clone();
            ss.push(item);
            // zero the trailing bit
            bitset &= bitset - 1;
        }
        v.push(ss);
    }
    v
}



/// Creates a hashset from a Vector. Not genetic. 
/// But can be tweaked to take other types carried by a vector.
fn hashset(data: &Vec<String>) -> HashSet<String> {
    HashSet::from_iter(data.iter().cloned())
}

fn get_possiable_hrd_phrases(candidate_vector_descriptions){
    let  p = powerset(candidate_vector_descriptions);

    for i in p {
        let vector_String: Vec<String> = i.iter().map(|s| s.to_string()).collect();
        println!("{:?}\n{:?}", &i, hashset(&vector_String));
}  
}




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
