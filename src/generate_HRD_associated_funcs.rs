use regex::Regex;
use std::collections::HashMap;
// use std::collections::HashSet;


/// Given filtered candidate descriptions it splits each word and returns a vector.
/// # Arguments
/// * `String of a candidate HRD description` & `Vector containing Regex strings`
// More robust function exist in the hit.rs file `description_words`
pub fn split_candidates_descripts(candidate_words:&str)-> Vec<&str> {    
    let mut word_vec : Vec<&str> = candidate_words.split(|c : char| c.is_whitespace()).collect();
    word_vec.push(candidate_words); // add the original candidate word
    word_vec

}

/// Creates a HashMap containing individual candidate descriptions(keys) and if its informative (value : bool).
/// The matches are converted by the if condition to match the logic. bool (True) is assigned to an informative word.
/// # Arguments
/// * `Vector of strings` & `Vector of Regex` 
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
/// # Arguments
/// *`String` - word & `Vec<Regex>` - Vector of Regex strings
pub fn matches_uninformative_list(word: &str, regex: &Vec<Regex>) -> bool {
    regex.iter().any(|x| x.is_match(&word.to_string()))
}


// pub struct PhrasesInfo {
//     phrases : HashSet<Vec<String>>,
//     word_classification: HashMap<String,bool>,
//     word_frequencies: HashMap<String,f32>
// }

// impl PhrasesInfo {
//     pub fn new(hdr : Vec<String>, input_word_classification : HashMap<String,bool>) -> PhraseInfo {
//         PhraseInfo {
//             phrases : hdr,
//             word_classification : input_word_classification,
//         }
//     }
//     pub fn get_word_universe(&self) {
//         let universe = (*self).word_classification.keys();
        
//     }
// }

// pub fn phrasesFromCandidateDescription( word_vec : Vec<String>, word_classif : HashMap<String,bool>) -> PhraseInfo{
//     let mut Phrase_info = PhraseInfo::new(word_vec, word_classif);
//     let universe = Phrase_info.get_word_universe();
//     Phrase_info   
// }


/// Powersets any slice
/// # Arguments
/// * `slice`- Given a slice vector get all possible powersets of the items.
pub fn powerset<T: Clone>(slice: &[T]) -> Vec<Vec<T>> {
    let mut v: Vec<Vec<T>> = Vec::new();

    for mask in 0..(1 << slice.len()) {
        let mut ss: Vec<T> = vec![];
        let mut bitset = mask;
        while bitset > 0 {
            let rightmost: u64 = bitset & !(bitset - 1);
            let idx = rightmost.trailing_zeros();
            let item = (*slice.get(idx as usize).unwrap()).clone();
            ss.push(item);
            bitset &= bitset - 1;
        }
        v.push(ss);
    }
    v
}


/// Given candidate descriptions gets all phrases from the candidate descriptions
/// Gets all phrases using the powerset function`powerset`
///  # Arguments
/// * `Vec<&str>`- Vector containing candidate descriptions as strings.
pub fn get_possible_hrd_phrases(candidate_vector_description: Vec<&str>)-> Vec<String> {
    let all_phrases = powerset(&candidate_vector_description);
    let mut return_vec = vec![];
    for vec in all_phrases {
        let vec_type_string: Vec<String> = vec.iter().map(|s| s.to_string()).collect(); // convert from Vec<&str> to Vec<String>
        return_vec.push(vec_type_string.join(" "));
    }
    return_vec.retain(|x| *x != ""); // remove any empty strings 
    return_vec
}

// TODO fix borrowing issues
/// Gets word frequencies of all candidate words
/// # Arguments 
/// * `Vector<String>` - vector of strings containing all candidate HRDS.
pub fn get_candidate_word_frequency(vec : Vec<String>) -> HashMap<String, f32> {
    let mut freq_map: HashMap<String, f32> = HashMap::new();
    let mut map: HashMap<String, i32> = HashMap::new();
    for i in &vec {
        *freq_map.entry(i.to_string()).or_default() += 1 as f32;
    }
    for i in &vec {
        *map.entry(i.to_string()).or_default() += 1 as i32;
    }
    let max_value = map.values().max().unwrap().clone();
    
    for key in freq_map.keys() {
        *freq_map.get_mut(key).unwrap() /= max_value as f32; //normalize
    };

freq_map
}


// fn get_possiable_hrd_phrases(candidate_vector_descriptions){
//     let  p = powerset(candidate_vector_descriptions);

//     for i in p {
//         let vector_String: Vec<String> = i.iter().map(|s| s.to_string()).collect();
//         println!("{:?}\n{:?}", &i, hashset(&vector_String));
// }  
// }




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
    // fn test_vec_to_hashset(){
    //     let test_vec = vec!["alcohol", "c", "alcohol dehydrogenase c terminal"];
    //     let v4: Vec<String> = test_vec.iter().map(|s| s.to_string()).collect();
    //     let result = {"dehydrogenase", "c", "terminal", "alcohol dehydrogenase c terminal", "alcohol"};
    //     assert_eq!(result, hashset(&v4));
    // }
    #[test]
    fn test_get_possible_hdr_phrases(){
        let candidate_descript_vec = vec!["alcohol", "dehydrogenase", "c", "terminal"];
        let result = vec!["alcohol", "dehydrogenase", "alcohol dehydrogenase", "c", "alcohol c", "dehydrogenase c", "alcohol dehydrogenase c", "terminal", "alcohol terminal", "dehydrogenase terminal", "alcohol dehydrogenase terminal", "c terminal", "alcohol c terminal", "dehydrogenase c terminal", "alcohol dehydrogenase c terminal"];
        assert_eq!(get_possible_hrd_phrases(candidate_descript_vec), result);
    }
//TODO:
//     #[test]
//     fn get_vector_word_frequency(){
//         let vec = vec!["alcohol".to_string(),"dehydrogenase".to_string(), "c".to_string(), "terminal".to_string()];
//         let result = {"alcohol": 1 , "dehydrogenase":1 , "c" : 1, "terminal":1};
//         assert_eq!(result, get_candidate_word_frequency(&vec));
//     }
}
