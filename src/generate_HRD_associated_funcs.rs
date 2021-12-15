use regex::Regex;
use std::collections::HashMap;
use super::default::SPLIT_DESCRIPTION_REGEX;
// use std::collections::HashSet;


/// Given filtered candidate descriptions it splits each word and returns a universe of all descriptions.
/// 
/// # Arguments
/// 
/// * `Vec<&str>` vector of all candidate descriptions
pub  fn candidate_hdrs_word_universe(vec : Vec<&str>) -> Vec<&str> {
    let mut splitted_candidates_universe = vec![];
    for candidate in vec.into_iter(){
        let mut split_candidate : Vec<&str> = (*SPLIT_DESCRIPTION_REGEX).split(&candidate).into_iter().collect();
        splitted_candidates_universe.append(& mut split_candidate);
    }
    splitted_candidates_universe.retain(|x | *x != "");
    splitted_candidates_universe
}



// /// Given filtered candidate descriptions it splits each word and returns a vector.
// /// 
// /// # Arguments
// /// 
// /// * `String of a candidate HRD description` & `Vector containing Regex strings`

// pub fn split_candidates_descripts(candidate_words:&str)-> Vec<&str> {    
//     let mut word_vec : Vec<&str> = candidate_words.split(|c : char| c.is_whitespace()).collect();
//     word_vec.push(candidate_words); // add the original candidate word
//     word_vec

// }

/// Creates a HashMap containing individual candidate descriptions(keys) and if its informative (value : bool).
/// The matches are converted by the if condition to match the logic. bool (True) is assigned to an informative word.
/// 
/// # Arguments
/// 
/// * `Vector of strings` & `Vector of Regex` 
pub fn word_classification_map(universe_hrd_words:Vec<&str>,re:Vec<Regex>) -> HashMap<&str,bool> {
    let mut word_classif = HashMap::new();
    for w in universe_hrd_words{
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
/// 
/// # Arguments
/// 
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
//     pub fn new(candidate_hdr : Vec<String>, input_word_classification : HashMap<String,bool>) -> PhraseInfo {
//         PhraseInfo {
//             phrases : candidate_hdr,
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
/// 
/// # Arguments
/// 
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
/// 
///  # Arguments
/// 
/// * `Vec<&str>`- Vector containing candidate descriptions as strings.
pub fn phrases(candidate_vector_description: Vec<&str>)-> Vec<String> {
    let all_phrases = powerset(&candidate_vector_description);
    let mut return_vec = vec![];
    for vec in all_phrases {
        let vec_type_string: Vec<String> = vec.iter().map(|s| s.to_string()).collect(); // convert from Vec<&str> to Vec<String>
        return_vec.push(vec_type_string.join(" "));
    }
    return_vec.retain(|x| *x != ""); // remove any empty strings 
    return_vec
}
/// Gets word frequencies of all candidate words
/// 
/// # Arguments
///  
/// * `Vector<String>` - vector of strings containing all candidate HRDS.
pub fn frequencies(vec : Vec<String>) -> HashMap<String, f32> {
    let mut freq_map: HashMap<String, f32> = HashMap::new();
    let mut map: HashMap<String, i32> = HashMap::new();
    for i in &vec {
        *freq_map.entry(i.to_string().to_lowercase()).or_default() += 1 as f32;
    }
    for i in &vec {
        *map.entry(i.to_string()).or_default() += 1 as i32;
    }
    let max_value = map.values().max().unwrap().clone();
    
    for key in freq_map.to_owned().keys() {
        *freq_map.get_mut(key).unwrap() /= max_value as f32; //normalize
    };

freq_map
}
   /// Computes the score of a word using 'inverse information content' calculated
    /// as -log( (1-p.word) ). In order to avoid infinite values for a word that is
    /// the single element of the word-set, i.e. it has a frequency of one, the
    /// score of one is returned.
    /// 
    /// # Arguments
    /// 
    /// * `word` - A string representing the word
    /// * `wrd.frequencies` - An instance of dictionary of all words with their frequencies.

    fn inverse_information_content(word : &str, wrd_frequencies : HashMap<&str,f32>) -> f32{
        let sum_wrd_frequencies : f32 = wrd_frequencies.values().into_iter().sum(); 
        if wrd_frequencies.values().len() as f32 > 1. {
            let pw = wrd_frequencies[word]/sum_wrd_frequencies;
            let iic = f32::log10(1. / (1. -  pw));
            iic
        } else if wrd_frequencies.values().len() as f32 == 1. && wrd_frequencies.contains_key(word) {
                1.
            } else {
                panic!("Invalid or no word frequency parsed");
            }
    }


    fn mean(list: &Vec<f64>) -> f64 {
        let sum: f64 = Iterator::sum(list.iter());
        f64::from(sum) / (list.len() as f64)
    }


    /// Computes the score of a word using 'inverse information content' and centering
    /// 
    /// # Arguments
    /// 
    /// *`word_frequencies`  An instance of dictionary of all words with their frequencies.

    fn centered_word_scores(wrd_frequencies: HashMap<&str, f64>) -> f64 /*HashMap<&str, f64>*/ {
        let mut all_iic = vec![];
        for word in wrd_frequencies.keys() {
            let mut iic = inverse_information_content(word, wrd_frequencies);
            all_iic.append(iic);
        }
        
        let ave = mean(all_iic);
        ave

        // println!("{:?}", all_iic);
    }


#[cfg(test)]
mod tests {
    use super::*;
    // use crate::default::SPLIT_DESCRIPTION_REGEX;
    // use crate::default::*;
    // use super::default::SPLIT_DESCRIPTION_REGEX;
    #[test]
    fn candidate_word_universe() {
        let candidates_vector = vec!["Alcohol dehydrogenase", "manitol dehydrogenase", "Cinnamyl alcohol-dehydrogenase", 
                                                "Geraniol dehydrogenase", "Geraniol|terminal dehydrogenase", "alcohol dehydrogenase c terminal"];
        let result = vec!["Alcohol", "dehydrogenase", "manitol", "dehydrogenase", "Cinnamyl", "alcohol", "dehydrogenase", "Geraniol",
                                    "dehydrogenase", "Geraniol", "terminal", "dehydrogenase", "alcohol", "dehydrogenase", "c", "terminal"];
        assert_eq!(result, candidate_hdrs_word_universe(candidates_vector));

    }
    #[test]
    fn test_match_uninformative_regex_list(){
        let re = vec![Regex::new(r"(?i)\bterminal\b").unwrap()];
        let word = "alcohol dehydrogenase terminal";
        assert_eq!(true, matches_uninformative_list(&word, &re));

    }
    #[test]
    // fn test_split_candidates_words(){
    //     let candidate_words = "alcohol dehydrogenase c terminal";
    //     let result = vec!["alcohol", "dehydrogenase", "c", "terminal", "alcohol dehydrogenase c terminal"];
    //     assert_eq!(result, split_candidates_descripts(candidate_words));
         
    // }
    #[test]
    // fn test_word_classification_map(){
    //     let candidate_words = "alcohol dehydrogenase c terminal";
    //     let re = vec![Regex::new(r"(?i)\bterminal\b").unwrap(),Regex::new(r"(?i)\bc\b").unwrap() ];
    //     let result = HashMap::from([("alcohol", true),("c", false), ("dehydrogenase", true), 
    //                                                     ("terminal", false), ("alcohol dehydrogenase c terminal", false)]);
    //     assert_eq!(result, word_classification_map(split_candidates_descripts(candidate_words),re))

    // }
    
    #[test]
    fn test_phrases(){
        let candidate_descript_vec = vec!["alcohol", "dehydrogenase", "c", "terminal"];
        let result = vec!["alcohol", "dehydrogenase", "alcohol dehydrogenase", "c", "alcohol c", "dehydrogenase c",
                                     "alcohol dehydrogenase c", "terminal", "alcohol terminal", "dehydrogenase terminal",  
                                     "alcohol dehydrogenase terminal", "c terminal", "alcohol c terminal", 
                                     "dehydrogenase c terminal", "alcohol dehydrogenase c terminal"];
        assert_eq!(phrases(candidate_descript_vec), result);
    }

    #[test]
    fn test_frequencies(){
        let vec = vec!["Alcohol".to_string(),"dehydrogenase".to_string(), "c".to_string(), "terminal".to_string(),];
        let mut result = HashMap::new();
        result.insert("terminal".to_string(),1.0);
        result.insert("dehydrogenase".to_string(),1.0);
        result.insert("alcohol".to_string(),1.0);
        result.insert("c".to_string(),1.0);

        assert_eq!(result, frequencies(vec));
    }
}
