// use super::default::SPLIT_DESCRIPTION_REGEX;
use regex::Regex;
use std::collections::HashMap;
use std::collections::HashSet;

pub fn hrds(){

    // input vector of candidate hrds
    let candidates_vector = vec!["alcohol dehydrogenase", "manitol dehydrogenase", "cinnamyl alcohol-dehydrogenase", 
                                "geraniol dehydrogenase", "geraniol|terminal dehydrogenase", 
                                "alcohol dehydrogenase c terminal"];

    // chop into words 
    pub  fn candidate_hdrs_word_universe(vec : Vec<&str>) -> Vec<&str> {
        let mut splitted_candidates_universe = vec![];
        let SPLIT_DESCRIPTION_REGEX = Regex::new(r"([-/|/\\;,':.\s+]+)").unwrap(); // TODO: remember to optimize to get it from default // (*SPLIT_DESCRIPTION_REGEX)
        for candidate in vec.into_iter(){
            let mut split_candidate : Vec<&str> = SPLIT_DESCRIPTION_REGEX.split(&candidate).into_iter().collect();
            splitted_candidates_universe.append(& mut split_candidate);
        }
        splitted_candidates_universe

    }
    
    // Build universe word-set
    let universe = candidate_hdrs_word_universe(candidates_vector.clone());
    println!("## Universe \n {:?}", &universe);


    // Filter universe word-set using input vec regex
    // Results contained to a dictionary/map
    pub fn matches_uninformative_list(word: &str, regex: &Vec<Regex>) -> bool {
        regex.iter().any(|x| x.is_match(&word.to_string()))
    }

    pub fn word_classification_map(universe_hrd_words:Vec<&str>) -> HashMap<String,bool> {
        let re = vec![Regex::new(r"(?i)\bterminal\b").unwrap(),Regex::new(r"(?i)\bc\b").unwrap() ]; // TODO: remember to optimize to get it from default
        let mut word_classif = HashMap::new();
        for w in universe_hrd_words{
            if matches_uninformative_list(&w, &re) == true {
                let tag = false;
                word_classif.insert(w.to_lowercase(),  tag);
            }else{
                let tag = true;
                word_classif.insert(w.to_lowercase(),  tag);
            }
        };
        word_classif
    }
    let mut classified_universe = word_classification_map(universe.clone());
    println!("## Classification dictionary  \n {:?}", classified_universe);
    classified_universe.retain(|_, v| *v != true);
    println!("## Filtered Classification dictionary  \n {:?}", classified_universe);

    // word frequency map from the universe. Only the informative // Normalize
    pub fn candidates_word_frequency_map(vec : Vec<&str>) -> HashMap<&str, f32> {
        let mut freq_map: HashMap<&str, f32> = HashMap::new();
        let mut map: HashMap<String, i32> = HashMap::new();
        for i in &vec {
            *freq_map.entry(i).or_default() += 1 as f32;
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
    
    let mut word_frequencies_map = candidates_word_frequency_map(universe.clone());
    println!("{:?}", word_frequencies_map);
    // Retain only the informative words and their frequencies
    for key in classified_universe.keys(){
        // println!("{:?}", key);
        word_frequencies_map.retain(|k,_| *k != key);
    }
    println!("## Frequency map of only true values \n {:?}", word_frequencies_map);

    // Powersets
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
    
            // This should return a vector of vector of vectors (each phrase should be a vector)
    pub fn phrases(candidate_vector_description: Vec<&str>)-> Vec<String> {
        let all_phrases = powerset(&candidate_vector_description);
        let mut return_vec = vec![];
        for vec in all_phrases {
            let vec_type_string: Vec<String> = vec.iter().map(|s| s.to_string()).collect(); 
            return_vec.push(vec_type_string.join(" "));
        }
        return_vec.retain(|x| *x != ""); // remove any empty strings 
        return_vec
    }

    pub fn split_candidates_descripts(candidate_words:&str)-> Vec<&str> {
        let re = Regex::new(r"([-/|/\\;,':.\s+]+)").unwrap(); // TODO regex should be parsed via defaults
        let mut word_vec : Vec<&str> = re.split(candidate_words).into_iter().collect();
        // word_vec.push(candidate_words); // add the original candidate word
        word_vec
    }

    let mut phrases_universe = vec![];
    for candidate in candidates_vector.clone().iter(){ 
        let word_candidate_vec = split_candidates_descripts(candidate);
        let phrase = phrases(word_candidate_vec);
        phrases_universe.push(phrase);
    }
    // Hash set of phrases from each candidate
    // let mut phrases_set: HashSet<&Vec<String>> = HashSet::new();

    // for ph  in phrases_universe.iter() {
    //     // println!("{:?}",ph);
    //     phrases_set.insert(ph);
    // }

    // println!("{:?}",phrases_set);
    // // for i in phrases_set.iter() {
    // // println!("{:?}", i );
    // // }

    // Hashset of phrases from all the candidates, a universe of all phrases
    let mut phrases_universe_set: HashSet<Vec<&str>> = HashSet::new();
    for vec_ph in phrases_universe.iter(){
        for ph in vec_ph.iter() {
            // convert the phrases to a vector
            phrases_universe_set.insert(split_candidates_descripts(ph));
            // println!("{:?}", split_candidates_descripts(ph));
        }
    }

    // for i in phrases_universe_set.iter() {
    //     println!("{:?}", i );
    //     }
    println!("## Phrases to be scored \n {:?}", phrases_universe_set);

    fn inverse_information_content(word : &str, wrd_frequencies : HashMap<&str,f32>) -> f32{
        let sum_wrd_frequencies : f32 = wrd_frequencies.values().into_iter().sum(); 
        println!("sum of word frequencies --> {:?}", &sum_wrd_frequencies);
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

    for phrase in phrases_universe_set.iter()  {
        for word in phrase{
            let iic = inverse_information_content(&word, word_frequencies_map.to_owned());
            println!("IIC --> {:?} {:?}", &word, &iic);
        }
    }
}
