use super::default::BLACKLIST_DESCRIPTION_WORDS_REGEXS;
use super::default::SPLIT_DESCRIPTION_REGEX;
use super::default::UNKNOWN_PROTEIN_DESCRIPTION;
use std::collections::HashMap;
use std::collections::HashSet;

/// Main function for generating human-readable descriptions (hrds).
///
/// # Arguments
///
/// * `hit_hrds: &Vec<String>` - A vector of strings containing all Hit descriptions.
pub fn generate_human_readable_description(descriptions: &Vec<String>) -> String {
    // Split the descriptions into vectors of words:
    let description_words: Vec<Vec<String>> = descriptions
        .iter()
        .map(|dsc| split_descriptions(dsc))
        .collect();

    // Build universe of words:
    let word_universe: Vec<String> = description_words.into_iter().flatten().collect();

    // Calculate the frequency of the universe words:
    let word_frequencies = frequencies(&word_universe);

    // All possible phrases using the powerset function
    let mut phrases_set: HashSet<Vec<String>> = HashSet::new();
    for desc_words in description_words {
        phrases_set.insert(desc_words);
    }
    // Convert to vector to ease usage and save memory in scoring:
    let phrases: Vec<Vec<String>> = phrases_set.into_iter().collect();

    // Score phrases using centered inverse information content
    let phrase_scores: Vec<f32> = score_phrases(&word_frequencies, &phrases);

    // Hrd description
    let hrd_res = predicted_hrd(&phrases, &phrase_scores);

    match hrd_res {
        Ok(hrd_vec) => hrd_vec.join(" "),
        Err(()) => (*UNKNOWN_PROTEIN_DESCRIPTION).to_string(),
    }
}

/// Given filtered Hit descriptions it splits each word and returns a universe of all descriptions.
///
/// # Arguments
///
/// * `descriptions: &Vec<String>` vector of all Hit descriptions
pub fn word_universe(descriptions: &Vec<String>) -> Vec<String> {
    let mut word_universe = vec![];
    for hit_description in descriptions.into_iter() {
        let mut split_hit: Vec<String> = (*SPLIT_DESCRIPTION_REGEX)
            // Generate a set of words from the hit_description:
            .split(&hit_description)
            .collect::<Vec<&str>>()
            .iter()
            // Filter out words that are not informative:
            .filter(|x| !matches_uninformative_list(x))
            .map(|x| (*x).to_string())
            .collect();
        // Add the words to the universe, retaining their frequency - this is a vector, not a set:
        word_universe.append(&mut split_hit);
    }
    word_universe
}

/// Given filtered Hit descriptions it splits each word and returns a vector.
///
/// # Arguments
///
/// * `String of a Hit HRD description` & `Vector containing Regex strings`
pub fn split_descriptions(description: &String) -> Vec<String> {
    (*SPLIT_DESCRIPTION_REGEX)
        .split(description)
        .map(|x| (*x).to_string())
        .collect()
}

/// Creates a HashMap containing individual Hit descriptions(keys) and if its informative (value : bool).
/// The matches are converted by the if condition to match the logic. bool (True) is assigned to an informative word.
///
/// # Arguments
///
/// * `Vector of strings`
/// * `Vector of Regex`
pub fn word_classification_map(universe_hrd_words: Vec<&str>) -> HashMap<String, bool> {
    let mut word_classif = HashMap::new();
    for w in universe_hrd_words {
        if matches_uninformative_list(&w) == true {
            let tag = false;
            word_classif.insert(w.to_lowercase(), tag);
        } else {
            let tag = true;
            word_classif.insert(w.to_lowercase(), tag);
        }
    }
    word_classif
}

/// Generates a vector of uninformative words from the universe of Hit hrd words
///
/// Arguments
///
/// * `Vector of strings`
///
pub fn uninformative_words_vec(universe_hrd_words: Vec<&str>) -> Vec<String> {
    let mut uninformative_words = vec![];
    for w in universe_hrd_words {
        if matches_uninformative_list(&w) == true {
            uninformative_words.push(w.to_lowercase());
        }
    }
    uninformative_words
}

/// When parsed a set of regex strings in a vector it matches and returns true if match.
///
/// # Arguments
///
/// * `String` - word & `Vec<Regex>` - Vector of Regex strings
pub fn matches_uninformative_list(word: &str) -> bool {
    (*BLACKLIST_DESCRIPTION_WORDS_REGEXS)
        .iter()
        .any(|x| x.is_match(&word.to_string()))
}

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

/// Given Hit descriptions gets all phrases from the Hit descriptions
/// Gets all phrases using the powerset function`powerset`
///
///  # Arguments
///
/// * `Vec<&str>`- Vector containing Hit descriptions as strings.
pub fn phrases(hit_vector_description: Vec<&str>) -> Vec<String> {
    let all_phrases = powerset(&hit_vector_description);
    let mut return_vec = vec![];
    for vec in all_phrases {
        let vec_type_string: Vec<String> = vec.iter().map(|s| s.to_string()).collect(); // convert from Vec<&str> to Vec<String>
        return_vec.push(vec_type_string.join(" "));
    }
    return_vec.retain(|x| *x != ""); // remove any empty strings
    return_vec
}

/// Calculates the word frequencies for argument `universe_words` and returns a `HashMap<String,
/// f32>` mapping the words to their respective frequency.
///
/// # Arguments
///
/// * `universe_words: &Vector<String>` - vector of words
pub fn frequencies(universe_words: &Vec<String>) -> HashMap<String, f32> {
    let mut word_freqs: HashMap<String, f32> = HashMap::new();
    for word in universe_words.iter() {
        if !word_freqs.contains_key(word) {
            let n_appearances = universe_words.iter().filter(|x| (*x) == word).count() as f32
                / universe_words.len() as f32;
            word_freqs.insert((*word).clone(), n_appearances);
        }
    }
    word_freqs
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
pub fn inverse_information_content(word: &str, wrd_frequencies: HashMap<&str, f32>) -> f32 {
    let sum_wrd_frequencies: f32 = wrd_frequencies.values().into_iter().sum();
    if wrd_frequencies.values().len() as f32 > 1. {
        let pw = wrd_frequencies[word] / sum_wrd_frequencies;
        let iic = f32::log(1. / (1. - pw), 2.71828182846); //TODO refactor base exp(1)
        iic
    } else if wrd_frequencies.values().len() as f32 == 1. && wrd_frequencies.contains_key(word) {
        let iic = 1.;
        iic
    } else {
        panic!("Invalid or no word frequency parsed");
    }
}

/// Returns the mean of list or vector containing f32 values.
///
/// # Arguments
///
/// * `Vec<f32>` A list of f32 values
pub fn mean(list: Vec<f32>) -> f32 {
    let sum: f32 = Iterator::sum(list.iter());
    f32::from(sum) / (list.len() as f32)
}

/// Computes the mean of all word's inverse information content
///  and returns type f32
///
/// # Arguments
///
/// * ` HashMap<&str, f32>` a map containing all words and their frequencies
///
/// * `HashSet<Vec<&str>>` a set containing all the phrases
pub fn mean_inv_inf_cntn(
    word_frequencies_map: HashMap<&str, f32>,
    phrases_universe_set: HashSet<Vec<&str>>,
) -> f32 {
    let mut all_iic = vec![];
    for phrase in phrases_universe_set.clone().iter() {
        for word in phrase {
            if word_frequencies_map.contains_key(word) {
                let iic = inverse_information_content(&word, word_frequencies_map.clone());
                all_iic.push(iic.clone());
            }
        }
    }

    let p_w_i = mean(all_iic.clone());
    p_w_i
}

/// Computes the score of each word in a phrase using 'inverse information content' and centering.
/// Outputs a HashMap of phrases and the sum of centered iic scores.
///
/// # Arguments
///
/// * `word_frequencies`  An instance of dictionary of all words with their frequencies.
/// * `phrases_universe_set` HashSet of Vectors containing all possible phrases (universe phrases)
pub fn score_phrases(
    word_frequencies_map: &HashMap<String, f32>,
    phrases_universe_set: &Vec<Vec<String>>,
) -> Vec<f32> {
    let p_w_i = mean_inv_inf_cntn(word_frequencies_map.clone(), phrases_universe_set.clone());
    let mut phrases_score_map = HashMap::new();
    for phrase in phrases_universe_set.iter().cloned() {
        let mut phrase_word_scores = vec![];
        for word in phrase.clone() {
            if word_frequencies_map.contains_key(word) {
                let iic = inverse_information_content(word, word_frequencies_map.to_owned());
                let word_score = iic - p_w_i;
                phrase_word_scores.push(word_score);
            }
        }
        let sum_phrase: f32 = phrase_word_scores.iter().sum();
        phrases_score_map.insert(phrase.clone(), sum_phrase);
    }
    phrases_score_map
}

/// Selects the highest scoring phrase from a map of phrases (keys) and overall scores (values)
/// outputs vector of the highest scoring phrases and length. Returns the index (`Option<usize>`)
/// of the best performing argument phrase (`phrases`).
///
/// # Arguments
///
/// * `phrases: &Vec<Vec<String>>` - The phrases that have been scored and of which the best
/// performing shall be selected and returned.
/// * `phrase_scores: &Vec<f32>` - The scores of the phrases matching by their index (`usize`).
pub fn find_best_scoring_phrase(
    phrases: &Vec<Vec<String>>,
    phrase_scores: &Vec<f32>,
) -> Option<usize> {
    if phrases.len() > 0 {
        // get max score without changing the order of phrase_scores
        phrase_scores.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let max_score: f32 = phrase_scores[0];
        // Vector of tuples (phrase index, phrase length):
        let mut high_scoring_phrases_tuples: Vec<(usize, usize)> = vec![];
        for (i, phr_scr) in phrase_scores.iter().enumerate() {
            if (*phr_scr) == max_score {
                (high_scoring_phrases_tuples.push(i),);
            }
        }
        if high_scoring_phrases_tuples.len() > 0 {
            Some(0)
        } else {
            Some(
                phrase_scores
                    .iter()
                    .position(|phr_scr| (*phr_scr) == max_score)
                    .unwrap(),
            )
        }
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;
    #[test]
    fn hit_word_universe() {
        let hits_vector = vec![
            "Alcohol dehydrogenase",
            "manitol dehydrogenase",
            "Cinnamyl alcohol-dehydrogenase",
            "Geraniol dehydrogenase",
            "Geraniol|terminal dehydrogenase",
            "alcohol dehydrogenase c terminal",
        ];
        let result = vec![
            "Alcohol",
            "dehydrogenase",
            "manitol",
            "dehydrogenase",
            "Cinnamyl",
            "alcohol",
            "dehydrogenase",
            "Geraniol",
            "dehydrogenase",
            "Geraniol",
            "terminal",
            "dehydrogenase",
            "alcohol",
            "dehydrogenase",
            "c",
            "terminal",
        ];
        assert_eq!(result, word_universe(hits_vector));
    }

    #[test]
    fn test_hit_hdrs_word_universe() {
        let hit_hrds = vec![
            "alcohol dehydrogenase",
            "manitol dehydrogenase",
            "cinnamyl alcohol-dehydrogenase",
            "geraniol dehydrogenase",
            "geraniol|dehydrogenase terminal",
            "manitol dehydrogenase",
            "alcohol dehydrogenase c-terminal",
        ];

        let result = vec![
            "alcohol",
            "dehydrogenase",
            "manitol",
            "dehydrogenase",
            "cinnamyl",
            "alcohol",
            "dehydrogenase",
            "geraniol",
            "dehydrogenase",
            "geraniol",
            "dehydrogenase",
            "terminal",
            "manitol",
            "dehydrogenase",
            "alcohol",
            "dehydrogenase",
            "c",
            "terminal",
        ];
        assert_eq!(result, word_universe(hit_hrds));
    }

    #[test]
    fn test_split_hits_descripts() {
        let hit_words = "alcohol dehydrogenase c terminal";
        let result = vec!["alcohol", "dehydrogenase", "c", "terminal"];
        assert_eq!(result, split_descriptions(hit_words));
    }
    #[test]
    fn test_match_uninformative_list() {
        let word = "protein";
        let another_word = "alcohol";
        assert_eq!(true, matches_uninformative_list(&word));
        assert_eq!(false, matches_uninformative_list(&another_word));
    }

    #[test]
    fn test_word_classification_map() {
        let hit_words = vec!["alcohol", "dehydrogenase", "protein", "homolog"];

        let result = HashMap::from([
            ("alcohol".to_string(), true),
            ("protein".to_string(), false),
            ("dehydrogenase".to_string(), true),
            ("homolog".to_string(), false),
        ]);
        assert_eq!(result, word_classification_map(hit_words));
    }

    #[test]
    fn test_uninformative_words_vector() {
        let hit_words = vec!["alcohol", "dehydrogenase", "protein", "homolog"];
        let result = vec!["protein".to_string(), "homolog".to_string()];
        assert_eq!(result, uninformative_words_vec(hit_words));
    }

    #[test]
    fn test_powerset() {
        let test_vec = vec!["alcohol", "dehydrogenase", "c", "terminal"];
        let result = vec![
            vec![], // this removed in the function that puts the phrases together
            vec!["alcohol"],
            vec!["dehydrogenase"],
            vec!["alcohol", "dehydrogenase"],
            vec!["c"],
            vec!["alcohol", "c"],
            vec!["dehydrogenase", "c"],
            vec!["alcohol", "dehydrogenase", "c"],
            vec!["terminal"],
            vec!["alcohol", "terminal"],
            vec!["dehydrogenase", "terminal"],
            vec!["alcohol", "dehydrogenase", "terminal"],
            vec!["c", "terminal"],
            vec!["alcohol", "c", "terminal"],
            vec!["dehydrogenase", "c", "terminal"],
            vec!["alcohol", "dehydrogenase", "c", "terminal"],
        ];
        assert_eq!(result, powerset(&test_vec))
    }

    #[test]
    fn test_phrases() {
        let hit_descript_vec = vec!["alcohol", "dehydrogenase", "c", "terminal"];
        let result = vec![
            "alcohol".to_string(),
            "dehydrogenase".to_string(),
            "alcohol dehydrogenase".to_string(),
            "c".to_string(),
            "alcohol c".to_string(),
            "dehydrogenase c".to_string(),
            "alcohol dehydrogenase c".to_string(),
            "terminal".to_string(),
            "alcohol terminal".to_string(),
            "dehydrogenase terminal".to_string(),
            "alcohol dehydrogenase terminal".to_string(),
            "c terminal".to_string(),
            "alcohol c terminal".to_string(),
            "dehydrogenase c terminal".to_string(),
            "alcohol dehydrogenase c terminal".to_string(),
        ];
        assert_eq!(result, phrases(hit_descript_vec));
    }

    #[test]
    fn test_frequencies() {
        let vec = vec!["alcohol", "dehydrogenase", "c", "terminal"];
        let mut result = HashMap::new();
        result.insert("terminal", 1.0);
        result.insert("dehydrogenase", 1.0);
        result.insert("alcohol", 1.0);
        result.insert("c", 1.0);

        assert_eq!(result, frequencies(vec));
    }

    #[test]
    fn test_inverse_information_content() {
        let word = "a";
        let mut freq_map = HashMap::new();
        freq_map.insert("a", 3. as f32);
        freq_map.insert("b", 2. as f32);
        freq_map.insert("c", 2. as f32);
        freq_map.insert("d", 1. as f32);
        freq_map.insert("e", 1. as f32);
        freq_map.insert("f", 1. as f32);

        let pw: f32 = freq_map.values().into_iter().sum();

        let result = f32::log(1. / (1. - 3. / pw), 2.71828182846);

        assert_eq!(result, inverse_information_content(word, freq_map))
    }

    #[test]
    fn test_mean() {
        let numbers = vec![2., 3., 4., 5., 6., 7., 8.];
        let result = 5. as f32;
        assert_eq!(result, mean(numbers));
    }

    #[test]
    fn test_mean_inv_inf_cntn() {
        let mut word_set = HashSet::new();
        word_set.insert(vec!["a"]);
        word_set.insert(vec!["b"]);
        word_set.insert(vec!["c"]);
        word_set.insert(vec!["d"]);
        word_set.insert(vec!["e"]);
        word_set.insert(vec!["f"]);

        let mut freq_map = HashMap::new();
        freq_map.insert("a", 3. as f32);
        freq_map.insert("b", 2. as f32);
        freq_map.insert("c", 2. as f32);
        freq_map.insert("d", 1. as f32);
        freq_map.insert("e", 1. as f32);
        freq_map.insert("f", 1. as f32);

        let mut iic_vec = vec![];
        for phrase in &word_set {
            for word in phrase {
                let iic = inverse_information_content(word, freq_map.clone());
                iic_vec.push(iic);
            }
        }
        let result = mean(iic_vec);

        assert_eq!(result, mean_inv_inf_cntn(freq_map, word_set))
    }

    #[test]
    fn test_centered_word_scores_phrases() {
        let mut word_set = HashSet::new();
        word_set.insert(vec!["a"]);
        word_set.insert(vec!["b"]);
        word_set.insert(vec!["c"]);
        word_set.insert(vec!["d"]);
        word_set.insert(vec!["e"]);
        word_set.insert(vec!["f"]);

        let mut freq_map = HashMap::new();
        freq_map.insert("a", 3. as f32);
        freq_map.insert("b", 2. as f32);
        freq_map.insert("c", 2. as f32);
        freq_map.insert("d", 1. as f32);
        freq_map.insert("e", 1. as f32);
        freq_map.insert("f", 1. as f32);

        let mut result = HashMap::new();
        result.insert(vec!["a"], 0.1701677 as f32);
        result.insert(vec!["d"], -0.08114673 as f32);
        result.insert(vec!["f"], -0.08114673 as f32);
        result.insert(vec!["c"], 0.036636263 as f32);
        result.insert(vec!["b"], 0.036636263 as f32);
        result.insert(vec!["e"], -0.08114673 as f32);
        assert_eq!(result, score_phrases(freq_map, word_set));
    }

    #[test]
    fn test_predicted_hrd() {
        let result = vec!["dehydrogenase", "c", "terminal"];
        let mut scored_phrases = HashMap::new();
        scored_phrases.insert(vec!["dehydrogenase", "c"], 0.265823);
        scored_phrases.insert(vec!["dehydrogenase", "c", "terminal"], 0.265823);
        scored_phrases.insert(vec!["alcohol", "dehydrogenase"], 0.22149113);
        scored_phrases.insert(vec!["cinnamyl", "alcohol", "dehydrogenase"], -0.110522866);

        assert_eq!(result, predicted_hrd(scored_phrases).unwrap());
    }

    #[test]
    fn test_generate_human_readable_description() {
        let hit_hrds = vec![
            "manitol dehydrogenase".to_string(),
            "cinnamyl alcohol-dehydrogenase".to_string(),
            "geraniol dehydrogenase".to_string(),
            "geraniol|dehydrogenase terminal".to_string(),
            "manitol dehydrogenase".to_string(),
            "alcohol dehydrogenase c-terminal".to_string(),
        ];
        let expected = "dehydrogenase".to_string();
        let result = generate_human_readable_description(&hit_hrds);
        println!("\n\n{}\n\n", result);

        assert_eq!(expected, result);
    }
}
