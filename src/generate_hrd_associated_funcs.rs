use super::default::SPLIT_DESCRIPTION_REGEX;
use std::collections::HashMap;
use std::collections::HashSet;

/// Main function for generating human-readable descriptions (hrds).
///
/// # Arguments
///
/// * `hit_hrds: &Vec<String>` - A vector of strings containing all Hit descriptions.
pub fn generate_human_readable_description(descriptions: &Vec<String>) -> Option<String> {
    // Split the descriptions into vectors of words:
    let description_words: Vec<Vec<String>> = descriptions
        .iter()
        .map(|dsc| split_descriptions(dsc))
        .collect();

    // All possible phrases using the power-set function and build the universe of words, i.e. the
    // flattened Vec<String> of description_words:
    let mut phrases_set: HashSet<Vec<String>> = HashSet::new();
    let mut word_universe: Vec<String> = vec![];
    for desc_words in description_words.into_iter() {
        word_universe.extend_from_slice(&desc_words);
        // Get all possible subsets of the current description `desc_words` vector:
        let desc_phrases: Vec<Vec<String>> = powerset(&desc_words);
        // Add those, that are not already present in the `phrases_set`:
        for phrase in desc_phrases.into_iter() {
            phrases_set.insert(phrase);
        }
    }

    // Calculate the frequency of the universe words:
    let word_frequencies = frequencies(&word_universe);

    // Convert to vector to ease usage and save memory in scoring. `phrases_set` is moved by the
    // `into_iter` and thus does not exist any more below the following expression:
    let phrases: Vec<Vec<String>> = phrases_set.into_iter().collect();

    // Score phrases using centered inverse information content:
    let phrase_scores: Vec<f32> = score_phrases(&word_frequencies, &phrases);

    // Find the index of the best performing phrase:
    let best_phrase_index = find_best_scoring_phrase(&phrases, &phrase_scores);

    // Return the result, if any:
    match best_phrase_index {
        Some(phrase_ind) => Some(phrases[phrase_ind].join(" ")),
        None => None,
    }
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

/// Computes the score of a word using 'inverse information content' calculated as `-1 * log(1 -
/// p.word)`. In order to avoid infinite values for a word that is the single element of the
/// word-set, i.e. it has a frequency of one, the score of one is returned.
///
/// # Arguments
///
/// * `word` - A string representing the word
/// * `wrd.frequencies` - An instance of dictionary of all words with their frequencies.
pub fn inverse_information_content(word: &String, wrd_frequencies: &HashMap<String, f32>) -> f32 {
    let sum_wrd_frequencies: f32 = wrd_frequencies.values().into_iter().sum();
    if wrd_frequencies.values().len() as f32 > 1. {
        let pw = wrd_frequencies[word] / sum_wrd_frequencies;
        f32::log(1. / (1. - pw), std::f32::consts::E)
    } else if wrd_frequencies.values().len() as f32 == 1. && wrd_frequencies.contains_key(word) {
        1.0
    } else {
        panic!("Invalid or no word frequency parsed");
    }
}

/// Returns the mean of list or vector containing f32 values.
///
/// # Arguments
///
/// * `score_map` A vector of f32 values
pub fn mean(score_map: &HashMap<String, f32>) -> f32 {
    let sum: f32 = Iterator::sum(score_map.values());
    f32::from(sum) / (score_map.len() as f32)
}

/// Computes the score of each word in a phrase using 'inverse information content' and centering.
/// Outputs a HashMap of phrases and the sum of centered iic scores.
///
/// # Arguments
///
/// * `word_frequencies` - An instance of dictionary of all words with their frequencies.
/// * `phrases` - HashSet of Vectors containing all possible phrases (universe phrases)
pub fn score_phrases(
    word_frequencies: &HashMap<String, f32>,
    phrases: &Vec<Vec<String>>,
) -> Vec<f32> {
    let mut word_inv_inf_cntnt: HashMap<String, f32> = HashMap::new();
    for word in word_frequencies.keys() {
        word_inv_inf_cntnt.insert(
            (*word).clone(),
            inverse_information_content(word, word_frequencies),
        );
    }
    let mean_word_inv_inf_cntnt = mean(&word_inv_inf_cntnt);
    let mut phrase_scores: Vec<f32> = vec![];
    for phrase in phrases.iter() {
        let mut phr_score = 0.0;
        for word in phrase.iter() {
            if word_inv_inf_cntnt.contains_key(word) {
                phr_score += word_inv_inf_cntnt.get(word).unwrap() / mean_word_inv_inf_cntnt;
            }
        }
        phrase_scores.push(phr_score);
    }
    phrase_scores
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
        let mut max_score: f32 = 0.0;
        for phr_scr in phrase_scores {
            if (*phr_scr) > max_score {
                max_score = *phr_scr
            }
        }
        // Vector of tuples (phrase index, phrase length):
        let mut high_scoring_phrases_tuples: Vec<(usize, usize)> = vec![];
        for (i, phr_scr) in phrase_scores.iter().enumerate() {
            if (*phr_scr) == max_score {
                high_scoring_phrases_tuples.push((i, phrases[i].len()));
            }
        }
        if high_scoring_phrases_tuples.len() > 0 {
            // Sort equally high scoring phrases by their length descending:
            high_scoring_phrases_tuples.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        }
        // TODO: If several high scoring phrases have the same length, apply some other measure to
        // select among these the best fitting one. For now, the one appearing first by chance is
        // returned.
        Some(high_scoring_phrases_tuples[0].0)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::vec;

    #[test]
    fn test_split_hits_descripts() {
        let hit_words = "alcohol dehydrogenase c terminal".to_string();
        let result = vec!["alcohol", "dehydrogenase", "c", "terminal"];
        assert_eq!(result, split_descriptions(&hit_words));
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
    fn test_frequencies() {
        let vec = vec![
            "alcohol".to_string(),
            "dehydrogenase".to_string(),
            "c".to_string(),
            "terminal".to_string(),
        ];
        let mut result = HashMap::new();
        result.insert("terminal".to_string(), 0.25);
        result.insert("dehydrogenase".to_string(), 0.25);
        result.insert("alcohol".to_string(), 0.25);
        result.insert("c".to_string(), 0.25);

        assert_eq!(result, frequencies(&vec));
    }

    #[test]
    fn test_inverse_information_content() {
        let word = "a".to_string();
        let mut freq_map = HashMap::new();
        freq_map.insert("a".to_string(), 3. as f32);
        freq_map.insert("b".to_string(), 2. as f32);
        freq_map.insert("c".to_string(), 2. as f32);
        freq_map.insert("d".to_string(), 1. as f32);
        freq_map.insert("e".to_string(), 1. as f32);
        freq_map.insert("f".to_string(), 1. as f32);

        let pw: f32 = freq_map.values().into_iter().sum();

        let result = f32::log(1. / (1. - 3. / pw), 2.71828182846);

        assert_eq!(result, inverse_information_content(&word, &freq_map))
    }

    #[test]
    fn test_mean() {
        let mut score_map = HashMap::new();
        score_map.insert("a".to_string(), 3. as f32);
        score_map.insert("b".to_string(), 2. as f32);
        score_map.insert("c".to_string(), 2. as f32);
        score_map.insert("d".to_string(), 3. as f32);
        score_map.insert("e".to_string(), 1. as f32);
        score_map.insert("f".to_string(), 1. as f32);
        let expected = 2. as f32;
        assert_eq!(expected, mean(&score_map));
    }

    #[test]
    fn test_scores_phrases() {
        let mut phrases: Vec<Vec<String>> = vec![];
        phrases.push(vec!["a".to_string()]);
        phrases.push(vec!["b".to_string()]);
        phrases.push(vec!["c".to_string()]);
        phrases.push(vec!["d".to_string()]);
        phrases.push(vec!["e".to_string()]);
        phrases.push(vec!["f".to_string()]);

        let mut freq_map: HashMap<String, f32> = HashMap::new();
        freq_map.insert("a".to_string(), 3. as f32);
        freq_map.insert("b".to_string(), 2. as f32);
        freq_map.insert("c".to_string(), 2. as f32);
        freq_map.insert("d".to_string(), 1. as f32);
        freq_map.insert("e".to_string(), 1. as f32);
        freq_map.insert("f".to_string(), 1. as f32);

        let expected: Vec<f32> = vec![
            0.1701677 as f32,
            -0.08114673 as f32,
            -0.08114673 as f32,
            0.036636263 as f32,
            0.036636263 as f32,
            -0.08114673 as f32,
        ];
        assert_eq!(expected, score_phrases(&freq_map, &phrases));
    }

    #[test]
    fn test_find_best_scoring_phrase() {
        let phrases: Vec<Vec<String>> = vec![
            vec!["dehydrogenase".to_string(), "c".to_string()],
            vec![
                "dehydrogenase".to_string(),
                "c".to_string(),
                "terminal".to_string(),
            ],
            vec!["alcohol".to_string(), "dehydrogenase".to_string()],
            vec![
                "cinnamyl".to_string(),
                "alcohol".to_string(),
                "dehydrogenase".to_string(),
            ],
        ];
        let phrase_scores: Vec<f32> = vec![0.265823, 0.265823, 0.22149113, -0.110522866];

        assert_eq!(
            1 as usize,
            find_best_scoring_phrase(&phrases, &phrase_scores).unwrap()
        );
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
        let expected = "alcohol dehydrogenase c terminal".to_string();
        let result = generate_human_readable_description(&hit_hrds).unwrap();

        assert_eq!(expected, result);
    }
}
