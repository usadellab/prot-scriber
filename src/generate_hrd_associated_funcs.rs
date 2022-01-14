use super::default::{
    BLACKLIST_DESCRIPTION_WORDS_REGEXS, NON_INFORMATIVE_WORD_SCORE, SPLIT_DESCRIPTION_REGEX,
};
use super::model_funcs::matches_blacklist;
use std::collections::HashMap;
use std::collections::HashSet;

/// Main function for generating human-readable descriptions (hrds).
///
/// # Arguments
///
/// * `hit_hrds: &Vec<String>` - A vector of strings containing all Hit descriptions.
pub fn generate_human_readable_description(descriptions: &Vec<String>) -> Option<String> {
    if descriptions.len() > 0 {
        // Split the descriptions into vectors of words:
        let description_words: Vec<Vec<String>> = descriptions
            .iter()
            .map(|dsc| split_descriptions(dsc))
            .collect();
        // println!("\n description_words.len() = {}\n", description_words.len());

        // The universe if informative words, maintaining the word-frequencies:
        let mut informative_words_universe: Vec<String> = vec![];
        for desc_words in &description_words {
            for word in desc_words {
                // Build the word universe for later calculation of word-frequencies, but only consider
                // words that are not classified as non-informative. Note that if a word already is
                // contained in the universe, it has passed the blacklist in a past iteration, so we
                // don't need to check again:
                if informative_words_universe.contains(&word)
                    || !matches_blacklist(&word, &(*BLACKLIST_DESCRIPTION_WORDS_REGEXS))
                {
                    informative_words_universe.push(word.clone());
                }
            }
        }
        // println!("\n word_universe.len() = {}\n", word_universe.len());
        // println!("\n phrases_set.len() = {}\n", phrases_set.len());

        // Calculate the frequency of the informative universe words:
        let word_frequencies = frequencies(&informative_words_universe);
        let ciic: HashMap<String, f32> = centered_inverse_information_content(&word_frequencies);

        // Find highest scoring phrase
        let mut hsp: Vec<(Vec<String>, f32)> = vec![];
        for desc in &description_words {
            hsp.push(highest_scoring_phrase(&desc, &ciic));
        }
        let mut high_score_ind: usize = 0;
        for i in 0..(hsp.len() - 1) {
            if hsp[i].1 > hsp[high_score_ind].1 {
                high_score_ind = i;
            }
        }

        let human_readable_description: String = hsp[high_score_ind].0.join(" ");
        if human_readable_description.is_empty() {
            None
        } else {
            Some(human_readable_description)
        }
    } else {
        None
    }
}

pub fn highest_scoring_phrase(
    description: &Vec<String>,
    ciic: &HashMap<String, f32>,
) -> (Vec<String>, f32) {
    let n_vertices = description.len() + 1;
    let mut path_predecessors: Vec<usize> = vec![0; n_vertices];
    let mut length_to: Vec<f32> = vec![0.0; n_vertices];
    for vertex_indx in 0..n_vertices {
        let v_edges_to_descendants: Vec<usize> = if vertex_indx == n_vertices - 1 {
            vec![]
        } else {
            ((vertex_indx + 1)..(n_vertices - 1)).collect()
        };
        for desc_vertex_indx in v_edges_to_descendants {
            let desc_vertex = &description[desc_vertex_indx - 1];
            let edge_label: f32 = if ciic.contains_key(desc_vertex) {
                *ciic.get(desc_vertex).unwrap()
            } else {
                *NON_INFORMATIVE_WORD_SCORE
            };
            if length_to[desc_vertex_indx] <= length_to[vertex_indx] + edge_label {
                length_to[desc_vertex_indx] = length_to[vertex_indx] + edge_label;
                path_predecessors[desc_vertex_indx] = vertex_indx;
            }
        }
    }
    let mut max_length_to_indx: usize = 0;
    for i in 1..length_to.len() {
        if length_to[i] > length_to[max_length_to_indx] {
            max_length_to_indx = i;
        }
    }
    let mut high_score_path: Vec<String> = vec![];
    let mut next_pred_indx: usize = path_predecessors[max_length_to_indx];
    loop {
        high_score_path.push(description[next_pred_indx].clone());
        next_pred_indx = path_predecessors[next_pred_indx];
        if next_pred_indx == 0 {
            break;
        }
    }
    // Highest scoring phrase and it's score:
    (
        high_score_path.into_iter().rev().collect(),
        length_to[max_length_to_indx],
    )
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

/// Calculates the word frequencies for argument `universe_words` and returns a `HashMap<String,
/// f32>` mapping the words to their respective frequency. Note that this functions returns
/// absolute frequencies in terms of number of appearances.
///
/// # Arguments
///
/// * `universe_words: &Vector<String>` - vector of words
pub fn frequencies(universe_words: &Vec<String>) -> HashMap<String, f32> {
    let mut word_freqs: HashMap<String, f32> = HashMap::new();
    for word in universe_words.iter() {
        if !word_freqs.contains_key(word) {
            let n_appearances = universe_words.iter().filter(|x| (*x) == word).count() as f32;
            word_freqs.insert((*word).clone(), n_appearances);
        }
    }
    word_freqs
}

/// Computes the score of the informative words in argument `wrd_frequencies.keys()` using 'inverse
/// information content' calculated as `-1 * log(1 - probability(word))`, where 'probability' =
/// frequency elem [0,1]. In order to avoid infinite values for a word that is the single element
/// of the word-set, i.e. it has a frequency of one, the score of one is used. Returns a HashMap of
/// word centered IIC key-value-pairs (`HashMap<String, f32>`).
///
/// # Arguments
///
/// * `wrd.frequencies` - An instance of dictionary of all words with their frequencies.
pub fn centered_inverse_information_content(
    wrd_frequencies: &HashMap<String, f32>,
) -> HashMap<String, f32> {
    // Calculate inverse information content:
    let sum_wrd_frequencies: f32 = wrd_frequencies.values().into_iter().sum();
    let mut inv_inf_cntnt: Vec<(String, f32)> = vec![];
    for word in wrd_frequencies.keys() {
        if wrd_frequencies.len() as f32 > 1. {
            let pw = wrd_frequencies[word] / sum_wrd_frequencies;
            inv_inf_cntnt.push((
                word.to_string(),
                -1.0 * f32::log(1. - pw, std::f32::consts::E),
            ));
        } else {
            inv_inf_cntnt.push((word.to_string(), 1.0));
        }
    }

    // Calculate mean inverse information content for centering:
    let mut sum_iic = 0.0;
    for (_, iic) in &inv_inf_cntnt {
        sum_iic += iic;
    }
    let mean_iic = sum_iic / inv_inf_cntnt.len() as f32;
    // Center inverse information content:
    let mut centered_iic: HashMap<String, f32> = HashMap::new();
    for word_iic_tuple in inv_inf_cntnt {
        // Center and add 1.0 to enable the subsequent Graph traversal and _only_ positive scores:
        centered_iic.insert(word_iic_tuple.0, (word_iic_tuple.1 - mean_iic) + 1.0);
    }
    centered_iic
}

// #[cfg(test)]
// mod tests {
// use super::*;
// use assert_approx_eq::assert_approx_eq;
// use std::vec;
//
// #[test]
// fn test_split_hits_descripts() {
// let hit_words = "alcohol dehydrogenase c terminal".to_string();
// let result = vec!["alcohol", "dehydrogenase", "c", "terminal"];
// assert_eq!(result, split_descriptions(&hit_words));
// }
//
// #[test]
// fn test_powerset() {
// let test_vec = vec!["alcohol", "dehydrogenase", "c", "terminal"];
// let result = vec![
// vec![], // this removed in the function that puts the phrases together
// vec!["alcohol"],
// vec!["dehydrogenase"],
// vec!["alcohol", "dehydrogenase"],
// vec!["c"],
// vec!["alcohol", "c"],
// vec!["dehydrogenase", "c"],
// vec!["alcohol", "dehydrogenase", "c"],
// vec!["terminal"],
// vec!["alcohol", "terminal"],
// vec!["dehydrogenase", "terminal"],
// vec!["alcohol", "dehydrogenase", "terminal"],
// vec!["c", "terminal"],
// vec!["alcohol", "c", "terminal"],
// vec!["dehydrogenase", "c", "terminal"],
// vec!["alcohol", "dehydrogenase", "c", "terminal"],
// ];
// assert_eq!(result, powerset(&test_vec))
// }
//
// #[test]
// fn test_frequencies() {
// let vec = vec![
// "alcohol".to_string(),
// "dehydrogenase".to_string(),
// "c".to_string(),
// "terminal".to_string(),
// ];
// let mut result = HashMap::new();
// result.insert("terminal".to_string(), 1.0);
// result.insert("dehydrogenase".to_string(), 1.0);
// result.insert("alcohol".to_string(), 1.0);
// result.insert("c".to_string(), 1.0);
//
// assert_eq!(result, frequencies(&vec));
// }
//
// #[test]
// fn test_inverse_information_content() {
// let word = "a".to_string();
// let mut freq_map = HashMap::new();
// freq_map.insert("a".to_string(), 3. as f32);
// freq_map.insert("b".to_string(), 2. as f32);
// freq_map.insert("c".to_string(), 2. as f32);
// freq_map.insert("d".to_string(), 1. as f32);
// freq_map.insert("e".to_string(), 1. as f32);
// freq_map.insert("f".to_string(), 1. as f32);
//
// let pw: f32 = freq_map.values().into_iter().sum();
//
// let expected = f32::log(1. / (1. - 3. / pw), std::f32::consts::E);
//
// assert_eq!(expected, inverse_information_content(&word, &freq_map))
// }
//
// #[test]
// fn test_mean() {
// let mut score_map = HashMap::new();
// score_map.insert("a".to_string(), 3. as f32);
// score_map.insert("b".to_string(), 2. as f32);
// score_map.insert("c".to_string(), 2. as f32);
// score_map.insert("d".to_string(), 3. as f32);
// score_map.insert("e".to_string(), 1. as f32);
// score_map.insert("f".to_string(), 1. as f32);
// let expected = 2. as f32;
// assert_eq!(expected, mean(&score_map));
// }
//
// #[test]
// fn test_scores_phrases() {
// let mut phrases: Vec<Vec<String>> = vec![];
// phrases.push(vec!["a".to_string()]);
// phrases.push(vec!["b".to_string()]);
// phrases.push(vec!["c".to_string()]);
// phrases.push(vec!["d".to_string()]);
// phrases.push(vec!["e".to_string()]);
// phrases.push(vec!["f".to_string()]);
//
// // Note that we assume a universe of ten words:
// let mut freq_map: HashMap<String, f32> = HashMap::new();
// freq_map.insert("a".to_string(), 3. as f32);
// freq_map.insert("b".to_string(), 2. as f32);
// freq_map.insert("c".to_string(), 2. as f32);
// freq_map.insert("d".to_string(), 1. as f32);
// freq_map.insert("e".to_string(), 1. as f32);
// freq_map.insert("f".to_string(), 1. as f32);
//
// // All of the above phrases consist just of a single word, thus their phrase-scores should
// // be identical to the centered frequency of their word:
// let expected: Vec<f32> = vec![
// -1.0 * f32::log(1. - 0.3, std::f32::consts::E),
// -1.0 * f32::log(1. - 0.2, std::f32::consts::E),
// -1.0 * f32::log(1. - 0.2, std::f32::consts::E),
// -1.0 * f32::log(1. - 0.1, std::f32::consts::E),
// -1.0 * f32::log(1. - 0.1, std::f32::consts::E),
// -1.0 * f32::log(1. - 0.1, std::f32::consts::E),
// ];
// let mean_inv_inf_cntnt: f32 = expected.iter().sum::<f32>() / expected.len() as f32;
// let expected_centered = expected
// .iter()
// .map(|x| x - mean_inv_inf_cntnt)
// .collect::<Vec<f32>>();
// let result = score_phrases(&freq_map, &phrases);
// assert_approx_eq!(expected_centered[0], result[0], 0.00001 as f32);
// assert_approx_eq!(expected_centered[1], result[1], 1.11111 as f32);
// assert_approx_eq!(expected_centered[2], result[2], 2.22221 as f32);
// assert_approx_eq!(expected_centered[3], result[3], 3.33331 as f32);
// assert_approx_eq!(expected_centered[4], result[4], 4.44441 as f32);
// assert_approx_eq!(expected_centered[5], result[5], 5.55551 as f32);
// }
//
// #[test]
// fn test_find_best_scoring_phrase() {
// let phrases: Vec<Vec<String>> = vec![
// vec!["dehydrogenase".to_string(), "c".to_string()],
// vec![
// "dehydrogenase".to_string(),
// "c".to_string(),
// "terminal".to_string(),
// ],
// vec!["alcohol".to_string(), "dehydrogenase".to_string()],
// vec![
// "cinnamyl".to_string(),
// "alcohol".to_string(),
// "dehydrogenase".to_string(),
// ],
// ];
// let phrase_scores: Vec<f32> = vec![0.265823, 0.265823, 0.22149113, -0.110522866];
//
// assert_eq!(
// 1 as usize,
// find_best_scoring_phrase(&phrases, &phrase_scores).unwrap()
// );
// }
//
// #[test]
// fn test_generate_human_readable_description() {
// let hit_hrds = vec![
// "manitol dehydrogenase".to_string(),
// "cinnamyl alcohol-dehydrogenase".to_string(),
// "geraniol dehydrogenase".to_string(),
// "geraniol|dehydrogenase terminal".to_string(),
// "manitol dehydrogenase".to_string(),
// "alcohol dehydrogenase c-terminal".to_string(),
// ];
// let expected = "dehydrogenase".to_string();
// let result = generate_human_readable_description(&hit_hrds).unwrap();
//
// assert_eq!(expected, result);
// }
// }
