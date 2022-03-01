use super::default::NON_INFORMATIVE_WORD_SCORE;
use super::model_funcs::matches_blacklist;
use regex::Regex;
use statrs::statistics::{Data, OrderStatistics};
use std::collections::HashMap;

/// Main function for generating human-readable descriptions (hrds).
///
/// # TODO
///
/// If one of the following 'manitol dehydrogenase' is removed the phrases 'manitol dehydrogenase'
/// *and* 'geraniol dehydrogenase' will receive identical scores. Currently, in such cases their
/// order of appearance decides, which is chosen as a result. Meaning in these cases the result is
/// chosen somewhat randomly. Find a solution for this!
///
/// # Arguments
///
/// * `hit_hrds: &Vec<String>` - A vector of strings containing all Hit descriptions.
/// * `split_regex` - The regular expression used to split descriptions (parsed `stitle`) into
/// vectors of words (`String`).
/// * `non_informative_words_regexs` - A reference to a vector holding regular expressions used to
/// identify non informative words, that receive only a minimum score.
pub fn generate_human_readable_description(
    descriptions: &Vec<String>,
    split_regex: &Regex,
    non_informative_words_regexs: &Vec<Regex>,
) -> Option<String> {
    // Initialize default result:
    let mut human_readable_rescription_result: Option<String> = None;

    if descriptions.len() > 0 {
        // Split the descriptions into vectors of words:
        let description_words: Vec<Vec<String>> = descriptions
            .iter()
            .map(|dsc| split_descriptions(dsc, split_regex))
            .collect();

        // The universe if informative words, maintaining the word-frequencies:
        let mut informative_words_universe: Vec<String> = vec![];
        for desc_words in &description_words {
            for word in desc_words {
                // Build the word universe for later calculation of word-frequencies, but only consider
                // words that are not classified as non-informative. Note that if a word already is
                // contained in the universe, it has passed the blacklist in a past iteration, so we
                // don't need to check again:
                if informative_words_universe.contains(&word)
                    || !matches_blacklist(&word, non_informative_words_regexs)
                {
                    informative_words_universe.push(word.clone());
                }
            }
        }

        // Calculate the frequency of the informative universe words:
        let word_frequencies = frequencies(&informative_words_universe);
        let ciic: HashMap<String, f64> = centered_inverse_information_content(&word_frequencies);

        println!("\n\nciic\n{:?}\n\n", ciic);

        // Find highest scoring phrase
        let mut phrases: Vec<(Vec<String>, f64)> = vec![];
        for desc in &description_words {
            let hsp_option = highest_scoring_phrase(&desc, &ciic);
            match hsp_option {
                Some(hsp) => {
                    phrases.push(hsp);
                }
                None => {}
            }
        }
        if phrases.len() > 0 {
            let mut high_score_ind: usize = 0;
            for i in 0..phrases.len() {
                if phrases[i].1 > phrases[high_score_ind].1 {
                    high_score_ind = i;
                }
            }

            let human_readable_description: String = phrases[high_score_ind].0.join(" ");
            human_readable_rescription_result = Some(human_readable_description);
        }
    }
    human_readable_rescription_result
}

/// Find the highest scoring "phrase" in argument `description`. A phrase is a subset of the
/// argument vector maintaining the order of elements. The highest scoring phrase is found using
/// the linear solution to the longest, or highest scoring, path problem in directed acyclic
/// graphs. The argument `description` is converted into a graph, in which each word has edges to
/// all words appearing after it in the `description`. An `Option<(Vec<String>, f64)>` is returned
/// holding the highest scoring phrase and that phrase's score.
///
/// # Arguments
///
/// * `description` - A vector of words representing the description for which to find the best
/// scoring phrase.
/// * `ciic` - A reference to a HashMap holding the centered inverse information content scores for
/// the informative words appearing in the argument `description`.
pub fn highest_scoring_phrase(
    description: &Vec<String>,
    ciic: &HashMap<String, f64>,
) -> Option<(Vec<String>, f64)> {
    // Initialize the default result:
    let mut result: Option<(Vec<String>, f64)> = None;
    // There's only work to do, if the argument `description` _has_ words:
    if description.len() > 0 {
        // Each word in argument `description` is a vertex in a directed acyclic graph (DAG). An
        // additional start vertex (index 0) is added that has edges to all words:
        let n_vertices = description.len() + 1;
        // Initialize backtracing for dynamic programming; that is the highest scoring path through the
        // word DAG:
        let mut path_predecessors: Vec<usize> = vec![0; n_vertices];
        // The score of the highest scoring path to each vertex _i_ is stored in this vector:
        let mut vertex_path_scores: Vec<f64> = vec![0.0; n_vertices];
        for vertex_indx in 0..n_vertices {
            let v_edges_to_descendants: Vec<usize> = if vertex_indx == n_vertices - 1 {
                // Last word in argument `description`
                vec![]
            } else {
                // Any word i in argument `description` has edges to all words k>i following it in
                // `description`:
                ((vertex_indx + 1)..n_vertices).collect()
            };
            for desc_vertex_indx in v_edges_to_descendants {
                // The word matching the vertex is index minus one, because we inserted a start vertex
                // at index zero:
                let desc_vertex = &description[desc_vertex_indx - 1];
                // Label edges with the score of the word (vertex) the respective edge leads to. If it
                // is an informative word, lookup its score, otherwise use the minimum default score
                // for non-informative words:
                let edge_label: f64 = if ciic.contains_key(desc_vertex) {
                    *ciic.get(desc_vertex).unwrap()
                } else {
                    *NON_INFORMATIVE_WORD_SCORE
                };
                // Set the score of the path to the currently processed vertex (word):
                if vertex_path_scores[desc_vertex_indx]
                    <= vertex_path_scores[vertex_indx] + edge_label
                {
                    vertex_path_scores[desc_vertex_indx] =
                        vertex_path_scores[vertex_indx] + edge_label;
                    path_predecessors[desc_vertex_indx] = vertex_indx;
                }
            }
        }
        // Find the path that yielded the highest score:
        let mut max_path_score_indx: usize = 0;
        for i in 1..vertex_path_scores.len() {
            if vertex_path_scores[i] > vertex_path_scores[max_path_score_indx] {
                max_path_score_indx = i;
            }
        }
        if max_path_score_indx > 0 {
            // Backtrace using dynamic programming the path with the highest score:
            let mut high_score_path: Vec<String> = vec![];
            let mut next_pred_indx: usize = max_path_score_indx;
            loop {
                // Get the word matching the vertex index `next_pred_indx` by subtracting one from it. This
                // needs to be done, because we inserted a start vertex with index zero:
                high_score_path.push(description[next_pred_indx - 1].clone());
                next_pred_indx = path_predecessors[next_pred_indx];
                if next_pred_indx == 0 {
                    break;
                }
            }
            // Highest scoring phrase and it's score:
            result = Some((
                high_score_path.into_iter().rev().collect(),
                vertex_path_scores[max_path_score_indx],
            ));
        }
    }
    result
}

/// Given filtered Hit descriptions it splits each word and returns a vector.
///
/// # Arguments
///
/// * `description` - A reference to the parsed `stitle` to be split into words
/// * `split_regex` - A reference to the regular expression to be used to split the argument
/// `description` into words.
pub fn split_descriptions(description: &String, split_regex: &Regex) -> Vec<String> {
    // Split the description using a simple regular expression:
    split_regex
        .split(&description.trim())
        .map(|wrd| wrd.to_string())
        .filter(|x| !x.is_empty())
        .collect()
}

/// Calculates the word frequencies for argument `universe_words` and returns a `HashMap<String,
/// f64>` mapping the words to their respective frequency. Note that this functions returns
/// absolute frequencies in terms of number of appearances.
///
/// # Arguments
///
/// * `universe_words: &Vector<String>` - vector of words
pub fn frequencies(universe_words: &Vec<String>) -> HashMap<String, f64> {
    let mut word_freqs: HashMap<String, f64> = HashMap::new();
    for word in universe_words.iter() {
        if !word_freqs.contains_key(word) {
            let n_appearances = universe_words.iter().filter(|x| (*x) == word).count() as f64;
            word_freqs.insert((*word).clone(), n_appearances);
        }
    }
    word_freqs
}

/// Computes the score of the informative words in argument `wrd_frequencies.keys()` using 'inverse
/// information content' calculated as `-1 * log(1 - probability(word))`, where 'probability' =
/// frequency tanges between zero and one. In order to avoid infinite values for a word that is the
/// single element of the word-set, i.e. it has a frequency of one, the score of one is used.
/// Returns a HashMap of word centered IIC key-value-pairs (`HashMap<String, f64>`).
///
/// # Arguments
///
/// * `wrd.frequencies` - An instance of dictionary of all words with their frequencies.
pub fn centered_inverse_information_content(
    wrd_frequencies: &HashMap<String, f64>,
) -> HashMap<String, f64> {
    // Initialize default result:
    let mut ciic_result: HashMap<String, f64> = HashMap::new();

    if wrd_frequencies.len() > 0 {
        // Calculate inverse information content (IIC):
        let sum_wrd_frequencies: f64 = wrd_frequencies.values().into_iter().sum();
        let mut inv_inf_cntnt: Vec<(String, f64)> = vec![];
        for word in wrd_frequencies.keys() {
            if wrd_frequencies.len() as f64 > 1. {
                let pw = wrd_frequencies[word] / sum_wrd_frequencies;
                let iic: f64 = -1.0 * f64::log(1. - pw, std::f64::consts::E);
                inv_inf_cntnt.push((word.to_string(), iic));
            } else {
                inv_inf_cntnt.push((word.to_string(), 1.0));
            }
        }

        // Center inverse information content (IIC) values, if and only if there is variation
        // between the calculated IIC values. Variation can only result from varying frequencies,
        // so find out if the argument `wrd_frequencies` contains such values:
        let wrd_frequency_vals: Vec<f64> = wrd_frequencies.values().map(|pw| *pw).collect();
        let mut iic_values_all_identical: (bool, f64) = (true, wrd_frequency_vals[0]);
        for i in 1..wrd_frequency_vals.len() {
            iic_values_all_identical.0 = iic_values_all_identical.1 == wrd_frequency_vals[i];
            if !iic_values_all_identical.0 {
                // Once a comparison was false, we _must not_ compare more pairs, because if the
                // last pair is in fact identical the boolean result would not be correct:
                break;
            }
            iic_values_all_identical.1 = wrd_frequency_vals[i];
        }

        // Calculate mean inverse information content for centering:
        let mut subtract_4_centering = 0.0;
        // Note that only in case of variance between IIC values, we calculate the indicated
        // quantile IIC to be subtracted from the actual IIC for centering. Otherwise the above
        // default zero will be subtracted:
        if !iic_values_all_identical.0 {
            subtract_4_centering = word_scores_quantile(&inv_inf_cntnt, 0.5);
        }
        // Center inverse information content:
        for word_iic_tuple in inv_inf_cntnt {
            ciic_result.insert(word_iic_tuple.0, word_iic_tuple.1 - subtract_4_centering);
        }
    }

    ciic_result
}

/// Computes and returns the argument `quantile` score of an argument word-score vector `values`.
/// The method used is explained here: https://www-users.york.ac.uk/~mb55/intro/quantile.htm
///
/// # Arguments
///
/// * `values` - A reference to a word-score vector
/// * `tau` - A value between 0.0 and 1.0 indicating the quantile to calculate
pub fn word_scores_quantile(values: &Vec<(String, f64)>, tau: f64) -> f64 {
    if tau < 0.0 || tau > 1.0 {
        panic!(
            "Cannot compute quantile {:?} because it is not a valid value between zero and one (inclusive).",
            &tau
        );
    }
    let scores: Vec<f64> = values.iter().map(|(_, s)| (*s)).collect();
    let mut scores_data = Data::new(scores);
    scores_data.quantile(tau)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::default::{NON_INFORMATIVE_WORDS_REGEXS, SPLIT_DESCRIPTION_REGEX};
    use assert_approx_eq::assert_approx_eq;
    use std::vec;

    #[test]
    fn test_word_scores_quantile() {
        let input = vec![
            ("".to_string(), 1.0),
            ("".to_string(), 2.0),
            ("".to_string(), 3.0),
            ("".to_string(), 4.0),
        ];
        assert_eq!(word_scores_quantile(&input, 0.5), 2.5);
        assert_approx_eq!(word_scores_quantile(&input, 1.0 / 3.0), 1.777777, 1e-5f64);
        assert_approx_eq!(
            word_scores_quantile(&input, 0.4),
            2.0 + 0.1 * 2.0 / 3.0,
            1e-5f64
        );
        assert_approx_eq!(
            word_scores_quantile(&input, 0.75),
            3.58 + 0.01 * 1.0 / 3.0,
            1e-5f64
        );
        assert_approx_eq!(word_scores_quantile(&input, 0.8), 3.8, 1e-5f64);
        assert_eq!(word_scores_quantile(&input, 1.0), 4.0);
        assert_eq!(word_scores_quantile(&input, 0.0), 1.0);
    }

    #[test]
    fn test_split_descriptions() {
        // Test 1:
        let hit_words = "alcohol dehydrogenase c terminal".to_string();
        let expected = vec!["alcohol", "dehydrogenase", "c", "terminal"];
        assert_eq!(
            expected,
            split_descriptions(&hit_words, &(*SPLIT_DESCRIPTION_REGEX))
        );
    }

    #[test]
    fn test_frequencies() {
        let mut words = vec![
            "alcohol".to_string(),
            "dehydrogenase".to_string(),
            "c".to_string(),
            "terminal".to_string(),
        ];
        let mut expected = HashMap::new();
        expected.insert("terminal".to_string(), 1.0);
        expected.insert("dehydrogenase".to_string(), 1.0);
        expected.insert("alcohol".to_string(), 1.0);
        expected.insert("c".to_string(), 1.0);
        assert_eq!(expected, frequencies(&words));

        words = vec![
            "importin".to_string(),
            "5".to_string(),
            "importin".to_string(),
            "5".to_string(),
            "ran".to_string(),
            "binding".to_string(),
            "6".to_string(),
            "ran".to_string(),
            "binding".to_string(),
            "6".to_string(),
            "importin".to_string(),
            "subunit".to_string(),
            "beta".to_string(),
            "3".to_string(),
            "importin".to_string(),
            "subunit".to_string(),
            "beta".to_string(),
            "3".to_string(),
        ];
        expected = HashMap::new();
        expected.insert("subunit".to_string(), 2.0);
        expected.insert("5".to_string(), 2.0);
        expected.insert("binding".to_string(), 2.0);
        expected.insert("3".to_string(), 2.0);
        expected.insert("ran".to_string(), 2.0);
        expected.insert("6".to_string(), 2.0);
        expected.insert("beta".to_string(), 2.0);
        expected.insert("importin".to_string(), 4.0);
        assert_eq!(expected, frequencies(&words));
    }

    #[test]
    fn test_centered_inverse_information_content() {
        let mut freq_map = HashMap::new();
        freq_map.insert("a".to_string(), 3. as f64);
        freq_map.insert("b".to_string(), 2. as f64);
        freq_map.insert("c".to_string(), 2. as f64);
        freq_map.insert("d".to_string(), 1. as f64);
        freq_map.insert("e".to_string(), 1. as f64);
        freq_map.insert("f".to_string(), 1. as f64);

        let mut freq_sum: f64 = freq_map.values().into_iter().sum();

        let mut expected: HashMap<String, f64> = HashMap::new();
        expected.insert(
            "a".to_string(),
            -1. * f64::log(1. - 3. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "b".to_string(),
            -1. * f64::log(1. - 2. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "c".to_string(),
            -1. * f64::log(1. - 2. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "d".to_string(),
            -1. * f64::log(1. - 1. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "e".to_string(),
            -1. * f64::log(1. - 1. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "f".to_string(),
            -1. * f64::log(1. - 1. / freq_sum, std::f64::consts::E),
        );
        let iic_scores: Vec<f64> = expected.iter().map(|(_, v)| *v).collect();
        let mut iic_data = Data::new(iic_scores);
        let mean_ciic: f64 = iic_data.quantile(0.5);
        // center the expected IIC:
        let mut centered_expected: HashMap<String, f64> = HashMap::new();
        for (word, iic) in &expected {
            centered_expected.insert((*word).clone(), iic - mean_ciic);
        }

        // test iteratively:
        let result: HashMap<String, f64> = centered_inverse_information_content(&freq_map);
        for (word, _) in &centered_expected {
            assert_approx_eq!(
                centered_expected.get(word).unwrap(),
                result.get(word).unwrap(),
                1e-6f64
            );
        }

        // Test 2:
        freq_map = HashMap::new();
        freq_map.insert("alcohol".to_string(), 2.0);
        freq_map.insert("terminal".to_string(), 2.0);
        freq_map.insert("geraniol".to_string(), 2.0);
        freq_map.insert("manitol".to_string(), 3.0);
        freq_map.insert("dehydrogenase".to_string(), 7.0);
        freq_map.insert("c".to_string(), 1.0);
        freq_map.insert("cinnamyl".to_string(), 1.0);
        freq_sum = freq_map.values().into_iter().sum();
        expected = HashMap::new();
        expected.insert(
            "alcohol".to_string(),
            -1. * f64::log(1. - 2. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "terminal".to_string(),
            -1. * f64::log(1. - 2. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "geraniol".to_string(),
            -1. * f64::log(1. - 2. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "manitol".to_string(),
            -1. * f64::log(1. - 3. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "dehydrogenase".to_string(),
            -1. * f64::log(1. - 7. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "c".to_string(),
            -1. * f64::log(1. - 1. / freq_sum, std::f64::consts::E),
        );
        expected.insert(
            "cinnamyl".to_string(),
            -1. * f64::log(1. - 1. / freq_sum, std::f64::consts::E),
        );
        let iic_scores: Vec<f64> = expected.iter().map(|(_, v)| *v).collect();
        let mut iic_data = Data::new(iic_scores);
        let mean_ciic: f64 = iic_data.quantile(0.5);
        // center the expected IIC:
        let mut centered_expected: HashMap<String, f64> = HashMap::new();
        for (word, iic) in &expected {
            centered_expected.insert((*word).clone(), iic - mean_ciic);
        }
        // test iteratively:
        let result: HashMap<String, f64> = centered_inverse_information_content(&freq_map);
        for (word, _) in &centered_expected {
            assert_approx_eq!(
                centered_expected.get(word).unwrap(),
                result.get(word).unwrap(),
                1e-6f64
            );
        }

        // test special case of all equally frequent words:
        freq_map = HashMap::new();
        freq_map.insert("foo".to_string(), 1.0);
        freq_map.insert("bar".to_string(), 1.0);
        freq_map.insert("baz".to_string(), 1.0);
        centered_expected = HashMap::new();
        // All words should have this NON CENTERED inverse information content:
        let iic: f64 = -1.0 * f64::log(1. - 1. / 3., std::f64::consts::E);
        centered_expected.insert("foo".to_string(), iic);
        centered_expected.insert("bar".to_string(), iic);
        centered_expected.insert("baz".to_string(), iic);
        assert_eq!(
            centered_expected,
            centered_inverse_information_content(&freq_map)
        );
    }

    #[test]
    fn test_highest_scoring_phrase() {
        let desc1: Vec<String> = vec!["importin".to_string(), "5".to_string()];
        let desc2: Vec<String> = vec![
            "ran".to_string(),
            "binding".to_string(),
            "protein".to_string(),
            "6".to_string(),
        ];
        let desc3: Vec<String> = vec!["ran".to_string(), "6".to_string()];
        let desc4: Vec<String> = vec![
            "protein".to_string(),
            "narrow".to_string(),
            "leaf".to_string(),
            "1".to_string(),
        ];

        let mut word_freqs: HashMap<String, f64> = HashMap::new();
        word_freqs.insert("6".to_string(), 2.0);
        word_freqs.insert("importin".to_string(), 5.0);
        word_freqs.insert("ran".to_string(), 2.0);
        word_freqs.insert("3".to_string(), 2.0);
        word_freqs.insert("subunit".to_string(), 2.0);
        word_freqs.insert("beta".to_string(), 2.0);
        word_freqs.insert("5".to_string(), 3.0);
        word_freqs.insert("binding".to_string(), 2.0);

        let mut ciic = centered_inverse_information_content(&word_freqs);

        let phrase1 = highest_scoring_phrase(&desc1, &ciic).unwrap();
        let expected1 = vec!["importin".to_string(), "5".to_string()];
        assert_eq!(expected1, phrase1.0);

        let phrase2 = highest_scoring_phrase(&desc2, &ciic).unwrap();
        let expected2 = vec![
            "ran".to_string(),
            "binding".to_string(),
            "protein".to_string(),
        ];
        assert_eq!(expected2, phrase2.0);

        let phrase3 = highest_scoring_phrase(&desc3, &ciic);
        assert!(phrase3.is_none());

        word_freqs = HashMap::new();
        for word in &desc4 {
            word_freqs.insert(word.clone(), 1.0);
        }
        ciic = centered_inverse_information_content(&word_freqs);
        let phrase4 = highest_scoring_phrase(&desc4, &ciic).unwrap();
        // Expect the full input description to be replicated:
        assert_eq!(desc4, phrase4.0);

        let desc5 = vec![
            "receptor".to_string(),
            "protein".to_string(),
            "eix2".to_string(),
        ];
        word_freqs = HashMap::new();
        word_freqs.insert("receptor".to_string(), 2.0);
        word_freqs.insert("eix1".to_string(), 1.0);
        word_freqs.insert("eix2".to_string(), 1.0);
        ciic = centered_inverse_information_content(&word_freqs);
        let phrase5 = highest_scoring_phrase(&desc5, &ciic).unwrap();
        assert_eq!(
            vec!["receptor".to_string(), "protein".to_string()],
            phrase5.0
        );
    }

    #[test]
    fn test_generate_human_readable_description() {
        // Test 1:
        // TODO: If one of the following 'manitol dehydrogenase' is removed the phrases 'manitol
        // dehydrogenase' *and* 'geraniol dehydrogenase' will receive identical scores. Currently,
        // in such cases their order of appearance decides, which is chosen as a result. Meaning
        // in these cases the result is chosen somewhat randomly. Find a solution for this!
        let mut hit_hrds = vec![
            "manitol dehydrogenase".to_string(),
            "cinnamyl alcohol-dehydrogenase".to_string(),
            "geraniol dehydrogenase".to_string(),
            "geraniol|dehydrogenase terminal".to_string(),
            "manitol dehydrogenase".to_string(),
            "manitol dehydrogenase".to_string(),
            "alcohol dehydrogenase c-terminal".to_string(),
        ];
        let mut expected = "manitol dehydrogenase".to_string();
        let mut result = generate_human_readable_description(
            &hit_hrds,
            &(*SPLIT_DESCRIPTION_REGEX),
            &(*NON_INFORMATIVE_WORDS_REGEXS),
        )
        .unwrap();
        assert_eq!(expected, result);

        // Test 2:
        hit_hrds = vec![
            "importin-5".to_string(),
            "importin-5".to_string(),
            "importin-5".to_string(),
            "ran-binding protein 6".to_string(),
            "ran-binding protein 6".to_string(),
            "importin subunit beta-3".to_string(),
            "importin subunit beta-3".to_string(),
        ];
        expected = "importin 5".to_string();
        result = generate_human_readable_description(
            &hit_hrds,
            &(*SPLIT_DESCRIPTION_REGEX),
            &(*NON_INFORMATIVE_WORDS_REGEXS),
        )
        .unwrap();
        assert_eq!(expected, result);

        // Test 3:
        hit_hrds = vec![
            "receptor protein eix1".to_string(),
            "receptor protein eix2".to_string(),
        ];
        expected = "receptor protein".to_string();
        result = generate_human_readable_description(
            &hit_hrds,
            &(*SPLIT_DESCRIPTION_REGEX),
            &(*NON_INFORMATIVE_WORDS_REGEXS),
        )
        .unwrap();
        assert_eq!(expected, result);
    }
}
