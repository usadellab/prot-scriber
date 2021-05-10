#' Calculates the prot.scriber score of a phrase as the sum of word scores.
#'
#' @param phrase - A character vector of words representing the phrase for
#' which to compute the score.
#' @param wrd.scores - A named numeric vector holding the scores of each word,
#' result of calling centerWordScores and inverseInformationContent on the
#' respective word-set.
#'
#' @return A numeric scalar value representing the score of the argument
#' phrase.
#' @export
phraseScore <- function(phrase, wrd.scores) {
    sum(wrd.scores[tolower(phrase)])
}

#' Computes the score of a word using 'inverse information content' calculated
#' as -log( (1-p.word) ).
#'
#' @param word - A string representing the word
#' @param wrd.frequencies - An instance of base::table holding the frequencies
#' of the dictionary (all words).
#' @param log.base - A scalar numeric value to be used as the base for the
#' logarithm. Default is getOption('inverseInformationContent.log.base',
#' exp(1))
#'
#' @return A scalar numeric value representing the score of the argument word.
#' @export
inverseInformationContent <- function(word, wrd.frequencies, 
    log.base = getOption("inverseInformationContent.log.base", 
        exp(1))) {
    x.s <- sum(wrd.frequencies)
    -log(1 - wrd.frequencies[[tolower(word)]]/x.s, base = log.base)
}

#' Computes the centered word scores using inverse information content and
#' centering.
#'
#' @param word.freqs - An instance of base::table holding the words'
#' frequencies. Words are expected to be in lower-case.
#' @param level.funk - The function to be used to calculate the value at which
#' to center the word scores. Default is
#' getOption('centerWordScores.level.funk', mean)
#'
#' @return A named numeric vector of same length as argument 'word.freqs'
#' holding the now centered scores. Return 'numeric(0)' if and only if
#' 'length(word.freqs) == 0'.
#' @export
centeredWordScores <- function(word.freqs, level.funk = getOption("centerWordScores.level.funk", 
    mean)) {
    if (length(word.freqs) > 0) {
        wrd.scores <- sapply(names(word.freqs), inverseInformationContent, 
            wrd.frequencies = word.freqs)
        wrd.scores - level.funk(wrd.scores)
    } else numeric(0)
}

#' Computes the centered word scores using the word frequency directly as
#' word-score. Centers the resulting scores.
#'
#' @param word.freqs - An instance of base::table holding the words'
#' frequencies. Words are expected to be in lower-case.
#' @param level.funk - The function to be used to calculate the value at which
#' to center the word scores. Default is
#' getOption('centerWordScores.level.funk', mean)
#'
#' @return A named numeric vector of same length as argument 'word.freqs'
#' holding the now centered scores. Return 'numeric(0)' if and only if
#' 'length(word.freqs) == 0'.
#' @export
centeredLinearWordScores <- function(word.freqs, level.funk = getOption("centerWordScores.level.funk", 
    mean)) {
    if (length(word.freqs) > 0) {
        wrd.scores <- word.freqs/sum(word.freqs)
        #' Just to be safe, convert to a named numeric vector:
        setNames(as.numeric(wrd.scores - level.funk(wrd.scores)), 
            names(wrd.scores))
    } else numeric(0)
}


#' A simple polynomic function used to compute centered scores for
#' probabilities p. If p < 0.5 the score will be negative, zero for p == 0, and
#' postive otherwise.
#'
#' @param p.x - A scalar or vector of numeric values (the probabilities or
#' frequencies) for which to compute a score.
#'
#' @return A numeric vector or scalar depending on the input, representing the
#' calculated score.
#' @export
polynomicScore <- function(p.x) {
    (2 * p.x - 1)^3
}

#' Computes the centered word scores using a polynomic function of the word
#' frequency f(p(x)) = (2*x-1)^3
#'
#' @param word.freqs - An instance of base::table holding the words'
#' frequencies. Words are expected to be in lower-case.
#' @param polynomic.funk - The function to be used to calculate the word
#' scores. Default is getOption('polynomicWordScores.polynomic.funk',
#' polynomicScore).
#'
#' @return A named numeric vector of same length as argument 'wrd.scores'
#' holding the polynomial word scores. Return 'numeric(0)' if and only if
#' 'length(word.freqs) == 0'.
#' @export
polynomicWordScores <- function(word.freqs, polynomic.funk = getOption("polynomicWordScores.polynomic.funk", 
    polynomicScore)) {
    if (length(word.freqs) > 0) {
        setNames(polynomic.funk(word.freqs), names(word.freqs))
    } else numeric(0)
}

#' Simple function calculating a phrase's score by summing up the scores of
#' contained words. Respective word scores a obtained from a named numeric
#' vector (see argument 'word.scores'). 
#'
#' @param phrase - A character vector representing a phrase
#' @param word.scores - A named numeric vector in which names represent words
#' and values their respective scores.
#'
#' @return A scalar numeric value, the score of the argument 'phrase' computed
#' as the sum of scores of its contained words.
#' @export
sumWordScores <- function(phrase, word.scores) {
    sum(word.scores[tolower(phrase)])
}
