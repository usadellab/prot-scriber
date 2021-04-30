#' Splits a source description into words, filters for informative ones (see
#' `wordSet`), and subsequently builds all possible phrases out of these
#' informative words. A phrase is a sub-set taken from all informative words
#' retaining the order of appearence in the source description.
#'
#' @param source.description - A string to be used as source for the generation
#' of sub-phrases.
#' @param split.regex - A regular expression used to split the argument HRD
#' (see base::strsplit). Default is
#' getOption("splitDescriptionIntoWordSet.split.regex", "\\s+|\\.|,|\\(|\\)").
#' See `wordSet` for more details.
#'
#' @return A list of named integer vectors. Each vector holds the indices of
#' the words selected from the argument `source.description`. The vector's
#' names are the actual words. So, in order to make a description out of such
#' an word-index-vector `wiv` execute `paste( names(wiv), collapse = " ")`.
#' @export
phrases <- function(source.description, split.regex = getOption("splitDescriptionIntoWordSet.split.regex", 
    "\\s+|\\.|,|\\(|\\)")) {
    source.words <- strsplit(source.description, split = split.regex)[[1]]
    sw.lower <- tolower(source.words)
    informative.words <- wordSet(source.description, split.regex = split.regex)
    #' Note that if a word appears multiple times, only the first appearance is
    #' considered (see below `[[1]]`):
    iw.inds <- sapply(informative.words, function(i.w) which(sw.lower == 
        i.w)[[1]])
    #' Generate all possible sub-phrases:
    res.phrs <- list()
    for (m in 1:length(iw.inds)) {
        m.word.inds <- combn(x = iw.inds, m = m, simplify = TRUE)
        for (k in 1:ncol(m.word.inds)) {
            select.i <- as.integer(m.word.inds[, k])
            res.phrs[[length(res.phrs) + 1]] <- setNames(select.i, 
                source.words[select.i])
        }
    }
    res.phrs
}

#' Tests its namesake `phrases`.
#'
#' @return TRUE if and only if all tests pass.
#' @export
testPhrases <- function() {
    split.regex <- "\\s+|\\.|,|\\(|\\)"
    source.description <- "Putative hypothetical Chlorophyll a-b binding protein 5, chloroplastic (Fragment)"
    expected.phrases <- c("Chlorophyll", "a-b", "binding", "5", 
        "chloroplastic", "Chlorophyll a-b", "Chlorophyll binding", 
        "Chlorophyll 5", "Chlorophyll chloroplastic", "a-b binding", 
        "a-b 5", "a-b chloroplastic", "binding 5", "binding chloroplastic", 
        "5 chloroplastic", "Chlorophyll a-b binding", "Chlorophyll a-b 5", 
        "Chlorophyll a-b chloroplastic", "Chlorophyll binding 5", 
        "Chlorophyll binding chloroplastic", "Chlorophyll 5 chloroplastic", 
        "a-b binding 5", "a-b binding chloroplastic", "a-b 5 chloroplastic", 
        "binding 5 chloroplastic", "Chlorophyll a-b binding 5", 
        "Chlorophyll a-b binding chloroplastic", "Chlorophyll a-b 5 chloroplastic", 
        "Chlorophyll binding 5 chloroplastic", "a-b binding 5 chloroplastic", 
        "Chlorophyll a-b binding 5 chloroplastic")
    res.phrases <- phrases(source.description, split.regex = split.regex)
    rp.strings <- c()
    for (x in res.phrases) {
        rp.strings[[length(rp.strings) + 1]] <- paste(names(x), 
            collapse = " ")
    }
    length(union(expected.phrases, rp.strings)) == length(expected.phrases)
}
