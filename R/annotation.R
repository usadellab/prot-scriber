#' Splits a source description into words, filters for informative ones (see
#' `wordSet`), and subsequently builds all possible phrases out of these
#' informative words. A phrase is a sub-set taken from all informative words
#' retaining the order of appearence in the source description.
#'
#' @param source.description - A string to be used as source for the generation
#' of sub-phrases.
#' @param split.regex - A regular expression used to split the argument HRD
#' (see base::strsplit). Default is
#' getOption('splitDescriptionIntoWordSet.split.regex',
#' '-|/|;|\\\\|,|:|\'|'|\\.|\\s+|\\||\\(|\\)'). See `wordSet` for more
#' details.
#'
#' @return A list with two entries: 'phrases' contains character vectors, each
#' of which represents a phrase composed of a subset of words extracted from
#' the source.description, maintaining their order of appearence. The second
#' entry 'source.informative.words' contains those words in the
#' source.description that passed the blacklist filtering (see `wordSet` for
#' more details). Returns NULL if and only if NO informative words are found in
#' the source.description.
#' @export
phrases <- function(source.description, split.regex = getOption("splitDescriptionIntoWordSet.split.regex", 
    "-|/|;|\\\\|,|:|\"|'|\\.|\\s+|\\||\\(|\\)")) {
    source.words <- strsplit(source.description, split = split.regex)[[1]]
    informative.words <- wordSet(source.description, split.regex = split.regex)
    if (length(informative.words) > 0) {
        sw.lower <- tolower(source.words)
        #' Note that if a word appears multiple times, only the first appearance is
        #' considered (see below `[[1]]`):
        iw.inds <- sapply(informative.words, function(i.w) which(sw.lower == 
            i.w)[[1]])
        #' Generate all possible sub-phrases:
        res.phrs <- list()
        for (m in 1:length(iw.inds)) {
            #' Need to convert iw.inds to character to avoid combn selecting
            #' out of seq(iw.inds) if iw.inds is a single integer:
            m.word.inds <- combn(x = as.character(iw.inds), m = m, 
                simplify = TRUE)
            for (k in 1:ncol(m.word.inds)) {
                #' Convert back to integer indices:
                select.i <- as.integer(m.word.inds[, k])
                res.phrs[[length(res.phrs) + 1]] <- source.words[select.i]
            }
        }
        list(phrases = res.phrs, source.informative.words = informative.words)
    } else NULL
}

#' Tests its namesake `phrases`.
#'
#' @return TRUE if and only if all tests pass.
#' @export
testPhrases <- function() {
    split.regex <- "\\s+|\\.|,|\\(|\\)"
    blacklist.words <- c("(?i)\\bunknown\\b", "(?i)\\bmember\\b", 
        "(?i)\\blike\\b", "(?i)\\bassociated\\b", "(?i)\\bcontaining\\b", 
        "(?i)\\bdomain\\b", "(?i)\\bactivated\\b", "(?i)\\bfamily\\b", 
        "(?i)\\bsubfamily\\b", "(?i)\\binteracting\\b", "(?i)\\bactivity\\b", 
        "(?i)\\bsimilar\\b", "(?i)\\bproduct\\b", "(?i)\\bexpressed\\b", 
        "(?i)\\bpredicted\\b", "(?i)\\bputative\\b", "(?i)\\bhypothetical\\b", 
        "(?i)\\buncharacterized\\b", "(?i)\\bprobable\\b", "(?i)\\bprotein\\b", 
        "(?i)\\bgene\\b", "(?i)\\btair\\b", "(?i)\\bfragment\\b", 
        "(?i)\\bhomolog\\b", "(?i)\\bcontig\\b", "(?i)\\brelated\\b", 
        "(?i)\\bremark\\b", "(?i)\\b\\w?orf(\\w?|\\d+)\\b", "(?i)\\bof\\b", 
        "(?i)\\bor\\b", "(?i)\\band\\b", "\\b\\d+\\b")
    mem.blacklist.words <- options("splitDescriptionIntoWordSet.blacklist.regexs")[[1]]
    options(splitDescriptionIntoWordSet.blacklist.regexs = blacklist.words)
    source.description <- "Putative hypothetical Chlorophyll a-b binding protein 5, chloroplastic (Fragment)"
    expected.phrases <- c("Chlorophyll", "a-b", "binding", "chloroplastic", 
        "Chlorophyll a-b", "Chlorophyll binding", "Chlorophyll chloroplastic", 
        "a-b binding", "a-b chloroplastic", "binding chloroplastic", 
        "Chlorophyll a-b binding", "Chlorophyll a-b chloroplastic", 
        "Chlorophyll binding chloroplastic", "a-b binding chloroplastic", 
        "Chlorophyll a-b binding chloroplastic")
    result <- phrases(source.description, split.regex = split.regex)
    rp.strings <- c()
    for (x in result$phrases) {
        rp.strings[[length(rp.strings) + 1]] <- paste(x, collapse = " ")
    }
    writeLines(rp.strings, "./tmp.txt")
    t1 <- length(union(expected.phrases, rp.strings)) == length(expected.phrases)
    t2 <- identical(c("chlorophyll", "a-b", "binding", "chloroplastic"), 
        result$source.informative.words)
    t3 <- is.null(phrases("putative hypothetical protein fragment", 
        split.regex = split.regex))
    #' Clean up:
    options(splitDescriptionIntoWordSet.blacklist.regexs = mem.blacklist.words)
    #' Test:
    all(t1, t2, t3)
}

#' Receives a phrases.lst, i.e. a list of depth two, in which the first level
#' names represent searched reference protein databases and the second level
#' hit-sseqid. Each hit-sseqid contains a list result from invoking 'phrases'.
#' From that lists the source.informative.words are extracted and joined in a
#' vector of all informative words found in all Hit human readable descriptions
#' (HRD). Finally, the frequency of all these words is assessed (see
#' base::table) for details.
#'
#' @param phrases.lst - A list of depth two. First level represents the
#' searched reference protein databases, and second level Hit sseqid mapped to
#' their respective results obtained from calling phrases. Use function
#' 'phrasesForQuery' to produce this argument.
#'
#' @return An instance of base::table holding the frequency of each informative
#' word in the respective Hit HRDs.
#' @export
wordFrequenciesFromPhrasesList <- function(phrases.lst) {
    source.informative.words <- c()
    if (length(phrases.lst) > 0 && !is.na(phrases.lst)) {
        for (ref.db.phrases in phrases.lst) {
            if (length(ref.db.phrases) > 0 && !is.na(ref.db.phrases)) {
                for (hit.phrases in ref.db.phrases) {
                  if (!is.null(hit.phrases) && "source.informative.words" %in% 
                    names(hit.phrases) && length(hit.phrases[["source.informative.words"]]) > 
                    0) {
                    source.informative.words <- c(source.informative.words, 
                      hit.phrases[["source.informative.words"]])
                  }
                }
            }
        }
    }
    table(source.informative.words)
}

#' Calculates the mean frequency of the words contained in the argument
#' `phrase`.
#'
#' @param phrase - A character vector representing a phrase.
#' @param wrd.frequencies - An instance of base::table holding the frequencies
#' of each word.
#'
#' @return A numeric value the computed mean frequency
#' @export
meanWordFreq <- function(phrase, wrd.frequencies) {
    n <- sum(wrd.frequencies)
    mean(sapply(phrase, function(wrd) wrd.frequencies[[tolower(wrd)]]/n))
}

#' Looks up the sequence similarity search result Hits for the argument query
#' 'prot.id'. From the Hits' human readable descriptions (HRD) phrases are
#' generated using function 'phrases'.
#'
#' @param prot.id - A string representing the query protein identifier (qseqid)
#' @param seq.sim.search.rslts - A named list of tabular sequence similarity
#' search results, e.g. using data provided with this package list(
#' 'Swissprot'=pc.sprot, 'trEMBL'=pc.trembl ). Note that the entries pc.sprot
#' and pc.trembl are generated by invoking parseMercator4Tblout on the
#' respective tabular output files.
#'
#' @return A list of depth two. First level represents the names of the
#' reference protein databases searched, each of which is assigned a list of
#' Hit (sseqid) and their respective phrases.
#' @export
phrasesForQuery <- function(prot.id, seq.sim.search.rslts) {
    if (is.null(names(seq.sim.search.rslts))) {
        stop("Argument 'seq.sim.search.rslts' MUST have the names of the searched reference protein databases.")
    }
    seq.sim.phrases <- setNames(lapply(names(seq.sim.search.rslts), 
        function(ref.db.name) {
            srch.tbl <- seq.sim.search.rslts[[ref.db.name]]
            query.ind <- which(srch.tbl$qseqid == prot.id)
            if (length(query.ind) > 0) {
                hits <- srch.tbl[query.ind, ]
                hits.uniq.sseqids <- unique(hits$sseqid)
                phrases.lst <- list()
                for (hit.sseqid in hits.uniq.sseqids) {
                  #' Each unique Hit should have a single distinct `HRD` (description)
                  hit.hrd <- unique(hits[which(hits$sseqid == 
                    hit.sseqid), "HRD"])[[1]]
                  hit.phrases <- phrases(hit.hrd)
                  #' Only add non NULL phrases:
                  if (!is.null(hit.phrases)) {
                    phrases.lst[[hit.sseqid]] <- hit.phrases
                  }
                }
                if (length(phrases.lst) > 0) {
                  phrases.lst
                } else NULL
            } else NULL
        }), names(seq.sim.search.rslts))
}

#' Computes the performance of all phrases generated for a query. Performance
#' is assessed in terms of F-Score, precision, and recall. Furthermore
#' identifies which phrase would be selected by selected prot.scriber scoring
#' schemes. Please note that all phrases are lower-cased in order to ease
#' identification of unique phrases and word indexing.
#'
#' @param prot.id - A string representing the query protein identifier (qseqid)
#' @param phrases.for.query - The result of calling function phrasesForQuery
#' @param hrd.references - The 'truth' i.e. a list of reference word sets for
#' query proteins (argument 'prot.id' should be among the names of this
#' argument list).
#' @param score.funks - A named list of two member lists (see default value
#' below for more details). Each first level member must contain two entries
#' 'word.score.funk' and 'phrase.score.funk'. 'word.score.funk' is a function
#' expected to receive at least a single first argument, an instance of
#' base::table holding the word frequencies. This function must return a named
#' numeric vector, whose names are the words and values are their respective
#' computed scores. This vector is used as basis to calculate prot.scriber
#' phrase-scores as the sum of the scores of the words a phrase contains.  The
#' second entry 'phrase.score.funk' is a function that receives at least two
#' arguments, the first being a character vector representing a phrase and the
#' second a named numeric vector result of invoking the respective
#' 'word.score.funk', i.e. this vector's names are words and values their
#' respective scores. The default of this argument is rather long:
#' getOption('statsOfPhrasesForQuery.score.funks',
#'  list(
#'    centered.inverse.inf.cntnt.mean = list(word.score.funk = function(x) {
#'          centeredWordScores(x, level.funk = mean)
#'      }, phrase.score.funk = sumWordScores),
#'    centered.inverse.inf.cntnt.median = list(word.score.funk = function(x) {
#'          centeredWordScores(x, level.funk = median)
#'      }, phrase.score.funk = sumWordScores),
#'    centered.inverse.inf.cntnt.quarterQuantile = list(word.score.funk = function(x) {
#'          centeredWordScores(x, level.funk = function(y) {
#'              quantile(y, probs = 0.25)
#'          })
#'      }, phrase.score.funk = sumWordScores),
#'    polynomial.word.scores = list(word.score.funk = polynomicWordScores, 
#'          phrase.score.funk = sumWordScores),
#'    centered.frequencies = list(word.score.funk = centeredLinearWordScores, 
#'          phrase.score.funk = sumWordScores)
#'  ))
#'
#' @return An instance of base::data.frame with the default columns Protein.ID,
#' precision, recall, F.Score, beta, Phrase, n.words, and one column for each
#' entry in argument 'word.score.funks' holding the respective phrase-score.
#' @export
statsOfPhrasesForQuery <- function(prot.id, phrases.for.query, 
    hrd.references, score.funks = getOption("statsOfPhrasesForQuery.score.funks", 
        list(centered.inverse.inf.cntnt.mean = list(word.score.funk = function(x) {
            centeredWordScores(x, level.funk = mean)
        }, phrase.score.funk = sumWordScores), centered.inverse.inf.cntnt.median = list(word.score.funk = function(x) {
            centeredWordScores(x, level.funk = median)
        }, phrase.score.funk = sumWordScores), centered.inverse.inf.cntnt.quarterQuantile = list(word.score.funk = function(x) {
            centeredWordScores(x, level.funk = function(y) {
                quantile(y, probs = 0.25)
            })
        }, phrase.score.funk = sumWordScores), polynomial.word.scores = list(word.score.funk = polynomicWordScores, 
            phrase.score.funk = sumWordScores), centered.frequencies = list(word.score.funk = centeredLinearWordScores, 
            phrase.score.funk = sumWordScores)))) {
    word.freqs <- wordFrequenciesFromPhrasesList(phrases.for.query)
    all.phrases <- list()
    #' Any informative words found?
    if (length(word.freqs) > 0) {
        for (ref.db.name in names(phrases.for.query)) {
            for (sseqid in names(phrases.for.query[[ref.db.name]])) {
                hit.phrases <- phrases.for.query[[ref.db.name]][[sseqid]]$phrases
                for (phrase in hit.phrases) {
                  all.phrases[[length(all.phrases) + 1]] <- phrase
                }
            }
        }
    }
    #' WARNING: Currently lowercasing all phrases in order to easily find
    #' unique phrases:
    uniq.phrases <- unique(lapply(all.phrases, tolower))
    rm(all.phrases)  #' Clean up
    if (length(uniq.phrases) > 0) {
        #' Prepare different prot.scriber phrase scoring as requested in
        #' argument 'score.funks':
        word.scores <- setNames(lapply(score.funks, function(score.funk.lst) {
            score.funk.lst[["word.score.funk"]](word.freqs)
        }), names(score.funks))
        do.call(rbind, lapply(uniq.phrases, function(phrase) {
            f.score.df <- fScore(phrase, hrd.references[[tolower(prot.id)]], 
                prot.id)
            #' Add prot.scriber phrase-scores as requested in argument
            #' 'word.score.funks':
            for (score.name in names(score.funks)) {
                wrd.scrs <- word.scores[[score.name]]
                f.score.df[, score.name] <- score.funks[[score.name]][["phrase.score.funk"]](phrase, 
                  wrd.scrs)
            }
            f.score.df
        }))
    } else {
        #' No phrases to process by Prot-Scriber:
        f.score.df <- fScore(character(0), hrd.references[[tolower(prot.id)]], 
            prot.id)
        f.score.df$HRD <- NA
        f.score.df$n.words <- NA
        #' Add prot.scriber phrase-scores as requested in argument
        #' 'word.score.funks':
        for (score.name in names(score.funks)) {
            f.score.df[, score.name] <- NA
        }
        f.score.df
    }
}

#' For multi domain query sequences this function creates a single prot-scriber
#' annotation for each disjoint alignment region and re-evaluates the
#' performance scores for the now concatonated human readable description.
#'
#' @param alignmnt.regions.phrases.stats - The result of invoking the following
#' code:
#' lapply(qseqid.sssr.per.alignment.regions, 
#'   function(sssr.i) {
#'     q.phrases <- phrasesForQuery(qseqid, sssr.i)
#'     q.phrases.stats <- statsOfPhrasesForQuery(qseqid, 
#'       q.phrases, ref.word.sets, score.funks = prot.scriber.score.funks)
#'     setNames(lapply(names(prot.scriber.score.funks), 
#'       function(ps.score.method) {
#'         bestProtScriberPhrases(qseqid, q.phrases.stats, 
#'           ps.score.method, ref.words, univ.words)
#'       }), names(prot.scriber.score.funks))
#'   }
#' )
#' @param hrd.ref - A vector holding the truth (references)
#' @param univ.words - A character vector holding all non blacklisted words
#' found in all considered candidate descriptions, the universe of words.
#' 
#' @return An instance of base::data.frame with the following columns:
#' Protein.ID, HRD, precision, recall, F.Score, beta, Method, Method.Score,
#' univ.words (number of words in the argument 'univ.words'), ref.words (number
#' of words in the argument 'ref'), univ.ref.words (number of words in the
#' intersection of arguments 'ref' and 'univ.words') and MCC
#' @export
joinMultiRegionStatsOfPhrasesForQuery <- function(alignmnt.regions.phrases.stats, 
    hrd.ref, univ.words) {
    if (length(alignmnt.regions.phrases.stats) == 1) {
        alignmnt.regions.phrases.stats[[1]]
    } else if (length(alignmnt.regions.phrases.stats) > 1) {
        prot.scriber.score.funks <- names(alignmnt.regions.phrases.stats[[1]])
        do.call(rbind, lapply(prot.scriber.score.funks, function(ps.sf) {
            hrd.lst <- list()
            method.score <- 0
            for (ar.ps in alignmnt.regions.phrases.stats) {
                ar.ps.i <- ar.ps[[ps.sf]]
                hrd.lst[[length(hrd.lst) + 1]] <- ar.ps.i$HRD
                #' print(ar.ps.i)
                method.score <- method.score + ar.ps.i$Method.Score
            }
            hrd.concat.words <- wordSet(paste(hrd.lst, collapse = " "), 
                blacklist.regexs = NULL)
            #' All tables holding the results for the current ps.sf have
            #' constant columns, so use the first for the constant data:
            ar.ps.tbl <- alignmnt.regions.phrases.stats[[1]][[ps.sf]]
            prot.id <- ar.ps.tbl$Protein.ID
            ps.best.f.score.df <- fScore(hrd.concat.words, hrd.ref, 
                prot.id)
            ps.best.mcc.df <- matthewsCorrelationCoefficient(hrd.concat.words, 
                hrd.ref, univ.words, prot.id)[, c("Protein.ID", 
                "MCC", "univ.words", "ref.words", "univ.ref.words")]
            ps.best.eval.df <- merge(ps.best.f.score.df, ps.best.mcc.df, 
                by = "Protein.ID")
            ps.best.eval.df$Method <- ar.ps.tbl$Method
            ps.best.eval.df$Method.Score <- method.score
            ps.best.eval.df
        }))
    }
}
