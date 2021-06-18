#' Computes Matthew's correlation coefficient to assess the performance of the
#' predicted Human Readable Description (HRD).
#'
#' @param pred - A vector of predictions
#' @param ref - A vector holding the truth (references)
#' @param univ.words - A character vector holding all non blacklisted words
#' @param prot.id - A string identifying the query protein for which the
#' prediction has been made.
#' @param pred.to.string.funk - A function receiving as input the argument
#' 'pred', i.e. a character vector and set of words, returning a single scalar
#' string the human readable description generated from the ordered set of
#' words ('pred'). Default is getOption('fScore.pred.to.string.funk',
#' function(wrds) {paste(wrds, collapse=' ')})
#'
#' @return An instance of base::data.frame with the following columns:
#' Protein.ID, HRD, n.words (number of words in the argument 'pred'),
#' univ.words (number of words in the argument 'univ.words'), ref.words (number
#' of words in the argument 'ref'), univ.ref.words (number of words in the
#' intersection of arguments 'ref' and 'univ.words') and MCC.
#' @export
matthewsCorrelationCoefficient <- function(pred, ref, univ.words, 
    prot.id, pred.to.string.funk = getOption("fScore.pred.to.string.funk", 
        function(wrds) {
            paste(wrds, collapse = " ")
        })) {
    hrd <- if (length(pred) > 0) {
        pred.to.string.funk(pred)
    } else {
        as.character(NA)
    }
    if (length(ref) == 0 || is.na(ref)) {
        #' No references means we cannot compute a MCC:
        data.frame(Protein.ID = prot.id, HRD = hrd, n.words = length(pred), 
            univ.words = length(univ.words), ref.words = length(ref), 
            univ.ref.words = length(intersect(univ.words, ref)), 
            MCC = NA, stringsAsFactors = FALSE)
    } else {
        true.pos <- length(intersect(pred, ref))
        true.neg <- length(setdiff(univ.words, union(pred, ref)))
        false.pos <- length(setdiff(pred, ref))
        false.neg <- length(setdiff(ref, pred))
        mcc.denominator <- sqrt((true.pos + false.pos) * (true.pos + 
            false.neg) * (true.neg + false.pos) * (true.neg + 
            false.neg))
        if (mcc.denominator == 0) {
            mcc.denominator <- 1
        }
        m.c.c <- if (identical(pred, ref)) {
            1
        } else {
            tp.tn <- if (length(setdiff(univ.words, ref)) == 
                0) {
                #' Do not punish not having any false negative candidate words in
                #' the current setting:
                true.pos
            } else {
                true.pos * true.neg
            }
            (tp.tn - false.pos * false.neg)/mcc.denominator
        }
        data.frame(Protein.ID = prot.id, HRD = hrd, n.words = length(pred), 
            univ.words = length(univ.words), ref.words = length(ref), 
            univ.ref.words = length(intersect(univ.words, ref)), 
            MCC = m.c.c, stringsAsFactors = FALSE)
    }
}

#' Function tests its namesake `matthewsCorrelationCoefficient`.
#'
#' @return TRUE if and only if all tests pass.
#' @export
testMatthewsCorrelationCoefficient <- function() {
    p1 <- c("Hello", "World")
    r1 <- p1
    u1 <- p1
    t1 <- identical(matthewsCorrelationCoefficient(p1, r1, u1, 
        "P1")$MCC, 1)
    p2 <- c("Hello", "World")
    r2 <- c("Not", "a", "single", "word")
    u2 <- union(p2, r2)
    t2 <- identical(matthewsCorrelationCoefficient(p2, r2, u2, 
        "P2")$MCC, -1)
    p3 <- c("Hello", "World", "Hola", "Mundo")
    r3 <- c("Some", "words", "Hello", "World")
    u3 <- union(p3, r3)
    t3 <- identical(matthewsCorrelationCoefficient(p3, r3, u3, 
        "P3")$MCC, -0.5)
    t4.res <- matthewsCorrelationCoefficient(p3, NULL, NULL, 
        "P3")
    t4 <- is.na(t4.res$MCC)
    all(t1, t2, t3, t4)
}

#' Computes the F-Score of a prediction given the truth (reference). For
#' details see Wikipedia, here:
#' https://en.wikipedia.org/wiki/F-score
#' 
#' @param pred - A vector of predictions
#' @param ref - A vector holding the truth (references)
#' @param prot.id - A string identifying the query protein for which the
#' prediction has been made.
#' @param beta - recall is considered beta times as important as precision.
#' Default is getOption('fScore.beta', 1) thus emphazising precision
#' and recall equally. 
#' @param pred.to.string.funk - A function receiving as input the argument
#' 'pred', i.e. a character vector and set of words, returning a single scalar
#' string the human readable description generated from the ordered set of
#' words ('pred'). Default is getOption('fScore.pred.to.string.funk',
#' function(wrds) {paste(wrds, collapse=' ')})
#'
#' @return An instance of `base::data.frame` with exactly one row and the
#' following columns: Protein.ID, HRD, n.words, precision, recall, F.Score,
#' beta
#' @export
fScore <- function(pred, ref, prot.id, beta = getOption("fScore.beta", 
    1), pred.to.string.funk = getOption("fScore.pred.to.string.funk", 
    function(wrds) {
        paste(wrds, collapse = " ")
    })) {
    hrd <- if (length(pred) > 0) {
        pred.to.string.funk(pred)
    } else {
        as.character(NA)
    }
    if (length(ref) == 0 || is.na(ref)) {
        #' No references means we cannot compute a F-Score:
        data.frame(Protein.ID = prot.id, HRD = hrd, n.words = length(pred), 
            precision = NA, recall = NA, F.Score = NA, fScore.beta = beta, 
            stringsAsFactors = FALSE)
    } else {
        true.pos <- intersect(pred, ref)
        precision <- if (length(pred) > 0) {
            length(true.pos)/length(pred)
        } else {
            #' Zero predicted words match any of the reference words,
            #' means precision is zero:
            0
        }
        recall <- length(true.pos)/length(ref)
        f.score <- if ((beta * precision + recall) == 0) {
            0
        } else {
            (1 + beta^2) * (precision * recall)/(beta * precision + 
                recall)
        }
        data.frame(Protein.ID = prot.id, HRD = hrd, n.words = length(pred), 
            precision = precision, recall = recall, F.Score = f.score, 
            fScore.beta = beta, stringsAsFactors = FALSE)
    }
}

#' Function tests its namesake `fScore`.
#'
#' @return TRUE if and only if all tests pass.
#' @export
testFScore <- function() {
    beta <- 1
    p1 <- c("Hello", "World")
    r1 <- p1
    t1 <- identical(fScore(p1, r1, "P1", beta = beta)$F.Score, 
        1)
    p2 <- c("Hello", "World")
    r2 <- c("Not", "a", "single", "word")
    t2 <- identical(fScore(p2, r2, "P2", beta = beta)$F.Score, 
        0)
    p3 <- c("Hello", "World", "Hola", "Mundo")
    r3 <- c("Some", "words", "Hello", "World")
    t3 <- identical(fScore(p3, r3, "P3", beta = beta)$F.Score, 
        0.5)
    t4.res <- fScore(p3, NULL, "P3", beta = beta)
    t4 <- is.na(t4.res$F.Score) && is.na(t4.res$precision) && 
        is.na(t4.res$recall)
    all(t1, t2, t3, t4)
}

#' Finds the best sequence similarity search result Hits (Best Blast Hit) for
#' the argument query sequence ('prot.id') in each search result table stored
#' in the argument 'sssr' list of seq-sim-search-result tables. The best Hit is
#' considered that of highest bitscore and, if tied, the longest human readable
#' description (HRD).
#'
#' @param prot.id - A string representing the query protein identifier (qseqid)
#' @param sssr - A named list, in which names represent searched reference
#' sequence databases and values the read in tabular output (see function
#' parseSeqSimSearchTable).
#' @param hrd.ref - A vector holding the truth (references)
#' @param univ.words - A character vector holding all non blacklisted words
#' found in all considered candidate descriptions, the universe of words.
#' @param lowercase.hrd - A boolean flag indicating whether to return the found
#' HRD in lowercase or not. Default is getOption('bestBlastHrds.lowercase.hrd',
#' TRUE)
#'
#' @return An instance of base::data.frame with the following columns:
#' Protein.ID, HRD, precision, recall, F.Score, beta, Method, Method.Score
#' @export
bestBlastHrds <- function(prot.id, sssr, hrd.ref, univ.words, 
    lowercase.hrd = getOption("bestBlastHrds.lowercase.hrd", 
        TRUE)) {
    bb.hrd.lst <- list()
    for (ref.db.name in names(sssr)) {
        seq.sim.search.tbl <- sssr[[ref.db.name]]
        if (prot.id %in% seq.sim.search.tbl$qseqid) {
            q.hits <- seq.sim.search.tbl[which(seq.sim.search.tbl$qseqid == 
                prot.id), ]
            q.hits$HRD.nchar <- as.integer(lapply(q.hits$HRD, 
                nchar))
            q.hits.best <- q.hits[order(q.hits$bitscore, q.hits$HRD.nchar, 
                decreasing = TRUE), ][1, ]
            bb.hrd <- q.hits.best$HRD
            if (lowercase.hrd) {
                bb.hrd <- tolower(bb.hrd)
            }
            bb.hrd.word.set <- wordSet(bb.hrd, blacklist.regexs = NULL)
            method.score <- q.hits.best$bitscore
        } else {
            bb.hrd.word.set <- character(0)
            method.score <- NA
        }
        bb.hrd.f.score.df <- fScore(bb.hrd.word.set, hrd.ref, 
            prot.id)
        bb.hrd.mcc.df <- matthewsCorrelationCoefficient(bb.hrd.word.set, 
            hrd.ref, univ.words, prot.id)[, c("Protein.ID", "MCC", 
            "univ.words", "ref.words", "univ.ref.words")]
        bb.hrd.eval.df <- merge(bb.hrd.f.score.df, bb.hrd.mcc.df, 
            by = "Protein.ID")
        bb.hrd.eval.df$Method <- paste0("BB.", ref.db.name)
        bb.hrd.eval.df$Method.Score <- method.score
        bb.hrd.lst[[length(bb.hrd.lst) + 1]] <- bb.hrd.eval.df
    }
    do.call(rbind, bb.hrd.lst)
}

#' Extracts the highest scoring prot.scriber phrase (candidate human readable
#' description) from the argument 'phrases.stats' for the argument prot.scriber
#' scoring-method 'prot.scriber.score.col'.
#'
#' @param prot.id - A string representing the query protein identifier (qseqid)
#' @param phrases.stats - An instance of base::data.frame result of calling
#' function 'statsOfPhrasesForQuery'.
#' @param prot.scriber.score.col - A string or integer identifying which column
#' of argument 'phrases.stats' to use to find the highest scoring phrase.
#' Refers to the different tested prot.scriber scoring methods.
#' @param hrd.ref - A vector holding the truth (references)
#' @param univ.words - A character vector holding all non blacklisted words
#' found in all considered candidate descriptions, the universe of words.
#'
#' @return An instance of base::data.frame with the following columns:
#' Protein.ID, HRD, precision, recall, F.Score, beta, Method, Method.Score
#' @export
bestProtScriberPhrases <- function(prot.id, phrases.stats, prot.scriber.score.col, 
    hrd.ref, univ.words) {
    if (nrow(phrases.stats) > 0) {
        ps.best <- phrases.stats[order(phrases.stats[, prot.scriber.score.col], 
            phrases.stats$n.words, decreasing = TRUE), ][1, ]
        #' No need to filter out blacklisted words again, set blacklist.regexs
        #' to 'NULL':
        ps.best.word.set <- wordSet(ps.best$HRD, blacklist.regexs = NULL)
        ps.best.f.score.df <- fScore(ps.best.word.set, hrd.ref, 
            prot.id)
        ps.best.mcc.df <- matthewsCorrelationCoefficient(ps.best.word.set, 
            hrd.ref, univ.words, prot.id)[, c("Protein.ID", "MCC", 
            "univ.words", "ref.words", "univ.ref.words")]
        ps.best.eval.df <- merge(ps.best.f.score.df, ps.best.mcc.df, 
            by = "Protein.ID")
        ps.best.eval.df$Method <- prot.scriber.score.col
        ps.best.eval.df$Method.Score <- ps.best[, prot.scriber.score.col]
        ps.best.eval.df
    }
}

#' Splits all human readable descriptions (HRD) of sequences similarity search
#' result Hits for the argument query 'prot.id' into words and returns them as
#' a set (unique). These words can be used to estimate e.g. the true negative
#' rate.
#'
#' @param prot.id - A string representing the query protein identifier (qseqid)
#' @param sssr - A named list, in which names represent searched reference
#' sequence databases and values the read in tabular output (see function
#' parseSeqSimSearchTable).
#' @param to.lower - A boolean flag indicating whether to return the found HRD
#' in lowercase or not. Default is getOption('wordUniverse.to.lower', TRUE)
#'
#' @return A character vector of unique words found in the HRDs of the argument
#' query sequence's ('prot.id') sequence similarity search Hits.
#' @export
wordUniverse <- function(prot.id, sssr, to.lower = getOption("wordUniverse.to.lower", 
    TRUE)) {
    univ.words <- unique(unlist(lapply(sssr, function(seq.sim.tbl) {
        hits.ind <- which(seq.sim.tbl$qseqid == prot.id)
        if (length(hits.ind) > 0) {
            unlist(lapply(seq.sim.tbl[hits.ind, ]$HRD, wordSet))
        } else {
            character(0)
        }
    })))
    if (to.lower) {
        tolower(univ.words)
    } else univ.words
}

#' Generate short human readable descriptions (HRD) for all query proteins in
#' the argument sequence similarity search results 'sssr' using the Best Blast
#' and registered argument 'prot.scriber.score.funks'. Evaluate the performance
#' (F-Score) of the HRDs if reference word sets are available for the
#' respective query proteins (qseqid).
#'
#' @param sssr - A named list, in which names represent searched reference
#' sequence databases and values the read in tabular output (see function
#' parseSeqSimSearchTable).
#' @param ref.word.sets - A named list in which the names are query protein
#' identifiers ('qseqid') in lowercase and the values are character vectors
#' holding reference words ('the truth').
#' @param alignmnt.regions - An instance of base::data.frame result of calling
#' function 'allQueriesAlignmentRegions'. The data.frame has three columns
#' 'Protein.ID', 'start.pos', and 'end.pos' indicating to which sequence region
#' in the query the sequence similarity searches generated local alignments
#' for. If non intersecting regions for the same query are found the
#' prot-scriber annotation is carried out independently for these and results
#' are concatonated, in the order of the regions.
#' @param prot.scriber.score.funks - A named list of two member lists (see
#' default value below for more details). Each first level member must contain
#' two entries 'word.score.funk' and 'phrase.score.funk'. 'word.score.funk' is
#' a function expected to receive at least a single first argument, an instance
#' of base::table holding the word frequencies. This function must return a
#' named numeric vector, whose names are the words and values are their
#' respective computed scores. This vector is used as basis to calculate
#' prot.scriber phrase-scores as the sum of the scores of the words a phrase
#' contains.  The second entry 'phrase.score.funk' is a function that receives
#' at least two arguments, the first being a character vector representing a
#' phrase and the second a named numeric vector result of invoking the
#' respective 'word.score.funk', i.e. this vector's names are words and values
#' their respective scores. The default of this argument is rather long:
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
#' @return An instance of base::data.table with the following columns:
#' Protein.ID, HRD, precision, recall, F.Score, beta, Method, Method.Score
#' @export
annotateProteinsAndEvaluatePerformance <- function(sssr, ref.word.sets, 
    alignmnt.regions, prot.scriber.score.funks = getOption("statsOfPhrasesForQuery.score.funks", 
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
    prot.ids <- unique(unlist(lapply(sssr, function(seq.sim.search.tbl) {
        unique(seq.sim.search.tbl$qseqid)
    })))
    message("Starting annotating human readable descriptions (HRD) for ", 
        length(prot.ids), " query sequences.")
    do.call(rbind, lapply(prot.ids, function(qseqid) {
        tryCatch({
            univ.words <- wordUniverse(qseqid, sssr)
            #' Best Blast:
            ref.words <- ref.word.sets[[tolower(qseqid)]]
            bb.hrds <- bestBlastHrds(qseqid, sssr, ref.words, 
                univ.words)
            #' Prot Scriber:
            qseqid.sssr.per.alignment.regions <- sssrForRegions(qseqid, 
                alignmnt.regions, sssr)
            alignmnt.regions.phrases.stats <- lapply(qseqid.sssr.per.alignment.regions, 
                function(sssr.i) {
                  q.phrases <- phrasesForQuery(qseqid, sssr.i)
                  q.phrases.stats <- statsOfPhrasesForQuery(qseqid, 
                    q.phrases, ref.word.sets, score.funks = prot.scriber.score.funks)
                  setNames(lapply(names(prot.scriber.score.funks), 
                    function(ps.score.method) {
                      bestProtScriberPhrases(qseqid, q.phrases.stats, 
                        ps.score.method, ref.words, univ.words)
                    }), names(prot.scriber.score.funks))
                })
            prot.scribr.hrds <- joinMultiRegionStatsOfPhrasesForQuery(alignmnt.regions.phrases.stats, 
                ref.words, univ.words)
            rbind(bb.hrds, prot.scribr.hrds)
        }, error = function(e) {
            #' browser()
            warning("Query '", qseqid, "' caused an error:\n", 
                e)
        })
    }))
}
