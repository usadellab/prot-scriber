#' Computes the Matthew's correlation coefficient.
#'
#' @param pred - A vector of predictions
#' @param ref - A vector holding the truth (references)
#' @param univ.words - A character vector holding all non blacklisted words
#'
#' @return A scalar numeric value between -1 and +1, the MCC.
#' @export
mcc <- function(pred, ref, univ.words) {
  true.pos <- length(intersect(pred, ref))
  true.neg <- length(setdiff(univ.words, union(
    pred,
    ref
  )))
  false.pos <- length(setdiff(pred, ref))
  false.neg <- length(setdiff(ref, pred))
  mcc.denominator <- sqrt((true.pos + false.pos) *
    (true.pos + false.neg) * (true.neg + false.pos) *
    (true.neg + false.neg))
  if ((true.pos + false.pos) == 0 || (true.pos +
    false.neg) == 0 || (true.neg + false.pos) ==
    0 || (true.neg + false.neg) == 0) {
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
    (tp.tn - false.pos * false.neg) / mcc.denominator
  }
  if (m.c.c > 1) {
    #' Did as well as possible:
    1
  } else {
    m.c.c
  }
}

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
#' intersection of arguments 'ref' and 'univ.words'), MCC (the Matthew's
#' correlation coefficient), MCC.best.possible (the theoretically highest
#' achievable MCC given the argument 'univ.words'), MCC.relative (the relative
#' fraction of the MCC.best.possible achieved by the argument 'pred').
#' @export
matthewsCorrelationCoefficient <- function(
    pred, ref,
    univ.words, prot.id, pred.to.string.funk = getOption(
      "fScore.pred.to.string.funk",
      function(wrds) {
        paste(wrds, collapse = " ")
      }
    )) {
  hrd <- if (length(pred) > 0 && all(!is.na(pred))) {
    pred.to.string.funk(pred)
  } else {
    NA
  }
  if (length(ref) == 0 || all(is.na(ref))) {
    #' No references means we cannot compute a MCC:
    data.frame(
      Protein.ID = prot.id, HRD = hrd,
      n.words = length(pred), univ.words = length(univ.words),
      ref.words = length(ref), univ.ref.words = length(intersect(
        univ.words,
        ref
      )), MCC = NA, MCC.relative = NA,
      MCC.best.possible = NA, stringsAsFactors = FALSE
    )
  } else {
    m.c.c <- mcc(pred, ref, univ.words)
    #' Best MCC achievable with the seq-sim-search results words:
    m.c.c.best <- mcc(
      intersect(ref, univ.words),
      ref, univ.words
    )
    if (m.c.c > m.c.c.best) {
      warning(
        "MCC.best.possible (", m.c.c.best,
        ") is smaller than MCC (", m.c.c, ")! This situation should NEVER arise."
      )
    }
    #' Fraction of the best MCC achievable actually reached:
    m.c.c.rel <- 1 - (m.c.c.best - m.c.c) / (m.c.c.best -
      -1)
    data.frame(
      Protein.ID = prot.id, HRD = hrd,
      n.words = length(pred), univ.words = length(univ.words),
      ref.words = length(ref), univ.ref.words = length(intersect(
        univ.words,
        ref
      )), MCC = m.c.c, MCC.relative = m.c.c.rel,
      MCC.best.possible = m.c.c.best, stringsAsFactors = FALSE
    )
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
  t1 <- identical(matthewsCorrelationCoefficient(
    p1,
    r1, u1, "P1"
  )$MCC, 1)
  p2 <- c("Hello", "World")
  r2 <- c("Not", "a", "single", "word")
  u2 <- union(p2, r2)
  t2 <- identical(matthewsCorrelationCoefficient(
    p2,
    r2, u2, "P2"
  )$MCC, -1)
  p3 <- c("Hello", "World", "Hola", "Mundo")
  r3 <- c("Some", "words", "Hello", "World")
  u3 <- union(p3, r3)
  t3 <- identical(matthewsCorrelationCoefficient(
    p3,
    r3, u3, "P3"
  )$MCC, -0.5)
  t4.res <- matthewsCorrelationCoefficient(
    p3, NULL,
    NULL, "P3"
  )
  t4 <- is.na(t4.res$MCC)
  p5 <- NA
  r5 <- c("kinase", "domain")
  t5.res <- matthewsCorrelationCoefficient(
    p5, r5,
    r5, "P5"
  )
  t5 <- identical(t5.res$MCC, -1)
  p6 <- character(0)
  r6 <- c("kinase", "domain")
  t6.res <- matthewsCorrelationCoefficient(
    p6, r6,
    r6, "P6"
  )
  t6 <- identical(t6.res$MCC, 0)
  #' Test what happens if we do not have any false negatives:
  p7 <- LETTERS[1:3]
  r7 <- LETTERS[1:3]
  t7.res <- matthewsCorrelationCoefficient(
    p7, r7,
    union(r7, LETTERS[4:5]), "P7"
  )
  t7 <- identical(t7.res$MCC, 1)
  all(t1, t2, t3, t4, t5, t6, t7)
}

#' Computes the F-Score from argument prediction and reference.
#'
#' @param pred - A vector of predictions
#' @param ref - A vector holding the truth (references)
#' @param beta - recall is considered beta times as important as precision.
#' Default is getOption('fScore.beta', 1) thus emphazising precision
#' and recall equally.
#'
#' @return An instance of base::data.frame with columns F.Score, precision,
#' recall, and fScore.beta.
#' @export
fScoreCalculator <- function(pred, ref,
                             beta = getOption(
                               "fScore.beta",
                               1
                             )) {
  true.pos <- intersect(pred, ref)
  # presision is also called specificity:
  precision <- if (length(pred) > 0) {
    length(true.pos) / length(pred)
  } else {
    #' Zero predicted words match any of the reference words,
    #' means precision is zero:
    0
  }
  # recall is also called sensitivity:
  recall <- length(true.pos) / length(ref)
  #' Calculate the F.Score with constant beta:
  f.c <- if ((beta * precision + recall) == 0) {
    0
  } else {
    (1 + beta^2) * (precision * recall) / (beta *
      precision + recall)
  }
  data.frame(
    F.Score = f.c, precision = precision,
    recall = recall, fScore.beta = beta, stringsAsFactors = FALSE
  )
}

#' Computes the F-Score of a prediction given the truth (reference). For
#' details see Wikipedia, here:
#' https://en.wikipedia.org/wiki/F-score
#'
#' @param pred - A vector of predictions
#' @param ref - A vector holding the truth (references)
#' @param univ.words - A character vector holding all non blacklisted words
#' found in all considered candidate descriptions, the universe of words.
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
#' F.Score.best.possible (the theoretically highest achievable F.Score given
#' the argument 'univ.words'), F.Score.relative (the relative fraction of the
#' highest possibly achievable F.Score given the argument 'pred'), and beta
#' @export
fScore <- function(
    pred, ref, univ.words, prot.id,
    beta = getOption("fScore.beta", 1), pred.to.string.funk = getOption(
      "fScore.pred.to.string.funk",
      function(wrds) {
        paste(wrds, collapse = " ")
      }
    )) {
  hrd <- if (length(pred) > 0) {
    pred.to.string.funk(pred)
  } else {
    as.character(NA)
  }
  if (length(ref) == 0 || all(is.na(ref))) {
    #' No references means we cannot compute a F-Score:
    data.frame(
      Protein.ID = prot.id, HRD = hrd,
      n.words = length(pred), precision = NA,
      recall = NA, F.Score = NA, F.Score.relative = NA,
      F.Score.best.possible = NA, fScore.beta = beta,
      stringsAsFactors = FALSE
    )
  } else {
    f.score.df <- fScoreCalculator(pred, ref, beta = beta)
    #' Best possibly achievable F.Score given the seq-sim-search results'
    #' words:
    f.score.best <- fScoreCalculator(intersect(
      ref,
      univ.words
    ), ref, beta = beta)$F.Score
    #' How much of the best possible was achieved?
    f.score.rel <- if (f.score.best == 0) {
      #' Any prediction made, cannot be correct, because the
      #' intersection between all words used in predictions and the
      #' reference is empty:
      0
    } else {
      f.score.df$F.Score / f.score.best
    }
    #' Return results:
    f.score.df$Protein.ID <- prot.id
    f.score.df$HRD <- hrd
    f.score.df$n.words <- length(pred)
    f.score.df$F.Score.relative <- f.score.rel
    f.score.df$F.Score.best.possible <- f.score.best
    f.score.df
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
  u1 <- p1
  t1 <- identical(
    fScore(p1, r1, u1, "P1", beta = beta)$F.Score,
    1
  )
  p2 <- c("Hello", "World")
  r2 <- c("Not", "a", "single", "word")
  u2 <- p2
  t2 <- identical(
    fScore(p2, r2, u2, "P2", beta = beta)$F.Score,
    0
  )
  p3 <- c("Hello", "World", "Hola", "Mundo")
  r3 <- c("Some", "words", "Hello", "World")
  u3 <- p3
  t3 <- identical(
    fScore(p3, r3, u3, "P3", beta = beta)$F.Score,
    0.5
  )
  t4.res <- fScore(p3, NULL, u3, "P3", beta = beta)
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
bestBlastHrds <- function(
    prot.id, sssr, hrd.ref, univ.words,
    lowercase.hrd = getOption(
      "bestBlastHrds.lowercase.hrd",
      TRUE
    )) {
  bb.hrd.lst <- list()
  for (ref.db.name in names(sssr)) {
    seq.sim.search.tbl <- sssr[[ref.db.name]]
    if (prot.id %in% seq.sim.search.tbl$qseqid) {
      q.hits <- seq.sim.search.tbl[which(seq.sim.search.tbl$qseqid ==
        prot.id), ]
      q.hits$HRD.nchar <- as.integer(lapply(
        q.hits$HRD,
        nchar
      ))
      q.hits.order <- if ("bitscore" %in% colnames(q.hits)) {
        order(q.hits$bitscore, q.hits$HRD.nchar,
          decreasing = TRUE
        )
      } else {
        order(q.hits$HRD.nchar, decreasing = TRUE)
      }
      q.hits.best <- q.hits[q.hits.order, ][1, ]
      bb.hrd <- q.hits.best$HRD
      if (lowercase.hrd) {
        bb.hrd <- tolower(bb.hrd)
      }
      bb.hrd.word.set <- wordSet(bb.hrd, blacklist.regexs = NULL)
      method.score <- if ("bitscore" %in% colnames(q.hits.best)) {
        q.hits.best$bitscore
      } else {
        NA
      }
    } else {
      bb.hrd.word.set <- character(0)
      method.score <- NA
    }
    bb.hrd.f.score.df <- fScore(
      bb.hrd.word.set,
      hrd.ref, univ.words, prot.id
    )
    bb.hrd.mcc.df <- matthewsCorrelationCoefficient(
      bb.hrd.word.set,
      hrd.ref, univ.words, prot.id
    )[, c(
      "Protein.ID",
      "MCC", "MCC.relative", "MCC.best.possible",
      "univ.words", "ref.words", "univ.ref.words"
    )]
    bb.hrd.eval.df <- merge(bb.hrd.f.score.df,
      bb.hrd.mcc.df,
      by = "Protein.ID"
    )
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
#' Protein.ID, HRD, precision, recall, F.Score, beta, Method, Method.Score,
#' univ.words (number of words in the argument 'univ.words'), ref.words (number
#' of words in the argument 'ref'), univ.ref.words (number of words in the
#' intersection of arguments 'ref' and 'univ.words') and MCC
#' @export
bestProtScriberPhrases <- function(
    prot.id, phrases.stats,
    prot.scriber.score.col, hrd.ref, univ.words) {
  if (nrow(phrases.stats) > 0) {
    ps.best <- phrases.stats[order(
      phrases.stats[
        ,
        prot.scriber.score.col
      ], phrases.stats$n.words,
      decreasing = TRUE
    ), ][1, ]
    #' No need to filter out blacklisted words again, set blacklist.regexs
    #' to 'NULL':
    ps.best.word.set <- wordSet(ps.best$HRD, blacklist.regexs = NULL)
    ps.best.f.score.df <- fScore(
      ps.best.word.set,
      hrd.ref, univ.words, prot.id
    )
    ps.best.mcc.df <- matthewsCorrelationCoefficient(
      ps.best.word.set,
      hrd.ref, univ.words, prot.id
    )[, c(
      "Protein.ID",
      "MCC", "MCC.relative", "MCC.best.possible",
      "univ.words", "ref.words", "univ.ref.words"
    )]
    ps.best.eval.df <- merge(ps.best.f.score.df,
      ps.best.mcc.df,
      by = "Protein.ID"
    )
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
wordUniverse <- function(prot.id, sssr, to.lower = getOption(
                           "wordUniverse.to.lower",
                           TRUE
                         )) {
  univ.words <- unique(unlist(lapply(sssr, function(seq.sim.tbl) {
    hits.ind <- which(seq.sim.tbl$qseqid == prot.id)
    if (length(hits.ind) > 0) {
      unlist(lapply(
        seq.sim.tbl[hits.ind, ]$HRD,
        wordSet
      ))
    } else {
      character(0)
    }
  })))
  if (to.lower) {
    tolower(univ.words)
  } else {
    univ.words
  }
}

#' Generate short human readable descriptions (HRD) for all query proteins in
#' the argument sequence similarity search results 'sssr' using the Best Blast
#' and registered argument 'prot.scriber.score.funks'. Evaluate the performance
#' (F-Score) of the HRDs if reference word sets are available for the
#' respective query proteins (qseqid). NOTE that this function does NOT
#' annotate each domain separately, but treats all Blast Hits in a single go.
#'
#' @param sssr - A named list, in which names represent searched reference
#' sequence databases and values the read in tabular output (see function
#' parseSeqSimSearchTable).
#' @param ref.word.sets - A named list in which the names are query protein
#' identifiers ('qseqid') in lowercase and the values are character vectors
#' holding reference words ('the truth').
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
annotateProteinsAndEvaluatePerformance <- function(
    sssr,
    ref.word.sets, prot.scriber.score.funks = getOption(
      "statsOfPhrasesForQuery.score.funks",
      list(centered.inverse.inf.cntnt.mean = list(word.score.funk = function(x) {
        centeredWordScores(x, level.funk = mean)
      }, phrase.score.funk = sumWordScores), centered.inverse.inf.cntnt.median = list(word.score.funk = function(x) {
        centeredWordScores(x, level.funk = median)
      }, phrase.score.funk = sumWordScores), centered.inverse.inf.cntnt.quarterQuantile = list(word.score.funk = function(x) {
        centeredWordScores(x, level.funk = function(y) {
          quantile(y, probs = 0.25)
        })
      }, phrase.score.funk = sumWordScores), polynomial.word.scores = list(
        word.score.funk = polynomicWordScores,
        phrase.score.funk = sumWordScores
      ), centered.frequencies = list(
        word.score.funk = centeredLinearWordScores,
        phrase.score.funk = sumWordScores
      ))
    )) {
  prot.ids <- unique(unlist(lapply(sssr, function(seq.sim.search.tbl) {
    unique(seq.sim.search.tbl$qseqid)
  })))
  message(
    "Starting annotating human readable descriptions (HRD) for ",
    length(prot.ids), " query sequences."
  )
  do.call(rbind, mclapply(prot.ids, function(qseqid) {
    tryCatch(
      {
        q.phrases <- phrasesForQuery(qseqid, sssr)
        q.phrases.stats <- statsOfPhrasesForQuery(qseqid,
          q.phrases, ref.word.sets,
          score.funks = prot.scriber.score.funks
        )
        ref.words <- ref.word.sets[[tolower(qseqid)]]
        univ.words <- wordUniverse(qseqid, sssr)
        rbind(bestBlastHrds(
          qseqid, sssr, ref.words,
          univ.words
        ), do.call(rbind, lapply(
          names(prot.scriber.score.funks),
          function(ps.score.method) {
            bestProtScriberPhrases(
              qseqid, q.phrases.stats,
              ps.score.method, ref.words, univ.words
            )
          }
        )))
      },
      error = function(e) {
        #' browser()
        warning(
          "Query '", qseqid, "' caused an error:\n",
          e
        )
      }
    )
  }))
}

#' Multiple alignment region version of the functions namesake
#' 'annotateProteinsAndEvaluatePerformance'. It does the same as its namesake
#' except prot scriber generates one human readable description per disjoint
#' region of query sequence to which sequence similarity searches generated
#' local alignments. Also Best Blast annotations are not generated, as this is
#' already covered by the namesake function. Read its documentation for more
#' details. Note that ONLY proteins with more than a single alignment region
#' AND reference word sets will be processed.
#'
#' @param sssr - A named list, in which names represent searched reference
#' sequence databases and values the read in tabular output (see function
#' parseSeqSimSearchTable).
#' @param ref.word.sets - A named list in which the names are query protein
#' identifiers ('qseqid') in lowercase and the values are character vectors
#' holding reference words ('the truth').
#' @param alignment.regions.df - An instance of base::data.frame the result of
#' calling function 'allQueriesAlignmentRegions'.
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
annotateProteinsAndEvaluatePerformance.MultiRegion <- function(
    sssr,
    ref.word.sets, alignment.regions.df, prot.scriber.score.funks = getOption(
      "statsOfPhrasesForQuery.score.funks",
      list(centered.inverse.inf.cntnt.mean = list(word.score.funk = function(x) {
        centeredWordScores(x, level.funk = mean)
      }, phrase.score.funk = sumWordScores), centered.inverse.inf.cntnt.median = list(word.score.funk = function(x) {
        centeredWordScores(x, level.funk = median)
      }, phrase.score.funk = sumWordScores), centered.inverse.inf.cntnt.quarterQuantile = list(word.score.funk = function(x) {
        centeredWordScores(x, level.funk = function(y) {
          quantile(y, probs = 0.25)
        })
      }, phrase.score.funk = sumWordScores), polynomial.word.scores = list(
        word.score.funk = polynomicWordScores,
        phrase.score.funk = sumWordScores
      ), centered.frequencies = list(
        word.score.funk = centeredLinearWordScores,
        phrase.score.funk = sumWordScores
      ))
    )) {
  prot.id.cands <- unique(alignment.regions.df$Protein.ID[duplicated(alignment.regions.df$Protein.ID)])
  prot.ids <- prot.id.cands[tolower(prot.id.cands) %in%
    names(ref.word.sets)]
  rm(prot.id.cands) #' Clean up
  message(
    "Starting annotating human readable descriptions (HRD) for ",
    length(prot.ids), " query sequences - one HRD per disjoint alignment region."
  )
  do.call(rbind, mclapply(prot.ids, function(qseqid) {
    tryCatch(
      {
        ref.words <- ref.word.sets[[tolower(qseqid)]]
        qseqid.sssr.regions <- sssrForRegions(
          qseqid,
          alignment.regions.df, sssr
        )
        qseqid.region.hrds <- do.call(rbind, lapply(
          qseqid.sssr.regions,
          function(sssr.i) {
            q.phrases <- phrasesForQuery(
              qseqid,
              sssr.i
            )
            q.phrases.stats <- statsOfPhrasesForQuery(qseqid,
              q.phrases, ref.word.sets,
              score.funks = prot.scriber.score.funks
            )
            univ.words <- wordUniverse(
              qseqid,
              sssr.i
            )
            do.call(rbind, lapply(
              names(prot.scriber.score.funks),
              function(ps.score.method) {
                bestProtScriberPhrases(
                  qseqid,
                  q.phrases.stats, ps.score.method,
                  ref.words, univ.words
                )
              }
            ))
          }
        ))
        univ.words <- wordUniverse(qseqid, sssr)
        do.call(rbind, lapply(
          names(prot.scriber.score.funks),
          function(ps.score.method) {
            qrh.method <- qseqid.region.hrds[which(qseqid.region.hrds$Method ==
              ps.score.method), ]
            concat.hrd <- wordSet(paste(qrh.method$HRD,
              collapse = " "
            ), blacklist.regexs = NULL)
            concat.hrd.f.score.df <- fScore(
              concat.hrd,
              ref.words, univ.words, qseqid
            )
            concat.hrd.mcc.df <- matthewsCorrelationCoefficient(
              concat.hrd,
              ref.words, univ.words, qseqid
            )[
              ,
              c(
                "Protein.ID", "MCC", "MCC.relative",
                "MCC.best.possible", "univ.words",
                "ref.words", "univ.ref.words"
              )
            ]
            concat.hrd.eval.df <- merge(concat.hrd.f.score.df,
              concat.hrd.mcc.df,
              by = "Protein.ID"
            )
            concat.hrd.eval.df$Method <- paste0(
              ps.score.method,
              ".concat.regions"
            )
            concat.hrd.eval.df$Method.Score <- sum(qrh.method$Method.Score)
            concat.hrd.eval.df
          }
        ))
      },
      error = function(e) {
        #' browser()
        warning(
          "Query '", qseqid, "' caused an error:\n",
          e
        )
      }
    )
  }))
}

#' Measures performance of query sequence function predictors a.k.a
#' 'annotators', comparing the 'Best Blast' method with 'prot-scriber'. Note
#' that performance is measured in terms of Matthew's Correlation Coefficient,
#' F-Score, and their 'relative' counterparts. See respectively named functions
#' in this package for details.
#'
#' @param query.id - A string representing the biological sequence identifier
#' for which the competing annotators have produced annotations, i.e. human
#' readable descriptions (HRDs)
#' @param sssr - A named list, in which names represent searched reference
#' sequence databases and values the read in tabular output (see function
#' parseSeqSimSearchTable).
#' @param ps.res - An instance of base::data.frame with two columns
#' 'Annotee.Identifier' and 'Human.Readable.Description'. The table holds the
#' results of prot-scriber for the respective queries. Note that this table is
#' expected to have been produced by the Rust implementation of prot-scriber.
#' @param refs - A named list in which for the argument query.id a vector of
#' strings is expected to hold the 'truth', i.e. the words expected to be
#' reproduced by the competing annotators. Currently generated by Mercator and
#' HMMER3 on Pfam.
#' @param ps.na.hrd - A string representing the default HRD prot-scriber
#' annotates if none of the input Blast Hits (of the respective query) survive
#' blacklisting. This default HRD must be interpreted as NA and not literally
#' in order to avoid scoring this default HRD more negatively than necessary.
#' Default is 'ps.na.hrd = getOption('ps.na.hrd', 'unknown protein')'
#' @param ahrd.res - An instance of base::data.frame with two columns
#' 'Annotee.Identifier' and 'Human.Readable.Description' holding the results of
#' AHRD for the queries. Note that this is an optional argument.
#' @param use.fair.references - A logical (boolean) switch indicating whether
#' the word set of the references should be reduced to the intersection of gold
#' standard (Pfam-A and Mercator/Mapman Bin descriptions) and the words
#' appearing in ALL Blast Hit Descriptions. This makes the performance
#' calculation fair, because the competitors "Best Blast" and "prot-scriber"
#' can only use words that appear in Blast Hit Descriptions, they cannot
#' annotate words that appear in Pfam or Mercator but not in the Blast Hit
#' Descriptions. Hence those words should not result in worse scores.
#'
#' @return An instance of base::data.frame with the following columns:
#' 'Protein.ID', 'F.Score', 'precision' 'recall', 'fScore.beta', 'HRD'
#' 'n.words', 'F.Score.relative', 'F.Score.best.possible' 'MCC',
#' 'MCC.relative', 'MCC.best.possible' 'univ.words', 'ref.words',
#' 'univ.ref.words' 'Method', 'Method.Score'
#' @export
measurePredictionsPerformance <- function(
    query.id,
    sssr, ps.res, refs, ps.na.hrd = getOption(
      "ps.na.hrd",
      "unknown protein"
    ), ahrd.res = NULL, use.fair.references = TRUE) {
  if (query.id %in% names(refs)) {
    #' Best Blast Hits' performance:
    query.ref <- refs[[query.id]]
    #' The Rust implementation of prot-scriber actually adds words to the
    #' universe, because it used regular expressions and capture-groups to
    #' split stitle (Blast Hit descriptions) into words. See Rust module
    #' 'generate_hrd_associated_funcs.rs' function 'split_descriptions' for
    #' details. Thus, we need to add words from the prot-scriber annotation to
    #' the word universe:
    ps.words <- if (query.id %in% ps.res$Annotee.Identifier) {
      wordSet(ps.res[[query.id, "Human.Readable.Description"]],
        blacklist.regexs = NULL
      )
    } else {
      c()
    }
    #' To Do: Word universe should be the union of references and all words
    #' that appear in Blast Hit Descriptions (all lowercase)
    univ.words <- union(ps.words, wordUniverse(
      query.id,
      sssr
    ))
    #' TODO:
    #' intersect references and blat hit description word-sets to be fair.
    #'
    #' If requested reduce reference words to the "fair" list of references, i.e.
    #' only those words that the predictors "Best Blast" and "prot-scriber" can
    #' choose from which is the intersection of our "gold standard" (Pfam-A and
    #' Mercator/Mapman Bin descriptions; queries.ref) and the words appearing in
    #' Blast Hit descriptions (queries.sssr):
    if (use.fair.references) {
      #' For each query identifier in 'query.ids' intersect reference (the query's
      #' "ref"; see queries.ref) and prediction (the query's "sssr"; see
      #' queries.sssr) word-sets:
      #' To Do: Faensern please write a function (or two) that does the above.
      #' intersectPredictionsAndReferences(query.ids, queries.sssr, queries.ref)
      #' gene-id is key in queries.ref, and in the table queries.sssr[["swissprot"]]
      #' sp.tbl <- queries.sssr[["swissprot"]]
      #'
      query.ref <- intersect(refs, univ.words)
    }
    best.blast.performance <- bestBlastHrds(
      query.id,
      sssr, query.ref, univ.words
    )

    #' prot-scriber performance:
    query.ps.hrd <- if (query.id %in% ps.res$Annotee.Identifier) {
      ps.hrd <- ps.res[[query.id, "Human.Readable.Description"]]
      if (ps.hrd == ps.na.hrd) {
        #' prot-scriber annotates queries with a default HRD if none of
        #' the Blast Hits survived black-listing. This default HRD must
        #' be interpreted as NA, not literally:
        NA
      } else {
        ps.hrd
      }
    } else {
      NA
    }
    ps.n.words <- if (is.na(query.ps.hrd)) {
      0
    } else {
      length(strsplit(query.ps.hrd, " ")[[1]])
    }
    #' "Fake" the original experiment that compared several prot-scriber
    #' algorithms. We are now in a version that already has a stable
    #' implementation in Rust, i.e. we know what the best solution for
    #' prot-scriber is:
    query.ps.tbl <- data.frame(
      HRD = query.ps.hrd,
      n.words = ps.n.words, stringsAsFactors = FALSE
    )
    prot.scriber.performance <- bestProtScriberPhrases(
      query.id,
      query.ps.tbl, 2, query.ref, univ.words
    )
    #' Correct problem arising from legacy code which causes prot-scriber
    #' to appear as '2' in column 'Method' of data.frame
    #' 'prot.scriber.performance':
    prot.scriber.performance$Method <- "prot-scriber"


    #' ahrd performance:
    if (!is.null(ahrd.res)) {
      query.ahrd.hrd <- if (query.id %in% ahrd.res$Annotee.Identifier) {
        ahrd.hrd <- ahrd.res[[query.id, "Human.Readable.Description"]]
        if (tolower(ahrd.hrd) == tolower(ps.na.hrd)) {
          #' ahrd annotates queries with a default HRD if none of
          #' the Blast Hits survived black-listing. This default HRD must
          #' be interpreted as NA, not literally:
          NA
        } else {
          ahrd.hrd
        }
      } else {
        NA
      }
      ahrd.n.words <- if (is.na(query.ahrd.hrd)) {
        0
      } else {
        #' AHRD generates its HRD by choosing the best fitting Blast
        #' Hit description ('stitle'). This still must be split into
        #' words to count the number of words:
        length(wordSet(query.ahrd.hrd, blacklist.regexs = NULL))
      }
      query.ahrd.tbl <- data.frame(
        HRD = query.ahrd.hrd,
        n.words = ahrd.n.words, stringsAsFactors = FALSE
      )
      ahrd.performance <- bestProtScriberPhrases(
        query.id,
        query.ahrd.tbl, 2, query.ref, univ.words
      )
      #' Correct problem arising from legacy code which causes ahrd
      #' to appear as '2' in column 'Method' of data.frame
      #' 'ahrd.performance':
      ahrd.performance$Method <- "AHRD"
    }


    #' Result:
    rslt <- rbind(best.blast.performance, prot.scriber.performance)
    if (exists("ahrd.performance")) {
      rbind(rslt, ahrd.performance)
    } else {
      rslt
    }
  } else {
    message(
      "WARNING: Could not find query.id '",
      query.id, "' in references."
    )
    NULL
  }
}
