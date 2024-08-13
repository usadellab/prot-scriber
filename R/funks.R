#' Wrapper to base::gsub accepting just two arguments. Is required to make
#' several other functions available as plug-ins into `applyRegexList`.
#'
#' @param pattern - A string representation of a regular expression to be
#' replaced by empty strings.
#' @param subject - A string to be subjected to the replacement procedure,
#' implemented here.
#'
#' @return The transformed argument in which _all_ matches of argument
#' `pattern` have been replaced by an empty string.
#' @export
subWithEmptyString <- function(pattern, subject) {
  gsub(
    pattern = pattern, replacement = "",
    x = subject
  )
}

#' Recursively applies a vector of regular expressions upon the argument
#' subject. The function that actually is executed for each pair of subject,
#' the result of the previous iteration, and regular expression is defined in
#' the argument `regex.executor.funk`. Returns the final result, stripped of
#' leading and trailing white-spaces (`base::trimws`).
#'
#' @param subject - A string to be subjected recursively to all regular
#' expressions in argument `regexs` using argument function
#' `regex.executor.funk`.
#' @param regexs - A vector of regular expressions. Default is ``.
#' @param regex.executor.funk - A function accepting exactly two arguments. The
#' first must be the regular expression ('pattern') and the second the subject.
#' Is expected to return a string. Default is `subWithEmptyString`.
#'
#' @return String - The final result of recursive regex application.
#' @export
applyRegexList <- function(subject, regexs = getOption(
                             "applyRegexList.regexs",
                             filter.descline.regexs
                           ), regex.executor.funk = subWithEmptyString) {
  for (regex.i in regexs) {
    subject <- regex.executor.funk(
      regex.i,
      subject
    )
  }
  trimws(subject, which = "both")
}

#' Reads in a sequence similarity search (`Blast` or `Diamond`) result table
#' and optionally applies on each read in Hit `stitle` a set of regular
#' expressions to delete unwanted parts and retain a human readable description
#' in column `HRD`. Uses data.table::fread.
#'
#' @param path.to.table - A valid path to a seq-sim-search result table.
#' Separator MUST be TAB.
#' @param col.names - A character vector with column names. Default is
#' getOption('parseSeqSimSearchTable.col.names', c('qseqid', 'sseqid', 'qlen',
#' 'qstart', 'qend', 'slen', 'sstart', 'send', 'bitscore', 'stitle', 'evalue'))
#' @param parse.hrd - A boolean flagging whether to add a column `HRD` holding
#' the parsed out human readable descriptions extracted from each Hit's
#' `stitle` field.
#'
#' @return An instance of `data.table::data.table`
#' @export
parseSeqSimSearchTable <- function(
    path.to.table,
    col.names = getOption(
      "parseSeqSimSearchTable.col.names",
      c(
        "qseqid", "sseqid", "qlen", "qstart",
        "qend", "slen", "sstart", "send",
        "bitscore", "stitle", "evalue"
      )
    ),
    parse.hrd = TRUE) {
  sss.tbl <- fread(path.to.table,
    sep = "\t",
    quote = "", header = TRUE, na.strings = "",
    stringsAsFactors = FALSE, col.names = col.names,
    blank.lines.skip = TRUE, data.table = TRUE
  )
  if (parse.hrd) {
    sss.tbl$HRD <- as.character(unlist(mclapply(
      sss.tbl$stitle,
      function(stitle) {
        tryCatch(
          {
            applyRegexList(stitle)
          },
          error = function(e) {
            message(
              "stitle '", stitle,
              "' caused an error in function 'applyRegexList':\n",
              e
            )
            #' leave the stitle unchanged:
            stitle
          }
        )
      }
    )))
  }
  sss.tbl
}

#' Reads in a mercator4 table and parses it into a data.table
#'
#' @param path.to.table - A valid path to a seq-sim-search result table.
#' Separator MUST be TAB.
#' @param col.names - A character vector with column names. Default is
#' c('BINCODE', 'NAME', 'IDENTIFIER', 'DESCRIPTION', 'TYPE'))
#' @param exclude.bin.regex - A regular expression to be applied to the Mapman
#' Bincodes. Where this regex matches Bins are excluded from the result.
#' Default is '^(35|50)".
#'
#' @return An instance of `data.table::data.table` containing the Mercator4
#' tabular output.
#' @export
parseMercator4Tblout <- function(
    path.to.table,
    col.names = c("BINCODE", "NAME", "IDENTIFIER", "TYPE"),
    exclude.bin.regex = "^(35|50)") {
  m.dt <- fread(path.to.table,
    drop = "DESCRIPTION", sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, na.strings = "",
    quote = ""
  )
  m.dt$TYPE <- !is.na(m.dt$TYPE)
  #' If requested, remove Mercator results for BINCODEs that match the argument
  #' `exclude.bin.regex`:
  if (is.null(exclude.bin.regex) &&
    !is.na(exclude.bin.regex) &&
    exclude.bin.regex != "") {
    m.dt <- m.dt[!grepl(exclude.bin.regex, m.dt$BINCODE), ]
  }
  #' Remove leading and trailing single quotes from the respective columns.
  #' Note usage of the argument `quote = '\''` to the above `fread` command
  #' will actually not solve this issue, because some of the Map Man Bin
  #' DESCRIPTIONS actually _conatain_ a single quote. Hence, the manual
  #' removal of leading and trailing single quotes.
  for (col.i in c("BINCODE", "NAME", "IDENTIFIER")) {
    m.dt[[col.i]] <- sub("^'", "", sub(
      "'$",
      "", m.dt[[col.i]]
    ))
  }
  m.dt
}

#' Parse HMMER3 `--tblout` tabular output. Because this is a fixed width table
#' and no quoting is done, parsing it is quite a challange. This function uses
#' `data.table::fread` with a sophisticated combination of `sed` and `awk` to
#' parse the table.
#'
#' @param path.to.hmmr3.tblout - The valid file path to the table to be read in
#' @param read.syscmd - The system command to be applied to the argument
#' `path.to.hmmr3.tblout` before reading the result in with
#' `data.table::fread`. Default is
#' 'sed -e '1,3d' @FILE@ | awk -F \'\' 'match($0,/^(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(.+)$/,a){print a[1] \'\\t\' a[2] \'\\t\' a[3] \'\\t\' a[4] \'\\t\' a[5] \'\\t\' a[6] \'\\t\' a[7] \'\\t\' a[8] \'\\t\' a[9] \'\\t\' a[10] \'\\t\' a[11] \'\\t\' a[12] \'\\t\' a[13] \'\\t\' a[14] \'\\t\' a[15] \'\\t\' a[16] \'\\t\' a[17] \'\\t\' a[18] \'\\t\' a[19]}''
#' @param col.names - A character vector of column names to be used for the
#' read in table.
#'
#' @return An instance of `data.table::data.table` containing the HMMER3
#' `--tblout` content.
#' @export
parseHmmer3Tblout <- function(
    path.to.hmmr3.tblout,
    read.syscmd = "sed -e '1,3d' @FILE@ | awk -F \"\" 'match($0,/^(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(.+)$/,a){print a[1] \"\\t\" a[2] \"\\t\" a[3] \"\\t\" a[4] \"\\t\" a[5] \"\\t\" a[6] \"\\t\" a[7] \"\\t\" a[8] \"\\t\" a[9] \"\\t\" a[10] \"\\t\" a[11] \"\\t\" a[12] \"\\t\" a[13] \"\\t\" a[14] \"\\t\" a[15] \"\\t\" a[16] \"\\t\" a[17] \"\\t\" a[18] \"\\t\" a[19]}'",
    col.names = c(
      "target.name", "accession",
      "query.name", "accession.2", "E-value",
      "score", "bias", "E-value", "score", "bias",
      "exp", "reg", "clu", "ov", "env", "dom",
      "rep", "inc", "description.of.target"
    )) {
  fread(
    cmd = sub("@FILE@", path.to.hmmr3.tblout,
      read.syscmd,
      fixed = TRUE
    ), sep = "\t",
    header = FALSE, quote = "", na.strings = "",
    stringsAsFactors = FALSE, strip.white = TRUE,
    blank.lines.skip = TRUE, col.names = col.names
  )
}

#' Splits the argument human readable protein description (HRD) into words,
#' filters them to retain only informative ones, and returns them as an
#' alphabetically sorted set. Optionally returns all words in lower case.
#'
#' @param in.desc - The input HRD
#' @param split.regex - A regular expression used to split the argument HRD
#' (see base::strsplit). Default is
#' getOption('splitDescriptionIntoWordSet.split.regex',
#' '-|/|;|\\\\|,|:|\'|'|\\.|\\s+|\\||\\(|\\)')
#' @param blacklist.regex - A character vector of regular expressions to be
#' used to identify non meaningful words, i.e. those to be excluded from the
#' result. Set to 'NULL', 'c()', or 'NA' if you do not want blacklist
#' filtering. Default is
#' getOption('splitDescriptionIntoWordSet.blacklist.regexs',
#' blacklist.word.regexs)
#' @param lowercase.words - A boolean flag indicating whether to return all
#' words in lower case. Default is
#' getOption('splitDescriptionIntoWordSet.lowercase.words', TRUE)
#'
#' @return A character vector of informative words extracted from the argument
#' `in.desc` HRD. Returns 'character(0)' if identical(character(0), in.desc) ||
#' is.null(in.desc) || is.na(in.desc) || !is.character(in.desc) ||
#' nchar(in.desc) == 0.
#' @export
wordSet <- function(
    in.desc, split.regex = getOption(
      "splitDescriptionIntoWordSet.split.regex",
      "-|/|;|\\\\|,|:|\"|'|\\.|\\s+|\\||\\(|\\)"
    ),
    blacklist.regexs = getOption(
      "splitDescriptionIntoWordSet.blacklist.regexs",
      blacklist.word.regexs
    ), lowercase.words = getOption(
      "splitDescriptionIntoWordSet.lowercase.words",
      TRUE
    )) {
  if (identical(character(0), in.desc) || is.null(in.desc) ||
    is.na(in.desc) || !is.character(in.desc) ||
    nchar(in.desc) == 0) {
    character(0)
  } else {
    desc.words <- strsplit(in.desc,
      split = split.regex,
      perl = TRUE
    )[[1]]
    dw.retain.bool <- if (length(blacklist.regexs) >
      0) {
      "" != lapply(desc.words, applyRegexList,
        regexs = blacklist.regexs
      )
    } else {
      rep(TRUE, length(desc.words))
    }
    dw.set <- desc.words[dw.retain.bool]
    if (lowercase.words) {
      dw.set <- tolower(dw.set)
    }
    unique(dw.set)
  }
}

#' Function tests its namesake `wordSet`.
#'
#' @return TRUE if and only if all tests pass.
#' @export
testWordSet <- function() {
  x.test <- c(
    "World Hello", "Phototransferase",
    "Alien contig similar product subfamily predicted"
  )
  t1 <- identical(
    wordSet(x.test[[1]], lowercase.words = FALSE),
    c("World", "Hello")
  )
  t2 <- identical(wordSet(x.test[[2]]), "phototransferase")
  t3 <- identical(wordSet(x.test[[3]]), "alien")
  t4 <- identical(
    wordSet("Protein REDOX 1"),
    "redox"
  )
  t5 <- identical(
    wordSet("Probable mannitol dehydrogenase"),
    c("mannitol", "dehydrogenase")
  )
  all(t1, t2, t3, t4, t5)
}

#' The Map-Man4 Bin Ontology has a root Bin '50' spanning a sub-tree of Bins
#' whose classification has been extracted from KEGG's enzyme classification,
#' which in their DESCRIPTION field also refer to the human readable
#' description (HRD) of a Best Blast Hit (BBH). As we want to use both the
#' NAME and the DESCRIPTION fields as reference and want to outcompete the Best
#' Blast method, we need to exclude the HRD of the BBH from all Map-Man4 Bin
#' '50' DESCRIPTIONs. This is done by this function.
#'
#' @param mm4.anno.tbl - An instance of `data.table` and the result of e.g.
#' calling function `parseMercator4Tblout`.
#'
#' @return A modified version of argument `mm4.anno.tbl` in which HRD from BBH
#' have been removed from the respective MapMan4 Bins' DESCRIPTION fields.
#' @export
curateMercator4Annos <- function(mm4.anno.tbl) {
  bin.50.annos <- grepl("^50", mm4.anno.tbl$BINCODE) &
    mm4.anno.tbl$TYPE
  #' The DESCRIPTIONs contain also those of the Best Blast Hit, but the NAMEs do
  #' not:
  mm4.anno.tbl$DESCRIPTION[bin.50.annos] <- mm4.anno.tbl$NAME[bin.50.annos]
  mm4.anno.tbl
}

#' Extract all Map-Man Bin descriptions from the protein function annotations
#' in the argument table `mm4.anno.tbl`, split them into words and return them
#' as references for performance evaluation. See function `wordSet` for
#' details, particularly on how to provide its respective arguments as
#' environment options.
#'
#' @param mm4.anno.tbl - An instance of `data.table` result from using function
#' `parseMercator4Tblout`.
#' @param exclude.mm4.root.bins - A character vector of Map-Man 4 root Bins to
#' be excluded as reference annotations. Default is
#' getOption('referenceWordListFromMercator4Annos.exclude.mm4.root.bins',
#' c('35', '50'))
#' @param curate.annos.funk - A function receiving the single argument
#' `mm4.anno.tbl` and returning a subset of it. Can be used to preprocess e.g.
#' annotation of BIN 50 to filter out best Blast Hit descriptions. Default does
#' exactly that and is
#' getOption('referenceWordListFromMercator4Annos.curate.annos.funk',
#' curateMercator4Annos). Use `base::identity` to avoid anything being done to
#' the argument `mm4.anno.tbl`.
#'
#' @return A list with names being the protein identifiers and values character
#' vectors of reference words.
#' @export
referenceWordListFromMercator4Annos <- function(
    mm4.anno.tbl,
    exclude.mm4.root.bins = getOption(
      "referenceWordListFromMercator4Annos.exclude.mm4.root.bins",
      c("35", "50")
    )) {
  #' Exclude the Mercator results where BINCODE matches the following regular
  #' expression, according the argument `exclude.mm4.root.bins`:
  filter.i <- mm4.anno.tbl$TYPE & !grepl(
    paste0(
      "^", "(",
      paste(exclude.mm4.root.bins, collapse = "|"), ")"
    ),
    mm4.anno.tbl$BINCODE
  )
  mm4.fltrd.tbl <- mm4.anno.tbl[filter.i, ]
  uniq.prot.ids <- unique(mm4.fltrd.tbl$IDENTIFIER)
  setNames(mclapply(uniq.prot.ids, function(prot.id) {
    i <- which(mm4.fltrd.tbl$IDENTIFIER ==
      prot.id)
    mm4.descs <- unlist(mm4.fltrd.tbl[i, c("NAME")])
    unique(unlist(lapply(mm4.descs, wordSet)))
  }), uniq.prot.ids)
}


#' Extract all Pfam-A descriptions from the protein function annotations
#' in the argument table `pfamA.tbl`, split them into words and return them
#' as references for performance evaluation. See function `wordSet` for
#' details, particularly on how to provide its respective arguments as
#' environment options.
#'
#' @param pfamA.tbl - An instance of `data.table` result from using function
#' `parseMercator4Tblout`.
#'
#' @return A list with names being the protein identifiers and values character
#' vectors of reference words.
#' @export
referenceWordListFromPfamAAnnos <- function(pfamA.tbl) {
  uniq.prot.ids <- unique(pfamA.tbl$query.name)
  setNames(mclapply(uniq.prot.ids, function(prot.id) {
    i <- which(pfamA.tbl$query.name == prot.id)
    mm4.descs <- unlist(pfamA.tbl[i, "description.of.target"])
    unique(unlist(lapply(mm4.descs, wordSet)))
  }), uniq.prot.ids)
}

#' Merges the two references obtained from Mercator4 and PfamA annotations.
#' Note that, because Mercator4 lowercases protein identifiers the result of
#' this function will use lowercase protein identifiers, too. It's just easier
#' that way - Sorry!
#'
#' @param ref.mercator - A list and the result of calling function
#' `referenceWordListFromMercator4Annos`.
#' @param ref.pfamA - A list and the result of calling function
#' `referenceWordListFromPfamAAnnos`
#'
#' @return A merged list containing the set union of reference word sets found
#' for protein identifiers in both argument references. Note that protein IDs
#' will be always in lowercase (see above for the reason).
#' @export
mergeMercatorAndPfamAReferences <- function(
    ref.mercator,
    ref.pfamA) {
  names(ref.pfamA) <- tolower(names(ref.pfamA))
  prot.ids <- union(names(ref.mercator), names(ref.pfamA))
  setNames(mclapply(prot.ids, function(prot.id) {
    union(ref.mercator[[prot.id]], ref.pfamA[[prot.id]])
  }), prot.ids)
}

#' Analyzes the argument sequence similarity search (SSS) results 'sssr.tbl'
#' for the argument query 'prot.id' and identifies the regions in the query
#' sequence for which the SSS generated local alignments. Local alignments that
#' overlap are going to be merged.
#'
#' @param prot.id - A string representing the query protein identifier (qseqid)
#' @param sssr.tbl - An instance of data.table holding the tabular sequence
#' similarity search results to be analyzed. See function
#' 'parseSeqSimSearchTable' for details.
#'
#' @return An instance of base::data.frame with the following columns:
#' Protein.ID, start.pos, end.pos. One row for each local alignment region
#' covered by at least a single hit.
#' @export
findAlignmentRegions <- function(prot.id, sssr.tbl) {
  if (!is.null(sssr.tbl) && !is.na(sssr.tbl) &&
    nrow(sssr.tbl) > 0) {
    p.tbl <- sssr.tbl[which(sssr.tbl$qseqid ==
      prot.id), ]
    hit.ids <- unique(p.tbl$sseqid)
    p.seq <- rep(0, p.tbl$qlen[[1]])
    for (h.i in hit.ids) {
      hit.tbl <- p.tbl[which(p.tbl$sseqid ==
        h.i), ]
      if (nrow(hit.tbl) > 1) {
        hit.tbl <- hit.tbl[order(hit.tbl$bitscore,
          decreasing = TRUE
        ), ]
      }
      hit.algn.region <- hit.tbl[[1, "qstart"]]:hit.tbl[[
        1,
        "qend"
      ]]
      p.seq[hit.algn.region] <- p.seq[hit.algn.region] +
        1
    }
    in.gap <- TRUE
    reg.start <- -1
    prot.regions <- list()
    for (pos.i in 1:length(p.seq)) {
      n.hits.at.pos <- p.seq[[pos.i]]
      if (in.gap && n.hits.at.pos > 0) {
        reg.start <- pos.i
      }
      if (n.hits.at.pos == 0 || pos.i ==
        length(p.seq)) {
        if (!in.gap) {
          prot.regions[[length(prot.regions) +
            1]] <- data.frame(
            Protein.ID = prot.id,
            start.pos = reg.start, end.pos = pos.i,
            stringsAsFactors = FALSE
          )
        }
        in.gap <- TRUE
      } else {
        in.gap <- FALSE
      }
    }
    do.call(rbind, prot.regions)
  }
}

#' Identifies the regions sequence similarity search local alignments (Blast
#' Hits) cluster to. This is done for each query sequence in argument 'sssr'.
#'
#' @param sssr - A named list, in which names represent searched reference
#' sequence databases and values the read in tabular output (see function
#' parseSeqSimSearchTable).
#'
#' @return An instance of base::data.frame with columns 'Protein.ID',
#' 'start.pos', and 'end.pos' indicating for a given query each local alignment
#' region any hits aligned to.
#' @export
allQueriesAlignmentRegions <- function(sssr) {
  prot.ids <- unique(unlist(lapply(sssr, function(sssr.tbl) {
    unique(sssr.tbl$qseqid)
  })))
  all.sssr.tbl <- do.call(rbind, sssr)
  do.call(rbind, mclapply(prot.ids, function(p.id) {
    findAlignmentRegions(p.id, all.sssr.tbl)
  }))
}

#' For the argument 'qseqid' identifies whether the sequence similarity
#' searches have generated local alignments to disjoint regions of the query
#' sequence. If so, for each region the respective sssr tables are returned.
#'
#' @param qseqid - The query sequence's unique identifier in the form of a
#' scalar string.
#' @param alignmnt.regions - An instance of base::data.frame result of calling
#' function 'allQueriesAlignmentRegions'. The data.frame has three columns
#' 'Protein.ID', 'start.pos', and 'end.pos' indicating to which sequence region
#' in the query the sequence similarity searches generated local alignments
#' for. If non intersecting regions for the same query are found the
#' prot-scriber annotation is carried out independently for these and results
#' are concatonated, in the order of the regions.
#' @param sssr - A named list, in which names represent searched reference
#' sequence databases and values the read in tabular output (see function
#' parseSeqSimSearchTable).
#'
#' @return Either a list with a single entry, the unchanged argument 'sssr' or
#' a longer list, if and only if the argument query 'qseqid' has more than one
#' alignment region in argument 'alignmnt.regions'. In that case a list of
#' length equal to the number of disjoint alignment regions is returned. Each
#' list entry is a subset of argument 'sssr' only containing the search results
#' aligning to one region. Note that if the argument 'qseqid' is not found in
#' any sequence similarity search result table, instead of returning 'NULL' the
#' table is returned unchanged; so not to brake legacy code.
#' @export
sssrForRegions <- function(
    qseqid, alignmnt.regions,
    sssr) {
  q.align.regs.i <- which(alignmnt.regions$Protein.ID ==
    qseqid)
  if (length(q.align.regs.i) < 2) {
    list(sssr)
  } else {
    lapply(q.align.regs.i, function(a.r.i) {
      a.r <- alignmnt.regions[a.r.i, ]
      setNames(lapply(sssr, function(sssr.table) {
        if (qseqid %in% sssr.table$qseqid) {
          sssr.table[which(sssr.table$qseqid ==
            qseqid & sssr.table$qstart >=
            a.r$start.pos & sssr.table$qend <=
            a.r$end.pos), ]
        } else {
          sssr.table
        }
      }), names(sssr))
    })
  }
}

#' From a list of reference words generated for single proteins and a gene
#' family definition table generate a list of reference words for gene
#' families.
#'
#' @param protein.reference.words - A list with names being protein identifiers
#' and values character vectors of reference words.
#' @param gene.fam.tbl - An instance of data.table with two character columns:
#' "Annotee.Identifier", "Proteins". The former holds the gene family
#' identifiers and the latter character vectors of proteins that belong to the
#' respective gene families.
#' @param gene.family.prot.id.split.regex - A regular expression used to split
#' a string holding a list of protein identiers into a character vector of
#' those identifiers using base::strsplit. Default is "\\s*,\\s*|\\s+".
#'
#' @return A list with names being the gene family identifiers and values
#' character vectors of reference words. Only gene family identifiers appear as
#' names for families for which reference words could be found in the argument
#' 'protein.reference.words'.
#' @export
collectGeneFamilyReferenceWords <-
  function(protein.reference.words,
           gene.fam.tbl,
           gene.family.prot.id.split.regex = "[\\s,|]") {
    Reduce(append, mclapply(1:nrow(gene.fam.tbl), function(row.i) {
      fam.id <- gene.fam.tbl$Annotee.Identifier[[row.i]]
      fam.prot.refs <- list(unique(unlist(
        protein.reference.words[
          strsplit(
            gene.fam.tbl$Proteins[[row.i]],
            gene.family.prot.id.split.regex,
            perl = TRUE
          )[[1]]
        ]
      )))
      if (length(fam.prot.refs) > 0) {
        setNames(fam.prot.refs, fam.id)
      }
    }))
  }
