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
    gsub(pattern = pattern, replacement = "", x = subject)
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
applyRegexList <- function(subject, regexs = getOption("applyRegexList.regexs", 
    filter.descline.regexs), regex.executor.funk = subWithEmptyString) {
    for (regex.i in regexs) {
        subject <- regex.executor.funk(regex.i, subject)
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
#' `c('qseqid', 'sseqid', 'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send',
#' 'bitscore', 'stitle', 'evalue')`
#' @param parse.hrd - A boolean flagging whether to add a column `HRD` holding
#' the parsed out human readable descriptions extracted from each Hit's
#' `stitle` field.
#'
#' @return An instance of `data.table::data.table`
#' @export
parseSeqSimSearchTable <- function(path.to.table, col.names = c("qseqid", 
    "sseqid", "qlen", "qstart", "qend", "slen", "sstart", "send", 
    "bitscore", "stitle", "evalue"), parse.hrd = TRUE) {
    sss.tbl <- fread(path.to.table, sep = "\t", quote = "", header = TRUE, 
        na.strings = "", stringsAsFactors = FALSE, col.names = col.names, 
        blank.lines.skip = TRUE, data.table = TRUE)
    if (parse.hrd) {
        sss.tbl$HRD <- as.character(unlist(mclapply(sss.tbl$stitle, 
            applyRegexList)))
    }
    sss.tbl
}

#' Reads in a mercator4 table and parses it into a data.table
#'
#' @param path.to.table - A valid path to a seq-sim-search result table.
#' Separator MUST be TAB.
#' @param col.names - A character vector with column names. Default is
#' c('BINCODE', 'NAME', 'IDENTIFIER', 'DESCRIPTION', 'TYPE'))
#'
#' @return  An instance of `data.table::data.table` containing the Mercator4
#' `--tblout` content.
#' @export
parseMercator4Tblout <- function(path.to.table, col.names = c("BINCODE", 
    "NAME", "IDENTIFIER", "DESCRIPTION", "TYPE")) {
    m.dt <- fread(path.to.table, sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
        na.strings = "", quote = "")
    m.dt$TYPE <- !is.na(m.dt$TYPE)
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
parseHmmer3Tblout <- function(path.to.hmmr3.tblout, read.syscmd = "sed -e '1,3d' @FILE@ | awk -F \"\" 'match($0,/^(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(.+)$/,a){print a[1] \"\\t\" a[2] \"\\t\" a[3] \"\\t\" a[4] \"\\t\" a[5] \"\\t\" a[6] \"\\t\" a[7] \"\\t\" a[8] \"\\t\" a[9] \"\\t\" a[10] \"\\t\" a[11] \"\\t\" a[12] \"\\t\" a[13] \"\\t\" a[14] \"\\t\" a[15] \"\\t\" a[16] \"\\t\" a[17] \"\\t\" a[18] \"\\t\" a[19]}'", 
    col.names = c("target.name", "accession", "query.name", "accession", 
        "E-value", "score", "bias", "E-value", "score", "bias", 
        "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", 
        "description.of.target")) {
    fread(cmd = sub("@FILE@", path.to.hmmr3.tblout, read.syscmd, 
        fixed = TRUE), sep = "\t", header = FALSE, quote = "", 
        na.strings = "", stringsAsFactors = FALSE, strip.white = TRUE, 
        blank.lines.skip = TRUE, col.names = col.names)
}

#' Splits the argument human readable protein description (HRD) into words,
#' filters them to retain only informative ones, and returns them as an
#' alphabetically sorted set. Optionally returns all words in lower case.
#'
#' @param in.desc - The input HRD
#' @param split.regex - A regular expression used to split the argument HRD
#' (see base::strsplit). Default is
#' getOption('splitDescriptionIntoWordSet.spit.regex', '\\s+|\\.')
#' @param blacklist.regex - A character vector of regular expressions to be used to
#' identify non meaningful words, i.e. those to be excluded from the result.
#' Default is getOption('splitDescriptionIntoWordSet.spit.regex', '\\s+|\\.')
#' @param lowercase.words - A boolean flag indicating whether to return all
#' words in lower case. Default is
#' getOption('splitDescriptionIntoWordSet.lowercase.words', TRUE)
#'
#' @return A character vector of informative words extracted from the argument
#' `in.desc` HRD.
#' @export
wordSet <- function(in.desc, split.regex = getOption("splitDescriptionIntoWordSet.spit.regex", 
    "\\s+|\\."), blacklist.regexs = getOption("splitDescriptionIntoWordSet.blacklist.regexs", 
    blacklist.word.regexs), lowercase.words = getOption("splitDescriptionIntoWordSet.lowercase.words", 
    TRUE)) {
    desc.words <- strsplit(in.desc, split = split.regex, perl = TRUE)[[1]]
    dw.retain.bool <- "" != lapply(desc.words, applyRegexList, 
        regexs = blacklist.word.regexs)
    dw.set <- desc.words[dw.retain.bool]
    if (lowercase.words) {
        dw.set <- tolower(dw.set)
    }
    sort(unique(dw.set))
}

#' Function tests its namesake `wordSet`.
#'
#' @return TRUE if and only if all tests pass.
#' @export
testWordSet <- function() {
    x.test <- c("World Hello", "Phototransferase", "Alien contig similar product subfamily predicted")
    t1 <- identical(wordSet(x.test[[1]], lowercase.words = FALSE), 
        c("Hello", "World"))
    t2 <- identical(wordSet(x.test[[2]]), "phototransferase")
    t3 <- identical(wordSet(x.test[[3]]), "alien")
    all(t1, t2, t3)
}
