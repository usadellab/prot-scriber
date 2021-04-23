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
#' @param path.to.table - A valid path to a seq-sim-search result table. Separator MUST be TAB.
#' @param col.names - A character vector with column names. Default is
#' `c("qseqid", "sseqid", "qlen", "qstart", "qend", "slen", "sstart", "send",
#' "bitscore", "stitle", "evalue")`
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

#' Reads in a mercator4 table and parses its 
#'
#' @param path.to.table - A valid path to a seq-sim-search result table. Separator MUST be TAB.
#' @param col.names - A character vector with column names. Default is
#' c("BINCODE", "NAME", "IDENTIFIER", "DESCRIPTION", "TYPE"))
#
#'
#' @return
#' @export
parseMercator4Tbl <- function(path.to.table, col.names = c("BINCODE", "NAME", "IDENTIFIER", "DESCRIPTION", "TYPE")) {
   m.dt <- fread(path.to.table, sep = "\t", header = TRUE, stringsAsFactors = FALSE,na.strings = "", quote = "")
   m.dt$TYPE <- ifelse(is.na(dt$TYPE), FALSE, TRUE)
   m.dt
}
