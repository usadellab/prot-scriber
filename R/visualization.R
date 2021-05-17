#' Convert a human readable description (HRD) to a HTML snippet in which true
#' and false positives are marked with CSS classes. This is a helper used in
#' the brew template 'hrds_visualization_html_brew.txt'.
#'
#' @param hrd.str - A scalar string, the argument HRD to convert to HTML
#' @param ref.words - A character vector of lower case reference words used to
#' identify true and false positives.
#' @param false.pos.template - The HTML template used to mark false positives.
#' Default is '<span class=\'falsePos\'><%= w.i -%></span>'.
#' @param true.pos.template - The HTML template used to mark true positives.
#' Default is '<span class=\'truePos\'><%= w.i -%></span>'
#'
#' @return A scalar string containing the argument 'hrd.str' converted to HTML.
#' @export
hrdToHtmlHelper <- function(hrd.str, ref.words, false.pos.template = "<span class=\"falsePos\"><%= w.i -%></span>", 
    true.pos.template = "<span class=\"truePos\"><%= w.i -%></span>", 
    na.html = "<span class=\"naHrd\">NA</span>") {
    hrd.words <- wordSet(hrd.str, blacklist.regexs = NULL, lowercase.words = FALSE)
    if (length(hrd.words) > 0) {
        hrd.words.lc <- tolower(hrd.words)
        hrd.words.html <- lapply(1:length(hrd.words), function(i) {
            w.i <- hrd.words[[i]]
            w.i.lc <- hrd.words.lc[[i]]
            brew.template <- if (w.i.lc %in% ref.words) {
                true.pos.template
            } else {
                false.pos.template
            }
            capture.output(brew(text = brew.template))
        })
        paste(hrd.words.html, collapse = " ")
    } else na.html
}

#' Uses 'brew' to render a HTML document that visualizes the result of
#' competing methods to annotate query proteins with human readable
#' descriptions (HRD). See brew template 'hrds_visualization_html_brew.txt' for
#' more details.
#'
#' @param ps.hrd.df - An instance of base::data.frame result of calling
#' function 'annotateProteinsAndEvaluatePerformance'.
#' @param ps.hrd.ref - A list with lowercase protein identifiers as names and
#' character vectors of, also lowercase, reference words to identify true and
#' false positives in the respective annotated HRDs.
#' @param output.file - A connection or string referencing the HTML file to
#' generate.
#' @param brew.template - A connection or string referencing the file path to
#' the brew template. Default is file.path(path.package('prot.scriber'),
#' 'hrds_visualization_html_brew.txt').
#'
#' @return The resut of calling 'brew'.
#' @export
visualizeHrdResults <- function(ps.hrd.df, ps.hrd.ref, output.file, 
    brew.template = file.path(path.package("prot.scriber"), "hrds_visualization_html_brew.txt")) {
    brew(file = brew.template, output = output.file)
}
