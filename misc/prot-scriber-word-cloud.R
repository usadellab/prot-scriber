#' Import dependencies:
require(parallel)
require(RColorBrewer)
require(wordcloud)
require(wordcloud2)
require(htmlwidgets)
require(webshot)
webshot::install_phantomjs()

#' Usage:
message("\nUsage:\nRscript prot-scriber-word-cloud.R input-prot-scriber-table.txt output-files-name (the file-name without file-extension)\n")

#' Get arguments:
script.args <- commandArgs(trailingOnly = TRUE)
#' Validate:
if (length(script.args) != 2) {
    message("Please provide two command line arguments. See above 'Usage'.")
    stop()
}
prot.scriber.inp.tbl <- script.args[[1]]
output.files.name <- script.args[[2]]

#' Speed up computation:
options(mc.cores = (detectCores() - 1))

#' For reproducibility:
set.seed(1234)

#' Load prot-scriber result table:
ps.tbl <- read.table(prot.scriber.inp.tbl, sep = "\t",
    header = TRUE, stringsAsFactors = FALSE)


#' Exclude 'unknown protein' or 'unknown sequence family':
ps.tbl.clean <- ps.tbl[which(!ps.tbl$Human.Readable.Description %in%
    c("unknown protein", "unknown sequence family")),
    ]

#' Exclude the following words as non informative in terms of the word-cloud:
non.inf.words <- c("protein", "containing", "family",
    "domain", "isoform", "subunit", "associated", "duf",
    "factor", "dependent", "binding", "of", "and",
    "or")


#' Extract words and their frequencies:
all.words <- unlist(lapply(ps.tbl.clean$Human.Readable.Description,
    function(hrd) strsplit(hrd, split = " ")[[1]]))
unq.wrds <- setdiff(unique(all.words), non.inf.words)
#' Exclude only numbers or single letters:
unq.wrds <- unq.wrds[!grepl("^(\\d+|[a-z])$", unq.wrds,
    perl = TRUE)]
words.df <- data.frame(word = unq.wrds, freq = as.numeric(mclapply(unq.wrds,
    function(wrd) {
        length(which(all.words == wrd))
    })), stringsAsFactors = FALSE)


#' Generate word-cloud:
pdf(paste0(output.files.name, "_word_cloud_type_one.pdf"))
wordcloud(words = words.df$word, freq = words.df$freq,
    min.freq = quantile(words.df$freq, 1/3), max.words = 200,
    random.order = FALSE, rot.per = 0.35, colors = brewer.pal(8,
        "Dark2"))
dev.off()


#' Generate word-cloud 2:
wrd.cld.2 <- wordcloud2(data = words.df, size = 1.6,
    color = "random-dark")
#' save as HTML:
html.out <- paste0(output.files.name, "_word_cloud_type_two.html")
saveWidget(wrd.cld.2, html.out, selfcontained = F)
#' and as PDF:
webshot(html.out, paste0(output.files.name, "_word_cloud_type_two.pdf"),
    delay = 5, vwidth = 480, vheight = 480)


#' The End
message("DONE")
