#' FUNCTIONS

toWords <- function(x) {
    strsplit(x, " ")[[1]]
}

sharedWords <- function(all.descs, closest.desc) {
    uniq.words <- unique(unlist(all.descs))
    uniq.words[uniq.words %in% unlist(closest.desc)]
}


#' Input
candidates <- readLines("./misc/Soltu.DM.02G015700.1_AA_vs_trEMBL_blastpout_candidate_descriptions.txt")
blast.tbl <- read.table("./misc/Soltu.DM.02G015700.1_AA_vs_trEMBL_blastpout.txt", 
    sep = "\t", comment.char = "", na.strings = "", quote = "", 
    stringsAsFactors = FALSE)


#' Best match in terms of Blast bit-score (column V10)
closest.desc.i <- which(blast.tbl$V10 == max(blast.tbl$V10))


#' Process
cands.words <- lapply(candidates, toWords)
shared.words <- sharedWords(cands.words, cands.words[closest.desc.i])
