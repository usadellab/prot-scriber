#!/usr/bin/env Rscript

require(prot.scriber)


#' Parse command line options:
option_list <- list(
  make_option(
    c("-a", "--amino-acid-fasta"),
    type = "character",
    default = NULL,
    help = "The path to the Faba amino acid fasta file",
    metavar = "character"
  ),
  make_option(
    c("-o", "--longest-variants-fasta"),
    type = "character",
    default = NULL,
    help = "The directory into which to write the result binary 'p_coccineus_seq_sin_search.RData'",
    metavar = "character"
  ),
  make_option(
    c("-n", "--n-cores"),
    type = "numeric",
    default = detectCores(),
    help = "The number of cores to use in parallel processes.",
    metavar = "numeric"
  )
)

opt_parser <- OptionParser(option_list = option_list)
script.args <- parse_args(opt_parser)


#' Set up:
options(mc.cores = script.args[["n-cores"]])

#' Read the Faba amino acid sequences:
faba.aas <- read.fasta(script.args[["amino-acid-fasta"]], as.string = TRUE, 
    seqtype = "AA", strip.desc = TRUE)
faba.no.var.ids <- unique(sub("\\d+$", "", names(faba.aas)))


#' For genes with multiple splice variants find and use only the longest
#' variant:
faba.longest.ids <- unlist(mclapply(faba.no.var.ids, function(gene.id) {
    var.ids <- which(grepl(gene.id, names(faba.aas), fixed = TRUE))
    if (length(var.ids) > 1) {
        var.ids.len <- unlist(lapply(faba.aas[var.ids], nchar))
        names(faba.aas)[var.ids[which(var.ids.len == max(var.ids.len))]]
    } else names(faba.aas)[[var.ids]]
}))


#' Save results:
write.fasta(faba.aas[faba.longest.ids], faba.longest.ids, script.args[["longest-variants-fasta"]], 
    as.string = TRUE)


#' The end
message("DONE")
