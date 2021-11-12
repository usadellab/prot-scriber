#!/usr/bin/env Rscript

require(prot.scriber)

option_list <- list(
  make_option(
    c("-s", "--seq-sim-search-vs-swissprot-tbl"),
    type = "character",
    default = NULL,
    help = "Valid file path to the sequence similarity search results of P.coccineus proteome against UniProt Swissprot",
    metavar = "character"
  ),
  make_option(
    c("-t", "--seq-sim-search-vs-trembl-tbl"),
    type = "character",
    default = NULL,
    help = "Valid file path to the sequence similarity search results of P.coccineus proteome against UniProt trEMBL",
    metavar = "character"
  ),
  make_option(
    c("-d", "--data-dir"),
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

#' Set mc.cores:
options(mc.cores = script.args$`n-cores`)

#' Load P.coccineus vs Swissprot:
me.sprot <- parseSeqSimSearchTable(script.args$`seq-sim-search-vs-swissprot-tbl`)

#' Load P.coccineus vs trEMBL:
me.trembl <- parseSeqSimSearchTable(script.args$`seq-sim-search-vs-trembl-tbl`)

#' Save results:
save(me.sprot, me.trembl, file = file.path(script.args$`data-dir`, "meta_euk_seq_sim_search.RData"))

#' DONE
message("DONE")
