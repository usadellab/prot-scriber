#!/usr/bin/env Rscript

require(prot.scriber)

#' Load data:
data("p_coccineus_seq_sim_search")
data("p_coccineus_reference_words")


#' Parse command line options:
option_list <- list(
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

#' Ensure the LOCAL sequence similarity search results are read, and NOT the
#' ones installed previously with `R INSTALL`:
load('./data/p_coccineus_seq_sim_search.RData')

#' Generate human readable descriptions (HRDs) and evaluate their performance
#' for the P. coccineus query proteins:
pc.sssr <- list(Swissprot = pc.sprot, trEMBL = pc.trembl)

#' MOCK pc.ref - we do NOT want to evaluate performance, we JUST want to
#' annotate the query proteins:
pc.ref <- list()

#' Annotate the query proteins with short human readable descriptions:
pc.hrds <- annotateProteinsAndEvaluatePerformance(pc.sssr, pc.ref)

#' Save results:
save(pc.hrds, file = file.path(script.args$`data-dir`, "p_coccineus_HRDs.RData"))

#' DONE
message("DONE")
