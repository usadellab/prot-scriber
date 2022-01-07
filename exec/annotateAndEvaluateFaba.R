#!/usr/bin/env Rscript

require(prot.scriber)

#' Load data:
data("faba_seq_sim_search")
data("faba_reference_words")


#' Parse command line options:
option_list <- list(
  make_option(
    c("-d", "--data-dir"),
    type = "character",
    default = NULL,
    help = "The directory into which to write the result binary 'faba_HRDs.RData'",
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

#' Generate human readable descriptions (HRDs) and evaluate their performance
#' for the P. coccineus query proteins:
faba.sssr <- list(Swissprot = faba.sprot, trEMBL = faba.trembl)
faba.hrds <- annotateProteinsAndEvaluatePerformance(faba.sssr, faba.ref)

#' Save results:
save(faba.hrds, file = file.path(script.args$`data-dir`, "faba_HRDs.RData"))

#' DONE
message("DONE")
