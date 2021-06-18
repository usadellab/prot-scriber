#!/usr/bin/env Rscript

require(prot.scriber)

#' Load data:
data("p_coccineus_seq_sim_search")
data("p_coccineus_reference_words")
data("p_coccineus_alignment_regions")

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
args <- parse_args(opt_parser)

#' Set mc.cores:
options(mc.cores = args$`n-cores`)

#' Generate human readable descriptions (HRDs) and evaluate their performance
#' for the P. coccineus query proteins:
pc.sssr <- list(Swissprot = pc.sprot, trEMBL = pc.trembl)
pc.hrds <- annotateProteinsAndEvaluatePerformance(pc.sssr, pc.ref, pc.alignment.regions.df)

#' Save results:
save(pc.hrds, file = file.path(args$`data-dir`, "p_coccineus_HRDs.RData"))

#' DONE
message("DONE")
