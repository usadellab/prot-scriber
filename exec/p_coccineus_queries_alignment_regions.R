#!/usr/bin/env Rscript

require(prot.scriber)

#' Load data:
data("p_coccineus_seq_sim_search")

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

#' Find alignment regions for all P. coccineus query proteins:
pc.sssr <- list(Swissprot = pc.sprot, trEMBL = pc.trembl)
pc.alignment.regions.df <- allQueriesAlignmentRegions(pc.sssr)

#' Save results:
save(pc.alignment.regions.df, file = file.path(args$`data-dir`, 
    "p_coccineus_alignment_regions.RData"))

#' DONE
message("DONE")
