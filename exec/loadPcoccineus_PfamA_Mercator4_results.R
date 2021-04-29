#!/usr/bin/env Rscript

require(prot.scriber)

option_list <- list(
  make_option(
    c("-m", "--mercator-table"),
    type = "character",
    default = NULL,
    help = "Valid file path to the result of running Mercator4",
    metavar = "character"
  ),
  make_option(
    c("-p", "--hmmer3-vs-pfamA-tblout"),
    type = "character",
    default = NULL,
    help = "Valid file path to `--tblout` result of running HMMER3 PfamA",
    metavar = "character"
  ),
  make_option(
    c("-d", "--data-dir"),
    type = "character",
    default = NULL,
    help = "The directory into which to write the result binary 'p_coccineus_seq_sim_search.RData'",
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

#' Load P.coccineus vs PfamA:
pc.mercator <- parseMercator4Tblout(args$`mercator-table`)

#' Load Mercator4 result for P.coccineus:
pc.pfamA <- parseHmmer3Tblout(args$`hmmer3-vs-pfamA-tblout`)

#' Save results:
save(pc.mercator, pc.pfamA, file = file.path(args$`data-dir`, 
    "p_coccineus_pfamA_Mercator4.RData"))

#' Generate reference word list from P.coccineus PfamA:
pc.ref.pfamA <- referenceWordListFromPfamAAnnos(pc.pfamA)

#' Generate reference word list from P.coccineus Mercator4:
pc.ref.mercator <- referenceWordListFromMercator4Annos(pc.mercator)

#' Save the two reference sets for future use:
save(pc.ref.mercator, pc.ref.pfamA, file = file.path(args$`data-dir`, 
    "p_coccineus_mercator_pfamA_references.RData"))

#' Merge (union) these reference lists:
pc.ref <- mergeMercatorAndPfamAReferences(pc.ref.mercator, pc.ref.pfamA)

#' Save merged references:
save(pc.ref, file = file.path(args$`data-dir`, "p_coccineus_reference_words.RData"))

#' DONE:
message("DONE")
