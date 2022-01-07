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
    help = "The directory into which to write the resulting binary 'faba_reference_words.RData'",
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

#' Load Faba vs PfamA:
faba.mercator <- parseMercator4Tblout(args$`mercator-table`)

#' Load Mercator4 result for Faba:
faba.pfamA <- parseHmmer3Tblout(args$`hmmer3-vs-pfamA-tblout`)

#' Save results:
save(faba.mercator, faba.pfamA, file = file.path(args$`data-dir`, 
    "faba_pfamA_Mercator4.RData"))

#' Generate reference word list from Faba PfamA:
faba.ref.pfamA <- referenceWordListFromPfamAAnnos(faba.pfamA)

#' Generate reference word list from Faba Mercator4:
faba.ref.mercator <- referenceWordListFromMercator4Annos(faba.mercator)

#' Save the two reference sets for future use:
save(faba.ref.mercator, faba.ref.pfamA, file = file.path(args$`data-dir`, 
    "faba_mercator_pfamA_references.RData"))

#' Merge (union) these reference lists:
faba.ref <- mergeMercatorAndPfamAReferences(faba.ref.mercator, faba.ref.pfamA)

#' Save merged references:
save(faba.ref, file = file.path(args$`data-dir`, "faba_reference_words.RData"))

#' DONE:
message("DONE")
