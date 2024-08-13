#!/usr/bin/env Rscript

require(prot.scriber)

option_list <- list(
  make_option(
    c("-s", "--seq-sim-search-vs-swissprot-tbl"),
    type = "character",
    default = NULL,
    help = "Valid file path to the sequence similarity search results (Blast) against UniProt Swissprot",
    metavar = "character"
  ),
  make_option(
    c("-t", "--seq-sim-search-vs-trembl-tbl"),
    type = "character",
    default = NULL,
    help = "Valid file path to the sequence similarity search results (Blast) against UniProt trEMBL",
    metavar = "character"
  ),
  make_option(
    c("-f", "--gene-families"),
    type = "character",
    default = NULL,
    help = "Valid file path to the gene families result file",
    metavar = "character"
  ),
  make_option(
    c("-x", "--prot-scriber"),
    type = "character",
    default = NULL,
    help = "Valid file path to the prot-scriber result file",
    metavar = "character"
  ),
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
    c("-r", "--result-table"),
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

#' Read command line arguments:
opt_parser <- OptionParser(option_list = option_list)
script.args <- parse_args(opt_parser)


#' Set mc.cores:
options(mc.cores = script.args$`n-cores`)

#' Non default column names?
if (!is.null(script.args$`seq-sim-tbl-column-names`)) {
  seq.sim.tbl.col.names <- strsplit(
    script.args$`seq-sim-tbl-column-names`,
    " "
  )[[1]]
  message(
    "ALL input sequence similarity search result table are expected to have the following columns - but no actual header:\n",
    paste(seq.sim.tbl.col.names, collapse = ", ")
  )
  options(parseSeqSimSearchTable.col.names = seq.sim.tbl.col.names)
}

#' Load Blast Search results of queries vs Swissprot:
blast.sprot <- parseSeqSimSearchTable(script.args$`seq-sim-search-vs-swissprot-tbl`)
blast.sprot$qseqid <- tolower(blast.sprot$qseqid)

#' Load Blast Search results of queries vs trEMBL:
blast.trembl <- parseSeqSimSearchTable(script.args$`seq-sim-search-vs-trembl-tbl`)
blast.trembl$qseqid <- tolower(blast.trembl$qseqid)

#' Load Mercator4 results for queries:
queries.mercator <- parseMercator4Tblout(script.args$`mercator-table`)

#' Load results of HMMER3 searches of queries vs PfamA:
queries.pfamA <- parseHmmer3Tblout(script.args$`hmmer3-vs-pfamA-tblout`)
queries.pfamA$query.name <- tolower(queries.pfamA$query.name)

#' Load prot-scriber results for queries:
queries.prot.scriber <- read.table(script.args$`prot-scriber`,
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
#' Speed up lockup times:
rownames(queries.prot.scriber) <- queries.prot.scriber$Annotee.Identifier

#' Load gene families results for queries:
queries.gene.families <- read.table(script.args$`gene-families`,
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
#' Change header
colnames(queries.gene.families) <- c("Annotee.Identifier", "Proteins")

#' Speed up lockup times:
rownames(queries.gene.families) <- queries.gene.families$Annotee.Identifier

#' merge HRD result table with gene families table
#' - In this, ignore prot-scriber annotations for single genes (see
#' prot-scriber's -a command line flag); we achieve this with the below `all.y
#' = FALSE`.
gene.families.expanded <- merge(
  x = queries.gene.families,
  y = queries.prot.scriber,
  by = "Annotee.Identifier",
  all.x = TRUE, all.y = FALSE
)

#' Only process true prediction:
gene.families.expanded <- gene.families.expanded[
  which(
    gene.families.expanded$Human.Readable.Description !=
      "unknown sequence family"
  ),
]

#' lowercase Gene Family and Protein Identifiers:
gene.families.expanded$Proteins <- tolower(gene.families.expanded$Proteins)
gene.families.expanded$Annotee.Identifier <- tolower(gene.families.expanded$Annotee.Identifier)
#' Speed up lookup:
rownames(gene.families.expanded) <- gene.families.expanded$Annotee.Identifier

#' All query identifier that have data for performance evaluation:
queries.sssr <- list(swissprot = blast.sprot, trembl = blast.trembl)

#' Get the gold standard (reference) prepared:
#' - from PfamA annotations:
queries.ref.pfamA <- referenceWordListFromPfamAAnnos(queries.pfamA)

#' - from Mercator (MapMan Bin) annotations:
queries.ref.mercator <- referenceWordListFromMercator4Annos(queries.mercator)

#' - join both above references:
queries.prot.ref <- mergeMercatorAndPfamAReferences(
  queries.ref.mercator,
  queries.ref.pfamA
)

#' - set union all reference words found for each gene family's proteins:
queries.ref <- collectGeneFamilyReferenceWords(
  queries.prot.ref, gene.families.expanded
)

#' Prepare the "word universes" for each gene family, i.e. the mathematical
#' sets of words that appear in any Blast Hit Human Readable Description that
#' was found for proteins belonging to a given gene family. These universe word
#' sets are used to calculate the true negative rates:
word.univs <- geneFamilyWordUniverses(gene.families.expanded, queries.sssr)

#' Measure performance scores for the prot-scriber generated gene family HRDs:
gene.fam.prfrmnc.tbl <- do.call(rbind, lapply(
  gene.families.expanded$Annotee.Identifier, function(query.id) {
    measurePredictionPerformanceForGeneFamily(
      query.id, gene.families.expanded, queries.ref, word.univs
    )
  }
))

#' Write result:
write.table(
  gene.fam.prfrmnc.tbl, script.args$`result-table`,
  row.names = FALSE, sep = "\t", quote = FALSE
)

#' DONE
message("DONE")
