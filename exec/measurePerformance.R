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


opt_parser <- OptionParser(option_list = option_list)
script.args <- parse_args(opt_parser)


#' Set mc.cores:
options(mc.cores = script.args$`n-cores`)

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
    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
queries.prot.scriber$Annotee.Identifier <- tolower(queries.prot.scriber$Annotee.Identifier)
#' Speed up lockup times:
rownames(queries.prot.scriber) <- queries.prot.scriber$Annotee.Identifier
#' Only process true prediction:
queries.prot.scriber <- queries.prot.scriber[which(queries.prot.scriber$Human.Readable.Description != 
    "unknown protein"), ]

#' All query identifiers that actually have both references and at least one
#' prediction either by Best Blast or prot-scriber:
query.ids <- intersect(union(queries.pfamA$query.name, queries.mercator[queries.mercator$TYPE, 
    "IDENTIFIER"]), unique(c(blast.sprot$qseqid, blast.trembl$qseqid, 
    queries.prot.scriber$Annotee.Identifier)))

#' All query identifier that have data for performance evaluation:
queries.sssr <- list(swissprot = blast.sprot, trembl = blast.trembl)

#' Get the gold standard (reference) prepared:
#' - from PfamA annotations:
queries.ref.pfamA <- referenceWordListFromPfamAAnnos(queries.pfamA)

#' - from Mercator (MapMan Bin) annotations:
queries.ref.mercator <- referenceWordListFromMercator4Annos(queries.mercator)

#' - join both above references:
queries.ref <- mergeMercatorAndPfamAReferences(queries.ref.mercator, 
    queries.ref.pfamA)


#' Measure performance for each query that has predictions AND references:
queries.performance <- do.call(rbind, mclapply(query.ids, function(q.id) {
    measurePredictionsPerformance(q.id, queries.sssr, queries.prot.scriber, 
        queries.ref)
}))


#' Save result table:
write.table(queries.performance, script.args$`result-table`, 
    sep = "\t", row.names = FALSE, quote = TRUE)


#' Test significance of differences in performance score distributions:
compared.methods <- paste0("BB.", names(queries.sssr))
perf.diff.sign.df <- do.call(rbind, lapply(c("F.Score", "F.Score.relative", 
    "MCC", "MCC.relative"), function(p.measr) {
    do.call(rbind, lapply(compared.methods, function(c.meth) {
        ps.t.test <- t.test(queries.performance[which(queries.performance$Method == 
            c.meth), p.measr], queries.performance[which(queries.performance$Method == 
            "prot-scriber"), p.measr], alternative = "less")$p.value
        ps.wilcox.test <- wilcox.test(queries.performance[which(queries.performance$Method == 
            c.meth), p.measr], queries.performance[which(queries.performance$Method == 
            "prot-scriber"), p.measr], alternative = "less")$p.value
        data.frame(Null.Hypothesis = paste(c.meth, "has better", 
            p.measr, "than prot-scriber"), t.test.p.value = ps.t.test, 
            wilcox.test.p.value = ps.wilcox.test, stringsAsFactors = FALSE)
    }))
}))
#' Correct for multiple hypothesis testing:
perf.diff.sign.df$t.test.p.adj <- p.adjust(perf.diff.sign.df$t.test.p.value, 
    method = "BH")
perf.diff.sign.df$wilcox.test.p.adj <- p.adjust(perf.diff.sign.df$t.test.p.value, 
    method = "BH")
#' Save hypothesis test results:
hyp.test.res.file <- sub("\\..+$", "_hypothesis_tests.txt", script.args$`result-table`)
write.table(perf.diff.sign.df, hyp.test.res.file, sep = "\t", 
    quote = TRUE, row.names = FALSE)


#' Generate plots:
plot.filename <- sub("\\.\\S+$", "", script.args$`result-table`)
competitors <- unique(queries.performance$Method)
plot.scores <- c("F.Score", "F.Score.relative", "MCC", "MCC.relative")
clrs <- brewer.pal(length(competitors), "Set1")

#' Boxplots for each score measure:
for (scr in plot.scores) {
    plot.file <- paste0(plot.filename, "_", scr, ".pdf")
    pdf(plot.file)
    i <- which(!is.na(queries.performance[[scr]]) & !is.infinite(queries.performance[[scr]]))
    plot.df <- queries.performance[i, c("Method", scr)]
    boxplot(formula(paste0(scr, " ~ Method")), data = plot.df, 
        col = clrs, horizontal = TRUE, ylab = scr)
    dev.off()
}


#' DONE
message("DONE")
