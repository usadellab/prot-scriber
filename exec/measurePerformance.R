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
    c("-a", "--ahrd"),
    type = "character",
    default = NULL,
    help = "Valid file path to the AHRD result file",
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
  ),
  make_option(
    c("-k", "--seq-sim-tbl-column-names"),
    type = "character",
    default = NULL,
    help = "The column names in Diamond terminology that are used in ALL input sequence similarity search result tables (-s, -t). This is optional. Default is 'qseqid sseqid qlen qstart qend slen sstart send bitscore stitle evalue'!",
    metavar = "character"
  )
)

#' Read command line arguments:
opt_parser <- OptionParser(option_list = option_list)
script.args <- parse_args(opt_parser)


#' Set mc.cores:
options(mc.cores = script.args$`n-cores`)

#' Non default column names?
if (!is.null(script.args$`seq-sim-tbl-column-names`)) {
    seq.sim.tbl.col.names <- strsplit(script.args$`seq-sim-tbl-column-names`, 
        " ")[[1]]
    message("ALL input sequence similarity search result table are expected to have the following columns - but no actual header:\n", 
        paste(seq.sim.tbl.col.names, collapse = ", "))
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

#' If present load AHRD results for queries:
queries.ahrd <- if (!is.null(script.args$ahrd)) {
    ahrd.tbl <- read.table(script.args$ahrd, sep = "\t", 
        header = TRUE, comment.char = "", quote = "", 
        na.strings = "", skip = 2, stringsAsFactors = FALSE)[, 
        c("Protein.Accession", "Human.Readable.Description")]
    #' Make AHRD result table look like a prot-scriber result table:
    colnames(ahrd.tbl) <- c("Annotee.Identifier", "Human.Readable.Description")
    #' Mercator lowercases the protein identifiers:
    ahrd.tbl$Annotee.Identifier <- tolower(ahrd.tbl$Annotee.Identifier)
    #' Speed up lookup times:
    rownames(ahrd.tbl) <- ahrd.tbl$Annotee.Identifier
    #' Only process true predictions:
    ahrd.tbl[which(ahrd.tbl$Human.Readable.Description != 
        "Unknown protein"), ]
} else NULL

#' Load prot-scriber results for queries:
queries.prot.scriber <- read.table(script.args$`prot-scriber`, 
    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
queries.prot.scriber$Annotee.Identifier <- tolower(queries.prot.scriber$Annotee.Identifier)
#' Speed up lockup times:
rownames(queries.prot.scriber) <- queries.prot.scriber$Annotee.Identifier
#' Only process true prediction:
queries.prot.scriber <- queries.prot.scriber[which(queries.prot.scriber$Human.Readable.Description != 
    "unknown protein"), ]

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

#' All query identifiers that actually have both references and at least one
#' prediction either by Best Blast or prot-scriber:
query.ids <- intersect(names(queries.ref), unique(c(blast.sprot$qseqid, 
    blast.trembl$qseqid, queries.prot.scriber$Annotee.Identifier, 
    queries.ahrd$Annotee.Identifier)))

#' Validate:
if(length(query.ids) == 0) {
  stop("Could not find any query identifiers that were present in both the references and the Blast/Diamond and prot-scriber results. Please check your input.")
}

#' Measure performance for each query that has predictions AND references:
queries.performance <- do.call(rbind, mclapply(query.ids, 
    function(q.id) {
            measurePredictionsPerformance(q.id, queries.sssr, 
                queries.prot.scriber, queries.ref, 
                ahrd.res = queries.ahrd)
    }))


#' Save result table:
write.table(queries.performance, script.args$`result-table`, 
    sep = "\t", row.names = FALSE, quote = TRUE)
message("Written performance scores into '", script.args$`result-table`, 
    "'.")


#' Test significance of differences in performance score distributions - in
#' this, only consider performance measurements, where there was some agreement
#' between the reference and the prediction 'ontologies' (Blast stitles on the
#' one hand and Pfam-A and Mapman Bins on the other):
compared.methods <- paste0("BB.", names(queries.sssr))
if (!is.null(queries.ahrd)) {
    compared.methods <- c(compared.methods, "AHRD")
}
queries.performance.san <- queries.performance[which(queries.performance$F.Score.best.possible > 
    0), ]

perf.diff.sign.df <- do.call(rbind, lapply(c("F.Score", 
    "F.Score.relative", "MCC", "MCC.relative"), function(p.measr) {
    do.call(rbind, lapply(compared.methods, function(c.meth) {
        #' In case of MCC only compare performance where both prot-scriber and 'c.meth'
        #' have made a prediction. This is done, because the 'ontologies' of the
        #' references, i.e. Mapman Bin annotations and Pfam-A annotations are quite
        #' different from the UniProt based Blast Hit descriptions ('stitle'). In the
        #' case of MCC this leads to an 'NA' annotation outperforming a prediction,
        #' even though the prediction might have some true positives.
        perf.tbl <- if (grepl("^MCC", p.measr, perl = TRUE)) {
            i <- which(queries.performance.san$Method %in% 
                c("prot-scriber", c.meth) & !is.na(queries.performance.san$HRD))
            queries.performance.san[i, ]
        } else {
            queries.performance.san
        }
        ps.t.test <- t.test(perf.tbl[which(perf.tbl$Method == 
            c.meth), p.measr], perf.tbl[which(perf.tbl$Method == 
            "prot-scriber"), p.measr], alternative = "less")$p.value
        ps.wilcox.test <- wilcox.test(perf.tbl[which(perf.tbl$Method == 
            c.meth), p.measr], perf.tbl[which(perf.tbl$Method == 
            "prot-scriber"), p.measr], alternative = "less")$p.value
        data.frame(Null.Hypothesis = paste(c.meth, 
            "has better", p.measr, "than prot-scriber"), 
            t.test.p.value = ps.t.test, wilcox.test.p.value = ps.wilcox.test, 
            stringsAsFactors = FALSE)
    }))
}))
#' Correct for multiple hypothesis testing:
perf.diff.sign.df$t.test.p.adj <- p.adjust(perf.diff.sign.df$t.test.p.value, 
    method = "BH")
perf.diff.sign.df$wilcox.test.p.adj <- p.adjust(perf.diff.sign.df$t.test.p.value, 
    method = "BH")
#' Save hypothesis test results:
hyp.test.res.file <- sub("\\.[^.]+$", "_hypothesis_tests.txt", 
    script.args$`result-table`)
write.table(perf.diff.sign.df, hyp.test.res.file, sep = "\t", 
    quote = TRUE, row.names = FALSE)
message("Written significance tests of differences in performance scores between compared annotators into '", 
    hyp.test.res.file, "'.")


#' Generate plots:
plot.filename <- sub("\\.[^.]+$", "", script.args$`result-table`)
competitors <- unique(queries.performance.san$Method)
plot.scores <- c("F.Score", "F.Score.relative", "MCC", 
    "MCC.relative")
clrs <- brewer.pal(length(competitors), "Set1")

#' Boxplots for each score measure:
for (scr in plot.scores) {
    plot.file <- paste0(plot.filename, "_", scr, ".pdf")
    pdf(plot.file)
    i <- which(!is.na(queries.performance.san[[scr]]) & 
        !is.infinite(queries.performance.san[[scr]]))
    plot.df <- queries.performance.san[i, c("Method", 
        scr)]
    boxplot(formula(paste0(scr, " ~ Method")), data = plot.df, 
        col = clrs, horizontal = TRUE, ylab = scr)
    dev.off()
}
message("Generated boxplots of the score distributions in files '", 
    plot.filename, "_*.pdf'.")


#' DONE
message("DONE")
