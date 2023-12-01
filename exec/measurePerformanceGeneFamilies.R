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
queries.pfamA<- parseHmmer3Tblout(script.args$`hmmer3-vs-pfamA-tblout`)
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
                   "Unknown protein"), ] #NA in the datatable
} else NULL


#########################################################
#' New code for gene families:

#' Load prot-scriber results for queries:
queries.prot.scriber <- read.table(script.args$`prot-scriber`, 
                                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' Speed up lockup times:
rownames(queries.prot.scriber) <- queries.prot.scriber$Annotee.Identifier

#' Load gene families results for queries:
queries.gene.families <- read.table(script.args$`gene-families`, 
                                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' Change header 
colnames(queries.gene.families)<-c("Annotee.Identifier", "Proteins")

#' Speed up lockup times:
rownames(queries.gene.families) <- queries.gene.families$Annotee.Identifier

#' merge HRD result table with gene families table
gene.families.expanded<-merge(x=queries.gene.families, y=queries.prot.scriber, by="Annotee.Identifier", all.x=TRUE, all.y=TRUE)

#' Only process true prediction:
gene.families.expanded <- gene.families.expanded[which(gene.families.expanded$Human.Readable.Description != 
                                                         "unknown sequence family"), ]

#'prot-scriber -a flag: -annotate-non-family-queries. Copy query ids to $Proteins column (2)
substrFamily<-"Seq-Fam_"
for(row in 1:nrow(gene.families.expanded)){
  substrID<-substr(gene.families.expanded[row, 1], 0, 8)
  if(substrFamily!=substrID){
    gene.families.expanded[row, 2] <- paste(gene.families.expanded[row, 1])
  }
}

#'lowercase query.id and Annotee.Identifier
gene.families.expanded$Proteins <- tolower(gene.families.expanded$Proteins)
gene.families.expanded$Annotee.Identifier <- tolower(gene.families.expanded$Annotee.Identifier)

#' add some columns to check intermediate results. 
gene.families.expanded <- cbind(gene.families.expanded,wordSet_hrd=NA, Merged.Ref.Family.Description=NA, blastWordUniverse=NA,univ.words=NA, mcc.simple=NA)

#########################################################

#' All query identifier that have data for performance evaluation:
queries.sssr <- list(swissprot = blast.sprot, trembl = blast.trembl)

#' Get the gold standard (reference) prepared:
#' - from PfamA annotations:
queries.ref.pfamA <- referenceWordListFromPfamAAnnos(queries.pfamA)
#########################################################
#new code for mercator4 references:

#' remove "DESCRIPTION" column from mercator4 result table
q.m <-subset(queries.mercator, select = -4)

#' test function from funks.R
#' filter 35 and 50 BINs
filterMercator4Table <- function(mm4.tbl, exclude.bin.regexs = c("35.", "50."), bincode.col = "BINCODE") {
  i <- do.call(`&`, lapply(exclude.bin.regexs,
                           function(bin.rgx) {
                             !grepl(bin.rgx, mm4.tbl[[bincode.col]],
                                    perl = TRUE)
                           }))
  mm4.tbl[i, ]
}
#' apply robust 35&50 filter
q.m<-filterMercator4Table(q.m)

#' overwritten function from funks.R without "DESCRIPTION" and 35 & 50 Bins
#' removed curateMercatorAnnos 
referenceWordListFromMercator4Annos_geneFamilies <- function(mm4.anno.tbl,
                                                       exclude.mm4.root.bins = getOption("referenceWordListFromMercator4Annos.exclude.mm4.root.bins",
                                                                                         c("35"))) {
  mm4.copy.tbl <- mm4.anno.tbl
  filter.i <- mm4.copy.tbl$TYPE & !grepl(paste0("^", "(",
                                                paste(exclude.mm4.root.bins, collapse = "|"), ")"),
                                         mm4.anno.tbl$BINCODE)
  mm4.fltrd.tbl <- mm4.copy.tbl[filter.i,
  ]
  uniq.prot.ids <- unique(mm4.fltrd.tbl$IDENTIFIER)
  setNames(mclapply(uniq.prot.ids, function(prot.id) {
    i <- which(mm4.fltrd.tbl$IDENTIFIER ==
                 prot.id)
    mm4.descs <- unlist(mm4.fltrd.tbl[i, c("NAME")])
    unique(unlist(lapply(mm4.descs, wordSet)))
  }), uniq.prot.ids)
}

#' - from Mercator (MapMan Bin) annotations:
queries.ref.mercator <- referenceWordListFromMercator4Annos_geneFamilies(q.m)

#' - join both references:
queries.ref <- mergeMercatorAndPfamAReferences(queries.ref.mercator, 
                                               queries.ref.pfamA)

#########################################################
#' New code for gene families:

#'apply filter wordSet prot-scriber HRD gene family results
#'find reference words set for each gene family
#'find blast word universe set for each gene family
#'find univ.words set for each gene family

rowList<-split(gene.families.expanded, 1:nrow(gene.families.expanded))
wordSetHRD<-mclapply(rowList,
                     function(row) {   
                       split.prot <-unlist(strsplit(row[,2], split=","))
                       all.ref.words<-character()
                       blast.word.universe<-character()
                       for (protein in (split.prot)){
                         all.ref.words<-union(all.ref.words, queries.ref[protein])
                         blast.word.universe<-union(blast.word.universe, wordUniverse(protein, queries.sssr))
                       }
                       #wordSet prot-scriber for univ.words parameter
                       ps.words <- wordSet(row[,3], blacklist.regexs = NULL)
                       univ.words<-unlist(strsplit(union(ps.words, blast.word.universe), split=' '))
                       
                       #Merged reference words for each family
                       row[,5]=paste(unique(unlist(all.ref.words)), collapse = ' ')
                       #wordSet prot-scriber HRD with blacklist filter
                       row[,4]=paste(wordSet(row[,3]), collapse = ' ')
                       #Blast word universe
                       row[,6]=paste(unique(unlist(blast.word.universe)), collapse = ' ')
                       #univ.words 
                       row[,7]=paste(unique(unlist(univ.words)), collapse = ' ')
                       
                       #intermediate results: 
                       #seq-fam id
                       #HRD prot-scriber
                       #wordSet filtered HRD prot-scriber
                       #reference words
                       #blast word universe
                       #univ.words
                       #mcc.simple(empty)
                       columns<-row[,c(1,3,4,5,6,7,8)]
                     })

wordSetTable<-do.call(rbind, wordSetHRD)


#'calculate Matthews Correlation Coefficient

wordSetHRD_rowList<-split(wordSetTable, 1:nrow(wordSetTable))
mccDf<-mclapply(wordSetHRD_rowList, function(row) {
  
  #wordSet prot-scriber HRD
  pred<-unlist(strsplit(row[,3], split=" "))
  #merged reference words
  ref<-unlist(strsplit(row[,4], split=" "))
  
  univ.words<-unlist(strsplit(row[,6], split=" "))
  
  mcc.simple <- mcc(pred, ref, univ.words)
  #fill in mcc.simple column
  row[,7]=paste(mcc.simple, collapse = ' ')
  columns<-row[,c(1,2,3,4,5,6,7)]
})
mccTable<-do.call(rbind, mccDf)
mccTable<-subset(mccTable[mccTable$Merged.Ref.Family.Description!="" & mccTable$wordSet_hrd!="", ])
#write.table(mccTable,"preliminary_results.txt",
#            sep = "\t", row.names = FALSE, quote = TRUE)



#'calculate FScore

mccTable_rowList<-split(mccTable, 1:nrow(mccTable))
print("testing")
fscoreDf<-mclapply(mccTable_rowList,
                   function(row) {
                     pred<-unlist(strsplit(row[,3], split=" "))
                     ref<-unlist(strsplit(row[,4], split=" "))
                     
                     fscore<-fScoreCalculator(pred, ref)
                     columns<-fscore[,1:4]
                   })
fscoreTable<-do.call(rbind, fscoreDf)
fscoreTable<-cbind(mccTable, fscoreTable[,1:4]) #c(1,2,3,4)])

#' Save result table:
write.table(fscoreTable, "EndResults.txt", 
            sep = "\t", row.names = FALSE, quote = TRUE)

################################################################                            


#' DONE
message("DONE")
