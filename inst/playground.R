require(prot.scriber)
options(mc.cores = detectCores() - 1)

data("p_coccineus_seq_sim_search")
data("p_coccineus_reference_words")

#' OVERLAP
#' currently not in use
overlap <- function(qend, qstart, send, sstart, qlen, slen) {
    ((qend - qstart + 1) + (send - sstart + 1))/(qlen + slen)
}
overlapFromDf <- function(sssr.df) {
    with(sssr.df, overlap(qend, qstart, send, sstart, qlen, slen))
}
#' END of OVERLAP definitions

q.1 <- pc.sprot$qseqid[[1]]
sssr <- list(Swissprot = pc.sprot, trEMBL = pc.trembl)
q.1.phrases <- phrasesForQuery(q.1, sssr)
q.1.phrases.stats <- statsOfPhrasesForQuery(q.1, q.1.phrases, 
    pc.ref)

#' Of all possible phrase, which'd be the best annotation?
max.F.Score <- max(q.1.phrases.stats$F.Score)
q.1.best.phrase <- q.1.phrases.stats[with(q.1.phrases.stats, 
    which(F.Score == max.F.Score)), "Phrase"]

#' Best performing prot.scriber phrases:
polynomial.best.phrase.i <- which(with(q.1.phrases.stats, polynomial.word.scores == 
    max(polynomial.word.scores)))
polynomial.best.phrase.max.F.Score <- max(q.1.phrases.stats[polynomial.best.phrase.i, 
    "F.Score"])
inverse.inf.cont.best.phrase.i <- which(with(q.1.phrases.stats, 
    centered.inverse.inf.cntnt == max(centered.inverse.inf.cntnt)))
inverse.inf.cont.best.phrase.max.F.Score <- max(q.1.phrases.stats[inverse.inf.cont.best.phrase.i, 
    "F.Score"])
frequencies.best.phrase.i <- which(with(q.1.phrases.stats, centered.frequencies == 
    max(centered.frequencies)))
frequencies.best.phrase.max.F.Score <- max(q.1.phrases.stats[frequencies.best.phrase.i, 
    "F.Score"])


#' Best Blast:
q.1.sprot <- pc.sprot[which(pc.sprot$qseqid == q.1), ]
max.bitscore.sprot <- max(q.1.sprot$bitscore)
bb.sprot.sseqid <- q.1.sprot[which(q.1.sprot$bitscore == max.bitscore.sprot), 
    ]$sseqid
bb.sprot.hrd <- q.1.sprot[which(q.1.sprot$bitscore == max.bitscore.sprot), 
    ]$HRD
bb.sprot.hrd.fscore <- fScore(tolower(strsplit(bb.sprot.hrd, 
    "\\s+")[[1]]), pc.ref[[tolower(q.1)]], q.1)
q.1.trembl <- pc.trembl[which(pc.trembl$qseqid == q.1), ]
max.bitscore.trembl <- max(q.1.trembl$bitscore)
bb.trembl.sseqid <- q.1.trembl[which(q.1.trembl$bitscore == max.bitscore.trembl), 
    ]$sseqid
bb.trembl.hrd <- q.1.trembl[which(q.1.trembl$bitscore == max.bitscore.trembl), 
    ]$HRD
bb.trembl.hrd.fscore <- fScore(tolower(strsplit(bb.trembl.hrd, 
    "\\s+")[[1]]), pc.ref[[tolower(q.1)]], q.1)


#' Test another query 'Pc53_1'
q.2 <- 'Pc53_1'
sssr <- list(Swissprot = pc.sprot, trEMBL = pc.trembl)
q.2.phrases <- phrasesForQuery(q.2, sssr)
q.2.phrases.stats <- statsOfPhrasesForQuery(q.2, q.2.phrases, 
    pc.ref)


#' Test time requirement to annotate P.coccineus:
sssr.orig <- sssr
sssr <- setNames(lapply(sssr.orig, function(x) x[1:1000,]), names(sssr))
t.start <- Sys.time()
pc.annos <- annotateProteinsAndEvaluatePerformance(sssr, pc.ref)
t.duration <- as.numeric(Sys.time()) * 1000 - as.numeric(t.start) * 1000
sssr <- sssr.orig
