require(prot.scriber)

#' Load data:
data("p_coccineus_seq_sim_search")
data("p_coccineus_reference_words")
data("p_coccineus_alignment_regions")
data("p_coccineus_HRDs")

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
script.args <- parse_args(opt_parser)

#' Set mc.cores:
options(mc.cores = script.args$`n-cores`)


#' Find alignment regions for all P. coccineus query proteins:
pc.sssr <- list(Swissprot = pc.sprot, trEMBL = pc.trembl)
pc.hrds.alignment.regions <- annotateProteinsAndEvaluatePerformance.MultiRegion(pc.sssr, 
    pc.ref, pc.alignment.regions.df)


#' Join 'standard' Prot-Scriber HRDs with concatonated region HRDs for
#' comparison:
prot.ids <- unique(pc.hrds.alignment.regions$Protein.ID)
pc.hrds.alignment.regions.df <- rbind(pc.hrds.alignment.regions, 
    pc.hrds[which(pc.hrds$Protein.ID %in% prot.ids), ])
pc.hrds.alignment.regions.score.diffs.df <- do.call(rbind, mclapply(prot.ids, 
    function(qseqid) {
        qseqid.df <- pc.hrds.alignment.regions[which(pc.hrds.alignment.regions$Protein.ID == 
            qseqid), ]
        ps.methods <- unique(qseqid.df$Method)
        ps.methods.no.concat <- sub(".concat.regions", "", ps.methods, 
            fixed = TRUE)
        do.call(rbind, lapply(1:length(ps.methods), function(i) {
            ps.m.i <- ps.methods.no.concat[[i]]
            ps.m.conc.i <- ps.methods[[i]]
            pc.hrd.i <- pc.hrds[which(pc.hrds$Protein.ID == qseqid & 
                pc.hrds$Method == ps.m.i), ]
            pc.hrd.conc.i <- pc.hrds.alignment.regions[which(pc.hrds.alignment.regions$Protein.ID == 
                qseqid & pc.hrds.alignment.regions$Method == 
                ps.m.conc.i), ]
            qseqid.df$MCC.diff <- ((pc.hrd.i$MCC + 1) - (pc.hrd.conc.i$MCC + 1)) / 2
            qseqid.df$F.Score.diff <- pc.hrd.i$F.Score - pc.hrd.conc.i$F.Score
            qseqid.df
        }))
    }))


#' Save results:
save(pc.hrds.alignment.regions.df, pc.hrds.alignment.regions.score.diffs.df, 
    file = file.path(script.args[["data-dir"]], "p_coccineus_HRDs_alignment_regios.RData"))


#' The End
message("DONE")
