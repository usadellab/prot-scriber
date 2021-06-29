#!/usr/bin/env Rscript

require(prot.scriber)

#' Load data:
data("p_coccineus_HRDs")
data("p_coccineus_HRDs_alignment_regios")


#' Parse command line options:
option_list <- list(
  make_option(
    c("-p", "--plot-dir"),
    type = "character",
    default = NULL,
    help = "The directory into which to write the resulting plots.",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
script.args <- parse_args(opt_parser)


#' Factorize the applied Methods:
pc.hrds$Method <- factor(pc.hrds$Method)


#' Plot F-Scores:
g.p <- ggplot(data = pc.hrds, aes(x = Method, y = F.Score, fill = Method, 
    color = Method)) + geom_point(color = "grey", position = position_jitter(width = 0.1), 
    size = 0.1, alpha = 0.8) + geom_boxplot(outlier.shape = NA, 
    alpha = 0.5, size = 1) + coord_flip() + theme_bw() + theme(legend.position = "none")
ggsave(plot = g.p, filename = file.path(script.args[["plot-dir"]], 
    "p_coccineus_eval_f-scores_dists.pdf"))


g.p <- ggplot(data = pc.hrds, aes(x = Method, y = F.Score.relative, 
    fill = Method, color = Method)) + geom_point(color = "grey", 
    position = position_jitter(width = 0.2), size = 0.1, alpha = 0.8) + 
    geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 1) + 
    coord_flip() + theme_bw() + theme(legend.position = "none")
ggsave(plot = g.p, filename = file.path(script.args[["plot-dir"]], 
    "p_coccineus_eval_relative_f-scores_dists.pdf"))


#' Plot MCC:
g.p <- ggplot(data = pc.hrds, aes(x = Method, y = MCC, fill = Method, 
    color = Method)) + geom_point(color = "grey", position = position_jitter(width = 0.1), 
    size = 0.1, alpha = 0.8) + geom_boxplot(outlier.shape = NA, 
    alpha = 0.5, size = 1) + coord_flip() + theme_bw() + theme(legend.position = "none")
ggsave(plot = g.p, filename = file.path(script.args[["plot-dir"]], 
    "p_coccineus_eval_MCC_dists.pdf"))


g.p <- ggplot(data = pc.hrds, aes(x = Method, y = MCC.relative, 
    fill = Method, color = Method)) + geom_point(color = "grey", 
    position = position_jitter(width = 0.1), size = 0.1, alpha = 0.8) + 
    geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 1) + 
    coord_flip() + theme_bw() + theme(legend.position = "none")
ggsave(plot = g.p, filename = file.path(script.args[["plot-dir"]], 
    "p_coccineus_eval_relative_MCC_dists.pdf"))


#' Plot differences in performance scores when concatonating HRDs generated for
#' disjoint alignment-regions of sequence similarity search results:
g.p <- ggplot(data = pc.hrds.alignment.regions.score.diffs.df, 
    aes(x = Method, y = F.Score.diff, fill = Method, color = Method)) + 
    geom_point(color = "grey", position = position_jitter(width = 0.1), 
        size = 0.1, alpha = 0.8) + geom_boxplot(outlier.shape = NA, 
    alpha = 0.5, size = 1) + coord_flip() + theme_bw() + theme(legend.position = "none")
ggsave(plot = g.p, filename = file.path(script.args[["plot-dir"]], 
    "p_coccineus_multi_region_HRDs_F_Score_diffs.pdf"))


#' Generate a HTML document of the Human Readable Descriptions assigned to the
#' P. coccineus query proteins by the various methods:
pc.hrds.html <- file.path(script.args[["plot-dir"]], "p_coccineus_HRDs.html")
visualizeHrdResults(pc.hrds, pc.ref, pc.hrds.html)
#' Compress the document, as it is rather large:
system(paste("bzip2", pc.hrds.html))


#' Generate a HTML document as above but for the multi region HRDs:
pc.hrds.mult.reg.html <- file.path(script.args[["plot-dir"]], 
    "p_coccineus_HRDs_multi_region.html")
visualizeHrdResults(pc.hrds.alignment.regions.df, pc.ref, pc.hrds.mult.reg.html)
#' Compress the document, as it is rather large:
system(paste("bzip2", pc.hrds.mult.reg.html))


#' The end
message("DONE")
