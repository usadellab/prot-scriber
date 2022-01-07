#!/usr/bin/env Rscript

require(prot.scriber)

#' Load data:
data("faba_HRDs")
data("faba_reference_words")


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
faba.hrds$Method <- factor(faba.hrds$Method)


#' Plot F-Scores:
g.p <- ggplot(data = faba.hrds, aes(x = Method, y = F.Score, fill = Method, 
    color = Method)) + geom_point(color = "grey", position = position_jitter(width = 0.1), 
    size = 0.1, alpha = 0.8) + geom_boxplot(outlier.shape = NA, 
    alpha = 0.5, size = 1) + coord_flip() + theme_bw() + theme(legend.position = "none")
ggsave(plot = g.p, filename = file.path(script.args[["plot-dir"]], 
    "faba_eval_f-scores_dists.pdf"))


g.p <- ggplot(data = faba.hrds, aes(x = Method, y = F.Score.relative, 
    fill = Method, color = Method)) + geom_point(color = "grey", 
    position = position_jitter(width = 0.2), size = 0.1, alpha = 0.8) + 
    geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 1) + 
    coord_flip() + theme_bw() + theme(legend.position = "none")
ggsave(plot = g.p, filename = file.path(script.args[["plot-dir"]], 
    "faba_eval_relative_f-scores_dists.pdf"))


#' Plot MCC:
g.p <- ggplot(data = faba.hrds, aes(x = Method, y = MCC, fill = Method, 
    color = Method)) + geom_point(color = "grey", position = position_jitter(width = 0.1), 
    size = 0.1, alpha = 0.8) + geom_boxplot(outlier.shape = NA, 
    alpha = 0.5, size = 1) + coord_flip() + theme_bw() + theme(legend.position = "none")
ggsave(plot = g.p, filename = file.path(script.args[["plot-dir"]], 
    "faba_eval_MCC_dists.pdf"))


g.p <- ggplot(data = faba.hrds, aes(x = Method, y = MCC.relative, 
    fill = Method, color = Method)) + geom_point(color = "grey", 
    position = position_jitter(width = 0.1), size = 0.1, alpha = 0.8) + 
    geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 1) + 
    coord_flip() + theme_bw() + theme(legend.position = "none")
ggsave(plot = g.p, filename = file.path(script.args[["plot-dir"]], 
    "faba_eval_relative_MCC_dists.pdf"))


#' Generate a HTML document of the Human Readable Descriptions assigned to the
#' Faba query proteins by the various methods:
faba.hrds.html <- file.path(script.args[["plot-dir"]], "faba_HRDs.html")
visualizeHrdResults(faba.hrds, faba.ref, faba.hrds.html)
#' Compress the document, as it is rather large:
system(paste("rm -f", paste0(faba.hrds.html, ".bz2"), "&& bzip2", 
    faba.hrds.html))


#' The end
message("DONE")
