#' Function ensured data-sets to be loaded, whenever this package is loaded.
.onLoad <- function(libname = find.package("prot.scriber"), pkgname = "prot.scriber") {
    data("regexLists", package = "prot.scriber")
    data("p_coccineus_reference_words", package = "prot.scriber")
    message("To load data from prot.scriber that is not loaded automatically, please use e.g.:\n", 
        "data( \"p_coccineus_seq_sim_search\", package = \"prot.scriber\" )")
    message("Available datasets in package 'prot.scriber' are:\n", 
        paste(data(package = "prot.scriber")$results[, "Item"], 
            collapse = ", "))
}
