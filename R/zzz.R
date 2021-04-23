# Function ensured data-sets to be loaded, whenever this package is loaded.
.onLoad <- function( libname = find.package( "prot.scriber" ), pkgname = "prot.scriber" ) {
    data( "regexLists", package = "prot.scriber" )
    data( "p_coccineus_seq_sin_search", package = "prot.scriber" )
}
