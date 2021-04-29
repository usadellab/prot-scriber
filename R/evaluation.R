#' Computes the F-Score of a prediction given the truth (reference). For
#' details see Wikipedia, here:
#' https://en.wikipedia.org/wiki/F-score
#' 
#' @param pred - A vector of predictions
#' @param ref - A vector holding the truth (references)
#' @param prot.id - A string identifying the query protein for which the
#' prediction has been made.
#' @param beta - recall is considered beta times as important as precision.
#' Default is getOption("fScore.beta", 0.25) thus emphazising precision four times
#' more than recall. 
#'
#' @return An instance of `base::data.frame` with exactly one row and the
#' following columns: Protein.ID, precision, recall, F.Score, beta
#' @export
fScore <- function(pred, ref, prot.id, beta = getOption("fScore.beta", 
    0.25)) {
    true.pos <- intersect(pred, ref)
    precision <- length(true.pos)/length(pred)
    recall <- length(true.pos)/length(ref)
    f.score <- if ((beta * precision + recall) == 0) {
        0
    } else {
        (1 + beta^2) * (precision * recall)/(beta * precision + 
            recall)
    }
    data.frame(Protein.ID = prot.id, precision = precision, recall = recall, 
        F.Score = f.score, fScore.beta = beta, stringsAsFactors = FALSE)
}

#' Function tests its namesake `fScore`.
#'
#' @return TRUE if and only if all tests pass.
#' @export
testFScore <- function() {
    beta <- 1
    p1 <- c("Hello", "World")
    r1 <- p1
    t1 <- identical(fScore(p1, r1, "P1", beta = beta)$F.Score, 
        1)
    p2 <- c("Hello", "World")
    r2 <- c("Not", "a", "single", "word")
    t2 <- identical(fScore(p2, r2, "P2", beta = beta)$F.Score, 
        0)
    p3 <- c("Hello", "World", "Hola", "Mundo")
    r3 <- c("Some", "words", "Hello", "World")
    t3 <- identical(fScore(p3, r3, "P3", beta = beta)$F.Score, 
        0.5)
    all(t1, t2, t3)
}
