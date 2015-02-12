#' @param x an \linkS4class{obkData} object.
#' @author Emmanuel Paradis

alignment_explorer <- function(x)
{
    require(pegas)

    DNA <- x@dna@dna
    LOCI <- names(DNA)
    cat("okbData object has", length(LOCI), "genetic data set(s)\n")
    for (i in seq_along(LOCI)) {
        cat("Processing data from", LOCI[i], "\n")
        y <- DNA[[i]]
        s <- seg.sites(y)
        ## op <- options(warn = -1)
        ss <- site.spectrum(y)
        ## msg <- names(last.warning)
        ## options(op)
        cat("Type enter to continue\n")
        readLines(n = 1)
        layout(matrix(c(1, 2, 1, 3), 2, 2))
        image(y, xlab = LOCI[i])
        rug(s, -0.05, 3, 1)
        mtext("Tick marks show segregating sites", at = 0, adj = 0, line = 2, font = 3)
        plot(ss)
        ## legend("topright", msg, bty = "n")
    }
}

## load("ToyOutbreak.RData")
## alignment_explorer(ToyOutbreak)
