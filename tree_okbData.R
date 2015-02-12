#' Function to plot phylogenies showing frequencies of the unique sequences
#'
#' @param x an \linkS4class{obkData} object.
#' @param ncatmax an integer giving the number of categories of haplotype frequency
#' @param colfun a function to define the colours
#' @author Emmanuel Paradis

tree_okbData <- function(x, ncatmax = 10, colfun = topo.colors)
{
    require(pegas)
    DNA <- x@dna@dna
    LOCI <- names(DNA)
    cat("okbData object has", length(LOCI), "genetic data set(s)\n")
    for (i in seq_along(LOCI)) {
        cat("Processing data from", LOCI[i], "\n")

        BF <- base.freq(dna, TRUE, TRUE)
        nb <- sum(BF)
        nambi <- sum(BF[-(1:4)])/nb
        warning(paste0(round(nambi, 3), "% of ambiguous or unknown sites (or alignment gaps) out of ", nb))

        h <- haplotype(dna)
        nh <- nrow(h)
        attr(h, "index") <- lapply(attr(h, "index"), function(x) rownames(dna)[x])
        names(attr(h, "index")) <- rownames(h) <- paste0("uniqseqID", seq_len(nh))
        phy <- nj(dist.dna(h, "N", pairwise.deletion = TRUE))

        hfreq <- sapply(attr(h, "index"), length)
        rangehfreq <- range(hfreq)
        if (rangehfreq[2] < 10) ncatmax <- rangehfreq[2]
        sq <- seq(rangehfreq[1] - 1, rangehfreq[2], length.out = ncatmax)
        cat <- cut(hfreq, sq)
        co <- colfun(ncatmax)
        n <- Ntip(phy)
        cat("Type enter to continue\n")
        readLines(n = 1)
        plot(phy, show.tip.label = FALSE, x.lim = sum(phy$edge.length)/Nedge(phy) + max(hfreq))
        phydataplot(hfreq, phy, "b")
        tiplabels(text = "   ", bg = co[cat], adj = 0)
        title(sub = LOCI[i])

        psr <- par("usr")
        xx <- psr[2]/2
        yy <- psr[4] * (0.5 + 0.5/par("plt")[4])
        legend(xx, yy, legend = round(sq[-1], 1), pch = 22, pt.bg = co,
               pt.cex = 2, bty = "n", xjust = 0.5, yjust = 0.5,
               horiz = TRUE, xpd = TRUE)
    }
}

