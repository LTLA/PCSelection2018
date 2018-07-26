# This explores the effect of removing technical and biological noise.

set.seed(1010101)

#' PCA with the Gavish-Donoho method to remove technical noise.
GV <- function(observed, SVD, SD) {
    m <- min(dim(observed))
    n <- max(dim(observed))
    beta <- m/n
    lambda <- sqrt( 2 * (beta + 1) + (8 * beta) / ( beta + 1 + sqrt(beta^2 + 14 * beta + 1) ) )
    gv <- sum(SVD$d > lambda * sqrt(n) * SD)
    max(1L, gv)
}

#' Comparing the number of neighbors from the same or different population.
library(cluster)
COMPARATOR <- function(mat, subpop) {
    dmat <- as.matrix(dist(mat))
    silhouette(dist=dmat, subpop)[,3]
}

# Creating the simulation.
ngenes <- 1000
subpop.means <- cbind(integer(ngenes), rep(0.5, ngenes))

ncells <- 100    
subpop.ids <- rep(1:2, each=ncells/2)
true.means <- subpop.means[,subpop.ids]

bio.sd <- 1
with.bio.noise <- true.means + rnorm(length(true.means), sd=bio.sd)
tech.sd <- 1
final <- with.bio.noise + rnorm(length(with.bio.noise), sd=tech.sd)

# Choosing the number of PCs to retain.
SVD <- svd(t(final - rowMeans(final)))
PCs <- sweep(SVD$u, 2, SVD$d, "*")

remove.tech <- GV(final, SVD, tech.sd)
remove.all <- GV(final, SVD, sqrt(bio.sd^2 + tech.sd^2))

# Computing distances between cells of the same and other subpopulation.
all.d <- COMPARATOR(PCs[, seq_len(remove.all), drop=FALSE], subpop.ids)    
tech.d <- COMPARATOR(PCs[, seq_len(remove.tech), drop=FALSE], subpop.ids)
true.d <- COMPARATOR(t(with.bio.noise), subpop.ids)

# Creating a plot of the PCs.
pdf("pics/bio_noise.pdf")
prog.exp <- SVD$d^2/(ncells-1)
prog.exp <- prog.exp/sum(prog.exp) * 100 
plot(PCs[,1], PCs[,2], 
    xlab=sprintf("PC1 (%.1f%%)", prog.exp[1]),
    ylab=sprintf("PC2 (%.1f%%)", prog.exp[2]),
    pch=c(4, 21)[subpop.ids])

# Creating the boxplot.
boxplot(list(All=all.d, Technical=tech.d, True=true.d), range=0, 
    ylab="Average silhouette width", col="grey80")
dev.off()
