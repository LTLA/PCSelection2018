# Setting up a function that decides how many PCs to remove.

computeMSE <- function(svd.out, center, truth, ncomponents=100) {
    running <- 0
    collected <- numeric(ncomponents)
    for (i in seq_along(collected)) { 
        running <- running + svd.out$u[,i,drop=FALSE] %*% (svd.out$d[i] * t(svd.out$v[,i,drop=FALSE]))
        collected[i] <- mean((t(running) + center - truth)^2)
    }
    return(collected)
}

chooseNumber <- function(observed, truth, max.rank=50)
# Assessing each strategy to choose the number of PCs.
{ 
    center <- rowMeans(observed)
    SVD <- svd(t(observed - center), nu=max.rank, nv=max.rank)
    prog.var <- SVD$d^2 / (ncol(observed) - 1) 
    total.var <- sum(prog.var)

    # Using our denoising approach.
    tech.comp <- apply(observed - truth, 1, var)
    tech.var <- sum(tech.comp)
    denoised <- scran:::.get_npcs_to_keep(prog.var, tech.var)
    
    # Applying parallel analysis.
    require(scran)
    parallel <- parallelPCA(observed, BPPARAM=MulticoreParam(3), value="n", threshold=0.05, approximate=TRUE, min.rank=1, max.rank=max.rank)

    # Using the Marchenko-Pastur limit on the _eigenvalues_.
    # (See https://www.wolfram.com/language/11/random-matrices/marchenko-pastur-distribution.html?product=mathematica)
    # Multiplying the limit by the mean technical component of the variance to account for scaling.
    library(RMTstat)
    ndf <- nrow(observed) - 1
    limit <- qmp(1, ndf=ndf, pdim=ncol(observed)) * mean(tech.comp)
    marchenko <- sum(SVD$d^2/ndf > limit)

    # Applying the Gavish-Donoho method, setting the noise 'sigma' at the square root of the mean technical component.
    m <- min(dim(observed))
    n <- max(dim(observed))
    beta <- m/n
    lambda <- sqrt( 2 * (beta + 1) + (8 * beta) / ( beta + 1 + sqrt(beta^2 + 14 * beta + 1) ) )
    gv <- sum(SVD$d > lambda * sqrt(n) * sqrt(mean(tech.comp)))

    # Detecting the elbow in a scree plot, based on distance from the line.
    v2last <- c(max.rank - 1L, prog.var[max.rank] - prog.var[1])
    v2last <- v2last/sqrt(sum(v2last^2))
    v2other <- rbind(seq_len(max.rank) - 1L, prog.var[1:max.rank] - prog.var[1])
    dist2point <- sqrt(colSums((v2other - outer(v2last, colSums(v2last * v2other)))^2))
    elbow <- which.max(dist2point)

    # Using Seurat's Jackstraw method, keeping up to the first PC with zero genes at a FDR of 5%.
    require(Seurat)
    rownames(observed) <- paste0("Gene", seq_len(nrow(observed)))
    colnames(observed) <- paste0("Cell", seq_len(ncol(observed)))
    Seu <- CreateSeuratObject(observed)
    Seu@scale.data <- observed
    Seu <- RunPCA(Seu, pc.genes=rownames(observed))
    Seu <- JackStraw(Seu)
    nsig <- colSums(apply(Seu@dr$pca@jackstraw@emperical.p.value, 2, p.adjust, method="BH") <= 0.05)
    jackstraw <- min(c(max.rank, which(nsig==0)-1L))

    # Determining the MSE at each number of components. 
    mse <- computeMSE(SVD, center, truth, ncomponents=max.rank)
    optimal <- which.min(mse)

    num.pcs <- c(elbow=elbow, parallel=parallel, marchenko=marchenko, gavish=gv, jackstraw=jackstraw, denoised=denoised, optimal=optimal)
    cur.mse <- mse[num.pcs]
    names(cur.mse) <- names(num.pcs)
    return(list(number=num.pcs, mse=cur.mse))
}

runSimulation <- function(prefix, truth.FUN, iters=10, observed.FUN=NULL)
# A convenience function to run simulations based on a function that generates a matrix of true signal.
{
    scn.fname <- paste0(prefix, "_scenarios.txt")
    npc.fname <- paste0(prefix, "_numbers.txt")
    mse.fname <- paste0(prefix, "_mse.txt")
    counter <- 1L

    if (is.null(observed.FUN)) {
        observed.FUN <- function(truth) {
            truth + rnorm(length(truth))
        }
    }

    for (ncells in c(200, 1000)) {
        for (ngenes in c(1000, 5000)) {
            for (affected in c(0.2, 0.5, 1)) { 
                cur.mse <- cur.retained <- NULL 

                for (it in seq_len(iters)) {
                    truth <- truth.FUN(ngenes*affected, ncells)
                    truth <- rbind(truth, matrix(0, ncol=ncells, nrow=(1-affected)*ngenes))
                    y <- observed.FUN(truth)
                    out <- chooseNumber(y, truth)

                    is.first <- counter==1L && it==1L
                    write.table(data.frame(Ncells=ncells, Ngenes=ngenes, Prop.DE=affected),
                        file=scn.fname, append=!is.first, col.names=is.first, row.names=FALSE, quote=FALSE, sep="\t")
                    write.table(data.frame(rbind(out$number)), file=npc.fname, append=!is.first, col.names=is.first, row.names=FALSE, quote=FALSE, sep="\t")
                    write.table(data.frame(rbind(out$mse)), file=mse.fname, append=!is.first, col.names=is.first, row.names=FALSE, quote=FALSE, sep="\t")
                }

                counter <- counter + 1L
            }
        }
    }

    return(NULL)
}

addNoise <- function(variance) 
# Adding noise with different strategies, with different variation of technical noise across genes.
{
    if (variance=="none") {
        function(truth) truth + rnorm(length(truth))
    } else if (variance=="moderate") {
        function(truth) {
            sd <- sqrt(rgamma(nrow(truth), 2, 2)) 
            truth + rnorm(length(truth), sd=sd)
        }
    } else if (variance=="high") {
        function(truth) {
            sd <- sqrt(runif(nrow(truth), 0, 6))
            truth + rnorm(length(truth), sd=sd)
        }
    } else {
        stop("unknown variance mode")
    }
}
