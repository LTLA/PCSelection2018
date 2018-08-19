#' Central function for SVD.
decomposer <- function(y, nrank, approximate=FALSE) {
    center <- rowMeans(y)
    if (approximate) {
        SVD <- irlba::irlba(t(y), center=center, nv=nrank)
    } else {
        SVD <- svd(as.matrix(t(y - center)), nu=nrank, nv=nrank)
        SVD$d <- SVD$d[seq_len(nrank)]            
    }
    list(SVD=SVD, center=center)
}

#' Computing the MSE for any number of components for a given SVD.
computeMSE <- function(svd.out, center, truth, ncomponents=100) {
    running <- 0
    collected <- numeric(ncomponents)
    for (i in seq_along(collected)) { 
        running <- running + svd.out$u[,i,drop=FALSE] %*% (svd.out$d[i] * t(svd.out$v[,i,drop=FALSE]))
        collected[i] <- mean((t(running) + center - truth)^2)
        gc()
    }
    return(collected)
}

#' Choosing the number of PCs via a variety of strategies.
chooseNumber <- function(observed, SVD, tech.comp) { 
    # Using our denoising approach.
    prog.var <- SVD$d^2 / (ncol(observed) - 1) 
    tech.var <- sum(tech.comp)
    denoised <- scran:::.get_npcs_to_keep(prog.var, tech.var, total=sum(DelayedMatrixStats::rowVars(DelayedArray::DelayedArray(observed))))
    
    # Applying parallel analysis.
    require(scran)
    max.rank <- ncol(SVD$v)
    parallel <- parallelPCA(observed, BPPARAM=MulticoreParam(3), value="n", threshold=0.05, approximate=TRUE, min.rank=1, max.rank=max.rank)
    gc()

    # Using the Marchenko-Pastur limit on the _eigenvalues_.
    # (See https://www.wolfram.com/language/11/random-matrices/marchenko-pastur-distribution.html?product=mathematica)
    # Multiplying the limit by the mean technical component of the variance to account for scaling.
    library(RMTstat)
    ndf <- nrow(observed) - 1
    limit <- qmp(1, ndf=ndf, pdim=ncol(observed)) * mean(tech.comp)
    marchenko <- sum(SVD$d^2/ndf > limit)
    marchenko <- max(1L, marchenko)

    # Applying the Gavish-Donoho method, setting the noise 'sigma' at the square root of the mean technical component.
    m <- min(dim(observed))
    n <- max(dim(observed))
    beta <- m/n
    lambda <- sqrt( 2 * (beta + 1) + (8 * beta) / ( beta + 1 + sqrt(beta^2 + 14 * beta + 1) ) )
    gv <- sum(SVD$d > lambda * sqrt(n) * sqrt(mean(tech.comp)))
    gv <- max(1L, gv)

    # Detecting the elbow in a scree plot, based on distance from the line.
    # We '-1' to get the element _before_ the elbow, which is what we actually want to keep.
    v2last <- c(max.rank - 1L, prog.var[max.rank] - prog.var[1])
    v2last <- v2last/sqrt(sum(v2last^2))
    v2other <- rbind(seq_along(prog.var) - 1L, prog.var - prog.var[1])
    dist2point <- sqrt(colSums((v2other - outer(v2last, colSums(v2last * v2other)))^2))
    elbow <- which.max(dist2point) - 1L

    # Using Seurat's Jackstraw method, keeping up to the first PC with zero genes at a FDR of 5%.
    require(Seurat)
    rownames(observed) <- paste0("Gene", seq_len(nrow(observed)))
    colnames(observed) <- paste0("Cell", seq_len(ncol(observed)))
    Seu <- CreateSeuratObject(observed)
    Seu <- ScaleData(Seu, do.scale=FALSE, display.progress=FALSE)
    Seu <- RunPCA(Seu, pc.genes=rownames(observed), seed.use=NULL, pcs.compute=max.rank)
    Seu <- JackStraw(Seu, display.progress=FALSE)

    nsig <- colSums(apply(Seu@dr$pca@jackstraw@emperical.p.value, 2, p.adjust, method="BH") <= 0.05)
    jackstraw <- min(c(max.rank, which(nsig==0)-1L))
    jackstraw <- max(1L, jackstraw)
    gc()

    return(c(elbow=elbow, parallel=parallel, marchenko=marchenko, gavish=gv, jackstraw=jackstraw, denoised=denoised))
}

#' Assessing each strategy to choose the number of PCs in simulations.
assessChoices <- function(observed, truth, max.rank=50, approximate=FALSE) {
    dec.out <- decomposer(observed, nrank=max.rank, approximate=approximate)
    SVD <- dec.out$SVD
    center <- dec.out$center

    # Making choices.
    tech.comp <- rowMeans((observed - truth)^2)
    choices <- chooseNumber(observed, SVD, tech.comp)

    # Determining the MSE at each number of components. 
    mse <- computeMSE(SVD, center, truth, ncomponents=max.rank)
    optimal <- which.min(mse)

    num.pcs <- c(choices, optimal=optimal)
    cur.mse <- mse[num.pcs]
    names(cur.mse) <- names(num.pcs)
    return(list(number=num.pcs, mse=cur.mse))
}

#' A convenience function to run simulations based on a function that generates a matrix of true signal.
runSimulation <- function(prefix, truth.FUN, iters=10, observed.FUN=NULL) {
    scn.fname <- paste0(prefix, "_scenarios.txt")
    npc.fname <- paste0(prefix, "_numbers.txt")
    mse.fname <- paste0(prefix, "_mse.txt")
    counter <- 1L

    if (is.null(observed.FUN)) {
        observed.FUN <- function(truth) {
            truth + rnorm(length(truth))
        }
    }

    for (ncells in c(200, 1000, 5000)) {
        for (ngenes in c(1000, 5000)) {
            for (affected in c(0.2, 0.5, 1)) { 
                cur.mse <- cur.retained <- NULL 

                for (it in seq_len(iters)) {
                    truth <- truth.FUN(ngenes*affected, ncells)
                    truth <- rbind(truth, matrix(0, ncol=ncells, nrow=(1-affected)*ngenes))
                    y <- observed.FUN(truth)
                    out <- assessChoices(y, truth, approximate=any(dim(y)>1000))

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

#' Adding noise with different strategies, with different variation of technical noise across genes.
addNoise <- function(variance) {
    if (variance=="none") {
        function(truth) truth + rnorm(length(truth))
    } else if (variance=="moderate") {
        function(truth) {
            sd <- sqrt(rgamma(nrow(truth), 2, 2)) 
            truth + rnorm(length(truth), sd=sd)
        }
    } else if (variance=="high") {
        function(truth) {
            sd <- sqrt(rgamma(nrow(truth), 0.2, 0.2)) 
            truth + rnorm(length(truth), sd=sd)
        }
    } else {
        stop("unknown variance mode")
    }
}

#' Generate a low-rank "truth" and flip the sign of the residuals to create the "observed" matrix.
generateReal <- function(original, SVD, center, noise, nrank=20) {
    i <- seq_len(nrank)
    recon <- t(SVD$u[,i,drop=FALSE] %*% (SVD$d[i] * t(SVD$v[,i,drop=FALSE]))) + center
    resim <- recon + rnorm(length(recon), sd=sqrt(noise))
    list(truth=recon, observed=resim) 
}

#' Execute PC selection methods on the simulated data based on real data.
simulateReal <- function(original, noise, prefix, iters=10) {
    sim.vals <- decomposer(original, approximate=TRUE, nrank=50)
    scn.fname <- paste0(prefix, "_scenarios.txt")
    npc.fname <- paste0(prefix, "_numbers.txt")
    mse.fname <- paste0(prefix, "_mse.txt")
    counter <- 1L

    for (recon.rank in c(10, 20, 30)) {
        for (it in seq_len(iters)) {
            sim <- generateReal(original, sim.vals$SVD, sim.vals$center, noise=noise, nrank=recon.rank)
            out <- assessChoices(sim$observed, sim$truth, approximate=TRUE)

            is.first <- counter==1L && it==1L
            write.table(data.frame(Rank=recon.rank), file=scn.fname, append=!is.first, col.names=is.first, row.names=FALSE, quote=FALSE, sep="\t")
            write.table(data.frame(rbind(out$number)), file=npc.fname, append=!is.first, col.names=is.first, row.names=FALSE, quote=FALSE, sep="\t")
            write.table(data.frame(rbind(out$mse)), file=mse.fname, append=!is.first, col.names=is.first, row.names=FALSE, quote=FALSE, sep="\t")

            gc()
        }
        counter <- counter + 1L
    }

    return(NULL)
}
