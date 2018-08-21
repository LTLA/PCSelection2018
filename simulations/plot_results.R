# Runs through the *.rds files in 'results/' and creates pictures.
# Each plot represents the average results of one simulation scenario.

dir.create("pics", showWarnings=FALSE)
all.files <- list.files("results", pattern="_scenarios.txt$", full=TRUE)
method.alias <- c(elbow="Elbow", gavish="GD", marchenko="MP", denoised="Summation", jackstraw="Jackstraw", parallel="Parallel")

################################################

collected.x <- collected.y <- collected.i <- NULL
counter <- 1L

for (res in all.files) { 
    prefix <- sub("_[^_]+.txt$", "", res)
    current.scen <- read.table(paste0(prefix, "_scenarios.txt"), header=TRUE)
    current.mse <- read.table(paste0(prefix, "_mse.txt"), header=TRUE)
    current.num <- read.table(paste0(prefix, "_numbers.txt"), header=TRUE)

    everything.else <- setdiff(colnames(current.num), "optimal")
    x.vals <- current.num[,everything.else] - current.num$optimal
    y.vals <- current.mse[,everything.else] / current.mse$optimal

    collected.y[[counter]] <- y.vals
    collected.x[[counter]] <- x.vals
    decode.name <- strsplit(basename(prefix), "-")[[1]]
    names(decode.name) <- c("Mode", "Npop", "Bio", "Noise")
    collected.i[[counter]] <- suppressWarnings(data.frame(rbind(decode.name), current.scen))
    counter <- counter + 1L

    current.scen <- lapply(current.scen, factor)
    pch <- ifelse(current.scen$Ngenes==1000, 1, 4)
    col <- viridis::plasma(3)[current.scen$Prop.DE]
    cex <- as.integer(current.scen$Ncells)

    pdf(file.path("pics", paste0(basename(prefix), ".pdf")))
    for (s in everything.else) {
        plot(x.vals[,s], y.vals[,s], xlim=range(x.vals), ylim=range(y.vals), log="y",
            xlab="Difference in PCs from optimal",
            ylab="Fold-change in MSE from optimal",
            main=method.alias[s], pch=pch, col=col, cex=cex)

        legend("bottomright", 
            col=c(viridis::plasma(3), rep("black", 4)),
            pch=c(16, 16, 16, 1, 4, 16, 16),
            pt.cex=c(1, 1, 1, 1, 1, 1, 2),
            legend=c(paste("Prop =", levels(current.scen$Prop.DE)),
                paste("Genes =", levels(current.scen$Ngenes)),
                paste("Cells =", levels(current.scen$Ncells))
            )
        )
    }
    dev.off()
}

################################################

all.mse <- do.call(rbind, collected.y)
all.num <- do.call(rbind, collected.x)
all.info <- do.call(rbind, collected.i)

FUN <- function(MSE, main=NULL) {
    MSE <- log2(MSE)
    breaks <- seq(0, max(MSE), length.out=50)
    zero <- colSums(MSE==0)

    out <- lapply(MSE, function(v) hist(v[v!=0], breaks=breaks, plot=FALSE))
    upper.y <- max(unlist(lapply(out, FUN=function(h) max(h$counts)))) * 1.1
    left.x <- -max(MSE) / 8

    par(mar=c(5.1, 5.2, 4.1, 2.0))
    plot(0,0,type="n", xlim=c(left.x, max(breaks)), ylim=c(0, upper.y*length(out)), yaxt="n", ylab="", xaxt="n", xlab="Fold-increase from optimal MSE", main=main)
    pbreaks <- pretty(breaks)
    axis(1, at=pbreaks, label=2^pbreaks)

    shift <- 0
    for (x in names(out)) {
        polygon(rep(out[[x]]$breaks, each=2), c(shift, rep(out[[x]]$counts, each=2) + shift, shift), col="grey80")
        prop <- zero[[x]] / nrow(MSE)
        rect(left.x * 0.8, shift, left.x * 0.5, shift + upper.y * 0.9 * prop, col="grey30")
        text(left.x * 0.9, shift + upper.y * 0.05, sprintf("%i%%", round(prop*100)), srt=90, pos=3)
        mtext(side=2, line=1, at=shift+ upper.y/3, text=method.alias[x], las=1)
        shift <- shift + upper.y
    }
}

# Creating a summary of all MSEs.

pdf("pics/overall.pdf")
FUN(all.mse)
dev.off()

# Summarizing by each metric in 'info'.
for (i in names(all.info)) {
    current <- factor(all.info[[i]])

    if (nlevels(current)==3) {
        pdf(sprintf("pics/%s.pdf", i), width=5, height=7)
    } else {
        pdf(sprintf("pics/%s.pdf", i))
    }

    for (lev in levels(current)) { 
        keep <- current==lev

        if (i=="Bio") {
            lev <- substitute(s^2~"="~X, list(X=lev))
        } else if (i=="Mode") {
            lev <- c(gaussclust="Cluster", trajectory="Trajectory")[lev]
        } else if (i=="Prop.DE") {
            lev <- substitute(P~"="~X, list(X=lev))
        } 

        FUN(all.mse[keep,], main=lev)
    }
    dev.off()
}

