# Runs through the *.rds files in 'results/' and creates pictures.
# Each plot represents the average results of one real data set. 

dir.create("pics", showWarnings=FALSE)
all.files <- list.files("results", pattern="_scenarios.txt$", full=TRUE)
decoder <- c('kolod'='mESC', 'pbmc4k'='PBMC')

for (res in all.files) { 
    prefix <- sub("_[^_]+.txt$", "", res)
    current.scen <- read.table(paste0(prefix, "_scenarios.txt"), header=TRUE)
    current.mse <- read.table(paste0(prefix, "_mse.txt"), header=TRUE)
    current.num <- read.table(paste0(prefix, "_numbers.txt"), header=TRUE)

    everything.else <- setdiff(colnames(current.num), "optimal")
    x.vals <- current.num[,everything.else] - current.num$optimal
    y.vals <- current.mse[,everything.else] / current.mse$optimal

    out <- split(y.vals, current.scen[,1])
    means <- do.call(rbind, lapply(out, colMeans))
    sds <- do.call(rbind, lapply(out, FUN=function(x) apply(x, 2, sd)/sqrt(nrow(x))))

    basefix <- basename(prefix)
    pdf(file.path("pics", paste0(basefix, ".pdf")))
    X <- barplot(means, beside=TRUE, ylab="Fold increase from optimal MSE", main=decoder[[basefix]])
    upper <- means + sds
    segments(X, means, X, upper)
    segments(X-0.1, upper, X+0.1, upper)
    legend("topright", title="rank", legend=rownames(means), fill=grey.colors(nrow(means)))
    dev.off()
}


