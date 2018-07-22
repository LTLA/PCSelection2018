library(scran)
library(scater)
library(DropletUtils)

# Pre-processing the 4kK PBMC dataset.
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)    
path.4k <- bfcrpath(bfc, "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz")
tmp.4k <- tempfile()
untar(path.4k, exdir=tmp.4k)
sce.4k <- read10xCounts(file.path(tmp.4k, "filtered_gene_bc_matrices/GRCh38/")) 

# Adding locational annotation (using a slightly off-version ensembl, but chromosome assignment shouldn't change).
library(EnsDb.Hsapiens.v86)
loc <- mapIds(EnsDb.Hsapiens.v86, keys=rownames(sce.4k), keytype="GENEID", column="SEQNAME")
rowData(sce.4k)$Chr <- loc

# Brief quality control.
sce.4k <- calculateQCMetrics(sce.4k, compact=TRUE, feature_controls=list(Mt=which(loc=="MT")))
lowlib <- isOutlier(sce.4k$scater_qc$all$log10_total_counts, type="lower", nmads=3)
lowfeat <- isOutlier(sce.4k$scater_qc$all$log10_total_features_by_counts, type="lower", nmads=3)
highmito <- isOutlier(sce.4k$scater_qc$feature_control_Mt$pct_counts, type="higher", nmads=3)
discard <- lowlib | lowfeat | highmito
##summary(discard)
sce.4k <- sce.4k[,!discard]

# Performing normalization.
set.seed(1000)
clusters <- quickCluster(sce.4k, min.mean=0.1, method="igraph")
##table(clusters)
sce.4k <- computeSumFactors(sce.4k, clusters=clusters, min.mean=0.1, BPPARAM=MulticoreParam(2))
##plot(sce.4k$scater_qc$all$total_counts, sizeFactors(sce.4k), log="xy")
sce.4k <- normalize(sce.4k)

# Modelling the mean-variance trend.
fit.4k <- trendVar(sce.4k, use.spikes=FALSE, loess.args=list(span=0.1))
##plot(fit.4k$mean, fit.4k$vars)
##curve(fit.4k$trend(x), add=TRUE, col="red")
dec.4k <- decomposeVar(fit=fit.4k)

keep <- dec.4k$bio > 0
dir.create("Processed", showWarnings=FALSE)
saveRDS(file=file.path("Processed", "pbmc4k.rds"), as.matrix(logcounts(sce.4k)[keep,]))
