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
fit.4k <- makeTechTrend(x=sce.4k)
dec.4k <- decomposeVar(sce.4k, fit=list(trend=fit.4k))
##plot(dec.4k$mean, dec.4k$total)
##curve(fit.4k(x), add=TRUE, col="red")

keep <- dec.4k$bio > 0
dir.create("Processed", showWarnings=FALSE)
saveRDS(file=file.path("Processed", "pbmc4k.rds"), list(exprs=logcounts(sce.4k)[keep,], noise=dec.4k$tech[keep]))
