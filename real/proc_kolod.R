library(scran)
library(scater)

# Pre-processing the Kolodziejczyk dataset.
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)    
path <- bfcrpath(bfc, "https://www.ebi.ac.uk/teichmann-srv/espresso/static/counttable_es.csv")
incoming <- read.table(path, header=TRUE, row.names=1)
cell.type <- sub("^ola_mES_(.*)_[0-9]+_[0-9]+.counts$", "\\1", colnames(incoming))
batch.num <- sub("^ola_mES_.*_([0-9]+)_[0-9]+.counts$", "\\1", colnames(incoming))

keep <- batch.num == "3"
count.data <- incoming[,keep]
metadata <- data.frame(Sample=colnames(incoming), Culture=cell.type, Batch=batch.num)[keep,]

count.data <- count.data[!grepl("^_", rownames(count.data)),]
sce <- SingleCellExperiment(list(counts=as.matrix(count.data)), colData=metadata)

# Adding locational annotation.
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rownames(sce), column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
##table(anno$CHR)

is.spike <- grepl("^ERCC", rownames(sce))
isSpike(sce, "ERCC") <- is.spike

# Brief quality control.
sce <- calculateQCMetrics(sce, compact=TRUE, feature_controls=list(Mt=which(rowData(sce)$CHR=="MT")))
lowlib <- isOutlier(sce$scater_qc$all$log10_total_counts, type="lower", nmads=3)
lowfeat <- isOutlier(sce$scater_qc$all$log10_total_features_by_counts, type="lower", nmads=3)
highspike <- isOutlier(sce$scater_qc$feature_control_ERCC$pct_counts, type="higher", nmads=3)
discard <- lowlib | lowfeat | highspike
##summary(discard)
sce <- sce[,!discard]

# Performing normalization.
sce <- computeSumFactors(sce, min.mean=1)
sce <- computeSpikeFactors(sce)
##plot(sce.4k$scater_qc$all$total_counts, sizeFactors(sce.4k), log="xy")
sce <- normalize(sce)

# Modelling the mean-variance trend.
fit <- trendVar(sce)
##plot(fit$mean, fit$vars)
##curve(fit$trend(x), add=TRUE, col="red")
dec <- decomposeVar(sce, fit)

keep <- dec$bio > 0 & !isSpike(sce)
dir.create("Processed", showWarnings=FALSE)
saveRDS(file=file.path("Processed", "kolod.rds"), as.matrix(logcounts(sce)[keep,]))
