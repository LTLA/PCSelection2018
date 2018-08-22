# Running PC number choice on the 4K data set.

source("../simulations/functions.R")
full <- readRDS("Processed/pbmc4k.rds")
dir.create("results", showWarnings=FALSE)
simulateReal(as.matrix(full$exprs), full$noise, prefix="results/pbmc4k")

# Seeing the values chosen in practice.
dec.out <- decomposer(as.matrix(full$exprs), nrank=50, approximate=TRUE)
chooseNumber(full$exprs, dec.out$SVD, full$noise)
