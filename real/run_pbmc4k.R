# Running PC number choice on the 4K data set.

source("../simulations/functions.R")
full <- readRDS("Processed/pbmc4k.rds")
dir.create("results", showWarnings=FALSE)
simulateReal(full$exprs, prefix="results/pbmc4k")
