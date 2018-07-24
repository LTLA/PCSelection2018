# Running PC number choice on the 4K data set.

source("../simulations/functions.R")
full <- readRDS("Processed/kolod.rds")
dir.create("results", showWarnings=FALSE)
simulateReal(full, prefix="results/kolod")
