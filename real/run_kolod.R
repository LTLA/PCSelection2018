# Running PC number choice on the Kolod data set.

source("../simulations/functions.R")
full <- readRDS("Processed/kolod.rds")
dir.create("results", showWarnings=FALSE)
simulateReal(full$exprs, full$noise, prefix="results/kolod")

# Seeing the values in practice.
dec.out <- decomposer(full$exprs, nrank=50)
chooseNumber(full$exprs, dec.out$SVD, full$noise)
