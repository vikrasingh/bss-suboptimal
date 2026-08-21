#!/usr/bin/env Rscript
options(warn = 2)

# ---- ABSOLUTE PATH TO FUNCTION FILE ----
source("/home/x-vsingh7/lls_package/utility/R_files/abessUbxR.R")

# ---- READ COMMAND-LINE INPUTS ----
args <- commandArgs(trailingOnly = TRUE)
inputFile  <- args[1]
outputFile <- args[2]

# ---- LOAD DATA ----
library(R.matlab)

data <- readMat(inputFile)

xMatrix <- data$xMatrix
yArray <- as.vector(data$yArray)
tmax1 <- as.numeric(data$tmax1)

# ---- VALIDATION (CRITICAL) ----
stopifnot(is.matrix(xMatrix))
stopifnot(is.numeric(xMatrix))
stopifnot(is.numeric(yArray))
stopifnot(nrow(xMatrix) == length(yArray))


# ---- CALL MODEL FUNCTION ----
res <- abessUbxR(xMatrix, yArray, tmax1)

writeMat(outputFile, beta = res[[1]], ctime= res[[2]])

# ---- OUTPUT (MATLAB READS THIS) ----
cat(sprintf("%.10f", res$out1))
  