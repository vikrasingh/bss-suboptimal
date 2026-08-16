.libPaths(c('/home/x-vsingh7/R', .libPaths()))
#!/usr/bin/env Rscript
options(warn = 2)

# ---- READ COMMAND-LINE INPUTS ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Missing inputFile or outputFile arguments!")
}
inputFile  <- args[1]
outputFile <- args[2]

# ---- ABSOLUTE PATH TO FUNCTION FILE ----
source("/home/x-vsingh7/lls_package/utility/R_files/bssR.R")

# ---- LOAD DATA ----
library(R.matlab)
library(gurobi)
library(bestsubset)

data <- readMat(inputFile)

xMatrix <- data$xMatrix
yArray <- as.vector(data$yArray)
tmax <- as.numeric(data$tmaxValuesBss)
cpulimit <- as.numeric(data$cpulimit)

# ---- VALIDATION (CRITICAL) ----
stopifnot(is.matrix(xMatrix))
stopifnot(is.numeric(xMatrix))
stopifnot(is.numeric(yArray))
stopifnot(nrow(xMatrix) == length(yArray))


# ---- CALL MODEL FUNCTION ----
result_bss <- bssR(xMatrix, yArray, tmax, cpulimit)
writeMat(outputFile, beta = as.matrix(result_bss$beta), ctime=as.numeric(result_bss$cputime), status=as.character(result_bss$status))

# ---- OUTPUT (MATLAB READS THIS) ----
cat("R_SCRIPT_FINISHED_SUCCESSFULLY")
cat(paste(result_bss$beta, collapse=" "))
