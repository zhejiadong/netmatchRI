#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[[1]] else "../TRIP/data_clear.RData"
out_dir <- if (length(args) >= 2) args[[2]] else "trip-application-output"

library(netmatchRI)

res <- trip_application(data_file = data_file, output_dir = out_dir)
print(res$tests)
print(res$diagnostics)
