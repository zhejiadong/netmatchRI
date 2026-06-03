#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
task <- if (length(args) >= 1) args[[1]] else "table1"
out_dir <- if (length(args) >= 2) args[[2]] else "paper-output"
n_rep <- if (length(args) >= 3) as.integer(args[[3]]) else 500L

library(netmatchRI)

res <- reproduce_paper(task = task, n_rep = n_rep, output_dir = out_dir)
print(names(res))
