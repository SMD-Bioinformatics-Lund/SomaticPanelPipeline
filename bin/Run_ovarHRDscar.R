args = commandArgs(trailingOnly=TRUE)
library("ovaHRDscar")
a <-read.table(args[1], header=T)
get.ovaHRDscars(a, chrominfo ="grch38")