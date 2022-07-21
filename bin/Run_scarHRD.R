args = commandArgs(trailingOnly=TRUE)
library("scarHRD")
scar_score(args[1],reference = "grch38", seqz=FALSE)