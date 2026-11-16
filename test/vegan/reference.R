# Reference solutions from vegan::metaMDS for test/vegan/compare.jl.
# Usage: Rscript reference.R OUTDIR
suppressMessages(library(vegan))
out <- commandArgs(TRUE)[1]
data(dune)
D <- vegdist(dune, "bray")
set.seed(3)
m <- metaMDS(D, k = 2, trymax = 100, autotransform = FALSE, trace = 0)
write.csv(as.matrix(D), file.path(out, "D.csv"), row.names = FALSE)
write.csv(scores(m), file.path(out, "points.csv"), row.names = FALSE)
write.csv(data.frame(stress = m$stress), file.path(out, "stress.csv"), row.names = FALSE)
cat("vegan", as.character(packageVersion("vegan")), "stress", m$stress, "\n")
