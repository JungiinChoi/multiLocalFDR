# Raw p-values for pathways -------------------------

pathways <- read.csv("data-raw/pathways.csv")[,-1]

## -Inf fix
pathways[238,4] = 9.99999e-01

save(pathways, file = "data/pathways.rdata")
