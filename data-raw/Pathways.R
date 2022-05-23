# Raw p-values for pathways -------------------------

pathways <- read.csv("data-raw/pathways.csv")[,-1]

save(pathways, file = "data/pathways.rdata")
