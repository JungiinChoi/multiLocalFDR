# P-values for microarray data -------------------------

microarrays <- read.csv("data-raw/microarrays.csv")[,-1]

save(microarrays, file = "data/microarrays.rdata")
