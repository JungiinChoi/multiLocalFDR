# Microarrays Dataset -------------------------

Microarrays <- read.table("data-raw/Microarrays.txt",skip=1)

save(Microarrays, file = "data/Microarrays.rdata")
