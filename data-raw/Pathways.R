# Microarrays Dataset -------------------------

Pathways <- read.table("data-raw/Pathways.txt",skip=1)
colnames(Pathways) <- c("z_value_L","z_value_O","z_value_S")

save(Pathways, file = "data/Pathways.rdata")
