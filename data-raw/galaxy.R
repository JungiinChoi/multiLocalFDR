# Velocity of Galaxies -------------------------

galaxy <- read.table("data-raw/galaxy.txt", skip=1)
colnames(galaxy) <- c("radius","velocity")

save(galaxy, file = "data/galaxy.rdata")
