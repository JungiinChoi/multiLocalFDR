# Galaxy Dataset -------------------------

Galaxy <- read.table("data-raw/Galaxy.txt", skip=1)
colnames(Galaxy) <- c("radius","velocity")

save(Galaxy, file = "data/Galaxy.rdata")
