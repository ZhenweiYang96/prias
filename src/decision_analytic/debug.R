remove.packages("JMbayes")
setwd("/home/a_tomer/Google Drive/PhD/src/JMBayes/")
devtools::install("Dimitris")

library(JMbayes)

pbc2.id$Time <- pbc2.id$years
pbc2.id$event <- as.numeric(pbc2.id$status != "alive")

pbc2.id = pbc2.id[1:30,]
pbc2.id$id = droplevels(pbc2.id$id)

pbc2 = pbc2[pbc2$id %in% pbc2.id$id,]
pbc2$id = droplevels(pbc2$id)
