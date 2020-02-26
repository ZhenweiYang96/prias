library(ggplot2)
library(ggpubr)

load(file = "Rdata/lastpaper/simulation/sim_res.Rdata")

levels(sim_res$schedule)

missing_sims = sim_res[0,]

seeds = 2001:2500

all_res=data.frame(expand.grid(seed=2001:2500, P_ID=751:1000, schedule=levels(sim_res$schedule)))
all_res = all_res[order(all_res$seed, all_res$P_ID, all_res$schedule),]
all_res$index = as.character(apply(all_res, 1, paste, collapse='-'))

sim_res$index = apply(sim_res[,c("seed", "P_ID", "schedule")], 1, paste, collapse='-')

all_res = all_res[!(all_res$index %in% sim_res$index),]
all_res = all_res[,-4]

all_res = droplevels(all_res)
levels(all_res$schedule)
