library(JMbayes)
library(MASS)
library(splines)
library(survival)

load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
source("src/lastpaper/simulation/simCommon_psa.R")

nSub = 1200
bNames = paste0("b", 0:4)
for(seed in 2020:2030){
  set.seed(seed)
  sim_data= generateSimulationData(nSub, bNames = bNames)
  sim_res = fitJointModelOnNewData(sim_data$simDs, sim_data$simDs.id, nSub * 0.75)
  save(sim_res, file=paste0("Rdata/lastpaper/sims/sim_seed_", seed, ".Rdata"))
}