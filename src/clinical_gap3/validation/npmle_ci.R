args = commandArgs(trailingOnly = T)

library(Icens)
library(interval)
library(survival)

load("Rdata/gap3/PRIAS_2019/motherdata.Rdata")

iteration = as.numeric(args[1])
seed = 2019 + iteration

set.seed(seed)
npmle_all = lapply(cohortnames, FUN = function(cohort){
  print(paste('doing for cohort: ', cohort))
  
  longdata.id = get(paste0("longdata_", cohort,".id"))
  
  if(iteration > 1){
    longdata.id = longdata.id[sample(1:nrow(longdata.id), replace = T),]
    longdata.id$P_ID = 1:nrow(longdata.id)
  }
  
  npmleres = icfit(Surv(latest_survival_time,earliest_failure_time,type="interval2")~1, 
                      data=longdata.id, control=icfitControl(maxit = 50000))
  
  return(npmleres)
})

save(npmle_all, file="Rdata/gap3/PRIAS_2019/npmle_all.Rdata")
