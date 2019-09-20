library(JMbayes)
library(splines)

load("Rdata/lastpaper/sims/sim_seedbig_light_2019.Rdata")
source("src/lastpaper/goodness_of_fit.R")
source("src/lastpaper/prediction.R")

MAX_FOLLOW_UP = 10

PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP, 0.5))
DRE_CHECK_UP_TIME = seq(0, MAX_FOLLOW_UP, 0.5)
BIOPSY_TEST_TIMES = PSA_CHECK_UP_TIME

roc_results_cache = vector("list", length(BIOPSY_TEST_TIMES))
names(roc_results_cache) = BIOPSY_TEST_TIMES
for(i in 1:length(BIOPSY_TEST_TIMES)){
  roc_results_cache[[i]] = vector("list", length(BIOPSY_TEST_TIMES))
  names(roc_results_cache[[i]]) = BIOPSY_TEST_TIMES
}

possible_latest_survival_times = c(0, BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES>=1 & BIOPSY_TEST_TIMES<=9])

for(latest_survival_time in possible_latest_survival_times){
  print(paste('latest survival time', latest_survival_time))
  for(cur_visit_time in BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES >= latest_survival_time + 1]){
    print(paste('current visit time', cur_visit_time))
    gof = goodness_of_fit(sim_res$mvJoint_psa_simDs, sim_res$trainingData$trainingDs,
                          T_start = latest_survival_time, T_horiz = cur_visit_time, horizon = MAX_FOLLOW_UP, M = 500)
    roc_results_cache[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]] = gof
    
    save(roc_results_cache, file = "Rdata/lastpaper/rocresultscache_forward.Rdata")
  }
}

pdf(file = "Rdata/lastpaper/roc_results.pdf", onefile = T)
for(latest_survival_time in possible_latest_survival_times){
  print(paste('latest survival time', latest_survival_time))
  for(cur_visit_time in BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES >= latest_survival_time + 1]){
    gof = sim_res$roc_results[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]]
    
    print(ggplot() + geom_line(data=gof$roc_results, aes(x=fpr, y=tpr, color='CUM_RISK')) + 
      geom_line(data=gof$roc_results_auc_T_start_T_horiz, aes(x=fpr, y=tpr, color='AUC')) + 
      geom_abline(slope=1, intercept = 0) + 
      theme_bw() + theme(legend.position = "bottom", 
                         legend.title = element_blank(), 
                         text = element_text(size=14)) + 
      ylab("TPR") + xlab("FPR") + xlim(0,1) + ylim(0,1) + 
        ggtitle(paste(latest_survival_time, "--", cur_visit_time)))
  }
}
dev.off()
