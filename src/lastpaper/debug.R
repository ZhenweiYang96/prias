files = list.files("Rdata/lastpaper/simulation/combined_results_exp_tests/", full.names = T)
file_normal_names = list.files("Rdata/lastpaper/simulation/combined_results_exp_tests/", full.names = F)

for(i in 1:length(files)){
  load(files[i])
  biopsyDf_summary=biopsyDf_summary[biopsyDf_summary$schedule %in% c("Annual", "Biennial", "PRIAS", "Risk: 5%", "Risk: 10%", "Risk: 15%"),]
  save(biopsyDf_summary,
       file=paste0("Rdata/lastpaper/simulation/combined_results_non_auto/", file_normal_names[i]))
}

