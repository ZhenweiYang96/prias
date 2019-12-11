seeds = 2021:2030

biopsyDf_summary = do.call('rbind', lapply(seeds, FUN = function(seed){
  load(paste0("Rdata/lastpaper/simulation/combined_results/seed_", seed, ".Rdata"))
  
  return(biopsyDf_summary)
}))


biopsyDf_summary$nb[biopsyDf_summary$delay<0] = biopsyDf_summary$nb[biopsyDf_summary$delay<0] + 1
biopsyDf_summary$delay[biopsyDf_summary$delay<0] = 10 - biopsyDf_summary$progression_time[biopsyDf_summary$delay<0]

a = ggplot(biopsyDf_summary[biopsyDf_summary$progressed==0,]) +
  geom_boxplot(aes(x=schedule, y=nb)) + ylim(0,10) +
  theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))

b = ggplot(biopsyDf_summary[biopsyDf_summary$progressed==1,]) +
  geom_boxplot(aes(x=schedule, y=nb)) + ylim(0,10)+
  theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))

c = ggplot(biopsyDf_summary[biopsyDf_summary$progressed==1,]) +
  geom_boxplot(aes(x=schedule, y=delay))+
  theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggpubr::ggarrange(a,b,c, nrow=1, ncol=3)
