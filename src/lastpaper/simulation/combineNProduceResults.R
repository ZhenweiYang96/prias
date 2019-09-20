# combined = do.call('rbind', lapply(c(2020, 2022:2024, 2026:2031), function(seed){
#   print(seed)
# 
#   load(paste0("Rdata/lastpaper/sims/sim_seed_", seed, ".Rdata"))
# 
#   load(paste0("Rdata/lastpaper/resfpr/res_seed_", seed, ".Rdata"))
#   fpr_res = res
#   load(paste0("Rdata/lastpaper/resuptothres/res_seed_", seed, ".Rdata"))
#   threshold_res = res
#   load(paste0("Rdata/lastpaper/res/res_seed_", seed, ".Rdata"))
#   res = cbind(res, threshold_res[, -c(1:12)], fpr_res[, -c(1:12)])
# 
#   res[, c("nb_PRIAS", "delay_PRIAS")] = runPRIASSchedule(sim_res$testData$testDs.id,
#                                                          sim_res$testData$testDs)
# 
#   for(gap in c(1, 1.5, 2, 3, 4)){
#     schedule = seq(gap, 10, gap)
#     if(max(schedule) < 10){
#       schedule = c(schedule, 10)
#     }
# 
#     delay = sapply(res$progression_time, function(x){
#       delays = schedule - x
#       delays[which(delays >= 0)[1]]
#     })
#     nb = sapply(res$progression_time, function(x){
#       which(schedule - x>=0)[1]
#     })
# 
#     res[, paste0("fixed_nb_", gap)] = nb
#     res[, paste0("fixed_delay_", gap)] = delay
#   }
# 
#   methodNames = c(paste0("MD_", c(seq(0,1, 0.1),seq(1.25,2, 0.25))),
#                   paste0("TH_", seq(0,1, 0.05)),
#                   paste0("FP_", seq(0,1, 0.05)),
#                   "PR",paste0("FI_", c(1,1.5,2,3,4)))
#   
#   res_long = reshape(res, varying = list(seq(13, ncol(res), 2),
#                                          seq(14, ncol(res), 2)),
#                      v.names = c("nb", "delay"), idvar = "P_ID",
#                      direction = "long",
#                      timevar = "methodName",
#                      times = methodNames)
# 
#   res_long$methodType = factor(c(rep("MD", 15*300),
#                                  rep("TH", 21*300),
#                                  rep("FP", 21*300),
#                                  rep("PR", 1*300),
#                                  rep("FI", 5*300)))
# 
#   res_long$methodName = factor(res_long$methodName)
#   res_long$seed = seed
# 
#   return(res_long)
# }))
# 
# temp_filter = is.na(combined$delay) | combined$delay < 0
# combined$nb[temp_filter] = combined$nb[temp_filter] + 1
# combined$delay[temp_filter] = 10 - combined$progression_time[temp_filter]
# 
# combined_summary = data.frame(methodName=levels(combined$methodName),
#                               median_nb = as.numeric(by(combined$methodName, data = combined$nb, median)),
#                               q1_nb = as.numeric(by(combined$methodName, data = combined$nb, quantile, probs=0.25)),
#                               q3_nb = as.numeric(by(combined$methodName, data = combined$nb, quantile, probs=0.75)),
#                               mean_nb = as.numeric(by(combined$methodName, data = combined$nb, mean)),
#                               median_nb_nonprog = as.numeric(by(combined$methodName[combined$progression_time==10], data = combined$nb[combined$progression_time==10], median)),
#                               q1_nb_nonprog = as.numeric(by(combined$methodName[combined$progression_time==10], data = combined$nb[combined$progression_time==10], quantile, probs=0.25)),
#                               q3_nb_nonprog = as.numeric(by(combined$methodName[combined$progression_time==10], data = combined$nb[combined$progression_time==10], quantile, probs=0.75)),
#                               mean_nb_nonprog = as.numeric(by(combined$methodName[combined$progression_time==10], data = combined$nb[combined$progression_time==10], mean)),
#                               median_nb_prog = as.numeric(by(combined$methodName[combined$progression_time<10], data = combined$nb[combined$progression_time<10], median)),
#                               q1_nb_prog = as.numeric(by(combined$methodName[combined$progression_time<10], data = combined$nb[combined$progression_time<10], quantile, probs=0.25)),
#                               q3_nb_prog = as.numeric(by(combined$methodName[combined$progression_time<10], data = combined$nb[combined$progression_time<10], quantile, probs=0.75)),
#                               mean_nb_prog = as.numeric(by(combined$methodName[combined$progression_time<10], data = combined$nb[combined$progression_time<10], mean)),
#                               median_delay_prog = as.numeric(by(combined$methodName[combined$progression_time<10], data = combined$delay[combined$progression_time<10], median)),
#                               q1_delay_prog = as.numeric(by(combined$methodName[combined$progression_time<10], data = combined$delay[combined$progression_time<10], quantile, probs=0.25)),
#                               q3_delay_prog = as.numeric(by(combined$methodName[combined$progression_time<10], data = combined$delay[combined$progression_time<10], quantile, probs=0.75)),
#                               mean_delay_prog = as.numeric(by(combined$methodName[combined$progression_time<10], data = combined$delay[combined$progression_time<10], mean)))
# combined_summary$methodType = rep(levels(combined$methodType), table(combined$methodType)/3000)

# nb all patients, but delay for progressions
label_methods = c("FI_1", "FI_2", 
                  "TH_0.05", "TH_0.1", "TH_0.15", "TH_0.5", "TH_1",
                  "MD_1", "MD_2", "MD_0.2", "MD_0.5",
                  "FP_0.1", "FP_0.4", "FP_0.8",
                  "PR")
combined_summary_forlabels = combined_summary[combined_summary$methodName %in% label_methods,]

pdf("Rdata/lastpaper/res_sim.pdf", width = 8, onefile = T, paper = "a4r")

pp = ggplot() +
  geom_point(data = combined_summary, 
             mapping = aes(x=mean_nb, y=mean_delay_prog, 
                 color=methodType, shape=methodType), size=3) +
  geom_label(data = combined_summary_forlabels,
        mapping = aes(x=mean_nb, y=mean_delay_prog, 
                 label=methodName, fill=methodType), color='white', size=3) +
  theme_bw() + theme(legend.position = "bottom") +
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(breaks=seq(0,7,1)) +
  facet_wrap(~methodType) +
  xlab("Mean #biopsies all") + ylab("Mean delay progressors (years)") 

print(pp)

pp = ggplot() +
  geom_point(data = combined_summary, 
             mapping = aes(x=mean_nb, y=mean_delay_prog, 
                           color=methodType, shape=methodType), size=3) +
  geom_label(data = combined_summary_forlabels,
             mapping = aes(x=mean_nb, y=mean_delay_prog, 
                           label=methodName, fill=methodType), color='white', size=3) +
  theme_bw() + theme(legend.position = "bottom") +
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(breaks=seq(0,7,1)) +
  xlab("Mean #biopsies all") + ylab("Mean delay progressors (years)") 

print(pp)

pp = ggplot() +
  geom_point(data = combined_summary, 
             mapping = aes(x=mean_nb_prog, y=mean_delay_prog, 
                           color=methodType, shape=methodType), size=3) +
  geom_label(data = combined_summary_forlabels,
             mapping = aes(x=mean_nb_prog, y=mean_delay_prog, 
                           label=methodName, fill=methodType), color='white', size=3) +
  theme_bw() + theme(legend.position = "bottom") +
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(breaks=seq(0,7,1)) +
  facet_wrap(~methodType) +
  xlab("Mean #biopsies progressors") + ylab("Mean delay progressors (years)") 

print(pp)

pp=ggplot(combined_summary) +
  geom_point(aes(x=median_nb, y=median_delay_prog, 
                 color=methodType, shape=methodType), size=3) +
  geom_label(data = combined_summary_forlabels,
             mapping = aes(x=median_nb, y=median_delay_prog, 
                           label=methodName, fill=methodType), color='white', size=3) +
  theme_bw() + theme(legend.position = "bottom") +
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(breaks=seq(0,7,1)) +
  facet_wrap(~methodType) +
  xlab("Median #biopsies all") + ylab("Median delay progressors (years)") 

print(pp)


pp=ggplot(combined_summary) +
  geom_point(aes(x=median_nb, y=median_delay_prog, 
                 color=methodType, shape=methodType), size=3) +
  geom_label(data = combined_summary_forlabels,
             mapping = aes(x=median_nb, y=median_delay_prog, 
                           label=methodName, fill=methodType), color='white', size=3) +
  theme_bw() + theme(legend.position = "bottom") +
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(breaks=seq(0,7,1)) +
  xlab("Median #biopsies all") + ylab("Median delay progressors (years)") 

print(pp)

pp=ggplot(combined_summary) +
  geom_point(aes(x=median_nb_prog, y=median_delay_prog, 
                 color=methodType, shape=methodType), size=3) +
  geom_label(data = combined_summary_forlabels,
             mapping = aes(x=median_nb_prog, y=median_delay_prog, 
                           label=methodName, fill=methodType), color='white', size=3) +
  theme_bw() + theme(legend.position = "bottom") +
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(breaks=seq(0,7,1)) +
  facet_wrap(~methodType) +
  xlab("Median #biopsies progressors") + ylab("Median delay progressors (years)") 

print(pp)

thres_filter = combined$methodName %in% c(paste0("MD_", c(seq(0,0.4,0.1),1)),
                                          paste0("TH_", seq(0,0.2, 0.05)),
                                          paste0("FP_", seq(0.5,1, 0.05)),
                                          "PR",paste0("FI_", c(1,1.5,2,3)))

#thres_filter = T
# #nb boxplot non progressions
ggplot(combined[combined$progression_time==10 & thres_filter,]) +
  geom_boxplot(aes(x=methodName, y=nb, color=methodType),outlier.shape = NA) +
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(breaks = seq(1,10,1)) +
  xlab("Method") + ylab("#Biopsies") + ggtitle("Non progressors")

# #nb boxplot progressions
ggplot(combined[combined$progression_time<10 & thres_filter,]) +
  geom_boxplot(aes(x=methodName, y=nb, color=methodType),outlier.shape = NA) +
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = seq(1,10,1)) +
  xlab("Method") + ylab("# Biopsies") + ggtitle("Progressors")

# #delay boxplot progressors
ggplot(combined[combined$progression_time<10 & thres_filter,]) +
  geom_boxplot(aes(x=methodName, y=delay, color=methodType),outlier.shape = NA) +
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(breaks = seq(0,7,1), limits = c(0,7)) +
  xlab("Method") + ylab("Delay (years)") + ggtitle("Progressors")

dev.off()
