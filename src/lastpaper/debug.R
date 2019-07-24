load("C:/Users/838035/Google Drive/PhD/src/prias/Rdata/lastpaper/sim_seed_2019_sched.Rdata")
t1 = sim_res$testData$testDs.id
rm(list=setdiff(ls(), "t1"))

load("C:/Users/838035/Google Drive/PhD/src/prias/Rdata/lastpaper/sim_seed_2019_gof.Rdata")

testDs.id$nb_0.05[testDs.id$delay_0.05 < 0] = 1 + testDs.id$nb_0.05[testDs.id$delay_0.05 < 0]
testDs.id$delay_0.05[testDs.id$delay_0.05 < 0] = 10 - testDs.id$progression_time[testDs.id$delay_0.05 < 0]

testDs.id$nb_0.1[testDs.id$delay_0.1 < 0] = 1 + testDs.id$nb_0.1[testDs.id$delay_0.1 < 0]
testDs.id$delay_0.1[testDs.id$delay_0.1 < 0] = 10 - testDs.id$progression_time[testDs.id$delay_0.1 < 0]

testDs.id$delay_0.15[195] = -15
testDs.id$nb_0.15[testDs.id$delay_0.15 < 0] = 1 + testDs.id$nb_0.15[testDs.id$delay_0.15 < 0]
testDs.id$delay_0.15[testDs.id$delay_0.15 < 0] = 10 - testDs.id$progression_time[testDs.id$delay_0.15 < 0]

#######Combining all results
testDs.id = cbind(testDs.id, t1[, 15:24])

testDs.id.long = reshape(testDs.id, varying = list(seq(15, 30, 2), seq(16, 30, 2)), timevar = "method",
             direction='long', idvar='P_ID', v.names = c('nb', 'delay'))

testDs.id.long$method_name = rep(c("5%", "10%", "15%",
                                   "Ratio\n0",
                                   "Ratio\n0.25",
                                   "Ratio\n0.5",
                                   "Ratio\n0.75",
                                   "Ratio\n1"), each=300)

testDs.id.long$method_type = rep(c("Risk", "Risk", "Risk",
                                   "Minimum Dist",
                                   "Minimum Dist",
                                   "Minimum Dist",
                                   "Minimum Dist",
                                   "Minimum Dist"), each=300)


a = ggplot(testDs.id.long[testDs.id.long$progression_time == 10,]) + 
  geom_boxplot(aes(y=nb, x=method_name, color=method_type), outlier.shape = NA) + 
  theme_bw() + xlab("Schedule") + ylab("Number of biopsies") +
  ggtitle("Non progressors")

b = ggplot(testDs.id.long[testDs.id.long$progression_time < 10,]) + 
  geom_boxplot(aes(y=nb, x=method_name, color=method_type), outlier.shape = NA) + 
  theme_bw() + xlab("Schedule") + ylab("Number of biopsies") +
  ggtitle("Progressors")

c = ggplot(testDs.id.long[testDs.id.long$progression_time < 10,]) + 
  geom_boxplot(aes(y=delay, x=method_name, color=method_type), outlier.shape = NA) +
  theme_bw() + xlab("Schedule") + ylab("Delay (years)")+
  ggtitle("Progressors")

ggpubr::ggarrange(a,b,c, nrow = 1, ncol=3, legend = "bottom", common.legend = T)
