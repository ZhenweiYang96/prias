load("Rdata/gap3/PRIAS_2019/km_all.Rdata")

kmdf = data.frame(cohort=unlist(lapply(names(km_all), FUN = function(name){rep(name, length(km_all[[name]]$time))})),
                  timePoints=unlist(lapply(km_all, "[[", "time")),
                  cumrisk = 1-unlist(lapply(km_all, "[[", "surv")))
kmdf = kmdf[order(kmdf$cohort, kmdf$timePoints),]

kmdf$time_10pat_risk_set = rep(reclassification_df$time_10pat_risk_set[order(reclassification_df$Cohort)], table(kmdf$cohort))
kmdf = kmdf[kmdf$timePoints <= kmdf$time_10pat_risk_set,]

cohort_names = levels(kmdf$cohort)
cohort_labpos_x = as.numeric(by(kmdf$cohort, data = kmdf$timePoints, max))
cohort_labpos_y = as.numeric(by(kmdf$cohort, data = kmdf$cumrisk, max))

FONT_SIZE=13
km_plot_all = ggplot() + 
  geom_line(aes(x=kmdf$timePoints, 
                y=kmdf$cumrisk, 
                group=kmdf$cohort, 
                color=kmdf$cohort)) +  
  geom_label(aes(x=cohort_labpos_x, 
                 y=cohort_labpos_y, 
                 label=cohort_names,
                 fill=cohort_names), color='white')+
  coord_cartesian(xlim=c(0,8)) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), 
        axis.text=element_text(size=FONT_SIZE),
        legend.position = "none",
        legend.text = element_text(size=FONT_SIZE-4),
        axis.line = element_line(), 
        plot.margin = margin(0, 0, 0, 0, "pt")) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"),
                     limits = c(0,1)) + 
  ylab("Cumulative risk of reclassification (%)") +
  xlab("Follow-up time (years)")

ggsave(filename = "report/clinical/images/km_plot.eps",
       plot=km_plot_all, device=cairo_ps, height=5.5, width=6, dpi = 500)


