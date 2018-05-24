########################################
# Graphs at baseline
########################################

#What age group do we have for our patients?
ggplot(data=prias, aes(x=Age)) + geom_histogram()

#Did the people who discontinued visited less number of times
ggplot(data=prias) + geom_histogram(aes(Nr_FUvisits)) + facet_grid(~DiscontinuedYesNo)

#Did the people who discontinued had low/high PSA scores at baseline
ggplot(data=prias) + geom_histogram(aes(PSA)) + facet_grid(~DiscontinuedYesNo)

#To less people to make much sense out of this graph
ggplot(data=prias) + geom_histogram(aes(PSA)) + facet_grid(~diGleason_sum)

#Did the people who discontinued had low/high Gleason scores at baseline
ggplot(data=prias) + geom_bar(aes(diGleason_sum)) + facet_grid(~DiscontinuedYesNo)

#PSA , Gleason and DRE score spaghetti plot
numSub = 1000
selectedSubPID = sample(prias$P_ID, numSub)
ggplot(data=prias_long[prias_long$P_ID %in% selectedSubPID,], 
       aes(x=visit_number, y=psa, group=P_ID)) + 
  geom_line() +
  stat_summary(aes(group=1), colour="#FD5F00", geom = "point", fun.y = median)+
  ylim(0, 75)+
  facet_grid(~DiscontinuedYesNo)

#Gleason should be plotted against visit_number or biopsy number?
ggplot(data=prias_long[prias_long$P_ID %in% selectedSubPID & !is.na(prias_long$gleason),]) + 
  geom_line(aes(x=biopsy_number, y=gleason, group=P_ID, color=P_ID)) + guides(color=FALSE) + 
  facet_grid(~DiscontinuedYesNo)

#Gleason and DRE percentage stacked bar charts
colorPalette = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                 "#cab2d6", "#6a3d9a")

ggplot(data = prias_long[!is.na(prias_long$biopsy_number),], aes(factor(biopsy_number))) + 
  geom_bar(aes(fill = digleason), position = "fill") + 
  xlab("Time") + ylab("Percent") + facet_grid(~DiscontinuedYesNo) + 
  scale_fill_manual(values=colorPalette[sample(1:10, 10)])

ggplot(data = prias_long[!is.na(prias_long$biopsy_number),], aes(factor(visit_number))) + 
  geom_bar(aes(fill = digleason), position = "fill") + 
  xlab("Time") + ylab("Percent") + facet_grid(~DiscontinuedYesNo) + 
  scale_fill_manual(values=colorPalette[sample(1:10, 10)])

ggplot(data = prias_long[!is.na(prias_long$dre),], aes(factor(visit_number))) + 
  geom_bar(aes(fill = dre), position = "fill") + 
  xlab("Time") + ylab("Percent") + facet_grid(~DiscontinuedYesNo) + 
  scale_fill_manual(values=colorPalette[sample(1:10, 10)])

#How are the PSA scores correlated with each other over time
gg

#Are the visits more or less periodic for every person. check for first 10 visits
ggplot(data=prias_long[prias_long$visit_number<10,], aes(x=diff_dom)) + geom_histogram() + facet_wrap(~visit_number)

prias_sub_dataset = prias_long[!is.na(prias_long$biopsy_number),]
prias_sub_dataset$gleason_risk=sapply(prias_sub_dataset$gleason, function(gleason){
  if(gleason<=6){
    "Low Risk"
  }else if(gleason==7){
    "Medium Risk"
  }else{
    "High Risk"
  }
}, simplify = T)

#Does the biopsy gleason score always increases
prias_biopsy = prias_long[!is.na(prias_long$biopsy_number), ]
is_unsorted_gleason = sapply(unique(prias_biopsy$P_ID), function(pid){
  is.unsorted(prias_biopsy[prias_biopsy$P_ID==pid, ]$gleason)
}, simplify = T)

pid_unsorted_gleason = unique(prias_biopsy$P_ID)[(1:length(is_unsorted_gleason))[is_unsorted_gleason]]
View(prias_biopsy[prias_biopsy$P_ID %in% pid_unsorted_gleason,])

##########
#Biopsies per person, impact of less number of biopsies on modeling
##########
table(prias$Gleason_sum)

table(prias_sub_dataset$biopsy_number, prias_sub_dataset$gleason_risk)

ggplot(data = prias_sub_dataset, aes(factor(biopsy_number))) + geom_bar(aes(fill = gleason_risk)) + 
  geom_text(stat='count',aes(label=..count..),vjust=-0.5) + xlab("Biopsy Number") + ylab("Count")

ggplot(data = prias_sub_dataset, aes(factor(biopsy_number))) + geom_bar(aes(fill = gleason_risk), position = "fill") + 
  xlab("Biopsy Number") + ylab("Percent")