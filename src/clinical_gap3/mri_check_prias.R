library(foreign)

load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")

prias=read.spss(file = file.choose(),
                to.data.frame = T)

#look for anything in colnames that says pirads...represents pirads score
# and the following colnumbers had pirads in them
prias_mri = prias[, c(1, which(grepl(colnames(prias), pattern = 'Mri_Date', fixed = F)==T),
                      which(grepl(colnames(prias), pattern = 'Pirad', fixed = F)==T))]

prias_mri_long = reshape(prias_mri, direction='long', idvar='P_ID', timevar = "visit_number",
                         varying=list(c(2:48), c(49:95), c(96:142), c(143:189)), 
                         v.names=c('date', 'pirads1', 'pirads2','pirads3'))

prias_mri_long = prias_mri_long[order(prias_mri_long$P_ID),]
prias_mri_long$any_pirads = sapply(1:nrow(prias_mri_long), function(i){
  sum(prias_mri_long$pirads1[i], prias_mri_long$pirads2[i], prias_mri_long$pirads3[i], na.rm = T)
})

prias_mri_long$any_pirads[prias_mri_long$any_pirads %in% 0]=NA
table(prias_mri_long$any_pirads, useNA = "always")
sum(!is.na(prias_mri_long$any_pirads))

YEAR_DIVIDER = 24 * 60 * 60 * 365
prias_mri_long$visitsecs = sapply(prias_mri_long$date, function(x){
  tt = try(difftime(x, "1582-10-14", units='secs'), silent = T)
  if(inherits(tt, "try-error")){
    return(NA)
  }else{
    return(as.numeric(tt))
  }
})

prias_mri_long = prias_mri_long[prias_mri_long$P_ID %in% prias_final.id$P_ID,]
prias_mri_long$dom_diagnosis = rep(prias_final.id$dom_diagnosis, 47)
prias_mri_long$reclassification = rep(prias_final.id$reclassification, 47)
prias_mri_long$latest_survival_time = rep(prias_final.id$latest_survival_time, 47)
prias_mri_long$earliest_failure_time = rep(prias_final.id$earliest_failure_time, 47)
prias_mri_long$year_max_followup = rep(prias_final.id$year_max_followup, 47)

prias_mri_long$year_visit = (prias_mri_long$visitsecs - prias_mri_long$dom_diagnosis)/YEAR_DIVIDER
prias_mri_long = prias_mri_long[!is.na(prias_mri_long$year_visit) & prias_mri_long$year_visit>0 & prias_mri_long$year_visit<10,]
prias_mri_long = prias_mri_long[!is.na(prias_mri_long$any_pirads),]
prias_mri_long = prias_mri_long[order(prias_mri_long$P_ID, prias_mri_long$year_visit),]


prias_mri_long = prias_mri_long[prias_mri_long$year_visit<=prias_mri_long$year_max_followup,]
View(prias_mri_long[abs(prias_mri_long$year_visit - prias_mri_long$earliest_failure_time)<=0.5,])


prias_mri_long = prias_mri_long[!is.na(prias_mri_long$any_pirads),]
prias_final.id_mri = prias_final.id[prias_final.id$P_ID %in% prias_mri_long$P_ID,]
#Only 115 could be detected due to MRI if that was the case.
sum(prias_final.id_mri$reclassification)

length(unique(prias_mri_long$P_ID[!is.na(prias_mri_long$any_pirads)]))

#How many patients had MRI within 6 months of their
