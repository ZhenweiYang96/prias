#The last visit date is 24 Dec 2016, that is 13701916800 seconds from 14 Oct 1582

prias.id$max_visit_years = c(by(prias_long$P_ID, data=prias_long$visitTimeYears, FUN = max))
prias.id$last_biopsy_time = c(by(prias_long$P_ID, data=prias_long, FUN = function(x){
  max(x$visitTimeYears[!is.na(x$gleason)])
}))
prias_long$max_visit_years = rep(prias.id$max_visit_years, prias.id$nr_visits)
prias_long$last_biopsy_time = rep(prias.id$last_biopsy_time, prias.id$nr_visits)

#########################################
# Question: how many people (%) are non-compliant for year 1 biopsy (plus minus 0.5 years)
#########################################
#We answer this by first taking all people who are starting at or before (24 Dec 2016 - 1.5 years)
#That is, people with firstVisit <= (13701916800 - 1.5 * 365 * 24 * 60 * 60)
year1 = prias_long[prias_long$firstVisitDom <= (13701916800 - 1.5 * 365 * 24 * 60 * 60),]

#Now we take all those patients who did not fail within 0.5 years of follow up
year1 = year1[year1$progression_time_end > 0.5,]

#Now we remove all those patients who have progression time end = Inf and PSA visit time < 1.5 years
year1 = year1[!(year1$progression_time_end==Inf & year1$max_visit_years<1.5 & year1$last_biopsy_time<0.5),]

year1$P_ID = droplevels(year1$P_ID)

sum(by(year1$P_ID, data = year1, FUN = function(x){
  biopsyTimes = x$visitTimeYears[!is.na(x$gleason)]
  any(biopsyTimes >= 0.5 & biopsyTimes <=1.5)
}))/length(unique(year1$P_ID))

# 84.38% are compliant for biopsy at year 1

#########################################
# Question: how many people (%) are non-compliant for biopsy at year 4
#########################################
#We answer this by first taking 
#That is, people with firstVisit <= (13701916800 - 4 * 365 * 24 * 60 * 60)
year4 = prias_long[prias_long$firstVisitDom <= (13701916800 - 4.5 * 365 * 24 * 60 * 60),]

#Now we take all those patients who did not fail within 3.5 years of follow up
year4 = year4[year4$progression_time_end > 3.5,]

#Now we remove all those patients who have progression time end = Inf and PSA visit time < 4.5 years
year4 = year4[!(year4$progression_time_end==Inf & year4$max_visit_years<4.5 & year4$last_biopsy_time<3.5),]

year4$P_ID = droplevels(year4$P_ID)

sum(by(year4$P_ID, data = year4, FUN = function(x){
  biopsyTimes = x$visitTimeYears[!is.na(x$gleason)]
  any(biopsyTimes >= 3.5 & biopsyTimes <=4.5)
}))/length(unique(year4$P_ID))

# 63.8% are compliant for biopsy at year 4
biopsiesBefore = c(by(year4$P_ID, data = year4, FUN = function(x){
  biopsyTimes = x$visitTimeYears[!is.na(x$gleason)]
  sum(biopsyTimes < 3.5)
}))
# I checked who were compliants and non compliants. It neither depends on Age, nor on previous number of biopsies

#########################################
# Question: how many people (%) are non-compliant for biopsy at year 7
#########################################
#We answer this by first taking 
#That is, people with firstVisit <= (13701916800 - 7.5 * 365 * 24 * 60 * 60)
year7 = prias_long[prias_long$firstVisitDom <= (13701916800 - 7.5 * 365 * 24 * 60 * 60),]

#Now we take all those patients who did not fail within 6.5 years of follow up
year7 = year7[year7$progression_time_end > 6.5,]

#Now we remove all those patients who have progression time end = Inf and PSA visit time < 6.5 years
year7 = year7[!(year7$progression_time_end==Inf & year7$max_visit_years<7.5 & year7$last_biopsy_time<6.5),]

year7$P_ID = droplevels(year7$P_ID)

sum(by(year7$P_ID, data = year7, FUN = function(x){
  biopsyTimes = x$visitTimeYears[!is.na(x$gleason)]
  any(biopsyTimes >= 6.5 & biopsyTimes <=7.5)
}))/length(unique(year7$P_ID))

# 60.5% are compliant for biopsy at year 7
biopsiesBefore = c(by(year7$P_ID, data = year7, FUN = function(x){
  biopsyTimes = x$visitTimeYears[!is.na(x$gleason)]
  sum(biopsyTimes < 6.5)
}))
# I checked who were compliants and non compliants. It neither depends on Age, nor on previous number of biopsies

#########################################
# Question: how many people (%) are non-compliant for biopsy at year 10
#########################################
#We answer this by first taking 
#That is, people with firstVisit <= (13701916800 - 10 * 365 * 24 * 60 * 60)
year10 = prias_long[prias_long$firstVisitDom <= (13701916800 - 10.5 * 365 * 24 * 60 * 60),]

#Now we take all those patients who did not fail within 9.5 years of follow up
year10 = year10[year10$progression_time_end > 9.5,]

#Now we remove all those patients who have progression time end = Inf and PSA visit time < 9.5 years
year10 = year10[!(year10$progression_time_end==Inf & year10$max_visit_years<10.5 & year10$last_biopsy_time<9.5),]

year10$P_ID = droplevels(year10$P_ID)

sum(by(year10$P_ID, data = year10, FUN = function(x){
  biopsyTimes = x$visitTimeYears[!is.na(x$gleason)]
  any(biopsyTimes >= 9.5 & biopsyTimes <=10.5)
}))/length(unique(year10$P_ID))

# 37.5% are compliant at year 10
biopsiesBefore = c(by(year10$P_ID, data = year10, FUN = function(x){
  biopsyTimes = x$visitTimeYears[!is.na(x$gleason)]
  sum(biopsyTimes < 9.5)
}))
# I checked who were compliants and non compliants. It neither depends on Age, nor on previous number of biopsies

######################################################################
# Question: How many people adhere to yearly biopsies due to PSA rise before year 2
######################################################################
switchToAnnual = function(ds){
  psaDt = 1/(lm(log2psa~visitTimeYears, data = ds)$coefficients[2])
  return(psaDt>=0 & psaDt<=10)
}

psa_ds = prias_long[!is.na(prias_long$psa) & prias_long$P_ID==101,]
gleason_ds = prias_long[!is.na(prias_long$gleason),]

by(psa_ds$P_ID, data=psa_ds, FUN = function(patientDs){
  lastBiopsyTime = 0
  proposedBiopsyTime = Inf
  
  proposedAndDoneList = list()
  
  fixedSchedule = c(1, 4, 7, 10)
  
  isAnnualBiopsyProposed = F
  
  gleason_ds_pat = gleason_ds[gleason_ds$P_ID == patientDs$P_ID[1],]
  usedIndices = c()
  
  if(nrow(patientDs)>=4){
    #Min number of measurements before which PSA-DT can't be used
    for(j in 4:nrow(patientDs)){
      curVisitTime = patientDs$visitTimeYears[j]
      
      if(curVisitTime > proposedBiopsyTime){
        proposedBiopsyTime = Inf
        minIndex = which.min(abs(gleason_ds_pat$visitTimeYears - proposedBiopsyTime))
        if(abs(gleason_ds_pat$visitTimeYears - proposedBiopsyTime)<=0.5 & !(minIndex %in% usedIndices)){
          usedIndices = c(usedIndices, minIndex)
          
          proposedAndDoneList[[length(proposedAndDoneList)+1]] = c("proposedBiopsyTime"=proposedBiopsyTime, 
                                                                   "conductedBiopsyTime" = gleason_ds_pat$visitTimeYears[minIndex],
                                                                   "isAnnualBiopsyProposed"=isAnnualBiopsyProposed)
        }else{
          proposedAndDoneList[[length(proposedAndDoneList)+1]] = c("proposedBiopsyTime"=proposedBiopsyTime, 
                                                                   "conductedBiopsyTime" = NA,
                                                                   "isAnnualBiopsyProposed"=isAnnualBiopsyProposed)
        }
      }
      
      lastBiopsyTime = max(gleason_ds_pat$visitTimeYears[gleason_ds_pat$visitTimeYears<=curVisitTime])
      
      if(!is.na(offset) & offset > 0){
        break
      }
      
      if(switchToAnnual(patientDs[1:j, ])){
        if((curVisitTime - lastBiopsyTime) >= 1){
          proposedBiopsyTime = curVisitTime
        }else{
          proposedBiopsyTime = lastBiopsyTime + 1
        }
        isAnnualBiopsyProposed = T
      }else{
        proposedBiopsyTime = fixedSchedule[which(fixedSchedule >= curVisitTime)[1]]
        
        if(proposedBiopsyTime - lastBiopsyTime < 1){
          proposedBiopsyTime = lastBiopsyTime + 1
        }
        isAnnualBiopsyProposed = F
      }
    }
  }else{
    return(NULL)
  }
})




