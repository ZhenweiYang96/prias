library(foreign)
library(doParallel)

############################################
#The list of patients to be ignored
#############################################
ignoreList = factor(read.csv(file.choose(), header = T)$P_ID)

##############################################
# Load the data set
##############################################
prias = read.spss(file.choose(), to.data.frame=TRUE)
prias = prias[order(prias$P_ID, decreasing = F),]
prias$P_ID = droplevels(as.factor(prias$P_ID))

prias_biopsy_cores = read.spss(file.choose(), to.data.frame = T)
prias_biopsy_cores = prias_biopsy_cores[order(prias_biopsy_cores$P_ID, decreasing = F),]
prias_biopsy_cores$P_ID = droplevels(as.factor(prias_biopsy_cores$P_ID))

prias_biopsy_dates = read.spss(file.choose(), to.data.frame = T)
prias_biopsy_dates = prias_biopsy_dates[order(prias_biopsy_dates$P_ID, decreasing = F),]
prias_biopsy_dates$P_ID = droplevels(as.factor(prias_biopsy_dates$P_ID))

#############################################
# Remove ignore list from each of the 3 data sets
#############################################
prias = prias[!(prias$P_ID %in% ignoreList),]
prias_biopsy_cores = prias_biopsy_cores[!(prias_biopsy_cores$P_ID %in% ignoreList),]
prias_biopsy_dates = prias_biopsy_dates[!(prias_biopsy_dates$P_ID %in% ignoreList),]

setdiff(prias$P_ID, prias_biopsy_cores$P_ID)
setdiff(prias_biopsy_cores$P_ID, prias$P_ID)
setdiff(prias$P_ID, prias_biopsy_dates$P_ID)
diff1 = setdiff(prias_biopsy_dates$P_ID, prias$P_ID)

prias_biopsy_dates = prias_biopsy_dates[!(prias_biopsy_dates$P_ID %in% diff1),]

prias$P_ID = droplevels(prias$P_ID)
prias_biopsy_cores$P_ID = droplevels(prias_biopsy_cores$P_ID)
prias_biopsy_dates$P_ID = droplevels(prias_biopsy_dates$P_ID)

any(prias$P_ID!=prias_biopsy_dates$P_ID)
any(prias$P_ID!=prias_biopsy_cores$P_ID)

save.image("Rdata/Gleason as event/cleandata.Rdata")

######################################################
#Combine the data sets into PRIAS
######################################################
prias = cbind(prias, prias_biopsy_cores[,-1], prias_biopsy_dates[,-1])
colnames(prias)

######################################################
#Now we start fixing the data set
######################################################
colnames(prias)[match("Discontinued", colnames(prias))] = "DiscontinuedType"
prias$REASO0 = NULL
prias$Dummy.0 = rep(0, nrow(prias))
prias$Dummy.0.repeat = rep(0, nrow(prias))

##############################################
# Remove subjects
##############################################
#Some patients are too young to be real.
prias = prias[prias$Age > 30,]
prias$P_ID = droplevels(prias$P_ID)

rm(prias_biopsy_cores)
rm(prias_biopsy_dates)
rm(ignoreList)
rm(diff1)

##############################################
# Repeat biopsy score update
##############################################
#First replace all original biopsy scores as NA if num of cores are 0 or 99, 999
prias[prias$Num_cores %in% c(0, 99, 999, NA) | prias$Num_Cores_PC %in% c(99, 999, NA),]$Gleason_sum = NA

#Replaces all repeat biopsy scores as NA if num of cores are 0 or 99, 999
prias[prias$Num_cores2 %in% c(0, 99, 999, NA) | prias$Num_Cores_PC2 %in% c(99, 999, NA),]$Gleason1_2 = NA
prias[prias$Num_cores2 %in% c(0, 99, 999, NA) | prias$Num_Cores_PC2 %in% c(99, 999, NA),]$Gleason2_2 = NA

#If the repeat gleason is not NA then replace original gleason score with NA
# Don't replace original gleason as we want to keep its date intact
prias$Gleason_sum_repeat = prias$Gleason1_2 + prias$Gleason2_2
prias[!is.na(prias$Gleason_sum_repeat),]$Gleason_sum = NA

#Remove those subjects whose repeat Gleason score or whose first Gleason score is more than 6
prias = prias[!(prias$Gleason_sum_repeat %in% c(7,8,9,10)) & !(prias$Gleason_sum %in% c(7,8,9,10)),]

#These two are added just to balance the data set before converting to long
prias$psa_repeat = rep(NA, nrow(prias))
prias$dre_repeat = rep(NA, nrow(prias))

##############################################
# Data type cleaning before conversion to long
##############################################

# DO NOT ATTEMPT THE COMMENTED. THIS GIVES ISSUE WHEN MERGING THIS COLUMN WITH THE OTHER 35 to make a long dataset
# prias$Gleason_sum = as.ordered(prias$Gleason_sum)
# prias$Gleason1_2 = as.ordered(prias$Gleason1_2)
# prias$Gleason2_2 = as.ordered(prias$Gleason2_2)

levels(prias$DRE) = trimws(levels(prias$DRE))
prias$DRE = as.ordered(droplevels(prias$DRE))

levels(prias$DiscontinuedType) = trimws(levels(prias$DiscontinuedType))
prias$DiscontinuedType[prias$DiscontinuedType==""] = NA
prias$DiscontinuedType = droplevels(prias$DiscontinuedType)

levels(prias$Reason_treatment) = trimws(levels(prias$Reason_treatment))
prias$Reason_treatment[prias$Reason_treatment %in% c("N/A", "")] = NA
prias$Reason_treatment = droplevels(prias$Reason_treatment)

prias$day_discontinued = (prias$Date_discontinued - prias$Date_diagnosis)/(24*60*60)
prias$year_discontinued = prias$day_discontinued/365

#################################################
#What to do with these patients (wide format data set)???
#################################################

#Patients who have Nr_FUvisits = NA. All of these have no scores for psa, dre or gleason on any of the follow up
condition1 = is.na(prias$Nr_FUvisits)
View(prias[condition1,])

#Temporary treatment: Remove them, their baseline information is useful for PSA but no idea why they got censored
#few of them are also censored because of high gleason on repeat biopsy
prias=prias[!condition1,]
prias$P_ID = droplevels(prias$P_ID)

###########
#Patients who are reported NOT discontinued yet have reasons for discontinuation. No date of discontinuation available
condition2 = prias$DiscontinuedYesNo==0 & 
  (!is.na(prias$DiscontinuedType) | !is.na(prias$Reason_treatment))

View(prias[condition2, ])

#Temporary treatment: None of them seems to be having abnormal levels of PSA, gleason, DRE on the various follow ups. 
#Hence do not remove them, keep them instead with their reasons set to NA
prias[condition2,]$DiscontinuedType = NA
prias[condition2,]$Reason_treatment = NA

###########
#3 Patients who are reported YES discontinued but don't have date of discontinuation
condition3 = prias$DiscontinuedYesNo==1 & is.na(prias$Date_discontinued)

View(prias[condition3, ])
View(prias[condition3, ])

#Temporary treatment: Remove them because we do not know if their data is before or after treatment
prias = prias[!condition3,]
prias$P_ID = droplevels(prias$P_ID)

rm(condition1)
rm(condition2)
rm(condition3)
save.image("Rdata/Gleason as event/cleandata.Rdata")
#################################################
#################################################
# Convert wide to long and order by patient id and time
#################################################
#################################################
prias_long=reshape(prias, direction='long', idvar='P_ID', timevar = "visit_number",
        varying=list(c(3, 11, 50:84), c(3, 11, 264:298), c(4, 372, 15:49), c(5, 373, 120:154), 
                     c(8, 371, 85:119), c(8, 12, 299:333), c(8, 13, 334:368), 
                     c(369, 370, 159:193), c(6, 9, 194:228), c(7, 10, 229:263)),
        v.names=c('dom', 'domgleason', 'psa', 'dre', 'gleason', 'gleason1', 'gleason2', 'dummy', 'ncores', 'ncorespc'))
prias_long = prias_long[order(prias_long$P_ID, prias_long$dom, prias_long$domgleason, na.last = T), ]

save.image("Rdata/Gleason as event/cleandata.Rdata")

#################################################
#What to do with these patients/measurements (long format data set)???
#################################################
#Some measurements were dummy. i.e. patient didnt turn up. 
condition_long_1 = prias_long$dummy %in% 1

#solution: remove these
prias_long = prias_long[!condition_long_1,]

#When number of cores are 0/99 the gleason should be NA
condition_long_2 = prias_long$ncores %in% c(0, 99, 999, NA) | prias_long$ncorespc %in% c(99, 999, NA)
prias_long$gleason[condition_long_2] = NA
prias_long$gleason1[condition_long_2] = NA
prias_long$gleason2[condition_long_2] = NA

#ncores>0, ncorespc =0 and gleason > 0: possible because someone entered wrong ncorespc: no action required

#ncores>0, ncorespc>0, and Gleason = 0. Should be removed, I checked for the patients and it seems the scores should've been > 0
View(prias_long[!condition_long_2 & prias_long$ncorespc>0 & prias_long$gleason %in% 0,])
prias_long$gleason[!condition_long_2 & prias_long$ncorespc>0 & prias_long$gleason %in% 0] = NA
prias_long$gleason1[!condition_long_2 & prias_long$ncorespc>0 & prias_long$gleason %in% 0] = NA
prias_long$gleason2[!condition_long_2 & prias_long$ncorespc>0 & prias_long$gleason %in% 0] = NA

#Set all those dom to be NA where psa is NA, and set all those domgleason to NA where both dre & gleason are NA
#The dom for gleason and DRE is shared in the original data set
prias_long$psa[prias_long$psa %in% c(0, 9999, 1931, 3048)] = NA
prias_long$dom[is.na(prias_long$psa)] = NA
prias_long$domgleason[is.na(prias_long$dre) & is.na(prias_long$gleason)]= NA

#Remove all observations where both dom and domgleason are NA
condition_long_3 = is.na(prias_long$dom) & is.na(prias_long$domgleason)
View(prias_long[condition_long_3,])

prias_long = prias_long[!condition_long_3,]

#Further I have tested that gleason1 + gleason2 leads to overall gleason correctly.
View(prias_long[(!is.na(prias_long$gleason) & (prias_long$gleason != (prias_long$gleason1 + prias_long$gleason2))) %in% TRUE,])

####################################################
# Now we gotta merge and have a common date variable
####################################################
idList = unique(prias_long$P_ID)

ct= makeCluster(detectCores())
registerDoParallel(ct)
prias_long_onetime=foreach(i=1:length(idList),.combine='rbind') %dopar%{
  prias_long_i = prias_long[prias_long$P_ID==idList[[i]],]
  
  dom_i = prias_long_i$dom
  dom_gleason_i = prias_long_i$domgleason
  
  dom_index_i = order(prias_long_i$dom, na.last = T, decreasing = F)
  dom_gleason_index_i = order(prias_long_i$domgleason, na.last = T, decreasing = F)
  
  ind1=1
  ind2=1
  
  #Column 12 is Gleason column
  GLEASONDOM_COL = match("domgleason", colnames(prias_long_i))
  newds = prias_long_i[FALSE, -GLEASONDOM_COL]
  
  while(!is.na(dom_i[dom_index_i[ind1]]) | !is.na(dom_gleason_i[dom_gleason_index_i[ind2]])){
    
    #If both are not NA, then choose the one with min value and increase corresponding indicator
    d1 = dom_i[dom_index_i[ind1]]
    d2 = dom_gleason_i[dom_gleason_index_i[ind2]]
    
    dcommon = NA
    indexcommon = NA
    psacommon = drecommon = gleasoncommon = NA
   
    if(!is.na(d1) & !is.na(d2)){
      if(d1<d2){
        dcommon = d1
        indexcommon = dom_index_i[ind1]
        ind1 = ind1 + 1
        psacommon = prias_long_i[indexcommon,]$psa
      }else if(d1>d2){
        dcommon = d2
        indexcommon = dom_gleason_index_i[ind2]
        ind2 = ind2 + 1
        drecommon = prias_long_i[indexcommon,]$dre
        gleasoncommon = prias_long_i[indexcommon,]$gleason
      }else{ #d1==d2
        dcommon = d1
        indexcommon = dom_index_i[ind1]
        ind1 = ind1 + 1
        ind2 = ind2 + 1
        psacommon = prias_long_i[indexcommon,]$psa
        drecommon = prias_long_i[indexcommon,]$dre
        gleasoncommon = prias_long_i[indexcommon,]$gleason
      }
    }else if(!is.na(d1)){
      dcommon = d1
      indexcommon = dom_index_i[ind1]
      ind1 = ind1 + 1
      psacommon = prias_long_i[indexcommon,]$psa
    }else{ #(!is.na(d2))
      dcommon = d2
      indexcommon = dom_gleason_index_i[ind2]
      ind2 = ind2 + 1
      drecommon = prias_long_i[indexcommon,]$dre
      gleasoncommon = prias_long_i[indexcommon,]$gleason
    }
    
    row = prias_long_i[indexcommon, -GLEASONDOM_COL]
    row$dom = dcommon
    row$psa = psacommon
    row$dre = drecommon
    row$gleason = gleasoncommon
    
    newds = rbind(newds, row)
  }
  
  newds
}
stopCluster(ct)

#check the above result and then execute the statement below
prias_long = prias_long_onetime

prias_long$firstVisitDom = unlist(by(prias_long$dom, INDICES=prias_long$P_ID, FUN=function(x){rep(x[1], length(x))}))
prias_long$visitTimeDays = unlist(by(prias_long$dom, INDICES=prias_long$P_ID, FUN=function(x){(x-x[1])/(24*60*60)}))
prias_long$visitTimeYears = prias_long$visitTimeDays/365

#Cases where first two measurements have same DOM but different PSA scores. There are two such people
p_id_dom1dom2same = unique(prias_long$P_ID)[unlist(tapply(prias_long$visitTimeDays, prias_long$P_ID, function(x){sum(x==0, na.rm = T)>1}))]
View(prias_long[prias_long$P_ID %in% p_id_dom1dom2same,])

#In this regard subject with P_ID = 2173 has same PSA on both times but P_ID = 452 has a new PSA.
# So for the P_ID = 2173 we just remove second obs and for P_ID = 452, we replace first psa with second
prias_long = prias_long[!(prias_long$P_ID == 2173 & prias_long$visit_number==3),]

prias_long[prias_long$P_ID==452,]$psa[1] = prias_long[prias_long$P_ID==452,]$psa[2]
prias_long = prias_long[!(prias_long$P_ID == 452 & prias_long$visit_number==3),]

# Check if there are observations after the last visit Time. Only applicable to people who discontinued.
# In some cases it is just a few days after
# Nevertheless, I remove all observations after the treatment
# So keep those who never discontinued and those who discontinued but obs time < date of discontinuation
prias_long = prias_long[is.na(prias_long$Date_discontinued) | (prias_long$dom <= prias_long$Date_discontinued),]

prias_long$P_ID = droplevels(prias_long$P_ID)

#Check now if visittimes are sorted for every person. yes they are sorted
#FALSE: should be the answer
any(unlist(tapply(prias_long$visitTimeDays, prias_long$P_ID, function(x){is.unsorted(x, na.rm = T)})))

#Give new visit numbers, and adjust number of visits
prias_long$visit_number = unlist(by(prias_long, INDICES=prias_long$P_ID, FUN=function(x){1:nrow(x)}))
prias_long$nr_visits = unlist(by(prias_long, INDICES=prias_long$P_ID, FUN=function(x){rep(nrow(x), nrow(x))}))

#Choose only those patients who have atleast one gleason score
condition_survival_gleason2pid = droplevels(unique(prias_long[!is.na(prias_long$gleason),]$P_ID))
prias_long = prias_long[prias_long$P_ID %in% condition_survival_gleason2pid,]
prias_long$P_ID = droplevels(prias_long$P_ID)

#Create prias.id 
prias.id = prias_long[!duplicated(prias_long$P_ID),]
pid_progressed = droplevels(unique(prias_long[prias_long$gleason %in% c(7,8,9,10),]$P_ID))

#0 corresponds to right censored, 3 corresponds to interval censored. Check survival::Surv()
prias.id$progressed = ifelse(prias.id$P_ID %in% pid_progressed, 3, 0)
prias.id$progression_time_start = c(by(prias_long, prias_long$P_ID, function(prias_long_i){
  if(any(prias_long_i$gleason %in% c(7,8,9,10))){
    lastTime = prias_long_i$visitTimeYears[which(prias_long_i$gleason %in% c(7,8,9,10))[1]]
    tail(prias_long_i$visitTimeYears[!is.na(prias_long_i$gleason) & prias_long_i$visitTimeYears<=lastTime], 2)[1]
  }else{
    max(prias_long_i$visitTimeYears, na.rm = T)
  }
}))

prias.id$progression_time_end = c(by(prias_long, prias_long$P_ID, function(prias_long_i){
  if(any(prias_long_i$gleason %in% c(7,8,9,10))){
    prias_long_i$visitTimeYears[which(prias_long_i$gleason %in% c(7,8,9,10))[1]]  
  }else{
    Inf
  }
}))

#For the following special case of subjects, it may happen that there was only biopsy done and that gave a score > 6
#such patients may not have a biopsy at induction. they were inducted probably thinking that they are in a good condition
#Their DRE and PSA look normal afterall. So we assume that Gleason to be 0. If we don't do so then the
#start and end progression time may become equal due to the code written above. Following is the list
#of such patients and the corresponding hacky solution
#View(prias[prias$Gleason_sum %in% NA & prias$Gleason_sum_repeat %in% NA,])
prias.id$progression_time_start[prias.id$progression_time_start == prias.id$progression_time_end] = 0

#nr_visits is not corrupted coz we dropped entire subjects but not individual observations
prias_long$progression_time_start = rep(prias.id$progression_time_start, prias.id$nr_visits)
prias_long$progression_time_end = rep(prias.id$progression_time_end, prias.id$nr_visits)
prias_long$progressed = rep(prias.id$progressed, prias.id$nr_visits)

prias_long = prias_long[prias_long$visitTimeYears<=prias_long$progression_time_end,]
prias_long$nr_visits = unlist(by(prias_long, INDICES=prias_long$P_ID, FUN=function(x){rep(nrow(x), nrow(x))}))

#Keep only the columns of interest
prias_long = prias_long[, c("P_ID", "Age", 
                            "DiscontinuedYesNo", "Date_discontinued", "DiscontinuedType", 
                            "Reason_treatment", "nr_visits", "visit_number", 
                            "dom", "visitTimeDays", "visitTimeYears",
                            "progression_time_start", "progression_time_end", 
                            "progressed" , "psa", "dre", "gleason")]

prias.id = prias.id[, c("P_ID", "Age",  "DiscontinuedYesNo", "Date_discontinued", 
                        "DiscontinuedType", "Reason_treatment", "nr_visits",
                        "progressed", "progression_time_start", "progression_time_end")]

prias_long$log2psa = log(prias_long$psa, 2)
#########################################################
# Data type cleaning for the long version of the data set
#########################################################
levels(prias_long$dre) = trimws(levels(prias_long$dre))
prias_long$dre = as.ordered(droplevels(prias_long$dre))

rm(list = setdiff(ls(), c("prias.id", "prias_long")))

save.image("Rdata/Gleason as event/cleandata.Rdata")

# #How balanced is the data set. as in how often measurements are taken
# prias_long$diff_dom = unlist(lapply(prias$P_ID, FUN=function(id){
#   tempds = prias_long[prias_long$P_ID==id,]
#   numvisits = tempds$Nr_FUvisits[1]
#   if(is.na(numvisits)){ #699 such people exist
#     rep(NA, 35)
#   }else{
#     c(0,diff(tempds$dom[1:numvisits]), rep(NA, 35-numvisits))/(24*60*60)
#     #the extra 10 because there was an extra zero
#   }
# }))
# 
# 
# # Every biopsy for a patient is given a biopsy number
# prias_long$biopsy_number = unlist(lapply(prias$P_ID, FUN = function(id){
#   biopsy_indicator= as.numeric(prias_long[prias_long$P_ID==id, ]$gleason>0)
#   count = 1
#   for(i in 1:35){
#     if(!is.na(biopsy_indicator[i]) & biopsy_indicator[i]>0){
#       biopsy_indicator[i] <- count
#       count = count + 1
#     }else{
#       biopsy_indicator[i] = NA
#     }
#   }
#   biopsy_indicator
# }))
