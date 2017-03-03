library(foreign)

##############################################
# Load the data set
##############################################
prias = read.spss(file.choose(), to.data.frame=TRUE)
prias_biopsy_cores = read.spss(file.choose(), to.data.frame = T)

prias$P_ID = droplevels(as.factor(prias$P_ID))
prias_biopsy_cores$P_ID = droplevels(as.factor(prias_biopsy_cores$P_ID))

setdiff(prias$P_ID, prias_biopsy_cores$P_ID)
#Patient 5987 is not in prias_biopsy_cores
prias = prias[!(prias$P_ID %in% c(5987)),]
prias$P_ID = droplevels(as.factor(prias$P_ID))

#Combine the two data sets, but skip the ID column
#check first if the two data sets are sorted, result should be FALSE
any(prias$P_ID != prias_biopsy_cores$P_ID)
prias = cbind(prias, prias_biopsy_cores[,-1])

######################################################
#Now we start fixing the data set
######################################################
colnames(prias)[match("Discontinued", colnames(prias))] = "DiscontinuedType"
prias$REASO0 = NULL
prias$Dummy.0 = rep(0, nrow(prias))
prias$Dummy.0.repeat = rep(0, nrow(prias))

long_columns = c(c(15:49), c(50:84), c(85:119), c(120:154), c(264, 265, 159:193), 
                 c(194:228), c(229:263))
psa_gleason_dre_columns = c(15:49, 85:119, 120:154)

##############################################
# Remove subjects
##############################################
#15 patients do not have any information and they all have Age=NA
empty_patients = prias[is.na(prias$Age),]

#Some patients are too young to be real. So basically remove both empty ones and these patients too
prias = prias[!is.na(prias$Age) & prias$Age > 5,]

##############################################
# Repeat biopsy score update
##############################################
#Original biopsy score is invalid
prias$Gleason_sum_repeat = prias$Gleason1_2 + prias$Gleason2_2
prias[!is.na(prias$Gleason_sum_repeat) & prias$Gleason_sum_repeat>0, ]$Gleason_sum = NA

#there are still some patients who had second biopsy but the date of their second biopsy was missing
#we cannot use their  biopsy scores either ar baseline or at second measurement. Don't remove them
# we will use rest of their information
View(prias[!is.na(prias$Gleason1_2 + prias$Gleason2_2) & 
             (prias$Gleason1_2 + prias$Gleason2_2)>0 & 
             is.na(prias$Date_dianosis2), ])

#These two are added just to balance the data set before converting to long
prias$psa_repeat = rep(NA, nrow(prias))
prias$dre_repeat = rep(NA, nrow(prias))

##############################################
# Data type cleaning before conversion to long
##############################################
prias$P_ID = droplevels(as.factor(prias$P_ID))

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
View(prias[condition1, psa_gleason_dre_columns])

#Temporary treatment: Remove them
prias=prias[!condition1,]

###########
#Patients who are reported NOT discontinued yet have reasons for discontinuation. No date of discontinuation available
condition2 = prias$DiscontinuedYesNo==0 & 
  (!is.na(prias$DiscontinuedType) | !is.na(prias$Reason_treatment))

View(prias[condition2, -long_columns])
View(prias[condition2, c(1,psa_gleason_dre_columns)])

#Temporary treatment: None of them seems to be having abnormal levels of PSA, gleason, DRE on the various follow ups. 
#Hence do not remove them, keep them instead with their reasons set to NA
prias[condition2,]$DiscontinuedType = NA
prias[condition2,]$Reason_treatment = NA

###########
#Patients who are reported YES discontinued but don't have date of discontinuation
condition3 = prias$DiscontinuedYesNo==1 & is.na(prias$Date_discontinued)

View(prias[condition3, -long_columns])
View(prias[condition3, c(1,psa_gleason_dre_columns)])

#Temporary treatment: Remove them because we do not know if their data is before or after treatment
prias = prias[!condition3,]

###########
#Patients who discontinued but neither have DiscontinuedType nor have reason of treatment
condition4 = prias$DiscontinuedYesNo==1 & 
  is.na(prias$DiscontinuedType) & is.na(prias$Reason_treatment)

View(prias[condition4, -long_columns])
View(prias[condition4, c(1,psa_gleason_dre_columns)])

#Temporary treatment: Do not remove them as event of interest is Gleason > 6
#prias=prias[!condition4,]

###########
#Patients who discontinued YES, don't have DiscontinuedType, yet they do have reason of treatment
#Among these those reasons which are "Based on protocol advice" and "Other"
condition5 = prias$DiscontinuedYesNo==1 & 
  is.na(prias$DiscontinuedType) & prias$Reason_treatment %in% c("Based on protocol advice", "Other")

View(prias[condition5, -long_columns])
View(prias[condition5, c(1,psa_gleason_dre_columns)])

#Temporary treatment: Do not remove them as event of interest is Gleason > 6
#prias=prias[!condition5,]

prias$P_ID = droplevels(prias$P_ID)

#################################################
#################################################
# Convert wide to long and order by patient id and time
#################################################
#################################################
prias_long=reshape(prias, direction='long', idvar='P_ID', timevar = "visit_number",
        varying=list(c(3, 11, 50:84), c(4, 197, 15:49), c(5, 198, 120:154), c(8, 196, 85:119), c(194, 195, 159:193)),
        v.names=c('dom', 'psa', 'dre', 'gleason', 'dummy'))
prias_long = prias_long[order(prias_long$P_ID, prias_long$dom, na.last = T), ]
prias_long$dummy = factor(prias_long$dummy, labels = c("No", "Yes"))

prias_long$firstVisitDom = unlist(by(prias_long$dom, INDICES=prias_long$P_ID, FUN=function(x){rep(x[1], length(x))}))
prias_long$visitTimeDays = unlist(by(prias_long$dom, INDICES=prias_long$P_ID, FUN=function(x){(x-x[1])/(24*60*60)}))
prias_long$visitTimeYears = prias_long$visitTimeDays/365

#################################################
#What to do with these patients/measurements (long format data set)???
#################################################
#Some measurements were dummy. i.e. patient didnt turn up. 
condition_long_1 = prias_long$dummy %in% "Yes"

#solution: remove these
prias_long = prias_long[!condition_long_1,]

#Gleason scores which are 0 should be counted as NA
condition_long_2 = prias_long$gleason %in% c(0)
prias_long$gleason[condition_long_2] = NA

#Cases where psa/gleason/dre is available but date of measurement is missing
condition_long_3 = is.na(prias_long$dom) & (!is.na(prias_long$psa) | !is.na(prias_long$dre) |
                                              !is.na(prias_long$gleason))
View(prias_long[condition_long_3,])
#Temporary solution: remove them for now, otherwise take the date it should have actually been measured
prias_long = prias_long[!condition_long_3,]

#Also remove every measurement with no date of measurement (these can be real + we introduced some of these above)
prias_long = prias_long[!is.na(prias_long$dom),]

#Also remove every measurement with non NA date of measurement, but all 3: psa, dre, gleason are NA
#Most probably these are introduced by us as all of them are second visit...i.e. visit of gleason repeat
condition_long_4 = !is.na(prias_long$dom) & is.na(prias_long$psa) & is.na(prias_long$dre) & is.na(prias_long$gleason)
View(prias_long[condition_long_4,])
prias_long = prias_long[!condition_long_4,]

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

#First: Choose only those patients who have Gleason <=6 at baseline
condition_survival_gleason1 = prias_long$visit_number==1 & prias_long$gleason %in% c(7,8,9,10)
View(prias_long[condition_survival_gleason1,])

prias_long = prias_long[!(prias_long$P_ID %in% droplevels(prias_long[condition_survival_gleason1, "P_ID"])) ,]
prias_long$P_ID = droplevels(prias_long$P_ID)

#Second: Choose only those patients who have atleast one gleason score
condition_survival_gleason2pid = droplevels(unique(prias_long[!is.na(prias_long$gleason),]$P_ID))
prias_long = prias_long[prias_long$P_ID %in% condition_survival_gleason2pid,]
prias_long$P_ID = droplevels(prias_long$P_ID)

#Third:  drop all observations after the first time Gleason score is > 7
prias.id = prias_long[!duplicated(prias_long$P_ID),]
pid_progressed = droplevels(unique(prias_long[prias_long$gleason %in% c(7,8,9,10),]$P_ID))

prias.id$progressed = ifelse(prias.id$P_ID %in% pid_progressed, 1, 0)
prias.id$progression_time = by(prias_long, prias_long$P_ID, function(prias_long_i){
  if(any(prias_long_i$gleason %in% c(7,8,9,10))){
    prias_long_i$visitTimeYears[which(prias_long_i$gleason %in% c(7,8,9,10))[1]]  
  }else{
    max(prias_long_i$visitTimeYears, na.rm = T)
  }
})

prias_long$progression_time = unlist(lapply(prias.id$P_ID, function(pid){
  rowCount = nrow(prias_long[prias_long$P_ID == pid,])
  rep(prias.id[prias.id$P_ID==pid,"progression_time"][1], rowCount)
}))

prias_long = prias_long[prias_long$visitTimeYears<=prias_long$progression_time,]
prias_long$nr_visits = unlist(by(prias_long, INDICES=prias_long$P_ID, FUN=function(x){rep(nrow(x), nrow(x))}))

#Keep only the columns of interest
prias_long = prias_long[, c("P_ID", "Age", 
                            "DiscontinuedYesNo", "Date_discontinued", "DiscontinuedType", 
                            "Reason_treatment", "nr_visits", "visit_number", 
                            "dom", "visitTimeDays", "visitTimeYears",
                            "progression_time", "psa", "dre", "gleason")]

prias.id = prias.id[, c("P_ID", "Age",  "DiscontinuedYesNo", "Date_discontinued", 
                        "DiscontinuedType", "Reason_treatment", "nr_visits",
                        "progressed", "progression_time")]

#########################################################
# Data type cleaning for the long version of the data set
#########################################################
levels(prias_long$dre) = trimws(levels(prias_long$dre))
prias_long$dre = as.ordered(droplevels(prias_long$dre))

#####################################################
# Create simplified versions of various columns
#####################################################
prias_long$didre = rep(NA, nrow(prias_long))
prias_long[!is.na(prias_long$dre), ]$didre = ifelse(prias_long[!is.na(prias_long$dre), ]$dre %in% c("T1c"), yes = "T1", no = ">T1")
prias_long$didre = ordered(prias_long$didre, levels=c("T1", ">T1"))

prias_long$digleason = rep(NA, nrow(prias_long))
prias_long[!is.na(prias_long$gleason), ]$digleason = ifelse(prias_long[!is.na(prias_long$gleason), ]$gleason<=6, yes="Low", no="High")
prias_long$digleason = ordered(prias_long$digleason, levels=c("Low", "High"))

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
