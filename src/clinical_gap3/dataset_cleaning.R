library(foreign)

YEAR_DIVIDER = 24 * 60 * 60 * 365

prias=read.spss("/home/a_tomer/Documents/ErasmusMC_datasets/PRIAS-2019/2019-04-10_dump.sav",
          to.data.frame = T)

usefulCols = c(1,2,3,4,6,7,8,14,15,16, 21,22,23,24,25,32,33,34,35,36,38,
                   63:200, 201:292, 385:430, 707:844,
                   1075:1120, 1213:1258, 1949:1994, 1998,1999, 2000, 2001,
                   2002, 2003, 2021, 2022, 2023)

colnames(prias[,usefulCols])

prias = prias[,usefulCols]
prias$Date_birth = (prias$Date_dianosis - prias$Date_birth)/YEAR_DIVIDER
colnames(prias)[2] = "age"
prias$Gleason_sum_2 = prias$Gleason1_2 + prias$Gleason2_2
prias$N_A = NA

#Remove patients with NA, and age less than 30 (includes negative ages)
prias = prias[!is.na(prias$age) & prias$age > 30,]

#16 patients had a repeat Gleason measurement
gl2_filter = !is.na(prias$Gleason_sum_2) & !is.na(prias$Date_dianosis2) & 
  prias$Gleason_sum_2 > 0
sum(gl2_filter)

#lets check these 16 patients Gleason first and second time
prias$Gleason_sum[gl2_filter]
prias$Gleason_sum_2[gl2_filter]
prias$Num_cores2[gl2_filter]
prias$Num_Cores_PC2[gl2_filter]
(prias$Date_dianosis2[gl2_filter] - prias$Date_dianosis[gl2_filter])/YEAR_DIVIDER
prias$P_ID[gl2_filter]

#No patient had a Gleason downgrade upon rebiopsy, but it seems that
#the 4th patient in this list should never have been there in the first place
prias = prias[prias$P_ID != 317,]

#No need to worry about Gleason_sum_2 again

#Lets remove all patients who have a Gleason > 6 at first visit
prias = prias[!(prias$Gleason_sum %in% c(7,8,9,10)),]

#The column 584 is a dummy column to handle Gleason2 which is 
#gleason taken again at the time of diagnosis
prias_long=reshape(prias, direction='long', idvar='P_ID', timevar = "visit_number",
                   varying=list(c(3, 584, 68:113), c(4, 584, 22:67), 
                                c(3, 18, 252:297), c(5, 584, 114:159), c(5, 584, 528:573),
                                c(8, 19, 298:343), c(9, 20, 344:389), c(10, 583, 390:435), 
                                c(6, 16, 160:205), c(7, 17, 206:251), 
                                c(584, 584, 436:481), c(21, 584, 482:527)),
                   v.names=c('dom_psa', 'psa', 
                             'dom_dre_gleason', 'dre', 'dre_recode',
                             'gleason_1', 'gleason_2', 'gleason_sum', 
                             'ncores', 'ncores_pc', 
                             'dummy', 'active'))

prias_long = prias_long[order(prias_long$P_ID), ]

#First lets check Gleason table
table(prias_long$gleason_sum, useNA = "always")

#Didn't find anything negative. Now lets check ncores and ncores_pc
table(prias_long$ncores, useNA = "always")
table(prias_long$ncores_pc, useNA = "always")

#Some cores are NA, 99, 999 or 0, and some cores_pc are 99 or NA. 
#Lets set those Gleason and cores as NA
temp_filter = prias_long$ncores %in% c(NA, 99, 999, 0) | prias_long$ncores_pc %in% c(NA, 99)
prias_long$gleason_1[temp_filter] = NA
prias_long$gleason_2[temp_filter] = NA
prias_long$gleason_sum[temp_filter] = NA
prias_long$ncores[temp_filter] = NA
prias_long$ncores_pc[temp_filter] = NA

#Sometimes ncores with pc are more than total number of cores taken.
#Lets set those Gleason and cores as NA
summary(prias_long$ncores - prias_long$ncores_pc)

temp_filter = is.na(prias_long$ncores - prias_long$ncores_pc) | (prias_long$ncores - prias_long$ncores_pc) < 0
prias_long$gleason_1[temp_filter] = NA
prias_long$gleason_2[temp_filter] = NA
prias_long$gleason_sum[temp_filter] = NA
prias_long$ncores[temp_filter] = NA
prias_long$ncores_pc[temp_filter] = NA

#Lets check first row data, we have 48 rows per patient
firstRowLongFilter = seq(1, nrow(prias_long), by = 48)

#Set Gleason = 6 for all patients on row 1
prias_long$gleason_sum[firstRowLongFilter] = 6

# add a date of diagnosis in long
prias_long$dom_diagnosis = rep(prias_long$dom_dre_gleason[firstRowLongFilter], each=48)

prias_long$year_dre_gleason = (prias_long$dom_dre_gleason - prias_long$dom_diagnosis)/YEAR_DIVIDER
prias_long$year_psa = (prias_long$dom_psa - prias_long$dom_diagnosis)/YEAR_DIVIDER
prias_long$year_surgery = (prias_long$Surgerydate - prias_long$dom_diagnosis)/YEAR_DIVIDER
prias_long$year_dead = (prias_long$Dead_Date - prias_long$dom_diagnosis)/YEAR_DIVIDER
prias_long$year_discontinued = (prias_long$Date_discontinued - prias_long$dom_diagnosis)/YEAR_DIVIDER

#Set dom as NA for patient data based on implausible dates of measurements
#2019 june is 13.5 years since PRIAS started in 2006
MAX_FOLLOW_UP_YEARS = 13.5

temp_filter = is.na(prias_long$year_psa) | prias_long$year_psa < 0 | prias_long$year_psa > MAX_FOLLOW_UP_YEARS
prias_long$psa[temp_filter] = NA
prias_long$dom_psa[temp_filter] = NA
prias_long$year_psa[temp_filter] = NA

#same we do for date of DRE and Gleason
temp_filter = is.na(prias_long$year_dre_gleason) | prias_long$year_dre_gleason < 0 | prias_long$year_dre_gleason > MAX_FOLLOW_UP_YEARS
prias_long$dre[temp_filter] = NA
prias_long$dre_recode[temp_filter] = NA
prias_long$dom_dre_gleason[temp_filter] = NA
prias_long$year_dre_gleason[temp_filter] = NA
prias_long$gleason_1[temp_filter] = NA
prias_long$gleason_2[temp_filter] = NA
prias_long$gleason_sum[temp_filter] = NA
prias_long$ncores[temp_filter] = NA
prias_long$ncores_pc[temp_filter] = NA

#Now some PSA are too big to be true, and some are 99, none are negative
summary(prias_long$psa)
sort(prias_long$psa, decreasing = T)[1:50]

temp_filter = is.na(prias_long$psa) | prias_long$psa > 95
prias_long$psa[temp_filter] = NA
prias_long$dom_psa[temp_filter] = NA
prias_long$year_psa[temp_filter] = NA

#Now lets focus on DRE
table(prias_long$dre, useNA = "always")
table(prias_long$dre_recode, useNA = "always")

#So do not use DRE recode as it has too many NA.
#Now remove whitespace from dre levels, and all empty should be set to NA
levels(prias_long$dre) = trimws(levels(prias_long$dre))
prias_long$dre[prias_long$dre==""] = NA
prias_long$dre = droplevels(prias_long$dre)
table(prias_long$dre, useNA = "always")

#Now lets focus on dummy indicator, dont think of active indicator
prias_long$dom_psa[prias_long$dummy==1] = NA
prias_long$psa[prias_long$dummy==1] = NA
prias_long$dom_psa[prias_long$dummy==1] = NA
prias_long$dom_dre_gleason[prias_long$dummy==1] = NA
prias_long$dre[prias_long$dummy==1] = NA
prias_long$dre_recode[prias_long$dummy==1] = NA
prias_long$gleason_1[prias_long$dummy==1] = NA
prias_long$gleason_2[prias_long$dummy==1] = NA
prias_long$gleason_sum[prias_long$dummy==1] = NA
prias_long$ncores[prias_long$dummy==1] = NA
prias_long$ncores_pc[prias_long$dummy==1] = NA
prias_long$year_dre_gleason[prias_long$dummy==1] = NA
prias_long$year_psa[prias_long$dummy==1] = NA

#Now lets focus on gleason reclassification and risk reclassification
#What do we do with Gleason equal to 0?
table(prias_long$gleason_sum)

table(prias_long$ncores[prias_long$gleason_sum %in% 0], useNA = 'always')
table(prias_long$ncores_pc[prias_long$gleason_sum %in% 0], useNA = 'always')
