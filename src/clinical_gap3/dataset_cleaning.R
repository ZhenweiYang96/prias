library(foreign)

prias=read.spss("/home/a_tomer/Documents/ErasmusMC_datasets/PRIAS-2019/2019-04-10_dump.sav",
          to.data.frame = T)

usefulCols = c(1,2,3,4,6,7,8,14,15,16, 21,22,23,24,25,32,33,34,35,36,38,
                   63:200, 201:292, 385:430, 707:844,
                   1075:1120, 1213:1258, 1949:1994, 1998,1999, 2000, 2001,
                   2002, 2003, 2021, 2022, 2023)

colnames(prias[,usefulCols])

prias = prias[,usefulCols]
prias$Gleason_sum_2 = prias$Gleason1_2 + prias$Gleason2_2
prias$N_A = NA
colnames(prias)

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

prias_long = prias_long[order(prias_long$P_ID, prias_long$dom_psa, prias_long$dom_dre_gleason, na.last = T), ]
