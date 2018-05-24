library(foreign)

prias = read.spss(file.choose(), to.data.frame=TRUE)
prias_biopsy_cores = read.spss(file.choose(), to.data.frame = T)

prias$P_ID = droplevels(as.factor(prias$P_ID))
prias_biopsy_cores$P_ID = droplevels(as.factor(prias_biopsy_cores$P_ID))

setdiff(prias$P_ID, prias_biopsy_cores$P_ID)

#subject 5987 gotta be removed
prias = prias[!(prias$P_ID %in% c(5987)),]


prias_long=reshape(prias, direction='long', idvar='P_ID', timevar = "visit_number",
                   varying=list(c(85:119), 160:194),
                   v.names=c('gleason', 'dummy'))
prias_long = prias_long[,c("P_ID","visit_number", "gleason", "dummy")]
prias_long = prias_long[order(prias_long$P_ID, prias_long$visit_number, na.last = T), ]


prias_biopsy_cores_long=reshape(prias_biopsy_cores, direction='long', idvar='P_ID', timevar = "visit_number",
                   varying=list(c(2:36), c(37:71)),
                   v.names=c('cores_total', 'cores_cancer'))
prias_biopsy_cores_long = prias_biopsy_cores_long[order(prias_biopsy_cores_long$P_ID, prias_biopsy_cores_long$visit_number, na.last = T), ]


prias_long[, c("cores_total", "cores_cancer")] = prias_biopsy_cores_long[, c("cores_total", "cores_cancer")] 


#zero or 99 cores means ignore it
falseBiopsyNumbers = prias_long$cores_total %in% c(0, 99, NA)
table(prias_long[falseBiopsyNumbers,]$gleason)
table(prias_long[!falseBiopsyNumbers & prias_long$dummy==0,]$gleason)

View(prias_long[falseBiopsyNumbers,])
View(prias_long[!falseBiopsyNumbers,])
