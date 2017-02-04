trt_patients = prias_long[prias_long$event_type %in% c("Anxiety/Other/Request"),]
trt_patients$P_ID = droplevels(trt_patients$P_ID) 

total = length(unique(trt_patients$P_ID))

nonempty_gleason = trt_patients[!is.na(trt_patients$gleason),]
nonempty_dre = trt_patients[!is.na(trt_patients$dre),]

first_gleason = c(unlist(by(nonempty_gleason$gleason, INDICES=list(nonempty_gleason$P_ID), FUN=function(x){x[1]}))) + 1
qplot(x=factor(first_gleason), geom="histogram")

last_gleason = c(unlist(by(trt_patients$gleason, INDICES=list(trt_patients$P_ID), FUN=function(x){x[length(x)]}))) + 1
qplot(x=last_gleason, geom="histogram") + scale_x_continuous(breaks = seq(2, 10, by = 1))

max_gleason = c(unlist(by(nonempty_gleason$gleason, INDICES=list(nonempty_gleason$P_ID), FUN=function(x){max(x)}))) + 1
qplot(x=max_gleason, geom="histogram")

last_dre = c(unlist(by(trt_patients$dre, INDICES=list(trt_patients$P_ID), FUN=function(x){x[1]})))
table((unlist(by(nonempty_dre$dre, INDICES=list(nonempty_dre$P_ID), FUN=function(x){x[length(x)]}))))

#####################################
pids = unique(trt_patients$P_ID)

ppcCheck=foreach(i=1:length(pids),.combine='c') %dopar%{
  x = trt_patients[trt_patients$P_ID %in% pids[i],]
  k = nrow(x)-1
  while(k > 0){
    if((x$visitTimeYears[nrow(x)] - x$visitTimeYears[k]) >= 3){
      if(x$psa[nrow(x)] > 2*x$psa[k]){
        return(1)
      }else{
        return(0)
      }
    }
    
    k = k-1
  }
  
  return(0)
}


psadt = by(trt_patients, INDICES=list(trt_patients$P_ID), FUN=function(x){
    k = nrow(x)-1
    while(k > 0){
      if((x$visitTimeYears[nrow(x)] - x$visitTimeYears[k]) >= 3){
        rr = range(x$psa[k:nrow(x)], na.rm = T)
        if(rr[2] > 2*rr[1]){
          return(1)
        }else{
          return(0)
        }
      }
      
      k = k-1
    }
    
    return(NA)
  })


