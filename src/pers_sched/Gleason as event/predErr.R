# training_psa_data_set$progressed = ifelse(training_psa_data_set$progressed>0, 1, 0)
# training_psa_data_set$progression_time = ifelse(training_psa_data_set$progression_time_end==Inf,
#                                                training_psa_data_set$progression_time_start,
#                                                training_psa_data_set$progression_time_end)


predErrList_tdboth_rc = vector("list", 3)

i=1
for(Tstart in c(1.0001, 2.0001, 3.0001)){
  predErrList_tdboth_rc[[i]] = prederrJM(mvJoint_psa_spline_pt1pt54_pt1_tdboth_rc, 
                training_psa_data_set, 
                Tstart = Tstart, Thoriz=Tstart + 1.0101, idVar="P_ID")  
  i = i + 1
}

predErrList_tdVal_rc = vector("list", 3)

i=1
for(Tstart in c(1.0001, 2.0001, 3.0001)){
  predErrList_tdVal_rc[[i]] = prederrJM_mod(mvJoint_psa_spline_pt1pt54_pt1_tdval_rc, 
                                          training_psa_data_set, 
                                          Tstart = Tstart, Thoriz=Tstart + 1.0101, idVar="P_ID")  
  i = i + 1
}

predErrList_rc = list(predErrList_tdboth_rc=predErrList_tdboth_rc, 
                   predErrList_tdVal_rc=predErrList_tdVal_rc)
save(predErrList_rc, file="Rdata/predErrList_rc.Rdata")
