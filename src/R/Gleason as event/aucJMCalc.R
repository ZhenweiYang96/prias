aucList_tdboth_rc = vector("list", 3)

i=1
for(Tstart in c(1.0001, 2.0001, 3.0001)){
  aucList_tdboth_rc[[i]] = aucJM(mvJoint_psa_spline_pt1pt54_pt1_tdboth_rc,
                                          training_psa_data_set[training_psa_data_set$P_ID %in% training.prias.id$P_ID[1:3500],],
                                          Tstart = Tstart, Thoriz=Tstart + 1.0101, idVar="P_ID")
  i = i + 1
}

aucList_tdval_rc = vector("list", 3)

i=1
for(Tstart in c(1.0001, 2.0001, 3.0001)){
  aucList_tdval_rc[[i]] = aucJM(mvJoint_psa_spline_pt1pt54_pt1_tdval_rc,
                                training_psa_data_set[training_psa_data_set$P_ID %in% training.prias.id$P_ID[1:3500],],
                                         Tstart = Tstart, Thoriz=Tstart + 1.0101, idVar="P_ID")
  i = i + 1
}

aucList = list(aucList_tdboth_rc=aucList_tdboth_rc, aucList_tdval_rc=aucList_tdval_rc)
save(aucList, file="Rdata/aucList.Rdata")
