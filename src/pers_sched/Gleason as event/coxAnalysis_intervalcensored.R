prias.id$progression_time_start[prias.id$progression_time_start == 0] = 10e-3

training.prias.id = prias.id[!(prias.id$P_ID %in% c(3174, 2340, 911)),]

###########################################################################
survModel = survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                       I(Age - 70) +  I((Age - 70)^2), data = prias.id, model = TRUE)
survModel.training= survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                              I(Age - 70) +  I((Age - 70)^2), data = training.prias.id, model = TRUE)

################################################################################
prias.id.rightCens = prias.id
prias.id.rightCens$progressed = ifelse(prias.id$progressed>0, 1, 0)
prias.id.rightCens$progression_time = ifelse(prias.id$progression_time_end==Inf, prias.id$progression_time_start, prias.id$progression_time_end)

survModel_rightCens = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2), 
                                    data=prias.id.rightCens, x = T, model = T)

training.prias.id.rightCens = training.prias.id
training.prias.id.rightCens$progressed = ifelse(training.prias.id$progressed>0, 1, 0)
training.prias.id.rightCens$progression_time = ifelse(training.prias.id$progression_time_end==Inf, 
                                                      training.prias.id$progression_time_start, 
                                                      training.prias.id$progression_time_end)

survModel.training_rightCens = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2), 
                            data=training.prias.id.rightCens, x = T, model = T)
