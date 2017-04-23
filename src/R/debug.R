last_biopsy_time = -Inf
for(j in minVisits:(timesPerSubject-1)){
  biopsy_gap_needed = (patientDs_i$visitTimeYears[j] - last_biopsy_time) < 1
  print(patientDs_i$survTimeYouden[j] - patientDs_i$visitTimeYears[j])
  if(biopsy_gap_needed==FALSE & (patientDs_i$survTimeYouden[j] - patientDs_i$visitTimeYears[j]) <= 1){
    nb2 = nb2 + 1
    biopsytimeOffset2 = patientDs_i$survTimeYouden[j] - trueProgressionTime
    last_biopsy_time = patientDs_i$survTimeYouden[j]
    
    if(biopsytimeOffset2 >= 0){
      break
    }
  }
}