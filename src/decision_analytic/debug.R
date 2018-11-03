by(scheduleResCombined, INDICES = scheduleResCombined$yearProgression, FUN = function(yearResCombined){
  by(yearResCombined, INDICES = yearResCombined$methodName, FUN = function(methodRes){
    c(mean(methodRes$nb), mean(methodRes$offset[methodRes$progression_time!=10]))
  })
})