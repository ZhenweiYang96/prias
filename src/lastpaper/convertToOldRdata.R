files = list.files("Rdata/lastpaper/simulation/light/", full.names = T)

sapply(files, function(file){
  load(file)
  save(jointModelData, file = file, version = 2)
})
