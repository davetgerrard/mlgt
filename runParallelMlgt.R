#library(mlgt)

library(seqinr)
source("C:/Users/dave/HalfStarted/mlgt/mlgt.R")









## NEED to check there are no temporary files with common names. 
# fIdFile.txt
# fRawSeqExtract.fasta


my.design.list <- list()

for(thisMarker in names(myMarkerList))  {
	my.design.list[[thisMarker]] <- my.mlgt.Design
	my.design.list[[thisMarker]]@markers <- myMarkerList[thisMarker]
	
}

system.time(
my.result.list <- lapply(my.design.list, FUN=function(x) mlgt(x))
)

#A_E3 
#B_E2 
#DPA1_E2 
#DQA1_E2 
#   user  system elapsed 
#   4.58    0.32   43.57 




# then need script to pull this back together. Or use snowfall package.
# Another reason for a mergeMlgt() function may be to combine results across runs?

library(snowfall)

sfInit(parallel=TRUE, cpus=4, type="SOCK")

#sfExport("blastAllPath","fastacmdPath","formatdbPath","musclePath")
sfExport(list=ls())
#sfLibrary(mlgt)
sfLibrary(seqinr)


system.time(
sf.result.list <- sfLapply(my.design.list, mlgt)
)

#   user  system elapsed 
#   0.02    0.01   18.13



