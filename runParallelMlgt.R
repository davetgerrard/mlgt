#library(mlgt)

library(seqinr)
source("C:/Users/dave/HalfStarted/mlgt/mlgt.R")




setwd("C:/Users/dave/HalfStarted/mlgt/testProject/runParallel")
# mlgt.Design as per mlgt_README
# Load MIDs used to mark samples
fTagList <- read.fasta(system.file("data/namedBarcodes.fasta", package="mlgt"), 
			as.string=T) 
# here we're using the same tags at both ends of the amplicons.
rTagList <- fTagList
#The names of the samples
sampleList <- names(fTagList)
# Load the marker sequences. 
myMarkerList <- read.fasta(system.file("data/HLA_namedMarkers.fasta", package="mlgt"),
			as.string=T)	
# The fasta file of sequence reads
inputDataFile <- system.file("data/sampleSequences.fasta", package="mlgt")

my.mlgt.Design <- prepareMlgtRun(projectName="myProject", 
				runName="myRun", samples=sampleList, 
				markers=myMarkerList, fTags=fTagList, 
				rTags=rTagList, inputFastaFile=inputDataFile, 
				overwrite="yes")



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

#http://journal.r-project.org/2009-1/RJournal_2009-1_Knaus+et+al.pdf
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

# then join the results back into one. 
project.mlgt.Results <- combineMlgtResults(sf.result.list)



########### another test using split samples.

my.design.list <- list()

	my.design.list[['A']] <- my.mlgt.Design
	my.design.list[['A']]@samples <- sampleList[1:5]

	my.design.list[['B']] <- my.mlgt.Design
	my.design.list[['B']]@samples <- sampleList[6:10]

system.time(
sample.result.list <- sfLapply(my.design.list, mlgt)
)


sample.mlgt.Result <- combineMlgtResults(sample.result.list)




