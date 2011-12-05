


###################### GET INFO
library(seqinr)

formatdbPath <- "C:/Users/Public/Apps/Blast/bin/formatdb.exe"
fastacmdPath <- "C:/Users/Public/Apps/Blast/bin/fastacmd.exe"
blastAllPath <- "C:/Users/Public/Apps/Blast/bin/blastall.exe"
musclePath <- "C:/Users/Public/Apps/Muscle/muscle3.8.31_i86win32.exe"





analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/testRun2"
setwd( analysisDir )

sampleFile <- "C:/Users/dave/HalfStarted/mlgt/samples_ALL.tab"
sampleTable <- read.delim(sampleFile, sep=";")

# markers
markerFile <- "C:/Users/dave/HalfStarted/mlgt/HLA_MARKERS.fasta"
markerList <- read.fasta(markerFile ,as.string=T)

# samples
sampleList <- as.character(unique(sampleTable$sample_name))

# MIDs/tags
tempTable <- unique(subset(sampleTable, select=c(tag_f, sample_name)))
fTagList <- as.list(as.character(tempTable$tag_f))
names(fTagList) <- tempTable$sample_name
tempTable <- unique(subset(sampleTable, select=c(tag_r, sample_name)))
rTagList <- as.list(as.character(tempTable$tag_r))
names(rTagList) <- tempTable$sample_name

# data in fasta file.
inputDataFile <- "C:/Users/dave/HalfStarted/mlgt/3.TCA.454Reads.fna"

##################### SET UP mlgt
source("C:/Users/dave/HalfStarted/mlgt/mlgt_classes.R")


thisDesign <- new("mlgtDesign", projectName="testProject", runName="testRun2", 
				samples=sampleList, markers=markerList,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(thisDesign)

prepareMlgtRun(thisDesign)



####################### RUN mlgt

thisResult <- mlgt(thisDesign)

#thisResult <- new("mlgtResult", thisDesign )



shortTest <- 

shortDesign <- new("mlgtDesign", projectName="testProject", runName="testRun2", 
				samples=sampleList[1:3], markers=markerList[7:9],
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

thisResult <- mlgt(shortDesign)


################Concat data
analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/testRun3"
setwd( analysisDir )

# data in fasta file.
inputDataFile <- "C:/Users/dave/HalfStarted/mlgt/catLiverpoolNewcastleData.fasta"
catDesign <- new("mlgtDesign", projectName="testProject", runName="ConcatData", 
				samples=sampleList, markers=markerList,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(catDesign )

prepareMlgtRun(catDesign )

catResult <- mlgt(catDesign )

save(catResult, file="catResult.RData")


##################

# Liverpool data alone
#C:/Users/dave/HLA/data/ID632_FM_preliminary_ALL.fasta

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/testRun4"
setwd( analysisDir )

# data in fasta file.
inputDataFile <- "C:/Users/dave/HLA/data/ID632_FM_preliminary_ALL.fasta"
liverDesign <- new("mlgtDesign", projectName="testProject", runName="LiverpoolData", 
				samples=sampleList, markers=markerList,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(liverDesign )

prepareMlgtRun(liverDesign )

liverResult <- mlgt(liverDesign )

save(liverResult, file="liverResult.RData")

##################

# NCL Run_760591

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/Run_760591"
setwd( analysisDir )

# data in fasta file.
#C:\Users\dave\HLA\data\NCL\Run_760591
inputDataFile <- "C:/Users/dave/HLA/data/NCL/Run_760591/4.TCA.454Reads.fna"
Run_760591.Design <- new("mlgtDesign", projectName="testProject", runName="Run_760591", 
				samples=sampleList, markers=markerList,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(Run_760591.Design )

prepareMlgtRun(Run_760591.Design)

Run_760591.Result <- mlgt(Run_760591.Design)

save(Run_760591.Result, file="Run_760591.RData")

##################

# NCL Run_756948

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/Run_756948"
setwd( analysisDir )

# data in fasta file.
#C:\Users\dave\HLA\data\NCL\Run_756948
inputDataFile <- "C:/Users/dave/HLA/data/NCL/Run_756948/3.TCA.454Reads.fna"
Run_756948.Design <- new("mlgtDesign", projectName="testProject", runName="Run_756948", 
				samples=sampleList, markers=markerList,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(Run_756948.Design )

prepareMlgtRun(Run_756948.Design)

Run_756948.Result <- mlgt(Run_756948.Design)

save(Run_756948.Result, file="Run_756948.RData")


#################

alleleDbObject <- new(alleleDb,...)