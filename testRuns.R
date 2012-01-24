#TODO
# Build in checking for dodgy alignments 
#	sub aligns longer than marker?
#	multiple '-' near ends of marker seq.
#	Compare to profile alignment?
#	Use other aligments software?
#	Replace ends of marker with NNN before alignment?
# Create exon-only marker sequences and re-run with this.  DONE
# Do another 'clean run' from an empty starting environment.
# Run new data through.
# Batch run all samples to date plus concat sample.





###################### GET INFO
library(seqinr)

formatdbPath <- "C:/Users/Public/Apps/Blast/bin/formatdb.exe"
fastacmdPath <- "C:/Users/Public/Apps/Blast/bin/fastacmd.exe"
blastAllPath <- "C:/Users/Public/Apps/Blast/bin/blastall.exe"
musclePath <- "C:/Users/Public/Apps/Muscle/muscle3.8.31_i86win32.exe"

source("C:/Users/dave/HalfStarted/mlgt/mlgt.R")



############## cleanRun


analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/cleanRun"
setwd( analysisDir )

sampleFile <- "C:/Users/dave/HalfStarted/mlgt/samples_ALL.tab"
sampleTable <- read.delim(sampleFile, sep=";")

intersectMarkerList <- read.fasta("C:/Users/dave/HLA/data/alleleSeqs/HLA_intersectMarkersDec11.fasta", as.string=T)

# samples
sampleList <- as.character(unique(sampleTable$sample_name))

# MIDs/tags
tempTable <- unique(subset(sampleTable, select=c(tag_f, sample_name)))
fTagList <- as.list(as.character(tempTable$tag_f))
names(fTagList) <- tempTable$sample_name
tempTable <- unique(subset(sampleTable, select=c(tag_r, sample_name)))
rTagList <- as.list(as.character(tempTable$tag_r))
names(rTagList) <- tempTable$sample_name

## Use smallest NCL data for a quick test run.
inputDataFile <- "C:/Users/dave/HLA/data/NCL/Run_763932_Modified/4.TCA.454Reads.fna"
intersect.cleanRun.Design <- new("mlgtDesign", projectName="testProject", runName="cleanRun", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(intersect.cleanRun.Design )

prepareMlgtRun(intersect.cleanRun.Design)		# about 2 mins

intersect.cleanRun.Result <- mlgt(intersect.cleanRun.Design)		# about 4 mins, normally 30mins+

str(intersect.cleanRun.Result, max.level=2)  # could use for print.mlgtResult()
save(intersect.cleanRun.Result, file="intersect.cleanRun.RData")

# need to create/load knownAlleleDb first.
#test.genotypes <- callGenotypes.mlgtResult(intersect.cleanRun.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)
#write.table(test.genotypes, file="intersect.cleanRun.genotypesAlleles.tab",  quote=F, sep="\t", row.names=F)








########## RUNS

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/testRun2"
setwd( analysisDir )

sampleFile <- "C:/Users/dave/HalfStarted/mlgt/samples_ALL.tab"
sampleTable <- read.delim(sampleFile, sep=";")

# markers
markerFile <- "C:/Users/dave/HalfStarted/mlgt/HLA_MARKERS.fasta"
markerList <- read.fasta(markerFile ,as.string=T)

intersectMarkerList <- read.fasta("C:/Users/dave/HLA/data/alleleSeqs/HLA_intersectMarkersDec11.fasta", as.string=T)

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
analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/testRun2"
setwd( analysisDir )




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
#catDesign <- prepareMlgtRun(projectName="testProject", runName="ConcatData", 
#				samples=sampleList, markers=intersectMarkerList,
#				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )


catDesign <- prepareMlgtRun(projectName="testProject", runName="cleanRun", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile, overwrite="yes" )

#print(catDesign )
catDesign 
#prepareMlgtRun(catDesign )

catResult <- mlgt(catDesign )

save(catResult, file="catResult.RData")

cat.genotypes <- callGenotypes(catResult)

plotGenotypeEvidence(callList=cat.genotypes, file="genotypeCallsAll.pdf")

writeGenotypeCallsToFile(cat.genotypes )

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

# New intersect markers - should align better with IMGT allleles. 
# Test with Liverpool dataset only

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/testRun5"
setwd( analysisDir )

#markerFile <- "C:/Users/dave/HalfStarted/mlgt/HLA_MARKERS.fasta"
#markerList <- read.fasta(markerFile ,as.string=T)
intersectMarkerList <- read.fasta("C:/Users/dave/HLA/data/alleleSeqs/HLA_intersectMarkersDec11.fasta", as.string=T)


# data in fasta file.
inputDataFile <- "C:/Users/dave/HLA/data/ID632_FM_preliminary_ALL.fasta"
intersectDesign.Liverpool <- new("mlgtDesign", projectName="testProject", runName="Intersect.Liverpool", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(intersectDesign.Liverpool )

prepareMlgtRun(intersectDesign.Liverpool )

intersectLiverpool.Result <- mlgt(intersectDesign.Liverpool)

save(intersectLiverpool.Result, file="intersectLiverResult.RData")



####################
analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/intersectCat"
setwd( analysisDir )

# data in fasta file.
inputDataFile <- "C:/Users/dave/HalfStarted/mlgt/catLiverpoolNewcastleData.fasta"
intersect.catDesign <- new("mlgtDesign", projectName="testProject", runName="intersectCat", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(intersect.catDesign )

prepareMlgtRun(intersect.catDesign )

intersect.catResult <- mlgt(intersect.catDesign )

save(intersect.catResult, file="intersect.catResult.RData")


##################

# NCL Run_760591

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/intersect760591"
setwd( analysisDir )

# data in fasta file.
#C:\Users\dave\HLA\data\NCL\Run_760591
inputDataFile <- "C:/Users/dave/HLA/data/NCL/Run_760591/4.TCA.454Reads.fna"
intersect.760591.Design <- new("mlgtDesign", projectName="testProject", runName="intersect760591", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(intersect.760591.Design )

prepareMlgtRun(intersect.760591.Design)

intersect.760591.Result <- mlgt(intersect.760591.Design)

save(intersect.760591.Result, file="intersect.760591.RData")

##################

# NCL Run_756948

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/intersect756948"
setwd( analysisDir )

# data in fasta file.
#C:\Users\dave\HLA\data\NCL\Run_756948
inputDataFile <- "C:/Users/dave/HLA/data/NCL/Run_756948/3.TCA.454Reads.fna"
intersect.756948.Design <- new("mlgtDesign", projectName="testProject", runName="intersect756948", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(intersect.756948.Design )

prepareMlgtRun(intersect.756948.Design)

intersect.756948.Result <- mlgt(intersect.756948.Design)

save(intersect.756948.Result, file="intersect.756948.RData")




################


# NCL Run_763932_modified

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/intersect763932m"
setwd( analysisDir )

# data in fasta file.
inputDataFile <- "C:/Users/dave/HLA/data/NCL/Run_763932_Modified/4.TCA.454Reads.fna"
intersect.763932m.Design <- new("mlgtDesign", projectName="testProject", runName="intersect763932m", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(intersect.763932m.Design )

prepareMlgtRun(intersect.763932m.Design)

intersect.763932m.Result <- mlgt(intersect.763932m.Design)

save(intersect.763932m.Result, file="intersect.763932m.RData")

##################


# NCL Run_763932_normal

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/intersect763932n"
setwd( analysisDir )

# data in fasta file.
inputDataFile <- "C:/Users/dave/HLA/data/NCL/Run_763932_Normal/4.TCA.454Reads.fna"
intersect.763932n.Design <- new("mlgtDesign", projectName="testProject", runName="intersect763932n", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(intersect.763932n.Design )

prepareMlgtRun(intersect.763932n.Design)

intersect.763932n.Result <- mlgt(intersect.763932n.Design)

save(intersect.763932n.Result, file="intersect.763932n.RData")

##################


# NCL Run_763932_normal

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/testRun6"
setwd( analysisDir )

# data in fasta file.
inputDataFile <- "C:/Users/dave/HLA/data/NCL/Run_763932_Normal/4.TCA.454Reads.fna"
test.763932n.Design <- new("mlgtDesign", projectName="testProject", runName="test763932n", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile )

print(test.763932n.Design)

prepareMlgtRun(test.763932n.Design)

test.763932n.Result <- mlgt(test.763932n.Design)

save(test.763932n.Result, file="test.763932n.RData")

test.genotypes <- callGenotypes.mlgtResult(test.763932n.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)
write.table(test.genotypes, file="test.763932n.genotypesAlleles.tab",  quote=F, sep="\t", row.names=F)

str(test.763932n.Result, max.level=2)  # could use for print.mlgtResult()

##################



# DNAseqLab Run
# Only 2 samples but Typically over 1000 reads per marker/sample pair. 
# Alignments were taking a long time. (over an hour for a pair) so added test to reduce MUSCLE iteration from default (16) to 2.
# Fast enough but alignments are rubbish. Not sure if due to MUSCLE, or the markers.
# RAM may be an issue too with this many seqs. 
# Would be nice to limit the number of sequences to align BUT, many of the raw variants will become the same variant once converted to a sub-alignment. 

# Variant counts are very impressive though. e.g. Var1 5712, Var 2 89 ! (~1%) 

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/dnaSeqLabTest"
setwd( analysisDir )

seqLabSamples <- c("MID1", "MID3")
seqLab.rTagList <-  seqLab.fTagList <- list(MID1="ACGAGTGCGT", MID3="AGACGCACTC")		# same MIDs front and back.


# data in fasta file. C:\Users\dave\NextGen\DNAseqLab
inputDataFile <- "C:/Users/dave/NextGen/DNAseqLab/1.TCA.454Reads.fna"
dnaSeqLabTest.Design <- prepareMlgtRun( projectName="testProject", runName="dnaSeqLabTest", 
				samples=seqLabSamples , markers=intersectMarkerList ,
				fTags=seqLab.fTagList, rTags=seqLab.rTagList, inputFastaFile=inputDataFile, overwrite="yes" )

printBlastResultGraphs(dnaSeqLabTest.Design)

#prepareMlgtRun(dnaSeqLabTest.Design)

dnaSeqLabTest.Result <- mlgt(dnaSeqLabTest.Design)
save(dnaSeqLabTest.Result, file="dnaSeqLabTest.RData")

dnaSeqLabTest.Result
dnaSeqLabTest.Result@markerSampleList[["A_E2"]]
test.genotypes <- callGenotypes(dnaSeqLabTest.Result)

test.genotypes[["A_E2"]]
test.genotypes[["C_E3"]]
test.genotypes[["DRB1_E2"]]

plotGenotypeEvidence(genotypeCall = test.genotypes[["A_E3"]])
plotGenotypeEvidence(genotypeCall = test.genotypes[["DRB1_E2"]])

#test.genotypes <- callGenotypes.mlgtResult(dnaSeqLabTest.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)
#write.table(test.genotypes, file="dnaSeqLabTest.genotypesAlleles.tab",  quote=F, sep="\t", row.names=F)

##################

#### USER problems
# not specifying the auxillary paths correctly and/or having spaces in the path to the working directory (WINDOWS!)
# TODO specify no spaces in file paths and test for this in scripts. Give parameter vals in single quotes?
library(mlgt)

formatdbPath <- "C:/Users/Public/Apps/Blast/bin/nothing.exe"
fastacmdPath <- "C:/Users/Public/Apps/Blast/bin/nothing.exe"
blastAllPath <- "C:/Users/Public/Apps/Blast/bin/nothing.exe"
musclePath <- "C:/Users/Public/Apps/Muscle/nothing.8.31_i86win32.exe"


analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/errors"
setwd( analysisDir )


intersectMarkerList <- read.fasta("C:/Users/dave/HLA/data/alleleSeqs/HLA_intersectMarkersDec11.fasta", as.string=T)

seqLabSamples <- c("MID1", "MID3")
seqLab.rTagList <-  seqLab.fTagList <- list(MID1="ACGAGTGCGT", MID3="AGACGCACTC")		# same MIDs front and back.


# data in fasta file. C:\Users\dave\NextGen\DNAseqLab
inputDataFile <- "C:/Users/dave/NextGen/DNAseqLab/1.TCA.454Reads.fna"
error.Design <- prepareMlgtRun( projectName="testProject", runName="dnaSeqLabTest", 
				samples=seqLabSamples , markers=intersectMarkerList ,
				fTags=seqLab.fTagList, rTags=seqLab.rTagList, inputFastaFile=inputDataFile, overwrite="yes" )







###########



alleleDbObject <- new(alleleDb,...)