




############################DEVEL FOR v0.12


### paths with spaces.
## Use this structure in the Vignette and general usage.
Sys.setenv(BLASTALL_PATH="C:/Users/Public/Apps/Blast/bin/blastall.exe",
		FORMATDB_PATH="C:/Users/Public/Apps/Blast/bin/formatdb.exe",
		FASTACMD_PATH="C:/Users/Public/Apps/Blast/bin/fastacmd.exe",
		MUSCLE_PATH="C:/Users/Public/Apps/Muscle/muscle3.8.31_i86win32.exe")


library(seqinr)
source("C:/Users/dave/HalfStarted/mlgt/mlgt.R")


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

setwd("C:/Users/dave/HalfStarted/mlgt/testProject/test space")

# test space in path to aux program
Sys.setenv(BLASTALL_PATH="C:/Users/dave/HalfStarted/mlgt/testProject/test space/blastall.exe")
my.mlgt.Design <- prepareMlgtRun(projectName="myProject", 
				runName="myRun", samples=sampleList, 
				markers=myMarkerList, fTags=fTagList, 
				rTags=rTagList, inputFastaFile=inputDataFile, 
				overwrite="yes")

Sys.setenv(BLASTALL_PATH="C:/Users/dave/HalfStarted/mlgt/testProject/test space/blastall.exe")
Sys.setenv(MUSCLE_PATH="C:/Users/dave/HalfStarted/mlgt/testProject/test space/muscle3.8.31_i86win32.exe")
my.mlgt.Result <- mlgt(my.mlgt.Design)

# test space in input fasta file location
inputDataFile <- "C:/Users/dave/HalfStarted/mlgt/testProject/test space/4.TCA.454Reads.fna"
my.mlgt.Design2 <- prepareMlgtRun(projectName="myProject", 
				runName="myRun", samples=sampleList, 
				markers=myMarkerList, fTags=fTagList, 
				rTags=rTagList, inputFastaFile=inputDataFile, 
				overwrite="yes")

#############



## Use this within the functions.
pathNames <- c("BLASTALL_PATH","FORMATDB_PATH","FASTACMD_PATH","MUSCLE_PATH")

for(thisPath in pathNames)  {
	if(nchar(Sys.getenv(thisPath)) < 1) {
		stop(paste(thisPath,"has not been set!"))
	}
	# shellQuote any paths containing spaces.
	if(length(grep(" ",Sys.getenv(thisPath), fixed=T))  > 0 )  Sys.setenv(thisPath, shQuote(thisPath))
}

#shQuote

#decided not to use because: doesn't work!
quoteIfSpaces <- function(pathString)  {
	if(length(grep(" ",pathString, fixed=T))  > 0 ) {
		  pathString <- shQuote(pathString)	
	}
	return(pathString)
}



#then need to test each of these and replace if necessary.

# The input file needs to be tested too. 




quotePath <- function(rawPathVarName) {
	assign(rawPathVarName,shQuote(get(rawPathVarName)))	

}


tempPath <- "C:/Users/Public/Apps/Blast/bin 2/fastacmd.exe"
quotePath("tempPath") 


grep(" ", tempPath)

grep(" ", tempPath) > 0




Sys.setenv(R_TEST="testit", "A+C"=123)




##############################################LATER TODO:


#Script to be used to run in parallel (unlikely to have true multicore for Windows)






############################################################


# Looking at alleleDb objects

str(intersect.cleanRun.Result@alleleDb[[7]], deparse.level=1)







#################################### New Clean Run

# prepare mlgt by setting up blast dbs and running blast searches.


library(seqinr)

formatdbPath <- "C:/Users/Public/Apps/Blast/bin/formatdb.exe"
fastacmdPath <- "C:/Users/Public/Apps/Blast/bin/fastacmd.exe"
blastAllPath <- "C:/Users/Public/Apps/Blast/bin/blastall.exe"
musclePath <- "C:/Users/Public/Apps/Muscle/muscle3.8.31_i86win32.exe"

source("C:/Users/dave/HalfStarted/mlgt/mlgt.R")



###################### GET INFO for this run


analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/cleanRun"
setwd( analysisDir )

#sampleFile <- "C:/Users/dave/HalfStarted/mlgt/samples_ALL.tab"
#sampleTable <- read.delim(sampleFile, sep=";")

#intersectMarkerList <- read.fasta("C:/Users/dave/HLA/data/alleleSeqs/HLA_intersectMarkersDec11.fasta", as.string=T)

# samples
#sampleList <- as.character(unique(sampleTable$sample_name))

## TODO: SIMPLIFY
# MIDs/tags
#tempTable <- unique(subset(sampleTable, select=c(tag_f, sample_name)))
#fTagList <- as.list(as.character(tempTable$tag_f))
#names(fTagList) <- tempTable$sample_name
#tempTable <- unique(subset(sampleTable, select=c(tag_r, sample_name)))
#rTagList <- as.list(as.character(tempTable$tag_r))
#names(rTagList) <- tempTable$sample_name


fTagList <- read.fasta("C:/Users/dave/HLA/data/fTags.fasta", as.string=T)
rTagList <- fTagList
sampleList <- names(fTagList)

intersectMarkerList <- read.fasta("C:/Users/dave/HLA/data/alleleSeqs/HLA_intersectMarkersDec11.fasta", 
				as.string=T)



########## START

inputDataFile <- "C:/Users/dave/HLA/data/NCL/Run_763932_Modified/4.TCA.454Reads.fna"

# Creates object to store run settings
intersect.cleanRun.Design <- prepareMlgtRun(projectName="testProject", runName="cleanRun", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile, overwrite="yes" )


# inspect BLAST results for a specific marker
thisMarker <- "DPA1_E2"
topHits <- getTopBlastHits(intersect.cleanRun.Design@markerBlastResults)
inspectBlastResults(topHits, thisMarker)
# automatic output to pdf of blast result graphs for a list of markers.
printBlastResultGraphs(intersect.cleanRun.Design)

# run mlgt
intersect.cleanRun.Result <- mlgt(intersect.cleanRun.Design)

# Call genotypes (standard, no allele mapping)
test.genotypes <- callGenotypes.mlgtResult(intersect.cleanRun.Result,  mapAlleles=FALSE)


# Call genotypes with marker specific parameter values for minimum number of reads (minTotalReads)
test.genotypes <- callGenotypes.mlgtResult(intersect.cleanRun.Result,  mapAlleles=FALSE, minTotalReads=seq(20,200, 10))

# Call genotypes on a single table (currently returns only a table, which is desirable as used internally by callGenotypes.mlgtResult
test.genotypes <- callGenotypes(table=intersect.cleanRun.Result@markerSampleList[["DPA1_E2"]])


# Write genotype call graphics to file
writeGenotypeCallsToFile(test.genotypes)								
writeGenotypeCallsToFile(test.genotypes, singleFile=F, file="genotypeTable.tab")
writeGenotypeCallsToFile(test.genotypes, singleFile=T, file="genotypeTable.tab")
writeGenotypeCallsToFile(genotypeCall=test.genotypes[[1]])

# Make known allele DB.
# Create a structured list of known alleles (e.g. from IMGT/HLA)

markerImgtFileTable <- read.delim("C:/Users/dave/HalfStarted/mlgt/marker.imgt.msf.list.tab", header=T)
imgtFileDir <- "C:/Users/dave/HLA/data/IMGT_manualDownload/"

# create a 'variantMap' object for each marker and store them all in a list.
knownAlleleDb <- list()
# takes about 2 minutes with 17 markers. 
for(thisMarker in names(intersectMarkerList)) {
	baseFile <- markerImgtFileTable$imgtAlignFile[markerImgtFileTable$marker==thisMarker]
	imgtAlignFile <- paste(imgtFileDir,baseFile , sep="") 
	knownAlleleDb[[thisMarker]] <- createKnownAlleleList(thisMarker,intersectMarkerList[[thisMarker]][1], imgtAlignFile)
}

# Call genotypes and map alleles to list of known alleles
#test.genotypes <- callGenotypes.mlgtResult(intersect.cleanRun.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)
#writeGenotypeCallsToFile(test.genotypes)	


##############################




#############################################################
#testing different overwrite values
intersect.cleanRun.Design <- prepareMlgtRun(projectName="testProject", runName="cleanRun", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile, overwrite="prompt" )


intersect.cleanRun.Design <- prepareMlgtRun(projectName="testProject", runName="cleanRun", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile, overwrite="no" )



intersect.cleanRun.Design <- prepareMlgtRun(projectName="testProject", runName="cleanRun", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile, overwrite="pants" )


intersect.cleanRun.Design <- prepareMlgtRun(projectName="testProject", runName="cleanRun", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile)