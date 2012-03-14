
sampleTable <- read.delim("C:/Users/dave/HalfStarted/mlgt/data/sampleMids.tab", header=T)

sampleTable <- read.delim("C:/Users/dave/HalfStarted/mlgt/data/sampleMids.tab", header=T, stringsAsFactors=F)
fTagList <- as.list(sampleTable$barcode.seq)
names(fTagList) <- sampleTable$sample.id
rTagList <- fTagList
#The names of the samples
sampleList <- names(fTagList)



######################### DEVEL for v0.14

## test if coping with spaces in the input file location.

Sys.setenv(BLASTALL_PATH="C:/Users/Public/Apps/Blast/bin/blastall.exe",
		FORMATDB_PATH="C:/Users/Public/Apps/Blast/bin/formatdb.exe",
		FASTACMD_PATH="C:/Users/Public/Apps/Blast/bin/fastacmd.exe",
		MUSCLE_PATH="C:/Users/Public/Apps/Muscle/muscle3.8.31_i86win32.exe")





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

setwd("C:/Users/Dave/HalfStarted/mlgt/testProject/linkTest/")
inputDataFile <- "C:/Temp/I am spaced/ID632_FM_preliminary_ALL.fasta"

my.mlgt.Design <- prepareMlgtRun(projectName="myProject", 
				runName="spaceNames", samples=sampleList, 
				markers=myMarkerList, fTags=fTagList, 
				rTags=rTagList, inputFastaFile=inputDataFile, 
				overwrite="yes")

######## below added to mlgt.R

alignReport.fix <- function(mlgtResultObject, markers=names(mlgtResultObject@markers), samples=mlgtResultObject@samples,
		correctThreshold = 0.01,  consThreshold = (1 - correctThreshold), profPlotWidth = 60, fileName=NULL, method="table", warn=TRUE)  {

	# need method for both plots (save processing time) but how to do without generating profiles twice.
	# need to tidy and label profile plots.

	reportList <- list()
	if(!is.null(fileName)) {
		if(length(grep(".pdf$", fileName))  < 1 & !is.null(fileName)) {
				fileName <- paste(fileName,"pdf", sep=".")
			}
		pdf(fileName)
	}

	colTable <- data.frame(profChar = c("-","a","c","g","t"),
					profCol = c("grey","green", "blue", "yellow", "red"))

	for(thisMarker in markers)  {
		thisTable <- mlgtResultObject@varCountTables[[thisMarker]]
		cat(paste(thisMarker,"\n"))
		if(nrow(thisTable) < 1)  {		# drops this marker if no data. Might be better to return and empty table?
			if(warn) { 
					warning(paste("No data for",  thisMarker)) 
				}
			#reportTable[thisSample,c("numbSeqs","numbVars")] <- 0
			#reportTable[thisSample,c("alignLength","invar.sites","mafBelowThreshold","mafAboveThreshold")] <- NA
			next;
		}
		reportTable <- data.frame()
		for(thisSample in samples) {
			if(is.na(match(thisSample,names(thisTable)))) {
				if(warn) { 
					warning(paste("No variants for", thisSample, "and", thisMarker)) 
				}
				#reportTable[thisSample,c("invar.sites","mafBelowThreshold","mafAboveThreshold")] <- NA
				reportTable[thisSample,"numbSeqs"] <- 0
				reportTable[thisSample,"numbVars"] <- 0
				reportTable[thisSample,"alignLength"] <- NA
				reportTable[thisSample,"invar.sites"] <- NA
				reportTable[thisSample,"mafAboveThreshold"] <- NA
				reportTable[thisSample,"mafBelowThreshold"] <- NA
				next;
			}
			#cat(paste(thisSample ," "))
			valueIndex <- !is.na(thisTable[,thisSample])
			seqCounts <- thisTable[valueIndex,thisSample]
			sampleSeqs <- rep(row.names(thisTable)[valueIndex ], seqCounts)
			thisAlign <- as.alignment(sum(seqCounts), sampleSeqs, sampleSeqs)

			if(thisAlign$nb < 2)  {	# too few sequences to plot.
				#reportTable[thisSample,c("invar.sites","mafBelowThreshold","mafAboveThreshold")] <- NA
				reportTable[thisSample,"numbSeqs"] <- 1
				reportTable[thisSample,"numbVars"] <- 1
				reportTable[thisSample,"alignLength"] <- nchar(thisAlign$seq[1])
				reportTable[thisSample,"invar.sites"] <- NA
				reportTable[thisSample,"mafAboveThreshold"] <- NA
				reportTable[thisSample,"mafBelowThreshold"] <- NA
				next; 
			}
			thisProfile <- consensus(thisAlign , method="profile")
			thisConsensus <- con(thisAlign, method="threshold", threshold=consThreshold)
			totalSeqs <- sum(seqCounts)
			totalLetters <- nrow(thisProfile)
			mafList <- apply(thisProfile,  2, FUN=function(x) (sort(x)[totalLetters-1] / totalSeqs)) 

			reportTable[thisSample, "numbSeqs"] <- totalSeqs 
			reportTable[thisSample, "numbVars"] <- length(seqCounts)
			reportTable[thisSample, "alignLength"] <- ncol(thisProfile)			
			reportTable[thisSample, "invar.sites"] <- sum(mafList == 0)
			## variable sites with minor allele > correction threshold.
			reportTable[thisSample, "mafAboveThreshold"] <- sum(mafList >= correctThreshold )
			## varaible sites with minor alleles < correction threshold.
			reportTable[thisSample, "mafBelowThreshold"] <- reportTable[thisSample, "alignLength"] - (reportTable[thisSample, "invar.sites"] + reportTable[thisSample, "mafAboveThreshold"])
	
			if(method=="profile")  {

				profColours <-  as.character(colTable$profCol[match(row.names(thisProfile),colTable$profChar)])
				## splits the plotting across so many lines. Buffers final plot to constant length.
				profLen <- ncol(thisProfile)
				n_plot <- ceiling(profLen / profPlotWidth )
				plotProfile <- thisProfile
				plotConsensus <-  toupper(thisConsensus) 
				remainder <- profLen %% profPlotWidth 

				if(remainder > 0) {
					extLen <- profLen + (profPlotWidth - remainder)
					profExtension <- matrix(0,nrow=nrow(thisProfile), ncol=(profPlotWidth - remainder), 
								dimnames=list(row.names(thisProfile),c((profLen+1): extLen)))
					plotProfile <- cbind(thisProfile,profExtension)
					plotConsensus <- c(toupper(thisConsensus), rep("",remainder)) 
				}

				old.o <- par(mfrow=c((n_plot+1),1), mar=c(2,4,2,2))
				plot.new()	
				title(paste(thisMarker, thisSample, sep=" : "), line=0)		
				title("", line=0,adj=0,	sub=paste(c("Total sites","Invariant sites","MAF above threshold","MAF below threshold"),reportTable[thisSample,3:6],sep=": ",collapse="\n") )
				
				legend("right", legend=toupper(colTable$profChar),fill=as.character(colTable$profCol), horiz=T)
				for(i in 1:n_plot) {
					start <- ((i-1)*profPlotWidth) + 1
					#index <- start:(min(start+profPlotWidth, profLen))

					index <- start:((start+profPlotWidth)-1)
					barplot(plotProfile[,index ], col=profColours , names.arg=toupper(plotConsensus[index]) )
				}
				par(old.o)
			}

			if(method=="hist")  {
				#if(!is.null(fileName)) pdf(fileName)
				if(sum(mafList > 0) > 0) {
					hist(mafList[mafList > 0], breaks=200, xlim=c(0,0.5), xlab="Site-specific minor allele frequency", sub="non-zero values only",main=paste(thisMarker, thisSample, sep=":")) 
					abline(v=correctThreshold, lty=2)
				}	
			}
			
		}			
		reportList[[thisMarker]] <- reportTable
	}
	if(!is.null(fileName)) {
		dev.off()
		cat(paste("Alignment figures(s) plotted to", fileName,"\n"))
	}
	#print(reportList)
	return(reportList)
}

alignReport.fix(my.mlgt.Result.dgM,markers="DPB1_E2", samples="MID1", method="hist")
alignReport.fix(my.mlgt.Result.dgM,markers="DPB1_E2", samples="MID1", method="profile")
alignReport.fix(my.mlgt.Result.dgM,markers="DPA1_E2", samples="MID1", method="profile")
alignReport.fix(my.mlgt.Result.dgM, method="profile",  file="alignReport.profile")



######################### DEVEL for v0.13

library(seqinr)
#source("C:/Users/dave/HalfStarted/mlgt/mlgt.R")
source("C:/Users/Dave/Work/mlgtProject/mlgt/mlgt.R")

setwd("C:/Users/Dave/Work/mlgtProject/testProject")
Sys.setenv(BLASTALL_PATH="C:/Users/Public/Apps/Blast/bin/blastall.exe",
		FORMATDB_PATH="C:/Users/Public/Apps/Blast/bin/formatdb.exe",
		FASTACMD_PATH="C:/Users/Public/Apps/Blast/bin/fastacmd.exe",
		MUSCLE_PATH="C:/Users/Public/Apps/Muscle/muscle3.8.31_i86win32.exe")
#fTagList <- read.fasta("C:/Users/Dave/Work/mlgtProject/data/namedBarcodes.fasta", 
#			as.string=T) 
#rTagList <- fTagList
#sampleList <- names(fTagList)
#myMarkerList <- read.fasta("C:/Users/Dave/Work/mlgtProject/data/HLA_namedMarkers.fasta",
#			as.string=T)	
		
#inputDataFile <- "C:/Users/Dave/Work/mlgtProject/data/sampleSequences.fasta"

#my.mlgt.Design <-  prepareMlgtRun(projectName="testProject", runName="cleanRun", 
#				samples=sampleList, markers=myMarkerList,
#				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile, overwrite="yes" )

#my.mlgt.Result <- mlgt(my.mlgt.Design)
#save(my.mlgt.Result, file= "my.mlgt.Result.RData")

load("my.mlgt.Result.RData")

thisMarker <- "DPA1_E2"
thisSample <- "MID-11"


# TODO: errorCorrection as option within mlgt()
# DONE: output fasta sequence of all variants in an alignment (useful for exploratory analysis)
# DONE: errorCorrect: strip all-gap columns from corrected alignment

#stripGapColumns(testAlign)


dumpVariants(my.mlgt.Result,fileSuffix="variantDump.fasta")

dumpVariants(my.corrected.mlgt.Result, markers="FLA_DRB", samples="cat348" )
dumpVariants(my.corrected.mlgt.Result,fileSuffix="corrected.0.05.variantDump.fasta")
dumpVariants(my.corrected.mlgt.Result,fileSuffix="corrected.0.05.variantDump.unique.fasta", uniqueOnly=T)

# having added varCountTables list to mlgtResult  
# Here is a good way to look at it:-
data.frame(my.mlgt.Result@varCountTables[["A_E3"]], row.names=NULL)


## DONE: function to errorCorrect. Could be used within mlgt or called afterwards.
## TODO: Incorporate into mlgt()

thisMarker <- "DPA1_E2"
thisSample <- "MID-22"

testCorrect <- errorCorrect(alignment=testAlign)

write.fasta(sequences=lapply(testAlign$seq,s2c), names=testAlign$nam, file.out="testUncorrection.fasta")
write.fasta(sequences=lapply(testCorrect$seq,s2c), names=testCorrect$nam, file.out="testCorrection.fasta")


my.corrected.mlgt.Result <- errorCorrect.mlgtResult(my.mlgt.Result)

my.corrected.mlgt.Result <- errorCorrect.mlgtResult(my.mlgt.Result,correctThreshold=0.05)

alignReport(my.mlgt.Result, fileName="alignReportOut.pdf", method="profile")
alignReport(my.corrected.mlgt.Result, fileName="alignReportOut.corrected.pdf", method="profile")

alignReport(my.mlgt.Result, fileName="alignReportOut.hist.pdf", method="hist")
alignReport(my.corrected.mlgt.Result, fileName="alignReportOut.hist.corrected.pdf")

my.genotypes <- callGenotypes(my.mlgt.Result)
my.corrected.genotypes <- callGenotypes(my.corrected.mlgt.Result)

my.genotypes[[thisMarker]]@genotypeTable
my.corrected.genotypes[[thisMarker]]@genotypeTable
## error correction has little effect on the the sample dataset even with high threshold of 0.05. Most samples with > 100 seqs are called correctly anyway.

 thisTable <- my.mlgt.Result@varCountTables[[thisMarker]]
 thisTable <- my.corrected.mlgt.Result@varCountTables[[thisMarker]]


alignReport(my.corrected.mlgt.Result,markers="DQA1_E2", samples="MID-5", method="profile")
alignReport(my.corrected.mlgt.Result,markers="DQA1_E2", samples="MID-8", method="profile")
alignReport(my.corrected.mlgt.Result,markers="DQA1_E2", samples="MID-20", method="profile")

alignReport(my.mlgt.Result,markers=thisMarker, samples=thisSample, method="profile")
alignReport(my.mlgt.Result,markers=thisMarker, samples=thisSample, method="hist", fileName="testOutHist")
alignReport(my.mlgt.Result, method="profile", fileName="testOutMultiProfile")

alignReport(my.mlgt.Result,markers="DPA1_E2", samples="MID-17", method="profile")
alignReport(my.mlgt.Result,markers="DPA1_E2", samples="MID-1", method="hist")
alignReport(my.mlgt.Result,markers="DPA1_E2", samples="MID-22", method="hist") # good example where change would be useful
alignReport(my.mlgt.Result,markers="DPA1_E2", samples="MID-1", method="profile", correctThreshold=0.02)
my.alignReport <- alignReport(my.mlgt.Result)




################# unused stuff

## This might be useful but is replicated by alignReport(..., method="hist"). Not added to mlgt.R
plotSiteMaf <- function(mlgtResultObject, markers=names(mlgtResultObject@markers), samples=mlgtResultObject@samples,
		correctThreshold = 0.01)  {
	for(thisMarker in markers)  {
		thisTable <- mlgtResultObject@varCountTables[[thisMarker]]
	
		for(thisSample in samples) {
			valueIndex <- !is.na(thisTable[,thisSample])

			if(length(valueIndex) <1) next;
			seqCounts <- thisTable[valueIndex,thisSample]
			sampleSeqs <- rep(row.names(thisTable)[valueIndex ], seqCounts)
		
			thisAlign <- as.alignment(sum(seqCounts), sampleSeqs, sampleSeqs)
			thisProfile <- consensus(thisAlign , method="profile")
			totalSeqs <- sum(seqCounts)
			totalLetters <- nrow(thisProfile)
			mafList <- apply(thisProfile,  2, FUN=function(x) (sort(x)[totalLetters-1] / totalSeqs)) 
			
			hist(mafList, breaks=200, xlab="Site-specific minor allele frequency", main=paste(thisMarker, thisSample, sep=":")) 
			abline(v=correctThreshold, lty=2)

		}

	}
	
}


## This might be useful but is replicated by alignReport(..., method="profile") Not added to mlgt.R
## currently fails to scale final (shorter) line to match previous lines
## also needs to display marker/sample info and summary info. 
plotAlignProfile <- function(mlgtResultObject, markers=names(mlgtResultObject@markers), samples=mlgtResultObject@samples,
		correctThreshold = 0.01, consThreshold = 0.95, profPlotWidth = 60)  {
	
	colTable <- data.frame(profChar = c("-","a","c","g","t"),
					profCol = c("grey","green", "blue", "yellow", "red"))
	for(thisMarker in markers)  {
		thisTable <- mlgtResultObject@varCountTables[[thisMarker]]
		cat(paste(thisMarker,"\n"))
		for(thisSample in samples) {
			if(is.na(match(thisSample,names(thisTable)))) {
				warning(paste("No variants for", thisSample, "and", thisMarker))
				next;
			}
			#cat(paste(thisSample ," "))
			valueIndex <- !is.na(thisTable[,thisSample])
			#cat(length(valueIndex))
			#if(length(valueIndex) <1) next;
			seqCounts <- thisTable[valueIndex,thisSample]
			sampleSeqs <- rep(row.names(thisTable)[valueIndex ], seqCounts)
		
			thisAlign <- as.alignment(sum(seqCounts), sampleSeqs, sampleSeqs)
			thisProfile <- consensus(thisAlign , method="profile")
			thisConsensus <- con(thisAlign, method="threshold", threshold=consThreshold)
			#totalSeqs <- sum(seqCounts)
			#totalLetters <- nrow(thisProfile)
			#mafList <- apply(thisProfile,  2, FUN=function(x) (sort(x)[totalLetters-1] / totalSeqs)) 

			profColours <-  as.character(colTable$profCol[match(row.names(thisProfile),colTable$profChar)])
			## splits the plotting across so many lines. Buffers final plot to constant length.
			profLen <- ncol(thisProfile)
			n_plot <- ceiling(profLen / profPlotWidth )
			plotProfile <- thisProfile
			plotConsensus <-  toupper(thisConsensus) 
			remainder <- profLen %% profPlotWidth 
			
			if(remainder > 0) {
				extLen <- profLen + (profPlotWidth - remainder)
				profExtension <- matrix(0,nrow=nrow(thisProfile), ncol=(profPlotWidth - remainder), 
							dimnames=list(row.names(thisProfile),c((profLen+1): extLen)))
				#row.names
				plotProfile <- cbind(thisProfile,profExtension)
				plotConsensus <- c(toupper(thisConsensus), rep("",remainder)) 
				#thisProfile[,(profLen+1): extLen ] <- 0
			
			}
			
			old.o <- par(mfrow=c((n_plot+1),1), mar=c(2,4,2,2))
			plot.new()	
			title(paste(thisMarker, thisSample, sep=" : "), line=0,	sub="no data")
			legend("center", legend=toupper(colTable$profChar),fill=as.character(colTable$profCol), horiz=T)
			for(i in 1:n_plot) {
				start <- ((i-1)*profPlotWidth) + 1
				#index <- start:(min(start+profPlotWidth, profLen))
				
				index <- start:((start+profPlotWidth)-1)
				barplot(plotProfile[,index ], col=profColours , names.arg=toupper(plotConsensus[index]) )
			}
			par(old.o)
			

		}

	}	
		
}




thisMarker <- "DPA1_E2"
thisSample <- "MID-11"

plotSiteMaf(my.mlgt.Result, markers=thisMarker, samples=thisSample)
plotAlignProfile(my.mlgt.Result, markers=thisMarker, samples=thisSample) 

#plotSiteMaf(my.mlgt.Result)		# tries to plot multiple to same window. 
#plotAlignProfile(my.mlgt.Result)	# tries to plot multiple to same window. 

my.alignReport <- alignReport(my.mlgt.Result,markers=thisMarker, samples=thisSample)

#my.alignReport <- alignReport(my.mlgt.Result)	# breaks with error about recursive indexing.



#################### old stuff

### NOT USING VARIANT COUNTS  # Is there any point?
correctThreshold <- 0.01
consThreshold <- 0.95
profPlotWidth <- 60
thisMarker <- "A_E3"
thisMarker <- "DPA1_E2"
thisSample <- "MID-8"
thisSample <- "MID-17"
thisTable <- my.mlgt.Result@varCountTables[[thisMarker]]
valueIndex <- !is.na(thisTable[,thisSample])
seqCounts <- thisTable[valueIndex,thisSample]
sampleSeqs <- row.names(thisTable)[valueIndex ]

thisAlign <- as.alignment(length(seqCounts), sampleSeqs, sampleSeqs)

#as.matrix.alignment()

thisProfile <- consensus(thisAlign , method="profile")
thisConsensus <- con(thisAlign, method="threshold", threshold=consThreshold)

#barplot(thisProfile, col=c("red", "green", "blue", "yellow"))

n_plot <- ceiling(ncol(thisProfile)/ profPlotWidth )
old.o <- par(mfrow=c(n_plot,1))
for(i in 1:n_plot) {
	start <- ((i-1)*profPlotWidth) + 1
	index <- start:(min(start+profPlotWidth, ncol(thisProfile)))
	barplot(thisProfile[,index ], col=c("red", "green", "blue", "yellow"), names.arg=toupper(thisConsensus[index ]) )
}
par(old.o)

##########USING VARIANT COUNTS
seqCounts <- thisTable[valueIndex,thisSample]
sampleSeqs <- rep(row.names(thisTable)[valueIndex ], seqCounts)


thisAlign <- as.alignment(sum(seqCounts), sampleSeqs, sampleSeqs)

thisProfile <- consensus(thisAlign , method="profile")
thisConsensus <- con(thisAlign, method="threshold", threshold=consThreshold)

#barplot(thisProfile, col=c("red", "green", "blue", "yellow"))

n_plot <- ceiling(ncol(thisProfile)/ profPlotWidth )
old.o <- par(mfrow=c(n_plot,1))
for(i in 1:n_plot) {
	start <- ((i-1)*profPlotWidth) + 1
	index <- start:(min(start+profPlotWidth, ncol(thisProfile)))
	barplot(thisProfile[,index ], col=c("red", "green", "blue", "yellow"), names.arg=toupper(thisConsensus[index ]) )
}
par(old.o)

## report on frequency spectra of variable sites.

## invariant sites
## variable sites with minor allele > correction threshold.
## varaible sites with minor alleles < correction threshold.

totalSeqs <- sum(seqCounts)
totalLetters <- nrow(thisProfile)
mafList <- apply(thisProfile,  2, FUN=function(x) (sort(x)[totalLetters-1] / totalSeqs)) 
## invariant sites
sum(mafList == 0)
## variable sites with minor allele > correction threshold.
sum(mafList > correctThreshold )
## varaible sites with minor alleles < correction threshold.
sum(mafList < correctThreshold )

hist(mafList, breaks=200, xlab="Site-specific minor allele frequency", main=paste(thisMarker, thisSample, sep=":")) ; abline(v=correctThreshold, lty=2)



############################DEVEL FOR v0.12


### paths with spaces.
## Use this structure in the Vignette and general usage.
Sys.setenv(BLASTALL_PATH="C:/Users/Public/Apps/Blast/bin/blastall.exe",
		FORMATDB_PATH="C:/Users/Public/Apps/Blast/bin/formatdb.exe",
		FASTACMD_PATH="C:/Users/Public/Apps/Blast/bin/fastacmd.exe",
		MUSCLE_PATH="C:/Users/Public/Apps/Muscle/muscle3.8.31_i86win32.exe")





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