
########################## DEV VERSION

# This alignment has ----- exactly at the primer limits.  Could be good. 
# Probably don't need to use primers at all and don't need to extract subseqs, only seqs, in correct orientation.
# How to extract ref, trim others  and keep track of counts for each. 
# Also, should alignment be across all seqs of a marker? Prob not. then becomes run dependent. Do within sample only.


# need to rework blastHit table to separate hits on forward and reverse strand.
# TODO: could use same method to improve MID hits and check both ends are marked. (or use different MIDs for each end).
	#topHits <- getTopBlastHits("blastOut.markers.tab")
	#topHits$strand <- ifelse(topHits$s.end > topHits$s.start, 1,2)
	#fMarkerMap <- split(as.character(topHits$query[topHits$strand==1]), topHits$subject[topHits$strand==1])
	#rMarkerMap <- split(as.character(topHits$query[topHits$strand==2]), topHits$subject[topHits$strand==2])

#


getSubSeqsTable <- function(thisMarker, thisSample, sampleMap, fMarkerMap,rmarkerMap, markerSeq)  {

	varCountTable <- data.frame()
	#thisMarker <- "DPA1_E2"
	#thisSample <- "MID-1"
	#intersect(sampleMap[[thisSample]], markerMap[[thisMarker]])
	fPairSeqList <- intersect(sampleMap[[thisSample]], fMarkerMap[[thisMarker]]) # fMarkerMap[[thisMarker]]
	rPairSeqList <- intersect(sampleMap[[thisSample]], rMarkerMap[[thisMarker]])
	#intersect(	fMarkerMap[[thisMarker]],	rMarkerMap[[thisMarker]])



# extract raw seqs from blastdb for forward hits. # THIS FUNCT CURRENTLY UNUSED BECAUSE NEED TO ASSESS IF ANYTHING TO QUERY BEFORE RUNNING fastacmd
	extractRawSeqsCommand <- function(idList,strand=1, fileName)  {
		queryList <- paste(idList , collapse=",")
		fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-S", strand, "-o", fileName,  "-s", queryList)	
		return(fastacmdCommand)
	}

	fRawSeqs <- rRawSeqs <- list()

	queryList <- paste(fPairSeqList , collapse=",")	
	if(length(queryList) > 0)  {
		fRawSeqFileName <- "fRawSeqExtract.fasta"
		strand <- 1
		fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-S", strand, "-o", fRawSeqFileName,  "-s", queryList)		
		system(fastacmdCommand)
		fRawSeqs <- read.fasta(fRawSeqFileName , as.string=T)
	}

	queryList <- paste(rPairSeqList , collapse=",")	
	if(length(queryList) > 0)  {
		rRawSeqFileName <- "rRawSeqExtract.fasta"
		strand <- 2
		fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-S", strand, "-o", rRawSeqFileName,  "-s", queryList)		
		system(fastacmdCommand)
		rRawSeqs <- read.fasta(rRawSeqFileName , as.string=T)
	}
#	system( extractRawSeqsCommand(idList=fPairSeqList ,strand=1, fileName=fRawSeqFileName ) )
#	rRawSeqFileName <- "rRawSeqExtract.fasta"
#	system( extractRawSeqsCommand(idList=rPairSeqList ,strand=2, fileName=rRawSeqFileName ) )



# Make file of unique seqs. Can name with sequence 

	### STRAND!!! - Done!

	rawSeqs <- c(fRawSeqs ,rRawSeqs )
	if(length(rawSeqs) < 1)  {
		return(varCountTable)
	}
	rawSeqCountTable <-  as.data.frame(table(unlist(rawSeqs)))
	names(rawSeqCountTable) <- c("rawSeq", "rawCount")
	rawVariantFile <- paste("test", thisMarker, thisSample, "raw.variants.fasta",sep=".")  #"localVariants.fasta"
	#rawVariantFile <- paste(runPath, rawVariantFileName, sep="/")
	write.fasta(as.list(c(markerSeq,as.character(rawSeqCountTable$rawSeq))) ,c(thisMarker,as.character(rawSeqCountTable$rawSeq)),file=rawVariantFile ,open="w")

# Align all seqs with reference
	rawAlignFile <- paste("test", thisMarker, thisSample, "raw.align.fasta",sep=".")  #"localAlign.fasta"
	#rawAlignFile <- paste(runPath, rawAlignFileName, sep="/")
	muscleCommand <- paste(musclePath, "-in", rawVariantFile , "-out", rawAlignFile , "-diags -quiet" )
	system(muscleCommand)


# Extract portion corresponding to reference. 

	rawAlignment <- read.fasta(rawAlignFile, as.string=T)		# do not use read.alignment() - broken
	alignedMarkerSeq <- s2c(rawAlignment[[thisMarker]])
	subStart <- min(grep("-",alignedMarkerSeq ,invert=T))
	subEnd <- max(grep("-",alignedMarkerSeq ,invert=T))
	alignedSubSeqs <- lapply(rawAlignment, FUN=function(x)	substr(x[1], subStart, subEnd))
	subAlignFile <- paste("test", thisMarker, thisSample, "sub.align.fasta",sep=".")  #"localAlign.fasta"
	#subAlignFile <- paste(runPath, subAlignFileName , sep="/")
	write.fasta(alignedSubSeqs , names(alignedSubSeqs ), file=subAlignFile )

	alignedSubTable <- data.frame(rawSeq =  names(alignedSubSeqs ) , subSeq= as.character(unlist(alignedSubSeqs )))

# R-apply count of each seq.  There may be some duplicated subSeqs.
	combTable <- merge(rawSeqCountTable ,alignedSubTable , by="rawSeq", all.x=T)
	varCount <- by(combTable, as.character(combTable$subSeq), FUN=function(x) sum(x$rawCount))
	varCountTable <- data.frame(alignedVar=names(varCount), count=as.numeric(varCount))	
	varCountTable$var <- gsub("-","",varCountTable$alignedVar)
	varCountTable <- varCountTable[order(varCountTable$count,decreasing=T),]
# Make unique list, summing counts where same seq found. (easier in table than list).  

	return(varCountTable)

}



#############################

# pre-specify marker reference sequence.
#markerList <- read.fasta("HLA_MARKERS.fasta",as.string=T)
	
################


mlgt <- function(object) attributes(object)
setGeneric("mlgt")



mlgt.mlgtDesign  <- function(designObject)  {
	topHits <- getTopBlastHits("blastOut.markers.tab")
	topHits$strand <- ifelse(topHits$s.end > topHits$s.start, 1,2)
	fMarkerMap <- split(as.character(topHits$query[topHits$strand==1]), topHits$subject[topHits$strand==1])
	rMarkerMap <- split(as.character(topHits$query[topHits$strand==2]), topHits$subject[topHits$strand==2])

## NEED TO MAKE SAMPLEMAP WITH HITS TO MID IN BOTH FORWARD AND REVERSE STRANDS like marker hits are split.
## Requires retention of 2 blast hits per sequence.
	topSampleHits <- read.delim("blastOut.rTags.tab", header=F)
	names(topSampleHits ) <- c("query", "subject", "percentId", "aliLength", "mismatches", "gapOpenings", "q.start","q.end", "s.start","s.end", "p_value", "e_value")
	topSampleHits$strand <- ifelse(topSampleHits$s.end > topSampleHits$s.start, 1,2)
	fSampleMap <- split(as.character(topSampleHits$query[topSampleHits$strand==1]), topSampleHits$subject[topSampleHits$strand==1])
	rSampleMap <- split(as.character(topSampleHits$query[topSampleHits$strand==2]), topSampleHits$subject[topSampleHits$strand==2])
	# combind sampleMaps to give sequences with MIDs in both orientations.
	pairedSampleMap <- lapply(names(fSampleMap), FUN=function(x) intersect(fSampleMap[[x]], rSampleMap[[x]]))
	names(pairedSampleMap) <- names(fSampleMap)

##########ITERATIONS

markerSampleList <- list()
runSummaryTable <- data.frame()
alleleDb <- list()

for(thisMarker in names(designObject@markers)) {
#for(thisMarker in names(markerMap)) {
#for(thisMarker in names(markerMap)[1:2]) {	# temp to finish off

cat(paste(thisMarker,"\n"))
#thisMarker <- "DQA1_E2"

## might need to combine all these to return a single item.
summaryList <- list()
summaryTable <- data.frame()
markerSequenceCount <- list("noSeq"=0)		#  BUG? requires some data otherwise won't sum properly with localSequenceCount.
alleleList <- list() 
variantList <- list()
alleleCount <- 1
markerSeq <- unlist(getSequence(markerList[[thisMarker]],as.string=T))

for(thisSample in designObject@samples) {
#for(thisSample in names(pairedSampleMap)[1:4]) {
	#print(thisSample)

	
	testPairSeqList <- intersect(pairedSampleMap[[thisSample]], markerMap[[thisMarker]])
	#testPairSeqList <- intersect(sampleMap[[thisSample]], markerMap[[thisMarker]])


seqTable <- data.frame()
localAlleleNames <- c("NA","NA","NA")
localAlleleFreqs <- c(0,0,0)

## go through all seq's mapped to this marker/sample pair.
## extract the corresponding sequence delimited by the top blast hits on the primers.  IS THIS THE BEST WAY?
##		Simple improvement: minimum blast hit length to primer to keep. 

## internal Function

recordNoSeqs <- function(summaryTable)  {		# to record no seqs before skipping out. 
		summaryRow <- data.frame(marker=thisMarker, sample=thisSample, numbSeqs=0,numbVars=0,
			varName.1="NA", varFreq.1= 0,
			varName.2="NA", varFreq.2= 0,
			varName.3="NA", varFreq.3= 0)
		summaryTable <- rbind(summaryTable, summaryRow)
		return(summaryTable)
}



if(length(testPairSeqList) < 1) {
	#summaryList[[thisMarker]][[thisSample]] <- NA	
	summaryTable  <- recordNoSeqs(summaryTable)
	next ;	# skip to next sample
} 


seqTable <- getSubSeqsTable(thisMarker, thisSample, pairedSampleMap, fMarkerMap,rmarkerMap, markerSeq)


# if no sequences returned, nothing to process. 
if(nrow(seqTable) < 1 )  {
	summaryTable  <- recordNoSeqs(summaryTable)
	#summaryList[[thisMarker]][[thisSample]] <- NA	
	next ;		# go to next sample.
}




#localSequenceMap <- split(seqTable[,2], seqTable[,1])

#localSequenceCount <- lapply(localSequenceMap , length)  # list named by sequence with counts.
#localSequenceCount <- localSequenceCount[order(as.numeric(localSequenceCount), decreasing=T)]

## test if variants are novel. 
## Give allele names?  
## Do with first three for now. 


alToRecord <- min(3,nrow(seqTable))
if(alToRecord > 0)  {
	for (a in 1:alToRecord )  {
		if(is.null(variantList[[seqTable$var[a]]]))  {    	# novel
			alleleName <- paste(thisMarker, alleleCount,sep=".")	
			variantList[[seqTable$var[a]]] <- alleleName
			localAlleleNames[a] <- alleleName 
			localAlleleFreqs[a] <- seqTable$count[a]
			alleleCount <- alleleCount + 1
		} else  {										# pre-existing alllele
			localAlleleNames[a] <- variantList[[seqTable$var[a]]]
			localAlleleFreqs[a] <- seqTable$count[a]		
		}
	}
}


# sequence correction?  


# compile stats

if(nrow(seqTable) >0 )  {	# cannot allow assignment from empty list as messes up class of list for remaining iterations
	summaryList[[thisMarker]] <- list()
	summaryList[[thisMarker]][[thisSample]] <- seqTable
}

summaryRow <- data.frame(marker=thisMarker, sample=thisSample, numbSeqs=sum(seqTable$count),numbVars=nrow(seqTable),
			varName.1=localAlleleNames[1], varFreq.1= localAlleleFreqs[1],
			varName.2=localAlleleNames[2], varFreq.2= localAlleleFreqs[2],
			varName.3=localAlleleNames[3], varFreq.3= localAlleleFreqs[3])
summaryTable <- rbind(summaryTable, summaryRow)

#sequence count across samples? 
# need to sum from summaryTable or from summaryList.
#markerSequenceCount <- 
#as.list(colSums(merge(m, n, all = TRUE), na.rm = TRUE))  # not working
localSequenceCount <- as.list(seqTable$count)
names(localSequenceCount) <- seqTable$var
markerSequenceCount   <- as.list(colSums(merge(markerSequenceCount  , localSequenceCount,  all = TRUE), na.rm = TRUE))
# might need to instantiate the markerSequenceCount if empty. 



}  # end of sample loop

markerSampleList[[thisMarker]] <- summaryTable
runSummaryRow <- data.frame(marker=thisMarker, assignedSeqs=sum(summaryTable$numbSeqs), assignedVariants=sum(summaryTable$numbVars), 
				minVariantLength=min(nchar(names(markerSequenceCount))), 
				maxVariantLength=max(nchar(names(markerSequenceCount))),
				minAlleleLength=min(nchar(names(variantList))), maxAlleleLength=max(nchar(names(variantList))))
runSummaryTable <- rbind(runSummaryTable, runSummaryRow)
if(length(variantList) > 0)  {
	alleleDb[[thisMarker]] <- variantList   # LATER: separate lists for alleles and variants?
}
}  # end of marker loop
	
	localMlgtResult <- new("mlgtResult", designObject,  runSummaryTable=runSummaryTable , alleleDb=alleleDb, markerSampleList=markerSampleList)
	return(localMlgtResult)

}  # end of mlgt function

###########################

setMethod("mlgt","mlgtDesign", definition=mlgt.mlgtDesign)

