

# seqGroup <- returnSubseqHits(testPairSeqList,fTopHits,rTopHits, minAmpliconLength=40, fastacmdPath)

returnSubseqHits <- function(testPairSeqList,fTopHits,rTopHits, minAmpliconLength=1, fastacmdPath  )	{
	seqGroup <-list()
	subseqFileName <- "returnSubseqHits.tempFile"
	for (i in 1:length(testPairSeqList)) {
		seqId <- testPairSeqList[i]
		primerCoords <- unlist((c(fTopHits[match(seqId,fTopHits$query),c("q.start","q.end")], rTopHits[match(seqId,rTopHits$query),c("q.start","q.end")])))
		if(length(na.omit(primerCoords)) != 4)  { 	# did not get positions for forward and reverse primers. 
			next ; 
		}
		strand <- ifelse(primerCoords[3] < primerCoords[2], 2,1)	# 2 -ve strand
		primerCoords <- sort(primerCoords)		# replace with tapply
		subCoords <- c(primerCoords[2]+1, primerCoords[3]-1)
		if((subCoords[2] - subCoords[1]) < minAmpliconLength ) {
			next ;
		}
		fastacmdCommand  <- paste(fastacmdPath, "-p F -d", "inputSeqs" , "-o", subseqFileName,  "-s", seqId, "-L", paste(subCoords[1],subCoords[2],sep=","), "-S", strand)
		system(fastacmdCommand)
		seqGroup[[seqId]] <- read.fasta(subseqFileName, as.string=T)
		
	}
	file.remove(subseqFileName)
	return(seqGroup)
}





#####################development
#getFrag


##mergedPrimerHits <- merge(fTopHits, rTopHits, by="query", 

# higher level specified.
bothPrimerHits <- intersect(fTopHits$query,rTopHits$query) # (almost all)
#testPairSeqList 


# local to marker sample.
testPairSeqList <- intersect(sampleMap[[thisSample]], markerMap[[thisMarker]])

listWithPrimerHits <- intersect(testPairSeqList , bothPrimerHits)
queryList <- paste(listWithPrimerHits, collapse=",")

fShort <- fTopHits[match(listWithPrimerHits ,fTopHits$query),]
rShort <- rTopHits[match(listWithPrimerHits ,rTopHits$query),]
# could put check in here that all subject == thisMarker.

primerHitTable <- merge(fShort, rShort, by=c("query", "subject"),suffixes = c(".f",".r"))
primerHitTable$strand <- ifelse(primerHitTable$q.start.f  > primerHitTable$q.start.r, 2,1)
primerHitTable$subStart <- ifelse(primerHitTable$strand == 1, primerHitTable$q.end.f,primerHitTable$q.end.r)
primerHitTable$subEnd <- ifelse(primerHitTable$strand == 2, primerHitTable$q.start.f,primerHitTable$q.start.r)

# pull all these seqs from the blastDb together.
#s
queryList <- paste(listWithPrimerHits, collapse=",")
	#fastacmdCommand  <- paste(fastacmdPath, "-p F -d", "inputSeqs" , "-o", subseqFileName,  "-s", seqId, "-L", paste(subCoords[1],subCoords[2],sep=","), "-S", strand)
	fastacmdCommand  <- paste(fastacmdPath, "-p F -d", "inputSeqs" , "-o", subseqFileName,  "-s", queryList)

	system(fastacmdCommand)
	#seqGroup[[seqId]] <- read.fasta(subseqFileName, as.string=T)
	tempSeqGroup <- read.fasta(subseqFileName, as.string=T)

# for each seq, extract the relevant portion, revers if strand =



#unusedSeqs <- setdiff(testPairSeqList, listWithPrimerHits )


## TRYING TO SOLVE TRIMMING PROBLEM
## How to get reliable amplicon that corresponds to the amplicon and only the amplicon?
### test gapped blast on primers?

#standara
fastaFile <- "fewReads.fasta"
dbName <- "fPrimers"
blastOutFileName <- "blastOut.fPrimers.test1.tab"
blastOutFile <- paste(runPath, blastOutFileName, sep="/")
blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile, "-W", 11 ,"-v", 4 , "-m 8 -b 1  -o", blastOutFile )
system.time(system(blastCommand ))

#gapped
fastaFile <- "fewReads.fasta"
dbName <- "fPrimers"
blastOutFileName <- "blastOut.fPrimers.test2.tab"
blastOutFile <- paste(runPath, blastOutFileName, sep="/")
blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile, "-W", 11 ,"-v", 4 , "-m 8 -g T -b 1  -o", blastOutFile )
system.time(system(blastCommand ))


#small word - getting back several sub-hits of 7 bp.
fastaFile <- "fewReads.fasta"
dbName <- "fPrimers"
blastOutFileName <- "blastOut.fPrimers.test3.tab"
blastOutFile <- paste(runPath, blastOutFileName, sep="/")
blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile, "-W", 7 ,"-v", 4 , "-m 8 -b 1  -o", blastOutFile )
system.time(system(blastCommand ))

#small word and gapped
fastaFile <- "fewReads.fasta"
dbName <- "fPrimers"
blastOutFileName <- "blastOut.fPrimers.test4.tab"
blastOutFile <- paste(runPath, blastOutFileName, sep="/")
blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile, "-W", 7 ,"-v", 4 , "-m 8 -g T -b 1  -o", blastOutFile )
system.time(system(blastCommand ))

#reduced mismatch penalty (less negative) - no change to alignments only worse p-values!
fastaFile <- "fewReads.fasta"
dbName <- "fPrimers"
blastOutFileName <- "blastOut.fPrimers.test5.tab"
blastOutFile <- paste(runPath, blastOutFileName, sep="/")
blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile, "-W", 11 ,"-v", 4 , "-m 8 -q -2 -b 1  -o", blastOutFile )
system.time(system(blastCommand ))

#################
# could try making alignment with reference sequences and trimming based ends of reference.


markerList <- read.fasta("HLA_MARKERS.fasta",as.string=T)

markerSeq <- unlist(getSequence(markerList[[thisMarker]],as.string=T))
# Condense number of sequence and map seqs to variants.  
#lapply(localSequenceMap , length)
localVariants <- names(localSequenceCount)

localVariantFileName <- paste("test", thisMarker, thisSample, "variants.fasta",sep=".")  #"localVariants.fasta"
localVariantFile <- paste(runPath, localVariantFileName, sep="/")
write.fasta(as.list(c(markerSeq,localVariants)) ,c(thisMarker,localVariants),file=localVariantFile ,open="w")

localAlignFileName <- paste("test", thisMarker, thisSample, "align.fasta",sep=".")  #"localAlign.fasta"
localAlignFile <- paste(runPath, localAlignFileName, sep="/")
muscleCommand <- paste(musclePath, "-in", localVariantFile , "-out", localAlignFile , "-diags -quiet" )
system(muscleCommand)



########################## DEV VERSION

# This alignment has ----- exactly at the primer limits.  Could be good. 
# Probably don't need to use primers at all and don't need to extract subseqs, only seqs, in correct orientation.
# How to extract ref, trim others  and keep track of counts for each. 
# Also, should alignment be across all seqs of a marker? Prob not. then becomes run dependent. Do within sample only.


# need to rework blastHit table to separate hits on forward and reverse strand.
# TODO: could use same method to improve MID hits and check both ends are marked. (or use different MIDs for each end).
	topHits <- getTopBlastHits("blastOut.markers.tab")
	topHits$strand <- ifelse(topHits$s.end > topHits$s.start, 1,2)
	fMarkerMap <- split(as.character(topHits$query[topHits$strand==1]), topHits$subject[topHits$strand==1])
	rMarkerMap <- split(as.character(topHits$query[topHits$strand==2]), topHits$subject[topHits$strand==2])

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
	rawVariantFileName <- paste("test", thisMarker, thisSample, "raw.variants.fasta",sep=".")  #"localVariants.fasta"
	rawVariantFile <- paste(runPath, rawVariantFileName, sep="/")
	write.fasta(as.list(c(markerSeq,as.character(rawSeqCountTable$rawSeq))) ,c(thisMarker,as.character(rawSeqCountTable$rawSeq)),file=rawVariantFile ,open="w")

# Align all seqs with reference
	rawAlignFileName <- paste("test", thisMarker, thisSample, "raw.align.fasta",sep=".")  #"localAlign.fasta"
	rawAlignFile <- paste(runPath, rawAlignFileName, sep="/")
	muscleCommand <- paste(musclePath, "-in", rawVariantFile , "-out", rawAlignFile , "-diags -quiet" )
	system(muscleCommand)


# Extract portion corresponding to reference. 

	rawAlignment <- read.fasta(rawAlignFile, as.string=T)		# do not use read.alignment() - broken
	alignedMarkerSeq <- s2c(rawAlignment[[thisMarker]])
	subStart <- min(grep("-",alignedMarkerSeq ,invert=T))
	subEnd <- max(grep("-",alignedMarkerSeq ,invert=T))
	alignedSubSeqs <- lapply(rawAlignment, FUN=function(x)	substr(x[1], subStart, subEnd))
	subAlignFileName <- paste("test", thisMarker, thisSample, "sub.align.fasta",sep=".")  #"localAlign.fasta"
	subAlignFile <- paste(runPath, subAlignFileName , sep="/")
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
markerList <- read.fasta("HLA_MARKERS.fasta",as.string=T)
	

##########ITERATIONS

markerSampleList <- list()
runSummaryTable <- data.frame()
alleleDb <- list()


for(thisMarker in names(markerMap)) {
#for(thisMarker in names(markerMap)[1:2]) {	# temp to finish off

print(thisMarker)
#thisMarker <- "DQA1_E2"

## might need to combine all these to return a single item.
summaryList <- list()
summaryTable <- data.frame()
markerSequenceCount <- list("noSeq"=0)		#  BUG? requires some data otherwise won't sum properly with localSequenceCount.
alleleList <- list() 
variantList <- list()
alleleCount <- 1
markerSeq <- unlist(getSequence(markerList[[thisMarker]],as.string=T))

for(thisSample in names(pairedSampleMap)) {
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

####################################






