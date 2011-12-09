

### TO DO

# Load reference alleles from DB/flatfile
# Construct 'comparable amplicon tree' from target markers and reference alleles. 
# 	problematic if differences between alleles are indels. 
# 	What about missing intronic sequence? Need to know source of reference allele. 
# 	Perhaps also need to mark projects amplicon as genomic or mRNA. 
# Compare with existing alleles
# checks for uniqueness of raw read IDs
# checks for uniqueness of marker IDs, MIDs
# Reduce I/O. Pipe or stdin/stdout from auxillary progs.
#	command <- "C:/Users/Public/Apps/Blast/bin/fastacmd.exe -p F -I  -d inputSeqs"
#	test <- readLines(pipe(command))
# Or set up proper MPI?
# Functions.
# Iterate across all loci/sample pairs.
# Handle/record lack of hits/sequences. 
# Check for really dodgy sequences. Differences in trimming?
# 	range of allele sizes. 
# set up directories to preserve results - fasta seqs dbs, alignments, summarytables.
# Log of files already processed as MD5_Checksum?
# Work RUN/PROJECT IDs. Could be used to test if a run has already been analysed. 
# Stats on variant frequencies
# Genotyping and genotyping parameters. (use "method" on a common structure so new methods can be updated)
# 	Assign alleles as iterate through marker/samples (i.e. ignore info on appearance of seqs across runs.)? 
# 	Or process all samples then go back and choose alleles and extract.
# allow use of BLAST+ 
# Test for hits to multiple MIDs (and remove). Limit valid hits to within X bp of start?
# Specify design of amplicons: number and type of MIDs, number and type of primers. Other seqs. 
# Sort out whether to use primer hits in assigning markers? It currently IS used indirectely by selection of subsequence with good hits to the correct primers.
# Fix warning from fastacmd caused by submission without accession list. MIGHT NOW BE FIXED
# make each marker run return an object. Good structure to run in parallel.
# SPEED: Do one massive alignment per marker? Might save time. Oringinal seqs with counts and aligned subseqs can be re-merged using orginal seq as name in sample/marker specific tables. 
# 	what if dodgy seqs in the mix? Will the alignment be any good?
# 	Might be worth throwing known alleles into the mix!
# MUSCLE alignments with dodgy ends e.g. A_E3.MID-13 - left a lone C at start of alignment because end of primer is a bit variable amogst seqs. 
# 	how to quickly test for and fix this problem?  (proportion of ---- before and after subSeq. 
# Decide on whether to do everything in a project/run directory. Blast DBs are easiest if located in working directory. 


### RECENTLY DONE



#################### Functions

getTopBlastHits <- function(blastTableFile)  {		# returns the first hit for each query in the table. May now be partially redundant if selecting for number of blast hits returned..
	blastResults <- read.delim(blastTableFile, header=F)
	names(blastResults) <- c("query", "subject", "percentId", "aliLength", "mismatches", "gapOpenings", "q.start","q.end", "s.start","s.end", "p_value", "e_value")
	topHits <- blastResults[match(unique(blastResults$query), blastResults$query),]
}


setUpBlastDb <- function(inputFastaFile, formatdbPath, blastdbName, indexDb="F")  {
	formatdbCommand <- paste(formatdbPath, "-i", inputFastaFile ,  "-p F -o", indexDb ,"-n", blastdbName)
	system(formatdbCommand)
}

getSubSeqs <- function(testPairSeqList, fTopHits, rTopHits, bothPrimerHits)  {
	seqGroup <- list()
	listWithPrimerHits <- intersect(testPairSeqList , bothPrimerHits)
	queryList <- paste(listWithPrimerHits, collapse=",")

	fShort <- fTopHits[match(listWithPrimerHits ,fTopHits$query),]
	rShort <- rTopHits[match(listWithPrimerHits ,rTopHits$query),]
	# could put check in here that all subject == thisMarker.

	primerHitTable <- merge(fShort, rShort, by=c("query", "subject"),suffixes = c(".f",".r"))
	primerHitTable$strand <- ifelse(primerHitTable$q.start.f  > primerHitTable$q.start.r, 2,1)
	primerHitTable$subStart <- ifelse(primerHitTable$strand == 1, primerHitTable$q.end.f,primerHitTable$q.end.r)
	primerHitTable$subEnd <- ifelse(primerHitTable$strand == 2, primerHitTable$q.start.f,primerHitTable$q.start.r)

	listWithBothGoodHits <- primerHitTable$query
	if(length(listWithBothGoodHits) < 1) {	# no remaining sequence IDs. Do not run fastacmd.
		return(seqGroup)
	}
 	# pull all these seqs from the blastDb together.
	queryList <- paste(listWithBothGoodHits , collapse=",")
	fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-o", subseqFileName,  "-s", queryList)
	system(fastacmdCommand)
	#seqGroup[[seqId]] <- read.fasta(subseqFileName, as.string=T)
	tempSeqGroup <- read.fasta(subseqFileName, as.string=T)
	# returned with "lcl|" prefix
	names(tempSeqGroup ) <-  sub("lcl\\|", "", names(tempSeqGroup ))


	# for each seq, extract the relevant portion, reverse if strand = 2
	for(seqId in listWithBothGoodHits)  {
		tableRow <- match(seqId, primerHitTable$query)
		strand <- primerHitTable$strand[tableRow]
		subStart <- primerHitTable$subStart[tableRow]
		subEnd <- primerHitTable$subEnd[tableRow]
		seq <- unlist(getSequence(tempSeqGroup[seqId], as.string=T))
		if(strand == 2)  {
			seq <- c2s(rev(comp(s2c((seq)))))
		}
		seqGroup[seqId] <- seq 
	}
	return(seqGroup)
}






####################

# seqinr http://cran.r-project.org/web/packages/seqinr/index.html
library(seqinr)


# system info
#formatdbPath
#blastAllPath
#musclePath
#params
#C:\Temp
# need to use correct system delimiter "\\" (Windows) or "/" (Windows or Unix)
formatdbPath <- "C:/Users/Public/Apps/Blast/bin/formatdb.exe"
fastacmdPath <- "C:/Users/Public/Apps/Blast/bin/fastacmd.exe"
blastAllPath <- "C:/Users/Public/Apps/Blast/bin/blastall.exe"
musclePath <- "C:/Users/Public/Apps/Muscle/muscle3.8.31_i86win32.exe"

#sysSep <- ifelse

###################

### set paramters. 

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt"
setwd( analysisDir )
projectName <- "testProject"
runName <- "testRun"

# load  marker sequences

# load sample barcodes/MIDs

sampleFile <- "samples_ALL.tab"
sampleTable <- read.delim(sampleFile, sep=";")

### test and create directories for this project/run.
overWriteBaseData <- FALSE
runPath <- paste(analysisDir, projectName, runName, sep="/")
if(file.exists(runPath))  {
	overWriteBaseData <- readline("A folder for this run already exists, do you want to write over this data? (")
} else {
	dir.create(runPath, recursive=T)
}




## could easy make this a function
## need to slim down to unique pairs (e.g. once for each MID)
fastaFile <- paste(runPath,"fTags.fasta", sep="/")
tempTable <- unique(subset(sampleTable, select=c(tag_f, sample_name)))
write.fasta(as.list(tempTable$tag_f), tempTable$sample_name, file=fastaFile )
#blastdbName <- "fTags"
#formatdbCommand <- paste(formatdbPath, "-i", fastaFile,  "-p F -n", blastdbName)
#system(formatdbCommand)
setUpBlastDb(fastaFile , formatdbPath, blastdbName="fTags")


fastaFile <- "rTags.fasta"
tempTable <- unique(subset(sampleTable, select=c(tag_r, sample_name)))
write.fasta(as.list(tempTable$tag_r), tempTable$sample_name, file=fastaFile )
blastdbName <- "rTags"
formatdbCommand <- paste(formatdbPath, "-i", fastaFile,  "-p F -n", blastdbName)
system(formatdbCommand)


### primers as blast DB
fastaFile <- "fPrimers.fasta"
tempTable <- unique(subset(sampleTable, select=c(primer_f, marker)))
write.fasta(as.list(tempTable$primer_f), tempTable$marker, file=fastaFile )
blastdbName <- "fPrimers"
formatdbCommand <- paste(formatdbPath, "-i", fastaFile,  "-p F -n", blastdbName)
system(formatdbCommand)

fastaFile <- "rPrimers.fasta"
tempTable <- unique(subset(sampleTable, select=c(primer_r, marker)))
write.fasta(as.list(tempTable$primer_r), tempTable$marker, file=fastaFile )
blastdbName <- "rPrimers"
formatdbCommand <- paste(formatdbPath, "-i", fastaFile,  "-p F -n", blastdbName)
system(formatdbCommand)

## markers as blast DB
fastaFile <- "HLA_MARKERS.fasta"
blastdbName <- "markerSeqs"
formatdbCommand <- paste(formatdbPath, "-i", fastaFile,  "-p F -n", blastdbName)
system(formatdbCommand)


## input as blast DB with index (useful for fast sub-sequence retrieval?)
inputFastaFile <- "3.TCA.454Reads.fna"		# NCL Run_756948
#blastdbName <- "inputSeqs"
#formatdbCommand <- paste(formatdbPath, "-i", inputFastaFile ,  "-p F -o T -n", blastdbName)	# indexed
#system(formatdbCommand)
setUpBlastDb(inputFastaFile, formatdbPath, blastdbName="inputSeqs", indexDb="T") 


#### blasting input seqs against markers, primers and MIDs.

#write.fasta(as.list(sampleTable$tag_r), sampleTable$sample_name, file="rTags.fasta")

fastaFile <- "fewReads.fasta"
fastaFile <- "moreTestSeqs.fasta"
fastaFile <- inputFastaFile 

## Blast against MIDs. reduce blast word size accordingly. Need exact match for full length or forward and reverse tags.

## If f and r tags are identical, don't need to do both?
#dbPath <- paste(getwd(), "rTags", sep)
rMinTagSize <- 10	# limit blast hits to perfect matches.
dbName <- "rTags"
blastOutFile <- "blastOut.rTags.tab"
blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile , "-W", rMinTagSize  ,"-m 8 -b 2 -S 3 -o", blastOutFile )
#blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile , "-W", rMinTagSize  ,"-m 8 -o", blastOutFile )	# working
system(blastCommand )

fMinTagSize <- 10	# limit blast hits to perfect matches.
dbName <- "fTags"
blastOutFile <- "blastOut.fTags.tab"
blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile , "-W", fMinTagSize ,"-m 8 -b 2 -S 3 -o", blastOutFile )
system(blastCommand )


## blast against primers.
#fastaFile <- "fewReads.fasta"
dbName <- "rPrimers"
blastOutFileName <- "blastOut.rPrimers.tab"
blastOutFile <- paste(runPath, blastOutFileName, sep="/")
blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile, "-W", 11 ,"-m 8 -b 1 -o", blastOutFile )
#blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile, "-W", 11 ,"-m 8 -o", blastOutFile )
system(blastCommand )

#fastaFile <- "fewReads.fasta"
dbName <- "fPrimers"
blastOutFileName <- "blastOut.fPrimers.tab"
blastOutFile <- paste(runPath, blastOutFileName, sep="/")
blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile, "-W", 11 ,"-v", 4 , "-m 8 -b 1  -o", blastOutFile )
system(blastCommand )

## blast against markers
#fastaFile <- "fewReads.fasta"
dbName <- "markerSeqs"
blastOutFileName <- "blastOut.markers.tab"
blastOutFile <- paste(runPath, blastOutFileName, sep="/")
# -b to restrict reported hits to 1. ### Could put hit quality filter in here or later. 
blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", fastaFile, "-W", 11 , "-m 8 -b 1 -o", blastOutFile )		
system(blastCommand )



## is blast 'v' useful. Trying to reduce number of reported hits. Only need top 2 or even less. 

## Need to make DB of marker sequences. Use primers too?  


# load run data (one seq at a time?). And store efficiently?
# there may be some identical sequences (even with random ends?). 



# assign sequences to samples
# 100% full length match to a single MID (or combination of MIDs if different ends are used)

##read blast result table.


# assign sequneces to markers. 
# Top BLAST hit with minimum level of identity
#topHits <- blastResults[match(unique(blastResults$query), blastResults$query),]
topHits <- getTopBlastHits("blastOut.markers.tab")		# redundant method if only storing top hit for each anyway.
markerMap <- split(as.character(topHits[,1]), topHits[,2])
## TO BE REPLACED BY:
	topHits <- getTopBlastHits("blastOut.markers.tab")
	topHits$strand <- ifelse(topHits$s.end > topHits$s.start, 1,2)
	fMarkerMap <- split(as.character(topHits$query[topHits$strand==1]), topHits$subject[topHits$strand==1])
	rMarkerMap <- split(as.character(topHits$query[topHits$strand==2]), topHits$subject[topHits$strand==2])



# read results from samples and make marker map. 

#blastResults <- read.delim("blastOut.rTags.tab", header=F)
#names(blastResults) <- c("query", "subject", "percentId", "aliLength", "mismatches", "gapOpenings", "q.start","q.end", "s.start","s.end", "p_value", "e_value")
#topHits <- blastResults[match(unique(blastResults$query), blastResults$query),]


topHits <- getTopBlastHits("blastOut.rTags.tab")
sampleMap <- split(as.character(topHits[,1]), topHits[,2])
## NEED TO MAKE SAMPLEMAP WITH HITS TO MID IN BOTH FORWARD AND REVERSE STRANDS like marker hits are split.
## Requires retention of 2 blast hits per sequence.
	topSampleHits <- getTopBlastHits("blastOut.rTags.tab")
	topSampleHits <- read.delim("blastOut.rTags.tab", header=F)
	names(topSampleHits ) <- c("query", "subject", "percentId", "aliLength", "mismatches", "gapOpenings", "q.start","q.end", "s.start","s.end", "p_value", "e_value")
	topSampleHits$strand <- ifelse(topSampleHits$s.end > topSampleHits$s.start, 1,2)
	fSampleMap <- split(as.character(topSampleHits$query[topSampleHits$strand==1]), topSampleHits$subject[topSampleHits$strand==1])
	rSampleMap <- split(as.character(topSampleHits$query[topSampleHits$strand==2]), topSampleHits$subject[topSampleHits$strand==2])
	# combind sampleMaps to give sequences with MIDs in both orientations.
	pairedSampleMap <- lapply(names(fSampleMap), FUN=function(x) intersect(fSampleMap[[x]], rSampleMap[[x]]))
	names(pairedSampleMap) <- names(fSampleMap)
	# test this worked correctly:  intersect(pairedSampleMap[[thisSample]], fSampleMap[[thisSample]])
	
#sampleMapClean <- intersect(fSampleMap, rSampleMap)


#integrate marker map and sample map to get fully assigned sequences. 
#testPairSeqList <- intersect(sampleMap[["MID-24"]], markerMap[["A_E2"]])

# Remove primer sequences and MIDs? Primers must be closest to amplicon. 
# In our study, primers are degenerate, the sequence corresponding to the primer is unreliable. 

# TO BE TAKEN OUT
fTopHits <- getTopBlastHits("blastOut.fPrimers.tab")	# is this a redundant method?
rTopHits <- getTopBlastHits("blastOut.rPrimers.tab")

bothPrimerHits <- intersect(fTopHits$query,rTopHits$query) # (almost all)

# subsequence required is marked by blast hits. But start/ end may either be in r or f hits. 
# sequences may fail if they don't have perfect hits in both.
# TO DO Might be best to do a merge on the two primer hit tables.
#primerCoords <- unlist((c(fTopHits[i,c("q.start","q.end")], rTopHits[i,c("q.start","q.end")])))
#strand <- ifelse(primerCoords[3] < primerCoords[2], 2,1)	# 2 -ve strand
#primerCoords <- sort(primerCoords)		# replace with tapply
#subCoords <- c(primerCoords[2]+1, primerCoords[3]-1)

#### is it best to extract full sequence then make subsequence or extract subsequence straighaway.

#retSeq <- as.character(fTopHits$query[i])
subseqFileName <- paste(runPath, "subseqFile.fasta",sep="/")

#fastacmdCommand  <- paste(fastacmdPath, "-p F -d", blastdbName , "-s", retSeq , "-L", paste(subCoords[1],subCoords[2],sep=","))
#fastacmdCommand  <- paste(fastacmdPath, "-p F -d", blastdbName , "-o", subseqFileName,  "-s", retSeq , "-L", paste(subCoords[1],subCoords[2],sep=","), "-S", strand)
#system(fastacmdCommand)
#retSeqObject <- read.fasta(file=stdin(system(fastacmdCommand)))



#seqGroup[i] <- read.fasta(subseqFileName, as.string=T)
#subseqFileName <- "subseqFile.fasta"
minAmpliconLength <- 40

#########################################
#### iterations
##########################################

markerSampleList <- list()
runSummaryTable <- data.frame()
alleleDb <- list()


for(thisMarker in names(markerMap)) {
#for(thisMarker in names(markerMap)[16:19]) {	# temp to finish off

print(thisMarker)
#thisMarker <- "DQA1_E2"

## might need to combine all these to return a single item.
summaryList <- list()
summaryTable <- data.frame()
markerSequenceCount <- list("noSeq"=0)		#  BUG? requires some data otherwise won't sum properly with localSequenceCount.
alleleList <- list() 
variantList <- list()
alleleCount <- 1


for(thisSample in names(sampleMap)) {
	#print(thisSample)

	
	#testPairSeqList <- intersect(pairedSampleMap[[thisSample]], markerMap[[thisMarker]])
	testPairSeqList <- intersect(sampleMap[[thisSample]], markerMap[[thisMarker]])


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



seqGroup <- getSubSeqs(testPairSeqList, fTopHits, rTopHits, bothPrimerHits)


# if no sequences returned, nothing to process. 
if(length(seqGroup) < 1 )  {
	summaryTable  <- recordNoSeqs(summaryTable)
	#summaryList[[thisMarker]][[thisSample]] <- NA	
	next ;		# go to next sample.
}

# looks pretty good, some very long seqs though. 
#str(seqGroup[[1]])

seqTable <- cbind(as.character(unlist(seqGroup)), names(unlist(seqGroup)))



localSequenceMap <- split(seqTable[,2], seqTable[,1])

localSequenceCount <- lapply(localSequenceMap , length)  # list named by sequence with counts.
localSequenceCount <- localSequenceCount[order(as.numeric(localSequenceCount), decreasing=T)]


# Condense number of sequence and map seqs to variants.  
#lapply(localSequenceMap , length)
localVariants <- names(localSequenceCount)

localVariantFileName <- paste("thisRun", thisMarker, thisSample, "variants.fasta",sep=".")  #"localVariants.fasta"
localVariantFile <- paste(runPath, localVariantFileName, sep="/")
write.fasta(as.list(localVariants) ,localVariants,file=localVariantFile ,open="w")

## test if variants are novel. 
## Give allele names?  
## Do with first three for now. 


alToRecord <- min(3,length(localSequenceCount))
if(alToRecord > 0)  {
	for (a in 1:alToRecord )  {
		if(is.null(variantList[[names(localSequenceCount[a])]]))  {    	# novel
			alleleName <- paste(thisMarker, alleleCount,sep=".")	
			variantList[[names(localSequenceCount[a])]] <- alleleName
			localAlleleNames[a] <- alleleName 
			localAlleleFreqs[a] <- localSequenceCount[[a]]
			alleleCount <- alleleCount + 1
		} else  {										# pre-existing alllele
			localAlleleNames[a] <- variantList[[names(localSequenceCount[a])]]
			localAlleleFreqs[a] <- localSequenceCount[[a]]		
		}
	}
}
# align sequences and store alignment.  # Is this needed for genotyping?  Maybe if some small differences in common seqs. Yes, if want to correct the sequences.
# Will need way to nicely visualise alignments. 

localAlignFileName <- paste("thisRun", thisMarker, thisSample, "align.fasta",sep=".")  #"localAlign.fasta"
localAlignFile <- paste(runPath, localAlignFileName, sep="/")
muscleCommand <- paste(musclePath, "-in", localVariantFile , "-out", localAlignFile , "-diags -quiet" )
system(muscleCommand)
# can read alignment in.
#localAlignment <- read.alignment(localAlignFile , format="fasta")


# sequence correction?  


# compile stats

if(length(localSequenceCount) >0 )  {	# cannot allow assignment from empty list as messes up class of list for remaining iterations
	summaryList[[thisMarker]] <- list()
	summaryList[[thisMarker]][[thisSample]] <- localSequenceCount 
}
summaryRow <- data.frame(marker=thisMarker, sample=thisSample, numbSeqs=nrow(seqTable),numbVars=length(localSequenceCount),
			varName.1=localAlleleNames[1], varFreq.1= localAlleleFreqs[1],
			varName.2=localAlleleNames[2], varFreq.2= localAlleleFreqs[2],
			varName.3=localAlleleNames[3], varFreq.3= localAlleleFreqs[3])
summaryTable <- rbind(summaryTable, summaryRow)

#sequence count across samples? 
# need to sum from summaryTable or from summaryList.
#markerSequenceCount <- 
#as.list(colSums(merge(m, n, all = TRUE), na.rm = TRUE))  # not working
markerSequenceCount   <- as.list(colSums(merge(markerSequenceCount  , localSequenceCount,  all = TRUE), na.rm = TRUE))
# might need to instantiate the markerSequenceCount if empty. 



}  # end of sample loop

markerSampleList[[thisMarker]] <- summaryTable
runSummaryRow <- data.frame(marker=thisMarker, assigneSeqs=sum(summaryTable$numbSeqs), assignedVariants=sum(summaryTable$numbVars), 
				minVariantLength=min(nchar(names(markerSequenceCount))), 
				maxVariantLength=max(nchar(names(markerSequenceCount))),
				minAlleleLength=min(nchar(names(variantList))), maxAlleleLength=max(nchar(names(variantList))))
runSummaryTable <- rbind(runSummaryTable, runSummaryRow)
if(length(variantList) > 0)  {
	alleleDb[[thisMarker]] <- variantList   # LATER: separate lists for alleles and variants?
}
}  # end of marker loop

####################################


# save.image("C:\\Users\\dave\\HalfStarted\\mlgt\\testProjectData.RData")

####################################



# http://stackoverflow.com/questions/5783241/variable-names-are-limited-to-256-bytes
# Running under R 2.13

length(variantList)
min(nchar(names(variantList)))		# shortest of variants retained as potential alleles
max(nchar(names(variantList)))		# longest of variants retained as potential alleles

min(nchar(names(markerSequenceCount)))	# shortest extracted variant
max(nchar(names(markerSequenceCount)))	# longest extracted variant

length(markerSequenceCount )
nrow(summaryTable)
sum(summaryTable$numbVars)
sum(summaryTable$numbSeqs)



addToSeqCount <- function(oldCount, newCount)  {
	

}


# How to combine results from one pair across all pairs or within a locus. Sum/merge the lists?






# choose alleles if possible


# output genotypes with alleles.



# match allleles to allele database. 



stopifnot(FALSE)
################DEVELOPMENT


### could be failing because merge collapses the results when the two list are identical to begin with.
markerSequenceCount <- list("noSeq"=0, "diff"=5)	# nope still not working. Adds the new names but doesn't summ the common values. 
markerSequenceCount   <- as.list(colSums(merge( localSequenceCount, markerSequenceCount  , all = TRUE), na.rm = TRUE))
newSequenceCount <- localSequenceCount
as.list(colSums(merge( localSequenceCount, markerSequenceCount  , all = TRUE), na.rm = TRUE))

markerSequenceCount   <- as.list(colSums(merge( localSequenceCount, markerSequenceCount  , all = TRUE), na.rm = TRUE))

markerSequenceCount <- lapply(markerSequenceCount, as.integer)
str(markerSequenceCount)
str(localSequenceCount)
localSequenceCount<- lapply(localSequenceCount, as.numeric)


markerSequenceCount   <- as.list(colSums(merge( localSequenceCount, markerSequenceCount  , all = TRUE), na.rm = TRUE))

localSequenceCount + markerSequenceCount  

localSequenceCounts <- markerSequenceCount   
#is this just going to be easier with a data.frame?

localSequenceCount["diff"] <- 2


list1 <- list("x"=16, "y"=2)
list2 <- list("x"=14, "y"=7)

