# require(seqinr)

setClass("mlgtDesign", 
	representation(
		projectName="character", 
		runName="character",
		markers="list",
		samples="character",
		fTags="list",
		rTags="list", 
		inputFastaFile="character", 
		markerBlastResults="character",
		fTagBlastResults="character",
		rTagBlastResults="character"
	)
)

setMethod("print", "mlgtDesign", definition= function(x, ...){
	cat("Design for mlgt run:\n")
	cat(paste("Project:",x@projectName,"\n"))
	cat(paste("Run:",x@runName,"\n"))
	cat(paste("Samples:",length(x@samples),"\n"))
	cat(paste("fTags:",length(x@fTags),"\n"))
	cat(paste("rTags:",length(x@rTags),"\n"))
	cat(paste("Markers:",length(x@markers),"\n"))
	#cat(paste(x[1:5]), "...\n")
	}
)


#setMethod("show", "mlgtDesign", definition= function(x, ...){
#	cat("Design for mlgt run:\n")
#	cat(paste("Project:",x@projectName,"\n"))
#	cat(paste("Run:",x@runName,"\n"))
#	cat(paste("Samples:",length(x@samples),"\n"))
#	cat(paste("fTags:",length(x@fTags),"\n"))
#	cat(paste("rTags:",length(x@rTags),"\n"))
#	cat(paste("Markers:",length(x@markers),"\n"))
#	#cat(paste(x[1:5]), "...\n")
#	}
#)



setClass("mlgtResult", 
	representation(
			runSummaryTable="data.frame",
			alleleDb="list" ,
			markerSampleList="list"
	),
	contains="mlgtDesign"
)





getTopBlastHits <- function(blastTableFile)  {		# returns the first hit for each query in the table. May now be partially redundant if selecting for number of blast hits returned..
	blastResults <- read.delim(blastTableFile, header=F)
	## Fields: 
	# Query id,Subject id,% identity,alignment length,mismatches,gap openings,q. start,q. end,s. start,s. end,e-value,bit score
	names(blastResults) <- c("query", "subject", "percent.id", "ali.length", "mismatches", "gap.openings", "q.start","q.end", "s.start","s.end", "e.value", "bit.score")
	topHits <- blastResults[match(unique(blastResults$query), blastResults$query),]
}



getSubSeqsTable <- function(thisMarker, thisSample, sampleMap, fMarkerMap,rMarkerMap, markerSeq)  {

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
	
	## 12th Dec 11. Edited following rawseq extraction because failed when using large dataset. Replaced '-s' (queryList) option with '-i' (inputfile)
	if(length(fPairSeqList ) > 0)  {
		#queryList <- paste(fPairSeqList , collapse=",")
		fIdFileName <- "fIdFile.txt"	
		write(fPairSeqList , file=fIdFileName )
		fRawSeqFileName <- "fRawSeqExtract.fasta"
		strand <- 1
		#fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-S", strand, "-o", fRawSeqFileName,  "-s", queryList)		
		fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-S", strand, "-o", fRawSeqFileName,  "-i", fIdFileName)	
		system(fastacmdCommand)
		fRawSeqs <- read.fasta(fRawSeqFileName , as.string=T)
	}


	if(length(rPairSeqList ) > 0)  {
		#queryList <- paste(rPairSeqList , collapse=",")	
		rIdFileName <- "rIdFile.txt"
		write(rPairSeqList , file=rIdFileName )		
		rRawSeqFileName <- "rRawSeqExtract.fasta"
		strand <- 2
		#fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-S", strand, "-o", rRawSeqFileName,  "-s", queryList)		
		fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-S", strand, "-o", rRawSeqFileName,  "-i", rIdFileName)
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
	rawVariantCount <- nrow(rawSeqCountTable)
	names(rawSeqCountTable) <- c("rawSeq", "rawCount")
	rawVariantFile <- paste("test", thisMarker, thisSample, "raw.variants.fasta",sep=".")  #"localVariants.fasta"
	#rawVariantFile <- paste(runPath, rawVariantFileName, sep="/")
	write.fasta(as.list(c(markerSeq,as.character(rawSeqCountTable$rawSeq))) ,c(thisMarker,as.character(rawSeqCountTable$rawSeq)),file=rawVariantFile ,open="w")

	# Align all seqs with reference
	rawAlignFile <- paste("test", thisMarker, thisSample, "raw.align.fasta",sep=".")  #"localAlign.fasta"
	#rawAlignFile <- paste(runPath, rawAlignFileName, sep="/")
	if(rawVariantCount > 800) {
		muscleCommand <- paste(musclePath, "-in", rawVariantFile , "-out", rawAlignFile , "-diags -quiet -maxiters 2" )	# faster alignment
		warning(paste("Using fast MUSCLE alignment for ", thisMarker, thisSample, rawVariantCount, "sequences\n"))
	} else {
		muscleCommand <- paste(musclePath, "-in", rawVariantFile , "-out", rawAlignFile , "-diags -quiet" )
	}
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

# Re-apply count of each seq.  There may be some duplicated subSeqs.
	combTable <- merge(rawSeqCountTable ,alignedSubTable , by="rawSeq", all.x=T)
	varCount <- by(combTable, as.character(combTable$subSeq), FUN=function(x) sum(x$rawCount))
	varCountTable <- data.frame(alignedVar=names(varCount), count=as.numeric(varCount))	
	varCountTable$var <- gsub("-","",varCountTable$alignedVar)
	varCountTable <- varCountTable[order(varCountTable$count,decreasing=T),]
# Make unique list, summing counts where same seq found. (easier in table than list).  
	# ?TODO?
	return(varCountTable)

}



#mlgt <- function(object) attributes(object)
#setGeneric("mlgt")

mlgt <- function(designObject) attributes(designObject)
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
	markerSeq <- unlist(getSequence(designObject@markers[[thisMarker]],as.string=T))

	for(thisSample in designObject@samples) {
		#for(thisSample in names(pairedSampleMap)[1:4]) {
			#print(thisSample)

			testPairSeqList <- intersect(pairedSampleMap[[thisSample]],union(fMarkerMap[[thisMarker]], rMarkerMap[[thisMarker]]))
			#testPairSeqList <- intersect(pairedSampleMap[[thisSample]], markerMap[[thisMarker]])
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


		seqTable <- getSubSeqsTable(thisMarker, thisSample, pairedSampleMap, fMarkerMap,rMarkerMap, markerSeq)


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
		# This line replaced. Not entirely tested the repurcussions. e.g. makeVarAlleleMap()?
		#alleleDb[[thisMarker]] <- variantList   # LATER: separate lists for alleles and variants? 
		#alleleDb[[thisMarker]] <- list(reference=as.SeqFastadna(markerSeq, thisMarker), alleleMap=variantList, inputAlleleCount = length(unlist(variantList)), uniqueSubAlleleCount=length(variantList))
		alleleDb[[thisMarker]] <- new("variantMap", reference=as.SeqFastadna(markerSeq, thisMarker), 
							variantSource=paste(designObject@projectName, designObject@runName,sep="."),
							variantMap=variantList, inputVariantCount = length(unlist(variantList)), uniqueSubVariantCount=length(variantList))
	}

}  # end of marker loop
	
	localMlgtResult <- new("mlgtResult", designObject,  runSummaryTable=runSummaryTable , alleleDb=alleleDb, markerSampleList=markerSampleList)
	return(localMlgtResult)

}  # end of mlgt function

###########################

setMethod("mlgt","mlgtDesign", definition=mlgt.mlgtDesign)


setUpBlastDb <- function(inputFastaFile, formatdbPath, blastdbName, indexDb="F")  {
	formatdbCommand <- paste(formatdbPath, "-i", inputFastaFile ,  "-p F -o", indexDb ,"-n", blastdbName)
	system(formatdbCommand)
}






#prepareMlgtRun <- function(object) attributes(object)
prepareMlgtRun <- function(designObject) attributes(designObject)
prepareMlgtRun <- function(designObject,projectName,runName, samples, markers,fTags,rTags, inputFastaFile, overwrite) attributes(designObject)

setGeneric("prepareMlgtRun")


prepareMlgtRun.listDesign <- function(projectName,runName, samples, markers,fTags,rTags, inputFastaFile,overwrite="prompt")  {
	designObject <- new("mlgtDesign", projectName=projectName, runName=runName, 
				samples=samples, markers=markers ,
				fTags=fTags, rTags=rTags, inputFastaFile=inputFastaFile)

	designObject <- prepareMlgtRun(designObject,overwrite=overwrite)
	
}

setMethod("prepareMlgtRun",
	signature(designObject="missing",projectName="character", runName="character", samples="character",markers="list", 
	fTags="list", rTags="list", inputFastaFile="character", overwrite="character"), 
	definition=prepareMlgtRun.listDesign)


prepareMlgtRun.mlgtDesign <- function(designObject, overwrite="prompt")  {
	cat(paste(designObject@projectName,"\n"))

	fTagsFastaFile <- paste(runPath,"fTags.fasta", sep="/")		
	rTagsFastaFile <- paste(runPath,"rTags.fasta", sep="/")
	markersFastaFile <- paste(runPath,"markers.fasta", sep="/")
	rTagsBlastOutFile <- paste(runPath,"blastOut.rTags.tab", sep="/")
	fTagsBlastOutFile <- paste(runPath,"blastOut.fTags.tab", sep="/")
	blastOutFileName <- "blastOut.markers.tab"
	markerBlastOutFile <- paste(runPath, blastOutFileName, sep="/")

	existingFiles <- (file.exists(fTagsFastaFile ) | file.exists(rTagsFastaFile ) | file.exists(markersFastaFile ) | 
					file.exists(rTagsBlastOutFile) | file.exists(fTagsBlastOutFile) | file.exists(markerBlastOutFile))
	if(existingFiles)  {
		overwrite <- tolower(overwrite)
		# set up new directories if required.
		if(overwrite == "prompt") overwrite <- readline("This folder already contains mlgt run files, do you want to write over this data? (yes/no)")

		if(!(overwrite=="y" | overwrite=="n" | overwrite=="yes" | overwrite=="no")) {stop(paste("Unrecognised value for overwrite:", overwrite))}
		overWriteBaseData <- switch(overwrite,
						"yes"=TRUE,
						"y"=TRUE,
						"no"=FALSE,
						"n"=FALSE)
		if(!overWriteBaseData) {stop("This folder already contains mlgt run files. Exiting")}
	}					

	runPath <- getwd()
	# set up blast DBs
	cat("Setting up BLAST DBs...\n")

	write.fasta(designObject@fTags, names(designObject@fTags), file=fastaFile )
	setUpBlastDb(fTagsFastaFile , formatdbPath, blastdbName="fTags")		# default is no index

	write.fasta(designObject@rTags, names(designObject@rTags), file=fastaFile )
	setUpBlastDb(rTagsFastaFile , formatdbPath, blastdbName="rTags")		# default is no index

	write.fasta(designObject@markers, names(designObject@markers), file=markersFastaFile )
	setUpBlastDb(fastaFile , formatdbPath, blastdbName="markerSeqs", indexDb="T") 

	## input as blast DB with index (useful for fast sub-sequence retrieval?)
	inputFastaFile <- designObject@inputFastaFile		#
	setUpBlastDb(inputFastaFile , formatdbPath, blastdbName="inputSeqs", indexDb="T") 

	# run preliminary blast
	cat("Running BLAST searches...\n")
	
	rMinTagSize <- 10	# limit blast hits to perfect matches.		### TODO: SET GLOBALLY
	dbName <- "rTags"

	blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", inputFastaFile , "-W", rMinTagSize  ,"-m 8 -b 2 -S 3 -o", rTagsBlastOutFile )
	system(blastCommand )

	fMinTagSize <- 10	# limit blast hits to perfect matches.
	dbName <- "fTags"

	blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", inputFastaFile , "-W", fMinTagSize ,"-m 8 -b 2 -S 3 -o", fTagsBlastOutFile )
	system(blastCommand )

	## blast against markers
	dbName <- "markerSeqs"

	blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", inputFastaFile , "-W", 11 , "-m 8 -b 1 -o", markerBlastOutFile )		
	system(blastCommand )

	designObject@markerBlastResults <- markerBlastOutFile
	designObject@fTagBlastResults <- fTagsBlastOutFile 
	designObject@rTagBlastResults <- rTagsBlastOutFile 

	return(designObject)

}	

##setMethod("prepareMlgtRun","mlgtDesign", definition=prepareMlgtRun.mlgtDesign)
#setMethod("prepareMlgtRun",signature(designObject="mlgtDesign"), definition=prepareMlgtRun.mlgtDesign)


#TODO: test if signature with overwrite="ANY" will work. First attempt, no.
setMethod("prepareMlgtRun",
	signature(designObject="mlgtDesign",projectName="missing", runName="missing", samples="missing",markers="missing", 
	fTags="missing", rTags="missing", inputFastaFile="missing", overwrite="character"),
	 definition=prepareMlgtRun.mlgtDesign)


#setGeneric("prepareMlgtRun","mlgtDesign", definition=prepareMlgtRun.mlgtDesign)

######################### genotyping & allele matching


setClass("genotypeCall", 
	representation(
		projectName="character", 
		runName="character",
		marker="character",
		genotypeTable="data.frame",
		callMethod="character",
		callParameters="list",
		mappedToAlleles="logical",
		alleleDbName="character"
	)
)


makeVarAlleleMap <- function(allele.variantMap, variant.variantMap)  {
			varAlleleMap <- data.frame()
			# to use stopifnot, need to ensure no empty maps are passed. Currently this is happening so taking out this check.
			#stopifnot(allele.variantMap@reference == variant.variantMap@reference)
			knownAlleleTable <- data.frame(alleleSeq=names(allele.variantMap@variantMap), knownAlleles=as.character(allele.variantMap@variantMap))


			dataAlleleTable <-  data.frame(alleleSeq=names(variant.variantMap@variantMap), varNames=as.character(variant.variantMap@variantMap))
			
			varAlleleMap <- merge(knownAlleleTable, dataAlleleTable , by="alleleSeq")
}

# deprecated and/or defunct
makeVarAlleleMap.list <- function(alleleDb, varDb,alleleMarkers=names(alleleDb),varMarkers=names(varDb))  {
			knownAlleleTable <- data.frame()
			for(thisMarker in alleleMarkers)  {
				knownAlleleTable <- rbind(knownAlleleTable , data.frame(alleleSeq=names(alleleDb[[thisMarker]]@alleleMap), knownAlleles=as.character(alleleDb[[thisMarker]]@alleleMap)))
			}

			dataAlleleTable <- data.frame()
			for(thisMarker in varMarkers)  {
				# this first line is what it SHOULD be like, once mlgtResult is updated to match the alleleDb format. Need new class: alleleDb
				dataAlleleTable <- rbind(dataAlleleTable , data.frame(alleleSeq=names(varDb[[thisMarker]]@alleleMap), varNames=as.character(varDb[[thisMarker]]@alleleMap)))
				#dataAlleleTable <- rbind(dataAlleleTable , data.frame(alleleSeq=names(varDb[[thisMarker]]), varNames=as.character(varDb[[thisMarker]])))
			}
			
			varAlleleMap <- merge(knownAlleleTable, dataAlleleTable , by="alleleSeq")
}


## how do I want to store these results?
callGenotypes <- function(table, alleleDb=NULL, method="custom", minTotalReads=50, maxPropUniqueVars=0.8, 
					minPropToCall=0.1, minDiffToVarThree=0.4,
					minPropDiffHomHetThreshold=0.3, mapAlleles=FALSE) {
	
	#table$genotype
	table$status <- "notCalled"
	enoughReads <- table$numbSeqs >= minTotalReads
	table$status[!enoughReads] <- "tooFewReads"
	
	# difference between sum of vars 1 + 2 and var 3, as proportion of total
	# > 0.5 is good. <0.3 not good. 0.4 as cut-off for now? 
	table$diffToVarThree <- with(table, ((varFreq.1+varFreq.2)-varFreq.3)/numbSeqs)

	distinctVars <- 	with(table, diffToVarThree  >= minDiffToVarThree)
	table$status[enoughReads & !distinctVars] <- "complexVars"

	# difference between var 1 and var2 as proportion of total
	# homozygote: >0.3, often > 0.5
	# heterozygote: <0.25, often < 0.1
	table$propDiffHomHet <- with(table, ((varFreq.1-varFreq.2)/numbSeqs))
	eligible <- (enoughReads & distinctVars) 
	table$status[eligible] <- ifelse((table$propDiffHomHet[eligible]  >= minPropDiffHomHetThreshold), "HOMOZYGOTE","HETEROZYGOTE")

	#table$homozygote <- 	with(table, ((varFreq.1-varFreq.2)/numbSeqs) >= minPropDiffHomHetThreshold)

	return(table)
}

vectorOrRepeat <- function(paramValue, requiredLength) {
	# Used by callGenotypes.mgltResult() to equalise lengths of parameter vectors when one has length > 1
	if(length(paramValue) == requiredLength) {
		return(paramValue)
	} else {
		if(length(paramValue) == 1) {
			return(rep(paramValue, requiredLength))
		} else {
			stop(paste("Parameter has length",length(paramValue) ,"but does not match required length", requiredLength))	
		}
	}

}

## call genotypes on an "mlgtResult" object. Can select specific markers/samples to return. 
callGenotypes.mlgtResult <- function(resultObject, alleleDb=NULL, method="custom", minTotalReads=50, maxPropUniqueVars=0.8, 
					minPropToCall=0.1, minDiffToVarThree=0.4,
					minPropDiffHomHetThreshold=0.3, markerList=names(resultObject@markers),
					sampleList=resultObject@samples, mapAlleles=FALSE	) {

	## FM requested marker specific parameters.
	## test for vectors in any of the calling parameters. 	
	## if all are single, apply genotyping to one large table of all results,	
	## if any are vectors, need to apply genotyping to each marker separately. 
	## Should mark if different thresholds used. 
	## need to test that the threshold vector matches the markerlist length.

	runTable <- data.frame()
	genotypeTable <- data.frame()
	callResults <- list()

#	if(length(minTotalReads) > 1 | length(minDiffToVarThree) > 1 | length(minPropDiffHomHetThreshold) > 1 )  {
		# find parameters as vectors and test the lengths. Set those which are only 1 to be length of markerList.
		minTotalReads <- vectorOrRepeat(minTotalReads, length(markerList))
		minDiffToVarThree <- vectorOrRepeat(minDiffToVarThree, length(markerList))
		minPropDiffHomHetThreshold<- vectorOrRepeat(minPropDiffHomHetThreshold, length(markerList))		

		# multiple parameter values set. Genotype each marker separately
		for(i in 1:length(markerList))  {	
			thisMarker <- markerList[i]	
			subTable <- resultObject@markerSampleList[[thisMarker]]
			subTable <- subTable[match(sampleList,subTable$sample),]	
			if(nrow(subTable) < 1)  {
				# do nothing with this marker
				warning(paste("No data for:",thisMarker))
			} else {
				genotypeTable <- callGenotypes(subTable , alleleDb=alleleDb, method=method,minTotalReads=minTotalReads[i], 
					minDiffToVarThree=minDiffToVarThree[i],
					minPropDiffHomHetThreshold=minPropDiffHomHetThreshold[i], mapAlleles=mapAlleles)
				#genotypeTable <- rbind(genotypeTable , subGenoTable)

				if(mapAlleles) {
					if(is.null(alleleDb)) {
						warning("No alleleDb specified\n")
					}	else {
						if(is.null(alleleDb[[thisMarker]])) {
							warning(paste("No known alleles for",thisMarker))
						} else  {
							if(is.null(resultObject@alleleDb[[thisMarker]])) {
								warning(paste("No variants for",thisMarker))
							} else  {
								#varAlleleMap <- makeVarAlleleMap(alleleDb, resultObject@alleleDb, alleleMarkers=markerList, varMarkers=markerList)
								varAlleleMap <- makeVarAlleleMap(allele.variantMap=alleleDb[[thisMarker]], variant.variantMap=resultObject@alleleDb[[thisMarker]])
								genotypeTable$allele.1 <- varAlleleMap$knownAlleles[match(genotypeTable$varName.1, varAlleleMap$varNames)]
								genotypeTable$allele.2 <- varAlleleMap$knownAlleles[match(genotypeTable$varName.2, varAlleleMap$varNames)]
							}
						}
					}
				}
				callResults[[thisMarker]] <- new("genotypeCall", 
						projectName=resultObject@projectName, 
						runName=resultObject@runName,
						marker=thisMarker ,
						genotypeTable=genotypeTable,
						callMethod=method,
						callParameters=list("minTotalReads"=as.integer(minTotalReads[i]),
							"minDiffToVarThree"=minDiffToVarThree[i],	"minPropDiffHomHetThreshold"=minPropDiffHomHetThreshold[i]),
						mappedToAlleles=mapAlleles,
						alleleDbName="NeedToSetThis" )
			}


		}
	
	#return(genotypeTable)
	return(callResults)
	#new("genotypeCall", 


}


# default for 'file' could be derived from project, run, marker attributes of genotypeCall.
writeGenotypeCallsToFile.genotypeCall <- function(genotypeCall, file=paste("genoCalls",genotypeCall@projectName,genotypeCall@runName,genotypeCall@marker,"tab",sep="."),
							writeParams=FALSE, appendValue=FALSE)  {
	if(writeParams)  {
		cat("# Genotype calls generated by R package 'mlgt' (X.X.X)", date(),"\n",file=file, append=appendValue)
		appendValue=TRUE
		cat("# Project:", genotypeCall@projectName,"\n", file=file, append=appendValue)
		cat("# Run", genotypeCall@runName,"\n", file=file, append=appendValue)
		cat("# Marker", genotypeCall@marker,"\n", file=file, append=appendValue)
		cat("# Call method:",  genotypeCall@callMethod,"\n", file=file, append=appendValue)		
		cat("# Call parameters:-\n", file=file, append=appendValue)
		for(i in 1:length(genotypeCall@callParameters)) {
			thisParam <- unlist(genotypeCall@callParameters[i])
			cat("#\t",names(thisParam),"=", thisParam,"\n", file=file, append=appendValue)
		}
		cat("# MappedToAlleles: ", genotypeCall@mappedToAlleles,"\n", file=file, append=appendValue)
		cat("# AlleleDb: ", genotypeCall@alleleDbName,"\n", file=file, append=appendValue)
	}	
	write.table(genotypeCall@genotypeTable, file=file, append=appendValue, row.names=F, quote=F, sep="\t")
	cat("Results written to",file,"\n")
	#return(invisible(1))
}


# Need function to write genotype calls as a single file.

writeGenotypeCallsToFile.list <- function(callList, file, singleFile=FALSE,writeParams=FALSE)  {
	if(singleFile)  {
		masterTable <- data.frame()
		for(i in 1:length(callList)) {
			masterTable <- rbind(masterTable , callList[[i]]@genotypeTable)
		}
		write.table(masterTable , file=file, row.names=F, quote=F, sep="\t")
		cat("Results written to",file,"\n")
	} else {
		#invisible(lapply(callList, FUN=function(x) writeGenotypeCallsToFile.genotypeCall(x,writeParams=writeParams)))
		cat(length(lapply(callList, FUN=function(x) writeGenotypeCallsToFile.genotypeCall(x,writeParams=writeParams))), "files written\n")
	}	

}

# generic method. Need to give defaults if defaults to be set in specific forms.
writeGenotypeCallsToFile <- function(callList, genotypeCall, file="", singleFile=FALSE, writeParams=FALSE, appendValue=FALSE) attributes(callList)
setGeneric("writeGenotypeCallsToFile")

setMethod("writeGenotypeCallsToFile", signature(callList="list", genotypeCall="missing",file="ANY", singleFile="ANY", 
									writeParams="ANY", appendValue="ANY"), 
				definition=writeGenotypeCallsToFile.list)

setMethod("writeGenotypeCallsToFile", signature(callList="missing", genotypeCall="genotypeCall",file="ANY", singleFile="ANY", 
									writeParams="ANY", appendValue="ANY"), 
			definition=writeGenotypeCallsToFile.genotypeCall )



########################## plotting

inspectBlastResults <- function(blastTable, subject)  {

	#topHits <- getTopBlastHits(resultFile)
	#topHits$strand <- ifelse(topHits$s.end > topHits$s.start, 1,2)
	hitCount <- length(which(blastTable$subject == subject))
	if(hitCount > 0) { 
		#subject <- "DPA1_E2"
		breakValue <- max(10, 10^(floor(log10(hitCount))))	# favour a large number of breaks. At least 10.
		par(mfrow=c(1,3))
		hist(blastTable$ali.length[blastTable$subject == subject], breaks=breakValue, xlab="Alignment Length", main=subject)
		hist(blastTable$bit.score[blastTable$subject == subject], breaks=breakValue, xlab="Bit Score", main=subject)
		hist(blastTable$percent.id[blastTable$subject == subject], breaks=breakValue, xlab="% identity", main=subject)
	}	else  {
		warning(paste("No data for ", subject))
	}
}

printBlastResultGraphs <- function(designObject, markerList=designObject@markers, fileName="blastResultGraphs.pdf") {
	topHits <- getTopBlastHits(designObject@markerBlastResults)
	pdf(fileName,height=4)
	for(thisMarker in names(markerList))  {
		inspectBlastResults(topHits, thisMarker )
	}
	dev.off()
}



########## Multiple methods for plotGenotypeEvidence
# 
plotGenotypeEvidence.genotypeCall <- function(genotypeCall)  {
	genotypeTable <- genotypeCall@genotypeTable
	thisMarker	<- genotypeCall@marker
	minTotalReads <- genotypeCall@callParameters['minTotalReads']
	minDiffToVarThree <- genotypeCall@callParameters['minDiffToVarThree']
	minPropDiffHomHetThreshold <- genotypeCall@callParameters['minPropDiffHomHetThreshold']

	if(sum(genotypeTable$numbSeqs) < 1)  {
		warning(paste("No seqs to plot for",thisMarker), call.=F)
		return()
	}
	
	statusList <- as.factor(genotypeTable$status)
	pchList <- statusList
	levels(pchList) <- (1:nlevels(pchList ))
	#levels(pchList) <- 20+(1:nlevels(pchList ))


	par(mfrow=c(2,3))
	hist( genotypeTable$numbSeqs, breaks=20, main=thisMarker, xlab="numbSeqs"); abline(v=minTotalReads , lty=2)
	hist( genotypeTable$diffToVarThree, breaks=20, main=thisMarker, xlab="diffToVarThree", xlim=c(0,1)); abline(v=minDiffToVarThree , lty=2)
	hist(genotypeTable$propDiffHomHet, breaks=20, main=thisMarker, xlab="propDiffHomHet", xlim=c(0,1)) ; abline(v=minPropDiffHomHetThreshold , lty=2)

	plot(genotypeTable$diffToVarThree,genotypeTable$propDiffHomHet, main=thisMarker, xlab="diffToVarThree", ylab="propDiffHomHet",xlim=c(0,1), ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minPropDiffHomHetThreshold , lty=2); abline(v=minDiffToVarThree , lty=2)
	legend("topleft", levels(as.factor(genotypeTable$status)), pch=as.numeric(levels(pchList)))
	plot(genotypeTable$numbSeqs,genotypeTable$diffToVarThree, main=thisMarker, xlab="numbSeqs", ylab="diffToVarThree", ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minDiffToVarThree , lty=2); abline(v=minTotalReads , lty=2)
	plot(genotypeTable$numbSeqs,genotypeTable$propDiffHomHet, main=thisMarker, xlab="numbSeqs", ylab="propDiffHomHet", ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minPropDiffHomHetThreshold , lty=2); abline(v=minTotalReads , lty=2)


}

plotGenotypeEvidence.genotypeCall.file <- function(genotypeCall, file)  {
	if(length(grep(".pdf$", file) ) < 1) {
		file <- paste(file,"pdf", sep=".")
	}
	pdf(file)
	plotGenotypeEvidence.genotypeCall(genotypeCall)
	dev.off()
	cat("Results output to", file, "\n")

}


# list must have a file specified for output
plotGenotypeEvidence.list <- function(callList, file) {
	if(length(grep(".pdf$", file) ) < 1) {
		file <- paste(file,"pdf", sep=".")
	}
	pdf(file)
	for(thisCall in callList) {
		#cat(thisCall@marker)
		plotGenotypeEvidence.genotypeCall(thisCall)
	}
	dev.off()
	cat("Results output to", file, "\n")	
	
}


plotGenotypeEvidence <- function(callList, genotypeCall, file) attributes(genotypeCall)
setGeneric("plotGenotypeEvidence")
setMethod("plotGenotypeEvidence", signature(genotypeCall="missing", callList="list", file="character"), definition=plotGenotypeEvidence.list)
setMethod("plotGenotypeEvidence", signature(genotypeCall="genotypeCall", callList="missing", file="character"), definition=plotGenotypeEvidence.genotypeCall.file)
setMethod("plotGenotypeEvidence", signature(genotypeCall="genotypeCall", callList="missing", file="missing"), definition=plotGenotypeEvidence.genotypeCall)








