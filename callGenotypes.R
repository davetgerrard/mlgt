# Sensible to call. Read/variant > 2. 
# Minreads per allele/sample 10 
# (have called some with just 3 if other info 
# e.g. it’s a common allele and only one with >1 reads)?
# Influenced by whether allele is present in other samples.
# Also if appears heterozygous, seems better to call 2nd allele if distinct from first 
# (e.g. >2% difference). 





## decide if genotype call can be made.
# total count, 100s is good, 10s is weak
# % of novel variants. 30% is very good, >80% is very poor.

## Call as homozygous or heterozygous or NA/NULL/?
# Ratio of first to second?
#Ratio of second to third, or first to third?

#var >20% for call.
#var2 >10% for heterozygotes
#var 3 < 15% ?  the highest for DPA1 was 12% (of 165 seqs)

# FUTURE: differences between alleles may be useful e.g. more confidence in two known but very different alleles than when two known, but very similar (1bp) alleles.


## Match up to known alleles if available. 

## Nice format




########################

### NB WILL BE FASTEST TO RUN GENOTYPE CALLS ON A MASTER TABLE OF ALL RESULTS FROM A RUN.
### E.G. cbind() the markerSampleList tables first.
### could integrate the cbind with a marker/sample list so that table is generated on the fly.
## e.g. markerList=names(table@markers)  as default for marker on object "mlgtResult" similar for samples.
## TODO:  Current merged varAlleleMap creates one line per sequence. If that sequence appears under more than one marker, only one will be retained. Or is this a match() problem?  

# mixed. Not great
#intersect.756948.Result@markerSampleList[["C_E3"]]
# very clean, not 100% perfect.
#intersect.756948.Result@markerSampleList[["DQA1_E2"]]

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

#methods("writeGenotypeCallsToFile")



######### examples of writeGenotypeCallsToFile methods. Use of generic below. 


writeGenotypeCallsToFile.genotypeCall(test.genotypes[[2]],file="test.genotypes.tab" )
writeGenotypeCallsToFile.genotypeCall(test.genotypes[[2]],file="test.genotypes.tab" ,writeParams=TRUE)
writeGenotypeCallsToFile.genotypeCall(test.genotypes[[2]],writeParams=TRUE)

# Write a whole list of genotype calls to separate files.
invisible(lapply(test.genotypes, FUN=function(x) writeGenotypeCallsToFile.genotypeCall(x,writeParams=TRUE)))
# is equivalent to 
writeGenotypeCallsToFile.list(test.genotypes)
# and with generic set also
writeGenotypeCallsToFile(test.genotypes)

writeGenotypeCallsToFile.list(test.genotypes, singleFile=T, file="genotypeTable.tab")

## once methods are set:-

writeGenotypeCallsToFile(test.genotypes, singleFile=T, file="genotypeTable.tab")	# list input, one file of concatenated table
writeGenotypeCallsToFile(test.genotypes)								# same
writeGenotypeCallsToFile(test.genotypes, singleFile=F, file="genotypeTable.tab")	# list input, one file per marker, ignores file name given
writeGenotypeCallsToFile(genotypeCall=test.genotypes[[1]],  file="genotypeTable.tab",writeParams=TRUE)		# genotypeCall input. 
writeGenotypeCallsToFile(genotypeCall=test.genotypes[[1]])

#combineTablesInList <- function(list, tableName)  {
#	resultTable <- data.frame()
#	for(x in 1:length(list)) {
#		resultTable <- rbind(resultTable[[i]][tableName])
#	lapply(list, function(x) rbind(resultTable, x@[tableName]))
#	}
#	return(resultTable)
#}#



stopifnot(FALSE)


######################### OUTPUT RESULTS TO DATE
setwd("C:/Users/dave/HalfStarted/mlgt/testProject/mapImgtAlleles/")

test.genotypes <- callGenotypes.mlgtResult(intersect.756948.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)
write.table(test.genotypes, file="intersect.756948.genotypesAlleles.tab",  quote=F, sep="\t", row.names=F)


test.genotypes <- callGenotypes.mlgtResult(intersect.760591.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)
write.table(test.genotypes, file="intersect.760591.genotypesAlleles.tab",  quote=F, sep="\t", row.names=F)

test.genotypes <- callGenotypes.mlgtResult(intersect.763932m.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)
write.table(test.genotypes, file="intersect.763932m.genotypesAlleles.tab",  quote=F, sep="\t", row.names=F)

test.genotypes <- callGenotypes.mlgtResult(intersect.763932n.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)
write.table(test.genotypes, file="intersect.763932n.genotypesAlleles.tab",  quote=F, sep="\t", row.names=F)

test.genotypes <- callGenotypes.mlgtResult(intersectLiverpool.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)
write.table(test.genotypes, file="intersectLiverpool.genotypesAlleles.tab",  quote=F, sep="\t", row.names=F)

test.genotypes <- callGenotypes.mlgtResult(intersect.catResult,  mapAlleles=TRUE, alleleDb=knownAlleleDb)
write.table(test.genotypes, file="intersect.cat.genotypesAlleles.tab",  quote=F, sep="\t", row.names=F)






########################## DEVELOPMENT


### this is broken. 
#setGeneric("callGenotypes")
#setMethod("callGenotypes",signature(resultObject="mlgtResult",alleleDb="list",minTotalReads="integer",maxPropUniqueVars="numeric",
#					minPropToCall="numeric",minDiffToVarThree="numeric",minPropDiffHomHetThreshold="numeric"), 
#		definition=callGenotypes.mlgtResult)
#
#method.skeleton("callGenotypes", "mlgtResult")

testTable <- intersect.756948.Result@markerSampleList[["C_E3"]]
testTable <- intersect.756948.Result@markerSampleList[["DQA1_E2"]]

(test.genotypes <- callGenotypes(testTable))

test.genotypes <- callGenotypes.mlgtResult(intersect.756948.Result)

test.genotypes <- callGenotypes.mlgtResult(intersect.756948.Result, markerList=c("DQA1_E2", "A_E2"), sampleList=c("MID-1", "MID-31"))
test.genotypes <- callGenotypes.mlgtResult(intersect.756948.Result, sampleList=c("MID-1"))
	


test.genotypes <- callGenotypes.mlgtResult(intersect.756948.Result, markerList=c("DQA1_E2", "A_E2"), sampleList=c("MID-1", "MID-31"), mapAlleles=TRUE)

(test.genotypes <- callGenotypes(testTable,mapAlleles=TRUE))



#mapTable <- mapVarsToKnownAlleles(knownAlleleDb, intersect.756948.Result@alleleDb)
#alleleDb <- knownAlleleDb
#varDb <- intersect.756948.Result@alleleDb
#rm(varDb, alleleDb, mapTable)


(test.genotypes <- callGenotypes.mlgtResult(intersect.756948.Result, markerList=c("DQA1_E2", "A_E2"), sampleList=c("MID-1", "MID-31"), mapAlleles=TRUE, alleleDb=knownAlleleDb))

test.genotypes <- callGenotypes.mlgtResult(intersect.756948.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)

write.table(test.genotypes, file="outTest.genotypesAlleles.tab",  quote=F, sep="\t", row.names=F)






