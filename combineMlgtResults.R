## combine mlgtResult objects.

		# @runSummaryTable
			# if above tests passed, should be able to rbind runSummaryTable
		# @markerSampleList. List (one element per marker) of data.frame, one row per sample
			# rbind within markers. or add if markers not in master@markerSampleList
		# @varCountTables. List (one element per marker) of data.frame, one column per sample, one row per variant.
		# @alleleDb. List (one element per marker) of class variantMap
			#variantMap @reference, @ variantSource, @variantMap, @inputVariantCount, @uniqueSubVariantCount
		# @projectName
		# @runName
		# @markers. List of class SeqFastadna
		# @samples
		# @fTags. List of class SeqFastadna
		# @rTags. List of class SeqFastadna
		# @inputFastaFile    
		# @markerBlastResults
		# @fTagBlastResults  
		# @rTagBlastResults  
		# combine variants/alleles.


updateSampleNames <- function(mlgtResultObject,newSampleNames)  {
	
}

updateMarkerNames <- function(mlgtResultObject,newMarkerNames)  {

}


mergeMlgtResults.complex <- function(result1, result2) {
	master <- result1

	# need to check if samples shared between results.
	# If not, then perform complex join. Combining results where approriate.
	# If they are, then test if marker/sample pairs shared. If yes, then stop(), otherwise perform complex join.

	master@markers <- c(master@markers,result2@markers)
	master@markers <- master@markers[!duplicated(master@markers)]

	newRunSummaryTable <- data.frame()
	for(thisMarker in names(master@markers)) {
		markerSampleOverlap <- intersect(master@markerSampleList[[thisMarker]]$sample,result2@markerSampleList[[thisMarker]]$sample)
		if(length(markerSampleOverlap) > 0)  {
			# STOP!
			stop(paste("Cannot join results sharing marker results for same samples\n",thisMarker,markerSampleOverlap,"\n"))
		}

		# need to combine alleleDBs and variantMaps first so that new consistent allele names can propagate to the markerSampleList
		#master@alleleDb <- as.list(merge(master@alleleDb, result2@alleleDb))
		#master@alleleDb <- mergeAlleleDbs(master@alleleDb)
		masterAlleleTable <- data.frame(seq=names(unlist(master@alleleDb[[thisMarker]]@variantMap)),
								masterName=as.character(unlist(master@alleleDb[[thisMarker]]@variantMap)), stringsAsFactors=F )
		alleleTable2 <- data.frame(seq=names(unlist(result2@alleleDb[[thisMarker]]@variantMap)),
								alleleName=as.character(unlist(result2@alleleDb[[thisMarker]]@variantMap)), stringsAsFactors=F )
		masterMatchTable <- merge(masterAlleleTable, alleleTable2, by="seq", all=T)
		# some alleles will be new and these need to be added to master table and given new allele names
		maxAllele <- max(na.omit(as.numeric(sub(paste(thisMarker,".",sep=""),"",masterMatchTable$masterName))))
		newAlleleCount <-sum(is.na(masterMatchTable$masterName))
		masterMatchTable$masterName[is.na(masterMatchTable$masterName)] <- paste(thisMarker, (1:newAlleleCount + maxAllele), sep=".")
	
		#create new variantMap from masterMatchTable
		newVarMap <- as.list(masterMatchTable$masterName)
		names(newVarMap) <- masterMatchTable$seq
		master@alleleDb[[thisMarker]]@variantMap <- newVarMap
	
		# update allele Names in result2 table
		nameCols <- grep("varName", names(result2@markerSampleList[[thisMarker]]))
		for(n in nameCols) {
			result2@markerSampleList[[thisMarker]][,n] <- masterMatchTable$masterName[match(result2@markerSampleList[[thisMarker]][,n] , masterMatchTable$alleleName)]
		}
		#master@markerSampleList <- as.list(merge(master@markerSampleList, result2@markerSampleList))
		master@markerSampleList[[thisMarker]] <- rbind(master@markerSampleList[[thisMarker]],result2@markerSampleList[[thisMarker]])
	
		# Combine result across varCountTables.
		#master@varCountTables <- as.list(merge(master@varCountTables, result2@varCountTables))
		# varCountTables[[thisMarker]][seqTable$alignedVar,thisSample] <- seqTable$count
		frame2 <- result2@varCountTables[[thisMarker]]
		for(thisCol in names(frame2)) {
			index <- which(!is.na(frame2[,thisCol] ))
			counts <- frame2[index,thisCol]
			names <- row.names(frame2)[index]
			master@varCountTables[[thisMarker]][names,thisCol] <- counts
		}

		index1 <- match(thisMarker, master@runSummaryTable$marker)
		index2 <- match(thisMarker, result2@runSummaryTable$marker)

		runSummaryRow <- data.frame(marker=thisMarker, assignedSeqs=sum(master@markerSampleList[[thisMarker]]$numbSeqs), assignedVariants=sum(master@markerSampleList[[thisMarker]]$numbVars), 
					minVariantLength=min(master@runSummaryTable$minVariantLength[index1], result2@runSummaryTable$minVariantLength[index2] ), 
					maxVariantLength=max(master@runSummaryTable$maxVariantLength[index1], result2@runSummaryTable$maxVariantLength[index2] ),
					minAlleleLength=min(master@runSummaryTable$minAlleleLength[index1], result2@runSummaryTable$minAlleleLength[index2] ), 
					maxAlleleLength=max(master@runSummaryTable$maxAlleleLength[index1], result2@runSummaryTable$maxAlleleLength[index2] ) )
		newRunSummaryTable <- rbind(newRunSummaryTable,runSummaryRow )
	}	
	master@samples <- union(master@samples,result2@samples)
	
	#master@runSummaryTable <- rbind(master@runSummaryTable, result2@runSummaryTable)  # these need to be rbinded then rows with same marker summed
	master@runSummaryTable <- newRunSummaryTable 
	# keep the following unless different between results.
	if(!identical(master@fTags,result2@fTags)) {	master@fTags <- list() }	# @fTags. List of class SeqFastadna
	if(!identical(master@rTags,result2@rTags)) {	master@rTags<- list() }		# @rTags. List of class SeqFastadna
	if(!identical(master@inputFastaFile,result2@inputFastaFile)) {	master@inputFastaFile <- '' }		# @inputFastaFile    
	if(!identical(master@markerBlastResults,result2@markerBlastResults)) {	master@markerBlastResults <- '' }		# @markerBlastResults
	if(!identical(master@fTagBlastResults,result2@fTagBlastResults)) {	master@fTagBlastResults <- '' }		# @fTagBlastResults  
	if(!identical(master@rTagBlastResults,result2@rTagBlastResults)) {	master@rTagBlastResults <- '' }		# @rTagBlastResults  
	if(!identical(master@runName,result2@runName)) {	master@runName<- 'CombinedResults' }			# @runName
	# @projectName is taken as the first result. May need to change this.
	return(master)
}

mergeMlgtResults.simple <- function(result1, result2) {
	master <- result1
	master@markers <- c(master@markers,result2@markers)
	master@markers <- master@markers[!duplicated(master@markers)]
	master@markerSampleList <- c(master@markerSampleList, result2@markerSampleList)
	master@markerSampleList <- master@markerSampleList[!duplicated(master@markerSampleList)]
	master@alleleDb <- c(master@alleleDb, result2@alleleDb)
	master@alleleDb <- master@alleleDb[!duplicated(master@alleleDb)]
	master@varCountTables <- c(master@varCountTables, result2@varCountTables)
	master@varCountTables <- master@varCountTables[!duplicated(master@varCountTables)]
	master@samples <- union(master@samples,result2@samples)
	master@runSummaryTable <- rbind(master@runSummaryTable, result2@runSummaryTable)
	# keep the following unless different between results.
	if(!identical(master@fTags,result2@fTags)) {	master@fTags <- list() }	# @fTags. List of class SeqFastadna
	if(!identical(master@rTags,result2@rTags)) {	master@rTags<- list() }		# @rTags. List of class SeqFastadna
	if(!identical(master@inputFastaFile,result2@inputFastaFile)) {	master@inputFastaFile <- '' }		# @inputFastaFile    
	if(!identical(master@markerBlastResults,result2@markerBlastResults)) {	master@markerBlastResults <- '' }		# @markerBlastResults
	if(!identical(master@fTagBlastResults,result2@fTagBlastResults)) {	master@fTagBlastResults <- '' }		# @fTagBlastResults  
	if(!identical(master@rTagBlastResults,result2@rTagBlastResults)) {	master@rTagBlastResults <- '' }		# @rTagBlastResults  
	if(!identical(master@runName,result2@runName)) {	master@runName<- 'CombinedResults' }			# @runName
	# @projectName is taken as the first result. May need to change this.
	return(master)
}

## I imagine most of the time, combining results can be done after mlgt is run by simple concatenation of genotype tables. 
## This might be the case if certain markers are rerun. 
## However, there are instances where combination within mlgt is desirable. e.g. after a parallel run. 
## Re-runs of certain marker samples could also be accomodated, but care would be have to be taken to update the correct results.

## 1. to combine and summarise results across runs.
## 2. to re-combine results after running parallel processing.
## Lots of the merging should be in common.
## Need checks that samples/markers are not duplicated (option to rename?)
## Needs to be fast enough for point 2 above to be worthwhile.
## Even if markers don't overlap, need to check that samples and mids are the same or different.
## Or do I? MIDs should be allowed to differ between runs. 

combineMlgtResults <- function(resultList,projectName=resultList[[1]]@projectName, runName="combinedMlgtResults")  {

	# set first result as master
	master <- resultList[[1]]

	for(i in 2:length(resultList)) {
		# cycle through remaining results adding each to master
		# check that sample/marker combinations do not overlap
		#Any overlap in markers?
		markerOverlap <- intersect(names(master@markers),names(resultList[[i]]@markers))
			

		if(length(markerOverlap) > 0) {	# overlapping markers, test for sample/marker overlap and maybe complex join.
			# need to check if overlapping marker sequences are identical
			cat("Complex join\n")
			if(!identical(master@markers[markerOverlap],resultList[[i]]@markers[markerOverlap]))  {
				#STOP
				stop("Cannot combine: markers with same name have different sequences.\n")
			}
			# need to check if samples shared between results. If not, then perform complex join. If they are, then test if marker/sample pairs shared.
			master <- mergeMlgtResults.complex(master, resultList[[i]])
		} else {					# no overlap in markers, 'simple' join
			cat("Simple join\n")
			master <- mergeMlgtResults.simple(master, resultList[[i]])

		}
		
	}
	return(master)
}


### Testing

resultList[[1]] <- my.mlgt.Result
resultList[[2]] <- my.mlgt.Result
resultList


project.mlgtResult <- combineMlgtResults(resultList)




project.mlgtResult <- combineMlgtResults(as.list(my.mlgt.Result.1, my.mlgt.Result.2))