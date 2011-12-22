
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


## semi-manual version. Good for single marker and/or in GUI
topHits <- getTopBlastHits(intersect.cleanRun.Design@markerBlastResults)

inspectBlastResults(topHits, thisMarker)

pdf("blastInspection.pdf",height=4)
for(thisMarker in names(intersectMarkerList))  {
	inspectBlastResults(topHits, thisMarker )
}
dev.off()

## automatic output to pdf of same for a list of markers.
printBlastResultGraphs(intersect.cleanRun.Design)


############ genotyping calls.

# need to distinguish sampleMarker table from same with genotype calls. Currently calss callGenotypes on the whole table.
plotGenotypeEvidence <- function(sampleMarkerTable,thisMarker,
			minDiffToVarThree = 0.4, minTotalReads = 50, minPropDiffHomHetThreshold = 0.3 )  {
	test.genotypes <- callGenotypes(sampleMarkerTable)

	statusList <- as.factor(test.genotypes$status)
	pchList <- statusList
	levels(pchList) <- (1:nlevels(pchList ))
	#levels(pchList) <- 20+(1:nlevels(pchList ))

	par(mfrow=c(2,3))
	hist( test.genotypes$numbSeqs, breaks=20, main=thisMarker, xlab="numbSeqs"); abline(v=minTotalReads , lty=2)
	hist( test.genotypes$diffToVarThree, breaks=20, main=thisMarker, xlab="diffToVarThree", xlim=c(0,1)); abline(v=minDiffToVarThree , lty=2)
	hist(test.genotypes$propDiffHomHet, breaks=20, main=thisMarker, xlab="propDiffHomHet", xlim=c(0,1)) ; abline(v=minPropDiffHomHetThreshold , lty=2)

	plot(test.genotypes$diffToVarThree,test.genotypes$propDiffHomHet, main=thisMarker, xlab="diffToVarThree", ylab="propDiffHomHet",xlim=c(0,1), ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minPropDiffHomHetThreshold , lty=2); abline(v=minDiffToVarThree , lty=2)
	legend("topleft", levels(as.factor(test.genotypes$status)), pch=as.numeric(levels(pchList)))
	plot(test.genotypes$numbSeqs,test.genotypes$diffToVarThree, main=thisMarker, xlab="numbSeqs", ylab="diffToVarThree", ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minDiffToVarThree , lty=2); abline(v=minTotalReads , lty=2)
	plot(test.genotypes$numbSeqs,test.genotypes$propDiffHomHet, main=thisMarker, xlab="numbSeqs", ylab="propDiffHomHet", ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minPropDiffHomHetThreshold , lty=2); abline(v=minTotalReads , lty=2)

			
}


plotGenotypeEvidence.mlgtResult <- function(mlgtResultObject, outFile="genotypeCallPlotsByMarker.pdf")  {
	pdf(outFile, width=10, height=6)
	for(thisMarker in names(mlgtResultObject@markers))  {
		testTable <- mlgtResultObject@markerSampleList[[thisMarker]]
		plotGenotypeEvidence(testTable, thisMarker)
	}
	dev.off()
}

#plotGenotypeEvidence(testTable, thisMarker)



analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/intersect756948"
setwd( analysisDir )
load("intersect.756948.RData")
plotGenotypeEvidence.mlgtResult(intersect.756948.Result, outFile="intersect756948.genotypeCallPlotsByMarker.pdf")

analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/intersectCat"
setwd( analysisDir )
load("intersect.catResult.RData")
plotGenotypeEvidence.mlgtResult(intersect.catResult, outFile="intersectCat.genotypeCallPlotsByMarker.pdf")



###########################  development


	topHits <- getTopBlastHits("blastOut.markers.tab")
	topHits$strand <- ifelse(topHits$s.end > topHits$s.start, 1,2)

	thisMarker <- "DPA1_E2"
	hist(topHits$e.value[topHits$subject == thisMarker], breaks=1000)


	hist(log(1/topHits$e.value[topHits$subject == thisMarker]), breaks=100) # inverse log. Good but hard to decipher

	# using plain bit.score is very similar
	hist(topHits$bit.score[topHits$subject == thisMarker], breaks=1000)
	hist(topHits$bit.score[topHits$subject == "DQA1_E2"], breaks=1000)
	hist(topHits$bit.score[topHits$subject == "DRB1_E2"], breaks=1000)

	# whole dataset
	hist(topHits$e.value, breaks=1000)
	hist(topHits$bit.score, breaks=1000)
	hist(topHits$percent.id, breaks=100)


	subject <- "DPA1_E2"
	breakValue <-  10^(floor(log10(length(topHits$subject == subject))))
	par(mfrow=c(1,3))
	hist(topHits$ali.length[topHits$subject == subject], breaks=1000, xlab="Alignment Length", main=subject)
	hist(topHits$bit.score[topHits$subject == subject], breaks=1000, xlab="Bit Score", main=subject)
	hist(topHits$percent.id[topHits$subject == subject], breaks=1000, xlab="% identity", main=subject)



inspectBlastResults(intersect.cleanRun.Design@markerBlastResults, "DPA1_E2")


inspectSampleBlastResults <- function(designObject, thisSample)  {
	# need to deal with separate results from forward and reverse strands.
	#inspectBlastResults(designObject@markerBlastResults, subject=thisSample)
}


inspectMarkerBlastResults <- function(designObject, thisMarker)  {
	inspectBlastResults(designObject@markerBlastResults, subject=thisMarker)
}

inspectMarkerBlastResults(intersect.cleanRun.Design, "DPA1_E2")


pdf("blastInspection.pdf")
for(thisMarker in names(intersectMarkerList))  {
	inspectMarkerBlastResults(intersect.cleanRun.Design, thisMarker )
}
dev.off()
















minDiffToVarThree <- 0.4
minTotalReads <- 50
minPropDiffHomHetThreshold <- 0.3


analysisDir <-  "C:/Users/dave/HalfStarted/mlgt/testProject/intersect756948"
setwd( analysisDir )
load("intersect.756948.RData")

pdf("genotypeCallPlotsByMarker.pdf", width=10, height=6)
for(thisMarker in names(intersect.756948.Result@markers))  {

testTable <- intersect.756948.Result@markerSampleList[[thisMarker]]
test.genotypes <- callGenotypes(testTable)

statusList <- as.factor(test.genotypes$status)
pchList <- statusList
levels(pchList) <- (1:nlevels(pchList ))
#levels(pchList) <- 20+(1:nlevels(pchList ))

par(mfrow=c(2,3))
hist( test.genotypes$numbSeqs, breaks=20, main=thisMarker, xlab="numbSeqs"); abline(v=minTotalReads , lty=2)
hist( test.genotypes$diffToVarThree, breaks=20, main=thisMarker, xlab="diffToVarThree", xlim=c(0,1)); abline(v=minDiffToVarThree , lty=2)
hist(test.genotypes$propDiffHomHet, breaks=20, main=thisMarker, xlab="propDiffHomHet") ; abline(v=minPropDiffHomHetThreshold , lty=2)

plot(test.genotypes$diffToVarThree,test.genotypes$propDiffHomHet, main=thisMarker, xlab="diffToVarThree", ylab="propDiffHomHet",xlim=c(0,1), ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minPropDiffHomHetThreshold , lty=2); abline(v=minDiffToVarThree , lty=2)
legend("topleft", levels(as.factor(test.genotypes$status)), pch=as.numeric(levels(pchList)))
plot(test.genotypes$numbSeqs,test.genotypes$diffToVarThree, main=thisMarker, xlab="numbSeqs", ylab="diffToVarThree", ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minDiffToVarThree , lty=2); abline(v=minTotalReads , lty=2)
plot(test.genotypes$numbSeqs,test.genotypes$propDiffHomHet, main=thisMarker, xlab="numbSeqs", ylab="propDiffHomHet", ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minPropDiffHomHetThreshold , lty=2); abline(v=minTotalReads , lty=2)


}

dev.off()



thisMarker <- "DPA1_E2"

testTable <- intersect.756948.Result@markerSampleList[[thisMarker]]
testTable$diffToVarThree <- with(testTable, ((varFreq.1+varFreq.2)-varFreq.3)/numbSeqs)
testTable$propDiffHomHet <- with(testTable, ((varFreq.1-varFreq.2)/numbSeqs))



hist( testTable$numbSeqs, breaks=20, main=thisMarker, xlab="numbSeqs"); abline(v=minTotalReads , lty=2)
hist( testTable$diffToVarThree, breaks=20, main=thisMarker, xlab="diffToVarThree", xlim=c(0,1)); abline(v=minDiffToVarThree , lty=2)
plot(testTable$numbSeqs,testTable$diffToVarThree, main=thisMarker, xlab="numbSeqs", ylab="diffToVarThree", ylim=c(0,1)); abline(h=minDiffToVarThree , lty=2); abline(v=minTotalReads , lty=2)

hist(testTable$propDiffHomHet, breaks=20, main=thisMarker, xlab="propDiffHomHet") ; abline(v=minPropDiffHomHetThreshold , lty=2)
plot(testTable$numbSeqs,testTable$propDiffHomHet, main=thisMarker, xlab="numbSeqs", ylab="propDiffHomHet", ylim=c(0,1)); abline(h=minPropDiffHomHetThreshold , lty=2); abline(v=minTotalReads , lty=2)

plot(testTable$diffToVarThree,testTable$propDiffHomHet, main=thisMarker, xlab="diffToVarThree", ylab="propDiffHomHet",xlim=c(0,1), ylim=c(0,1)); abline(h=minPropDiffHomHetThreshold , lty=2); abline(v=minDiffToVarThree , lty=2)


