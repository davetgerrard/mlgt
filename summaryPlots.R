	topHits <- getTopBlastHits("blastOut.markers.tab")
	topHits$strand <- ifelse(topHits$s.end > topHits$s.start, 1,2)


	hist(topHits$p_value[topHits$subject == thisMarker], breaks=1000)


	hist(log(1/topHits$p_value[topHits$subject == thisMarker]), breaks=100) # inverse log. Good but hard to decipher

	# using plain e-value is very similar
	hist(topHits$e_value[topHits$subject == thisMarker], breaks=1000)
	hist(topHits$e_value[topHits$subject == "DQA1_E2"], breaks=1000)
	hist(topHits$e_value[topHits$subject == "DRB1_E2"], breaks=1000)

	# whole dataset
	hist(topHits$p_value, breaks=1000)
	hist(topHits$e_value, breaks=1000)
	hist(topHits$percentId, breaks=100)





############ genotyping calls.

# need to distinguish sampleMarker table from same with genotype calls. Currently calss callGenotypes on the table.
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


