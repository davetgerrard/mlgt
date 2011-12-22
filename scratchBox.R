













#############################################################


intersect.cleanRun.Design <- prepareMlgtRun(projectName="testProject", runName="cleanRun", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile, overwrite="yes" )


intersect.cleanRun.Result <- mlgt(intersect.cleanRun.Design)


test.genotypes <- callGenotypes.mlgtResult(intersect.cleanRun.Result,  mapAlleles=FALSE)


# testing marker specific parameters. FAILED
test.genotypes <- callGenotypes.mlgtResult(intersect.cleanRun.Result,  mapAlleles=FALSE, minTotalReads=seq(20,200, 10))

writeGenotypeCallsToFile(test.genotypes)								
writeGenotypeCallsToFile(test.genotypes, singleFile=F, file="genotypeTable.tab")
writeGenotypeCallsToFile(genotypeCall=test.genotypes[[1]])



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