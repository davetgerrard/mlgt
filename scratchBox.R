





# Looking at alleleDb objects

str(intersect.cleanRun.Result@alleleDb[[7]], deparse.level=1)







#################################### New Clean Run

# prepare mlgt by setting up blast dbs and running blast searches.
# Creates object to store run settings
intersect.cleanRun.Design <- prepareMlgtRun(projectName="testProject", runName="cleanRun", 
				samples=sampleList, markers=intersectMarkerList ,
				fTags=fTagList, rTags=rTagList, inputFastaFile=inputDataFile, overwrite="yes" )


# inspect BLAST results for a specific marker
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

# Write genotype call graphics to file
writeGenotypeCallsToFile(test.genotypes)								
writeGenotypeCallsToFile(test.genotypes, singleFile=F, file="genotypeTable.tab")
writeGenotypeCallsToFile(genotypeCall=test.genotypes[[1]])

# Create a structured list of known alleles (e.g. from IMGT/HLA)

# Call genotypes and map alleles to list of known alleles
test.genotypes <- callGenotypes.mlgtResult(intersect.cleanRun.Result,  mapAlleles=TRUE, alleleDb=knownAlleleDb)

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