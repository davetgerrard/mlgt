## combine mlgtResult objects.
## 1. to combine and summarise results across runs.
## 2. to re-combine results after running parallel processing.
## Lots of the merging should be in common.
## Need checks that samples/markers are not duplicated (option to rename?)
## Needs to be fast enough for point 2 above to be worthwhile.


combineMlgtResults <- function(resultList,projectName=resultList[[1]]@projectName)  {


	# check that sample/marker combinations do not overlap


	# set first result as master
	master <- resultList[[1]]

	# cycle through remaining results adding each to master
	# runSummaryTable
		# if above tests passed, should be able to rbind runSummaryTable
	# markerSampleList 
		# rbind within markers. or add if markers not in master@markerSampleList
	# varCountTables
	# alleleDb.

	

	# combine variants/alleles.
		

	return(master)
}



project.mlgtResult <- combineMlgtResults(as.list(my.mlgt.Result.1, my.mlgt.Result.2))