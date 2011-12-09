

###############################################################################

## In Linux system, can get msf alignments for all loci using:-
#	 wget -r -l1 --no-parent -A gen.msf ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/

## Notes on nuc.msf files. (what are gen.msf files?)
# there are separate DRB, DRB1 and DRB345 alignments. DRB is about as large as DRB1 and DRB345 put together. 
# DRB contains 1150 alleles from 1-9 (most are DRB1)


#tempList <- c( dir("C:/Users/dave/HLA/data/IMGT_manualDownload/", pattern="nuc.msf"), names(markerList))
#write(tempList, file="marker.imgt.msf.list.tab")  # edited manually to pair up markers with imgt align files.

# Need checks if any alleles are empty sequence.
# TODO: knownAlleleDb and mlgtResult@alleleDb need to use same form (from the former). Update code in makeVarAlleleMap()


# FUTURE: Might be able to take alignment of known alleles and design primers to amplify regions which differentiate them. 


###########################################

library(seqinr)

# not needed anymore as fixed the issue with read.alignment. 
read.msf <- function(file, forceToLower = TRUE)  {
	if (file.access(file, mode = 4) != 0) 
		stop(paste("File", file, "is not readable"))
	ali  <- .Call("read_msf_align", file, PACKAGE = "seqinr")
	ali <- lapply(ali, as.character)
	ali[[3]] <- lapply(ali[[3]], tolower)
	ali <- list(nb = as.numeric(ali[[1]]), nam = ali[[2]], 
		seq = ali[[3]], com = NA)
	class(ali) <- "alignment"	
	return(ali)
}

#TODO: get this to return a proper classed object.
createKnownAlleleList <- function(markerName, markerSeq, alignedAlleleFile)  {
	## The aligned alleles must have unique names and the markerName must be different too. TODO: test for this.
	## TODO put default input file format (MSF). Make check for fasta format (can skip first part if already fasta).
	## clean up (remove) files. This function probably doesn't need to keep any files
	## Use defined class for return object giving marker sequence used as reference. 
	#alignedAlleles <- read.msf(alignedAlleleFile)
	alignedAlleles <- read.alignment(alignedAlleleFile, format="msf")
	alignedAlleleFastaFile <- paste(markerName, "alignedAlleles.fasta", sep=".")
	write.fasta(alignedAlleles[[3]], alignedAlleles[[2]], file=alignedAlleleFastaFile)
#}	# TEMP END
	markerSeqFile <- paste(markerName, "markerSeq.fasta", sep=".")
	write.fasta(markerSeq, markerName, file=markerSeqFile )
	#muscle -profile -in1 existing_aln.afa -in2 new_seq.fa -out combined.afa
	markerToAlleleDbAlign <- paste(markerName, "allignedToAlleles.fasta", sep=".")
	muscleCommand <- paste(musclePath, "-quiet -profile -in1", alignedAlleleFastaFile, "-in2", markerSeqFile, "-out" ,markerToAlleleDbAlign )
	system(muscleCommand)

#}	# TEMP END

	### this section copied from getSubSeqsTable()  Could be recoded as function?

	# Extract portion corresponding to reference. 

	rawAlignment <- read.fasta(markerToAlleleDbAlign , as.string=T)		# do not use read.alignment() - broken
	alignedMarkerSeq <- s2c(rawAlignment[[thisMarker]])
	subStart <- min(grep("-",alignedMarkerSeq ,invert=T))
	subEnd <- max(grep("-",alignedMarkerSeq ,invert=T))
	alignedSubSeqs <- lapply(rawAlignment, FUN=function(x)	substr(x[1], subStart, subEnd))
	#subAlignFile <- paste(thisMarker, "IMGT", "sub.align.fasta",sep=".")  #"localAlign.fasta"
	#write.fasta(alignedSubSeqs , names(alignedSubSeqs ), file=subAlignFile )

	alignedSubTable <- data.frame(name =  names(alignedSubSeqs ) , subSeq.aligned= as.character(unlist(alignedSubSeqs )))
	alignedSubTable$subSeq.stripped <-  gsub("-", "",alignedSubTable$subSeq.aligned )
	## TODO: remove marker sequence from allele subalignment list. DONE!
	alignedSubTable <- subset(alignedSubTable, name != markerName)

	alleleMap <- split(as.character(alignedSubTable$name), alignedSubTable$subSeq.stripped)
	return(list(reference=as.SeqFastadna(markerSeq, markerName), alleleMap=alleleMap, inputAlleleCount = length(unlist(alleleMap)), uniqueSubAlleleCount=length(alleleMap)))
}


setwd("C:/Users/dave/HalfStarted/mlgt/testProject/mapImgtAlleles/")

markerImgtFileTable <- read.delim("C:/Users/dave/HalfStarted/mlgt/marker.imgt.msf.list.tab", header=T)
imgtFileDir <- "C:/Users/dave/HLA/data/IMGT_manualDownload/"


## USING ORIGINAL MARKER LIST
knownAlleleDb <- list()
# takes about 2 minutes with 17 markers. 
for(thisMarker in names(markerList)) {
	baseFile <- markerImgtFileTable$imgtAlignFile[markerImgtFileTable$marker==thisMarker]
	imgtAlignFile <- paste(imgtFileDir,baseFile , sep="") 
	knownAlleleDb[[thisMarker]] <- createKnownAlleleList(thisMarker,markerList[[thisMarker]][1], imgtAlignFile)
}



str(knownAlleleDb)







## How to compare/assign sequence alleles with previously known alleles. 

## Positive assignment does depend on the two alleleDb constructions arriving at the exact same subsequence. Better than BLAST though. 

thisMarker <- "DQA1_E2"
thisMarker <- "A_E2"		#A_E2 allele alignment is pretty bad. 
thisMarker <- "DRB1_E2"

# Probably need to limit to exonic only. for most of the alleles. 
match(names(catResult@alleleDb[[thisMarker]]), names(knownAlleleDb[[thisMarker]]$alleleMap ))

match(names(catResult@alleleDb[[thisMarker]]), names(knownAlleleDb[[thisMarker]]$alleleMap ))

for(thisMarker in names(markerList)) {
	cat(paste(thisMarker, length(na.omit(match(names(catResult@alleleDb[[thisMarker]]), names(knownAlleleDb[[thisMarker]]$alleleMap )))), "\n"))
}


# use sequences to match sequenced alleles with known alleles
# works reasonably well with DRB1_E2
thisMarker <- "DRB1_E2"
thisMarker <- "DPA1_E2"

alleleMatchTable <- cbind(as.character(catResult@alleleDb[[thisMarker]]), as.character(knownAlleleDb[[thisMarker]]$alleleMap[ match(names(catResult@alleleDb[[thisMarker]]), names(knownAlleleDb[[thisMarker]]$alleleMap ))]))
catResult@markerSampleList[[thisMarker ]]



################################################
# TEST AGAINST RUN WITH INTERSECT MARKERS.



## USING INTERSECT MARKER LIST
knownAlleleDb <- list()
# takes about 2 minutes with 17 markers. 
for(thisMarker in names(intersectMarkerList)) {
	baseFile <- markerImgtFileTable$imgtAlignFile[markerImgtFileTable$marker==thisMarker]
	imgtAlignFile <- paste(imgtFileDir,baseFile , sep="") 
	knownAlleleDb[[thisMarker]] <- createKnownAlleleList(thisMarker,intersectMarkerList[[thisMarker]][1], imgtAlignFile)
}




# intersectLiverpool.Result

for(thisMarker in names(intersectMarkerList)) {
	cat(paste(thisMarker, length(na.omit(match(names(intersectLiverpool.Result@alleleDb[[thisMarker]]), names(knownAlleleDb[[thisMarker]]$alleleMap )))), "\n"))
}
# YES! Recognised alleles from all requested markers!



for(thisMarker in names(intersectMarkerList)) {
	cat(paste(thisMarker, length(na.omit(match(names( intersect.756948.Result@alleleDb[[thisMarker]]), names(knownAlleleDb[[thisMarker]]$alleleMap )))), "\n"))
}


for(thisMarker in names(intersectMarkerList)) {
	cat(paste(thisMarker, length(na.omit(match(names( intersect.760591.Result@alleleDb[[thisMarker]]), names(knownAlleleDb[[thisMarker]]$alleleMap )))), "\n"))
}


for(thisMarker in names(intersectMarkerList)) {
	cat(paste(thisMarker, length(na.omit(match(names( intersect.catResult@alleleDb[[thisMarker]]), names(knownAlleleDb[[thisMarker]]$alleleMap )))), "\n"))
}





thisMarker <- "DPA1_E2"

alleleMatchTable <- cbind(as.character(intersectLiverpool.Result@alleleDb[[thisMarker]]), as.character(knownAlleleDb[[thisMarker]]$alleleMap[ match(names(intersectLiverpool.Result@alleleDb[[thisMarker]]), names(knownAlleleDb[[thisMarker]]$alleleMap ))]))
alleleMatchTable 
intersectLiverpool.Result@markerSampleList[[thisMarker ]]


thisMarker <- "A_E2"

 names(knownAlleleDb[[thisMarker]]$alleleMap[1])

names(intersectLiverpool.Result@alleleDb[[thisMarker]][1])

getClass("mlgtResult")

length(intersectLiverpool.Result@alleleDb[[thisMarker]])


#########################################################################
#################################################### DEVELOPMENT 

testAlign <- read.msf("C:/Users/dave/HLA/data/IMGT_manualDownload/DQA1_nuc.msf")
#testAlign <- read.msf("C:/Users/dave/HLA/data/IMGT_manualDownload/DQA_gen.msf")

nchar(testAlign[[3]][1])		# length of one sequence
length(testAlign[[3]])			# number of sequences  (also in $nb)

thisMarker <- "DQA1_E2"
setwd("C:/Users/dave/HalfStarted/mlgt/")
write.fasta(testAlign[[3]], testAlign[[2]], file="DQA1_nuc.mlgt.fasta")
write.fasta(markerList[[thisMarker]][1], thisMarker, file="DQA1_E2_markerSeq.fasta")
#muscle -profile -in1 existing_aln.afa -in2 new_seq.fa -out combined.afa
markerToAlleleDbAlign <- "DQA1_E2_align_IMGT.fasta"
muscleCommand <- paste(musclePath, "-quiet -profile -in1", "DQA1_nuc.mlgt.fasta", "-in2", "DQA1_E2_markerSeq.fasta", "-out" ,markerToAlleleDbAlign )
system(muscleCommand)


### this section copied from getSubSeqsTable()  Could be recoded as function?

# Extract portion corresponding to reference. 

	rawAlignment <- read.fasta(markerToAlleleDbAlign , as.string=T)		# do not use read.alignment() - broken
	alignedMarkerSeq <- s2c(rawAlignment[[thisMarker]])
	subStart <- min(grep("-",alignedMarkerSeq ,invert=T))
	subEnd <- max(grep("-",alignedMarkerSeq ,invert=T))
	alignedSubSeqs <- lapply(rawAlignment, FUN=function(x)	substr(x[1], subStart, subEnd))
	subAlignFile <- paste("test", thisMarker, "IMGT", "sub.align.fasta",sep=".")  #"localAlign.fasta"
	#subAlignFile <- paste(runPath, subAlignFileName , sep="/")
	write.fasta(alignedSubSeqs , names(alignedSubSeqs ), file=subAlignFile )

	alignedSubTable <- data.frame(name =  names(alignedSubSeqs ) , subSeq.aligned= as.character(unlist(alignedSubSeqs )))
	alignedSubTable$subSeq.stripped <-  gsub("-", "",alignedSubTable$subSeq.aligned )

	
## remove markerSeq?


### what to do if markerSeq forms an overhang at one end?  Could extract just the same, then remove marker seq.
# All others will begin with "----" and these could be removed by removing all "-" columns.


#### TODO



### how many unique alleles remain. 
## Need to use aligned seq or gap-stripped seq?
## Will need gap stripped seq for comparison with results from mlgt.

nrow(alignedSubTable )

length(alignedSubTable$subSeq)

length(unique(alignedSubTable$subSeq))

alleleMap <- split(as.character(alignedSubTable$name), alignedSubTable$subSeq.stripped)




###
## e.g. of matching of allleles. 
catResult@alleleDb[[thisMarker]][1]
names(catResult@alleleDb[[thisMarker]][1])

match(names(catResult@alleleDb[[thisMarker]]), names(alleleMap) )

catResult@alleleDb[[thisMarker]][1]
alleleMap[11]


catResult@markerSampleList[[thisMarker]]

catResult@alleleDb[["DRB1_E2"]][1]



