
setwd("C:/Users/dave/HalfStarted/mlgt/testProject/intersect756948/")

### want to make 'nice' dataset for inclusion in package.
markerBlastTable <- read.delim("C:/Users/dave/HalfStarted/mlgt/testProject/intersect756948/blastOut.markers.tab")
names(markerBlastTable) <- c("query", "subject", "percent.id", "ali.length", "mismatches", "gap.openings", "q.start","q.end", "s.start","s.end", "e.value", "bit.score")
#head(markerBlastTable)

markerBlastTable <- subset(markerBlastTable, ( subject=="B_E2" | subject=="A_E3" | subject=="DPA1_E2" | subject=="DQA1_E2"  ))
markerBlastTable <- subset(markerBlastTable, ali.length > 100 )
nrow(markerBlastTable)

table(markerBlastTable$subject)


##samples. 
wantedSamples <- paste("MID", c(1,3,5,8,11,17,22,31),sep="-")
sampleBlastTable <- read.delim("C:/Users/dave/HalfStarted/mlgt/testProject/intersect756948/blastOut.fTags.tab")
names(sampleBlastTable ) <- c("query", "subject", "percent.id", "ali.length", "mismatches", "gap.openings", "q.start","q.end", "s.start","s.end", "e.value", "bit.score")
#head(sampleBlastTable )

wantedIndex <- !is.na(match(sampleBlastTable$subject,wantedSamples))
sampleBlastTable <- sampleBlastTable[wantedIndex ,]

wantedReads <- intersect(sampleBlastTable$query, markerBlastTable$query)
length(wantedReads)

# add some random reads as noise? 
randReads <- as.character(sample(markerBlastTable$query,200))
wantedReads <- unique(c(wantedReads,randReads))
length(wantedReads)

wantedReadIdFile <- "wantedReadIds.txt"
write(wantedReads, file=wantedReadIdFile )

wantedReadsOutputFile <- "sampleSequences.fasta"
fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" ,  "-o", wantedReadsOutputFile ,  "-i", wantedReadIdFile )
system(fastacmdCommand)





