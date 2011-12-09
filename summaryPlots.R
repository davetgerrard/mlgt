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


