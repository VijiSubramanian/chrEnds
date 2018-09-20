##################################################################################################
####                                                                                         #####
####         					           function: Plot centromeres                					         #####
####                                                                                         #####
####                           Viji Subramanian: June 21, 2017                               #####
#### 																						                                             #####
####                          																                               #####
####                                                                                         #####
##################################################################################################

# genome data SK1 Yue Nature Genetics 2017
plotCen = function(genome = 'SK1Yue') {
	 ptm <- proc.time()
	 
	 # Get centromeres
	if (genome == 'SK1Yue') {
    centromeres <- hwglabr2::get_chr_coordinates(genome)
	}
	
	if (genome == 'sacCer2') {
    centromeres <- get(load("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/sacCer2/sacCer2cen.RData"))
	}

	# Replace start and end centromere positions by midpoint
	library(GenomicRanges)

	#midpoint <- floor(width(centromeres) / 2)
	#start(centromeres) <- start(centromeres) + midpoint
	#end(centromeres) <- start(centromeres)
	
	positions <- centromeres@ranges@start + floor(centromeres@ranges@width / 2)
	df2 <- data.frame(chr=centromeres@seqnames, position=positions / 1000)
    df2$chrSize <- GenomeInfoDb::seqlengths(centromeres)[df2$chr]
    df2 <- dplyr::arrange(df2, df2$chrSize)
    
	 message(paste0('\n\nCompleted in ', round((proc.time()[3] - ptm[3]), 2), ' sec.\n'))
	 return(df2)
	}
