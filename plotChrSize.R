
##################################################################################################
####                                                                                         #####
####         		function: Chromosome size plot from ChIP-Seq data                    #####
####                                                                                         #####
####                        Viji Subramanian: July 7, 2017                                   #####
#### 											     #####										 #####
####                          								     #####							 #####
####                                                                                         #####
##################################################################################################

# genome data SK1 Yue Nature Genetics 2017
# scripts from hwglabr2, using Enriched Heat Maps
# modified March 21, 2018 to include sacCer2

# 
plotChrSize = function(Sample, genome = 'SK1Yue') {
    ptm <- proc.time()
    if (genome == 'SK1Yue'){
    	genome_info <- hwglabr2::get_chr_coordinates(genome)
    	genome_wide_mean <- hwglabr2::average_chr_signal(Sample)[[2]]
        
        # Get genome-wide mean and normalize Sample data
        Sample$score <- Sample$score/genome_wide_mean        
        # chromosome averaged data
        SampleChrSize <- hwglabr2::average_chr_signal(Sample, genome = 'SK1Yue')
        # Make data frame 
		df <- data.frame(SampleChrSize$seq_avrg)
	}
	if (genome == "sacCer2"){
        genome_info <- get(load("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/sacCer2/sacCer2cen.RData"))
        source("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/sacCer2/average_chr_signal_sacCer2.R")
        genome_wide_mean <- average_chr_signal_sacCer2(Sample)[[2]]
        
        # Get genome-wide mean and normalize Sample data
		Sample$score <- Sample$score/genome_wide_mean	
		# chromosome averaged data
		SampleChrSize <- average_chr_signal_sacCer2(Sample)	
		df <- data.frame(SampleChrSize$Sample_avrg)	
	}
	# Make data frame 
    df$chrSize <- GenomeInfoDb::seqlengths(genome_info)[df$chr]
    df <- dplyr::arrange(df, df$chrSize)
	
	message(paste0('\n\nCompleted in ', round((proc.time()[3] - ptm[3]), 2), ' sec.\n'))
	return(df)
	}


####### 
####### source above script
#source("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/SK1Yue_Bed/plotChrSize.R")

###### use bedgraph data as input for source script 
#AH6179Hop1T3chrSize <- plotChrSize(AH6179Hop1T3)
#AH6179Hop1T6chrSize <- plotChrSize(AH6179Hop1T6)

###### plot chromosome size
#plot(log2(AH6179Hop1T3chrSize[,3]), log2(AH6179Hop1T3chrSize[,2]), type="p", ylim=c(-0.4,1.2), xlim=c(17.7, 20.7), pch = 15, cex = 3, col="darkseagreen", xlab="Chromosome Size (Kb) log2 scale", ylab="Average Hop1 binding (log2 scale)", bty='n')
#points(log2(AH6179Hop1T6chrSize[,3]), log2(AH6179Hop1T6chrSize[,2]), type="p", pch = 16, cex = 3, col="forestgreen")
#abline(lm(log2(AH6179Hop1T3chrSize[,2])~ log2(AH6179Hop1T3chrSize[,3])), col="darkseagreen", lwd=3, lty=3)
#abline(lm(log2(AH6179Hop1T6chrSize[,2])~ log2(AH6179Hop1T6chrSize[,3])), col="forestgreen", lwd=3, lty=3)
#legend(19.5,1.19, c("T = 3hr ndt80", "T = 6hr ndt80"), text.col = c( "darkseagreen", "forestgreen"), bty = "n")

###### calculate fit of the regression line
#AH6179Hop1T3chrSizeSlope <- capture.output(summary(lm(log2(avrg_signal)~log2(chrSize), AH6179Hop1T3chrSize)))
#cat(AH6179Hop1T3chrSizeSlope, file= "AH6179Hop1T3chrSizeSlope", sep = '\n', append = TRUE)
#AH6179Hop1T6chrSizeSlope <- capture.output(summary(lm(log2(avrg_signal)~log2(chrSize), AH6179Hop1T6chrSize)))
#cat(AH6179Hop1T6chrSizeSlope, file= "AH6179Hop1T6chrSizeSlope", sep = '\n', append = TRUE)

