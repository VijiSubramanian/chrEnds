
##################################################################################################
####                                                                                         #####
####         		function: Chromosome size plot from ChIP-Seq data                		 #####
####                                                                                         #####
####                        Viji Subramanian: July 7, 2017                                   #####
#### 																						 #####
####                          																 #####
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
	
# par(las=1)
# plot(log2(AH574Hop1T2_avgChr[,3]), log2(AH574Hop1T2_avgChr[,2]), type="p", ylim=c(-0.4,1.3), xlim=c(17.7, 20.7), pch = 16, cex = 3, col="ivory3", xlab="Chromosome Size (lKb) log2 scale", ylab="Average Hop1 binding (log2 scale)", bty='n')
# points(log2(AH574Hop1T3_avgChr[,3]), log2(AH574Hop1T3_avgChr[,2]), type="p", pch = 16, cex = 3, col="ivory4")
# points(log2(AH574Hop1T4_avgChr[,3]), log2(AH574Hop1T4_avgChr[,2]), type="p", pch = 16, cex = 3, col="grey30")
#plot(log2(AH6179T3Hop1rep_avgChr[,3]), log2(AH6179T3Hop1rep_avgChr[,2]), type="p", ylim=c(-0.3,1.2), xlim=c(17.7, 20.7), pch = 16, cex = 3, col="lawngreen", xlab="Chromosome Size (Kb) log2 scale", ylab="Average binding log2 scale", bty='n')
#points(log2(AH6179T6Hop1rep_avgChr[,3]), log2(AH6179T6Hop1rep_avgChr[,2]), type="p", pch = 16, cex = 3, col="forestgreen")
#abline(lm(log2(AH6179T3Hop1rep_avgChr[,2])~log2(AH6179T3Hop1rep_avgChr[,3])), col = "lawngreen", lwd=3, lty=3)
#abline(lm(log2(AH6179T6Hop1rep_avgChr[,2])~log2(AH6179T6Hop1rep_avgChr[,3])), col = "forestgreen", lwd=3, lty=3)

#plot(log2(AH6639T3Hop1rep_avgChr[,3]), log2(AH6639T3Hop1rep_avgChr[,2]), type="p", ylim=c(-0.3,0.8), xlim=c(17.7, 20.7), pch = 16, cex = 3, col="darkolivegreen1", xlab="Chromosome Size (Kb) log2 scale", ylab="Average binding log2 scale", bty='n')
#points(log2(AH6639T6Hop1rep_avgChr[,3]), log2(AH6639T6Hop1rep_avgChr[,2]), type="p", pch = 16, cex = 3, col="darkolivegreen4")
#abline(lm(log2(AH6639T3Hop1rep_avgChr[,2])~log2(AH6639T3Hop1rep_avgChr[,3])), col = "darkolivegreen1", lwd=3, lty=3)
#abline(lm(log2(AH6639T6Hop1rep_avgChr[,2])~log2(AH6639T6Hop1rep_avgChr[,3])), col = "darkolivegreen4", lwd=3, lty=3)

#AH6179T3T6Hop1rep_avgChrAnova <- anova_seq_signal(AH6179T3Hop1rep,AH6179T6Hop1rep)
##log2(chrSize):timepointlate -0.29832    0.05441  -5.483 7.42e-06 ***
#AH6639Hop1anova <- anova_seq_signal(AH6639T3Hop1rep,AH6639T6Hop1rep)
##log2(chrSize):timepointlate   0.2041     0.0512   3.987 0.000435 ***