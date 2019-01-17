##################################################################################################
####                                                                                         #####
####         function: Centromere enrichment of signal from ChIP-Seq data                    #####
####                                                                                         #####
####                        Viji Subramanian: January 11, 2017                               #####
####               modified to include sacCer2, SK1Yue: March 25, 2018                       #####
##################################################################################################

cenSeqSignal = function(Sample, genome = 'SK1Yue', extend = 5000) {
  	ptm <- proc.time()
  	library(hwglabr2)
  	library(dplyr)
  	library(magrittr)
  	library(pbapply)
	library(GenomicRanges)
	library(EnrichedHeatmap)
    
    if (genome == 'SK1Yue') {
    genome_info <- hwglabr2::get_chr_coordinates(genome)
	}
	
	if (genome == 'sacCer2') {
    genome_info <- get(load("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/sacCer2/sacCer2cen.RData"))
	}
	
    # Get genome-wide mean and normalize Sample data
	message('Computing average signal ...')
	avrg <- function(x) (sum(GenomicRanges::width(x) * GenomicRanges::score(x), na.rm=TRUE) / sum(GenomicRanges::width(x)))
	genome_wide_mean <- avrg(Sample)	
	Sample$score <- Sample$score/genome_wide_mean
	
	#sampleName <- deparse(substitute(Sample))
	
	# Get centromeres
	midpoint <- floor(width(genome_info) / 2)
	start(genome_info) <- start(genome_info) + midpoint
	end(genome_info) <- start(genome_info)
	
	signal_at_cens <- EnrichedHeatmap::normalizeToMatrix(Sample, genome_info,
                                                     extend, w=1,
                                                     mean_mode="weighted",
                                                     value_column="score")
                                   
     signal_at_cens_avrg <- colMeans(signal_at_cens, na.rm = T)
     message(paste0('\n\nCompleted in ', round((proc.time()[3] - ptm[3]), 2), ' sec.\n'))
     return(signal_at_cens_avrg)
	}

####### to plot ChIP-Seq signal around centromere

####### source above script
#source("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/SK1Yue_Bed/cenSeqSignal_SK1Yue.R")

####### Collect signal around centromeres
#AH6179Hop1T3cen <- cenSeqSignal(AH6179Hop1T3)
#AH6179Hop1T6cen <- cenSeqSignal(AH6179Hop1T6)

####### plot the signal
par(las=1)
#plot(x=seq(-4999, 5000), y=AH6179Hop1T6cen, col="forestgreen", xlab='Distance from cen', ylab='Average Hop1 signal', type='l', bty='n', lwd=5, ylim = c(0.4,7))
#lines(x=seq(-4999, 5000), y=AH6179Hop1T3cen, type='l', lwd=5, col="darkseagreen")
#abline(h = 1, col = "gray60", lwd=3, lty=3)
#legend(2000,2.49, c("T = 3hr ndt80","T = 6hr ndt80"), text.col = c("darkseagreen", "forestgreen"), bty = "n")


