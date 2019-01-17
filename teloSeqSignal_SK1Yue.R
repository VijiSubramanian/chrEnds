##################################################################################################
####                                                                                         #####
####         function: Chromosome end enrichment of signal from ChIP-Seq data                #####
####                                                                                         #####
####                        Viji Subramanian: June 21, 2017                                  #####
#### 																						                                             #####
####                          																                               #####
####                                                                                         #####
##################################################################################################

# genome data SK1 Yue Nature Genetics 2017
# scripts from hwglabr2, using Enriched Heat Maps
# 
 

teloSeqSignal = function(Sample, remove_cen = FALSE, cen_region_length=50000, length_to_collect = 120000, Genome = 'SK1Yue') {
	ptm <- proc.time()
    library(hwglabr2)
    library(hwglabr)
    library(dplyr)
    library(magrittr)
    library(pbapply)
    library(readr)
	
	# Sample data is gRanges from either bedGraph or Wiggle data
	#AH6179T3Hop1_newBed <- hwglabr2::import_bedGraph("AH6179Hop1T3rep-118-214-Reps-SacCer3-B3W3-MACS2/AH6179Hop1T3rep-118-214-Reps-SacCer3-2mis_B3W3_MACS2_FE.bdg.gz")
	
	# load centromere data
	if (Genome == 'SK1Yue'){
		cen <- hwglabr2::get_chr_coordinates('SK1Yue')	
	}
	
	
	# if removing cen
	if (remove_cen){
		# get centromere midpoint and length to remove from either side of cen midpoint
		half_length <- floor(cen_region_length / 2) 
		offset <- floor(GenomicRanges::width(cen) / 2)
		GenomicRanges::start(cen) <- (GenomicRanges::start(cen) + offset - half_length)
		GenomicRanges::end(cen) <- (GenomicRanges::end(cen) - offset + half_length)
		# remove centromere regions
		Sample <- Sample[!IRanges::overlapsAny(Sample, cen)]
	}
	
	# Calculate genome average
	Sample_genomeAvg <- hwglabr2::average_chr_signal(Sample)[[2]]
	
	# Collect signal from chromosome ends
	Sample_telo2 <- hwglabr2::signal_from_telomeres2(Sample, length_to_collect, genome = Genome)
	
	# remove first three non-data columns
	Sample_telo1 <- Sample_telo2[,4:ncol(Sample_telo2)]
	
	# average data from different chromosomes
	Sample_telo0 <- colMeans(Sample_telo1, na.rm = T)
	
	# normalize to genome average
	Sample_teloGA <- Sample_telo0/Sample_genomeAvg
    
    #convert to dataframe
    Sample_teloGADF <- data.frame(position=seq(1, length_to_collect), value=Sample_teloGA)
	
    # function to compress data
    #source('/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/SK1Yue_Bed/Compress.R')
    Sample_teloCompressed <- hwglabr2::compress_signal_track(Sample_teloGADF, window_size = 200)
    
    # smooth data
    Sample_telo <- ksmooth(x=Sample_teloCompressed$position, y=Sample_teloCompressed$window_mean, bandwidth = 25000)
    
    message(paste0('\n\nCompleted in ', round((proc.time()[3] - ptm[3]), 2), ' sec.\n'))
    
	return(Sample_telo)
	}


########

######## source above script
#source("teloSeqSignal_SK1Yue.R")

######## calculate signal in telomeres
#AH6179Hop1T3telo <- teloSeqSignal(AH6179Hop1T3)
#AH6179Hop1T6telo <- teloSeqSignal(AH6179Hop1T6)

######## plot telo signal
#par(las=1)
#plot(AH6179Hop1T6telo, col="forestgreen", xlab='Distance to chr end', ylab='Average Hop1 signal', type='l', bty='n', lwd=5, ylim = c(0.5,2))
#lines(AH6179Hop1T3telo, type='l', lwd=5, col="darkseagreen")
#abline(h = 1, col = "gray60", lwd=3, lty=3)
#legend(80000,1.79, c("T = 3hr ndt80","T = 6hr ndt80"), text.col = c("darkseagreen", "forestgreen"), bty = "n")


