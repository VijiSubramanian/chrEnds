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

