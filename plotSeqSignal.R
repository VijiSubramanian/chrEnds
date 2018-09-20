##################################################################################################
####                                                                                         #####
####         		function: Plot signal from ChIP-Seq data                 	     #####
####                                                                                         #####
####                        Viji Subramanian: June 21, 2017                                  #####
####                                                                                         ##### 																						 #####
####                                                                                         #####                          																 #####
####                                                                                         #####
##################################################################################################

# genome data SK1 Yue Nature Genetics 2017
# scripts from hwglabr2, using Enriched Heat Maps
# modified to map to sacCer2 also March 25, 2018

plotSeqSignal = function(Sample, genome = 'SK1Yue', tw=300) {
    ptm <- proc.time()
    
    if (genome == 'SK1Yue') {
    genome_info <- hwglabr2::get_chr_coordinates(genome)
	}
	
	if (genome == 'sacCer2') {
    genome_info <- get(load("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/sacCer2/sacCer2cen.RData"))
            
        # remove rDNA information from Spo11-oligo data
        SampleXIIl <- Sample[GenomicRanges::seqnames(Sample) == "chrXII" & GenomicRanges::start(Sample) <= 451000 ] 
        SampleXIIr <- Sample[GenomicRanges::seqnames(Sample) == "chrXII" & GenomicRanges::start(Sample) >= 471000 ] 
        subSample <- Sample[GenomicRanges::seqnames(Sample) != "chrXII"]
        Sample <- c(SampleXIIl, SampleXIIr, subSample)

	}

	
	# Get genome-wide mean and normalize Sample data
	message('Computing average signal ...')
	avrg <- function(x) (sum(GenomicRanges::width(x) * GenomicRanges::score(x), na.rm=TRUE) / sum(GenomicRanges::width(x)))
	genome_wide_mean <- avrg(Sample)	
	Sample$score <- Sample$score/genome_wide_mean
	
	# Sort sequences and levels to make sure they match
	sortedSample <- sort(GenomeInfoDb::sortSeqlevels(Sample))
	genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

	# Add chromosome length info to signal object
	GenomeInfoDb::seqlengths(sortedSample) <- GenomeInfoDb::seqlengths(genome_info)

	# Compute 100-bp tiling windows
	bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(sortedSample),
                                  tilewidth=tw, cut.last.tile.in.chrom=TRUE)
                                  
	# Get signal as "RleList"; the signal is stored in the "score" metadata column
	score <- GenomicRanges::coverage(sortedSample, weight="score")

	# Compute average signal per tile
	bins <- GenomicRanges::binnedAverage(bins, score, "binned_score")
	
    # Get positions as the midpoints of the intervals
	positions <- bins@ranges@start + floor(bins@ranges@width / 2)

	# Make data frame (convert positions to Kb; signal is the binned score)
	df <- data.frame(chr=bins@seqnames, position=positions / 1000, signal=bins$binned_score)
    df$chrSize <- GenomeInfoDb::seqlengths(bins)[df$chr]
    df <- dplyr::arrange(df, chrSize)
	message(paste0('\n\nCompleted in ', round((proc.time()[3] - ptm[3]), 2), ' sec.\n'))
	return(df)
	}
#source("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/SK1Yue_Bed/plotSeqSignal.R")
#source("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/SK1Yue_Bed/plotCen.R")
#df2 <- plotCen()
# par(las=1)
# library(ggplot2)
# AH6179Hop1T6all <- plotSeqSignal(AH6179Hop1T6)
# p <- ggplot(AH6179Hop1T6all, aes(x=position, y=signal)) + geom_line(colour="forestgreen") + geom_point(data =df2, aes(x=position, y=-2), color='grey30', size=1.5, pch=17)
# p + facet_grid(chrSize ~ .) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text.y = element_blank(), panel.background = element_blank())
