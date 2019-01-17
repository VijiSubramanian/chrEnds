########################################################################
########################################################################
####                                                                ####
####                Created by Viji November 06, 2018               ####
####                                                                ####
####                                                                ####
########################################################################
########################################################################

# to collect centromere-adjacent signal values using granges

cenBedBarRanges = function(Sample, genome="sacCer2", lengthToCollect=5000) {
	library(magrittr)
	library("GenomicRanges")
	
	if (genome == "SK1Yue"){
        chrLen <- read.table('/Volumes/LabShare/GenomeSequences/SK1_Yue_et_al_2017/chr_len.txt',header=F)
        }
	
	if (genome == "sacCer2"){
        chrLen <- read.table('/Volumes/LabShare/GenomeSequences/SacCer2_chrLen.txt',header=F)
        
        # remove rDNA information from Spo11-oligo data
        SampleXIIl <- Sample[GenomicRanges::seqnames(Sample) == "chrXII" & GenomicRanges::start(Sample) <= 451000 ] 
        SampleXIIr <- Sample[GenomicRanges::seqnames(Sample) == "chrXII" & GenomicRanges::start(Sample) >= 471000 ] 
        subSample <- Sample[GenomicRanges::seqnames(Sample) != "chrXII"]
        Sample <- c(SampleXIIl, SampleXIIr, subSample)
		centromeres <- get(load("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/sacCer2/sacCer2cen.RData"))
		midpoint <- floor(GenomicRanges::width(centromeres) / 2)
		GenomicRanges::start(centromeres) <- GenomicRanges::start(centromeres) + midpoint
		GenomicRanges::end(centromeres) <- GenomicRanges::start(centromeres)
	}
	
	if (genome == "mouse"){
        chrLen <- read.table('/Volumes/LabShare/Viji/OthersData/langeKeeney/mouseChrSize.txt',header=T)
	}
	
	if (!genome %in% c('SK1Yue', 'sacCer2', 'mouse')){
		stop('"genome" must be either "SK1Yue" or "sacCer2"', call. = FALSE)
    }
    
	# sort Sample data
	Sample <- sort(GenomeInfoDb::sortSeqlevels(Sample))
	Sample <- sort(Sample)
	
    # creat empty GRanges objects
	cenAll <- GenomicRanges::GRanges()
    
    # Get genome-wide mean and normalize Sample data
	message('Computing average signal ...')
	avrg <- function(x) (sum(GenomicRanges::width(x) * GenomicRanges::score(x), na.rm=TRUE) / sum(GenomicRanges::width(x)))
	genome_wide_mean <- avrg(Sample)	
	Sample$score <- Sample$score/genome_wide_mean

    for (i in chrLen$V1){
        # assign length of the chromosome
        cenPos <- GenomicRanges::start(centromeres[GenomicRanges::seqnames(centromeres) == i]) 
        # subset data according to chromosome number
        # collect data for the right arm, left arm and internal regions while dropping 30 Kb (begin) from either end of the chromosome
        cen <- Sample[GenomicRanges::seqnames(Sample) == i & 
        	GenomicRanges::start(Sample) >= (cenPos  - lengthToCollect) & 
        	GenomicRanges::end(Sample) <= (cenPos  + lengthToCollect)]
        
        # append data for each chromsome to a new object
        cenAll <- c(cenAll, cen)
        }        
    meanCen <- (avrg(cenAll))
    
    #mean for each chromosome signal at ends
    listOfCen <- sapply(GenomicRanges::split(cenAll, GenomicRanges::seqnames(cenAll)),avrg)
    # convert to dataframe
    listCen <- data.frame(chr=names(listOfCen), avrg_signal=listOfCen, row.names=NULL, stringsAsFactors = F)
    
    # return data
    return(list(meanCen, listCen))
}


###### import data as bedgraph
#flag34Avg <- hwglabr2::import_bedGraph("flag34Avg.bedgraph")
#pch2Avg <- hwglabr2::import_bedGraph("pch2Avg.bedgraph")

###### sorce the script above
#source("cenBedBarRanges.R")
#wtCen_6Kb <- cenBedBarRanges(wtAvg, lengthToCollect = 3000)
#pch2Cen_6Kb <- cenBedBarRanges(pch2Avg, lengthToCollect = 3000)

###### plot data
#col <- c("forestgreen", "thistle4")
#barplot(c(flagCen_6Kb[[1]], pch2Cen_6Kb[[1]]), col = col)

###### statistical test
#wilcox.test(flagCen_6Kb[[2]]$avrg_signal, pch2Cen_6Kb[[2]]$avrg_signal)

