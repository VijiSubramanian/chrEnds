########################################################################
########################################################################
####                                                                ####
####    Created by Viji from hwglabr2 scripts September 15, 2017    ####
####                                                                ####
####                                                                ####
########################################################################
########################################################################

# to collect telomere-adjacent and internal signal values using granges
# script to calculate mean signal from both chromosome ends. 
# Can specify beginning and end positions of these data collectively for the ends.

teloBedBarRanges20 = function(Sample, genome, begin=20000, lengthToCollect=110000) {
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
	  rightArmAll <- GRanges()
    leftArmAll <- GRanges()
    intChrAll <- GRanges()
    
    for (i in chrLen[,1]){
        # assign length of the chromosome
        len <- chrLen[chrLen[,1]==i,2]
        # subset data according to chromosome number
        # collect data for the right arm, left arm and internal regions while dropping 30 Kb (begin) from either end of the chromosome
        rightArm <- Sample[GenomicRanges::seqnames(Sample) == i & 
			   GenomicRanges::start(Sample) >= (len  - lengthToCollect) & GenomicRanges::end(Sample) <= (len - begin)]
        leftArm <- Sample[GenomicRanges::seqnames(Sample) == i & 
			  GenomicRanges::start(Sample) >= (begin) & GenomicRanges::end(Sample) <= (lengthToCollect)]
        intChr <- Sample[GenomicRanges::seqnames(Sample) == i & 
			 GenomicRanges::start(Sample) >= (lengthToCollect + 1) & GenomicRanges::end(Sample) <= (len - lengthToCollect -1)]
        
        # append data for each chromsome to a new object
        rightArmAll <- c(rightArmAll, rightArm)
        leftArmAll <- c(leftArmAll, leftArm)
        intChrAll <- c(intChrAll, intChr)
    }
    
    # sort GRanges data
    rightArmAll <- sort(GenomeInfoDb::sortSeqlevels(rightArmAll))
    rightArmAll <- sort(rightArmAll)
    leftArmAll <- sort(GenomeInfoDb::sortSeqlevels(leftArmAll))
    leftArmAll <- sort(leftArmAll)
    intChrAll <- sort(GenomeInfoDb::sortSeqlevels(intChrAll))
    intChrAll <- sort(intChrAll)
    #genome <- c(leftArmAll, intChrAll, rightArmAll)
    #genome <- sort(GenomeInfoDb::sortSeqlevels(genome))
    genome <- sort(genome)
    ends <- c(rightArmAll, leftArmAll)
    ends <- sort(GenomeInfoDb::sortSeqlevels(ends))
    ends <- sort(ends)
    
    # calculate genome average excluding 30 Kb at either ends of all chromosomes
    avrg <- function(x) (sum(GenomicRanges::width(x) * GenomicRanges::score(x))/ sum(GenomicRanges::width(x)))
    genome_avrg <- avrg(Sample)
    meanEnds <- (avrg(ends))/genome_avrg
    meanInt <- (avrg(intChrAll))/genome_avrg
    
    #mean for each chromosome signal at ends
    listEnds <- sapply(GenomicRanges::split(ends, GenomicRanges::seqnames(ends)),avrg)
    # convert to dataframe
    listEnds <- data.frame(chr=names(listEnds), avrg_signal=listEnds, row.names=NULL, stringsAsFactors = F)
    listEnds$avrg_signal <- listEnds$avrg_signal/genome_avrg
    
    #mean for each chromosome signal at internal regions
    listInt <- sapply(GenomicRanges::split(intChrAll, GenomicRanges::seqnames(intChrAll)),avrg)
    # convert to dataframe
    listInt <- data.frame(chr=names(listInt), avrg_signal=listInt, row.names=NULL, stringsAsFactors = F)
    listInt$avrg_signal <- listInt$avrg_signal/genome_avrg
    
    # return data
    return(list(meanEnds, meanInt, listEnds, listInt))
}
