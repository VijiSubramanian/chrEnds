########################################################################
########################################################################
####                                                                ####
####               Created by Viji September 27, 2017               ####
####                                                                ####
####                                                                ####
########################################################################
########################################################################

# to collect mean signal/chromosome with or without the EARs (end-adjacent domains) using granges 

noTeloChrSize = function(Sample, genome, lengthToCollect=110000) {
	library(magrittr)
	library("GenomicRanges")
	
	ptm <- proc.time()
	if (genome == "SK1Yue"){
        chrLen <- read.table('/Volumes/LabShare/GenomeSequences/SK1_Yue_et_al_2017/chr_len.txt',header=F)
	}
	
	if (genome == "sacCer2"){
        chrLen <- read.table('/Volumes/LabShare/GenomeSequences/SacCer2_chrLen.txt',header=F)
	}
	
	if (!genome %in% c('SK1Yue', 'sacCer2')){
		stop('"genome" must be either "SK1Yue" or "sacCer2"', call. = FALSE)
    }
    
	# sort Sample data
	Sample <- sort(GenomeInfoDb::sortSeqlevels(Sample))
	Sample <- sort(Sample)
	avrg <- function(x) (sum(GenomicRanges::width(x) * GenomicRanges::score(x))/ sum(GenomicRanges::width(x)))	
	#normalize to each dataset to genome average: for sampleChrLarge and for intChrAll
    sample_avrg <- avrg(Sample)
    Sample$score <- Sample$score / sample_avrg
    message("normalized data to genome average, now capturing telomere-distal regions...")
    # creat empty GRanges objects
    intChrAll <- GRanges()
    intChrLarge <- GRanges()
    chrAll <- paste0("chr",as.roman(1:16))
    #chrLarge <- paste0("chr",as.roman(1:16)) [!paste0("chr",as.roman(1:16)) %in% c("chrI", "chrIII", "chrVI",  "chrXII")]
    # capture Sample data for only the large chromosomes
    #sampleChrLarge <- GRanges()
    #sampleChrLarge <- Sample[GenomicRanges::seqnames(Sample) %in% chrLarge]
    
    message("Chopping EARs ...")
    # chop EAR regions (110 Kb from telomeres) from the chromosomes
    for (i in chrAll){
        # assign length of the chromosome
        len <- chrLen[chrLen[,1]==i,2]
        # collect data for the internal regions while dropping 110 Kb (begin) from either end of the chromosome
        intChr <- Sample[GenomicRanges::seqnames(Sample) == i & GenomicRanges::start(Sample) >= (lengthToCollect + 1) & GenomicRanges::end(Sample) <= (len - lengthToCollect -1)]
        # append data for each chromsome to a new object
        intChrAll <- c(intChrAll, intChr)
    }
    
    #### chop 100 Kb from rDNA borders on chrXII
    ##xii_L <- intChrAll[GenomicRanges::seqnames(intChrAll) == "chrXII" & GenomicRanges::start(intChrAll) <= (447012 - lengthToCollect)]
    ##xii_R <- intChrAll[GenomicRanges::seqnames(intChrAll) == "chrXII" & GenomicRanges::start(intChrAll) >= (461699 + lengthToCollect)]
    ##GenomicRanges::start(xii_R) <- GenomicRanges::start(xii_R)-200000
	##GenomicRanges::end(xii_R) <- GenomicRanges::end(xii_R)-200000
    ##xii <- c(xii_L, xii_R)
    ##intChrAll <- dropSeqlevels(intChrAll, "chrXII")
    ##intChrAll <- c(intChrAll, xii)
    
    # sort GRanges data
    intChrAll <- sort(GenomeInfoDb::sortSeqlevels(intChrAll))
    intChrAll <- sort(intChrAll)
    
	message("calculating standard deviation values for telomere-distal regions...")
    # create 10 bins of equal size for each of the int regions from each chromosome
 	# to calculate median and standard deviation
    bins <- list()
    stdDev <- list()
    for (i in chrAll){
    	# create 10 equal bins. Cut function creates ranges (same last value for a range as the next first value, so kind of overlapping)
    	x <- cut(c(GenomicRanges::start(intChrAll[intChrAll@seqnames==i][1]),GenomicRanges::start(tail(intChrAll[intChrAll@seqnames==i],1))), 10)
  		# capture the first value of each of the 10 ranges
    	y <- vector()
    	for(n in 1:10) {y <- c(y, as.numeric(strsplit(strsplit(levels(x)[n],",")[[1]][1],"\\(")[[1]][2]))}
    	# capture the second value of each of the 10 ranges
    	z <- vector()
    	for(n in 1:10) {z <- c(z, as.numeric(strsplit((strsplit(levels(x)[n],",")[[1]][2]),"\\]")[[1]][1]))}
    	# chop the internal gRanges data into the 10 bins
    	for(n in 1:10) {bins[[i]] <- c(bins[[i]], avrg(intChrAll[GenomicRanges::seqnames(intChrAll) == i & GenomicRanges::start(intChrAll) >= y[n] & GenomicRanges::end(intChrAll) < z[n]]))}
    	# calculate std deviation and median for each bin
    	stdDev[[i]] <- sd(bins[[i]], na.rm=TRUE)
    	#median[[i]] <- median(bins[[i]], na.rm=TRUE)
    	}

	message("calculating standard deviation values for entire Sample...")
    # create 10 bins of equal size for each of the int regions from each chromosome
 	# to calculate median and standard deviation
    bins <- list()
    stdDevAll <- list()
    for (i in chrAll){
    	# create 10 equal bins. Cut function creates ranges (same last value for a range as the next first value, so kind of overlapping)
    	x <- cut(c(GenomicRanges::start(Sample[Sample@seqnames==i][1]),GenomicRanges::start(tail(Sample[Sample@seqnames==i],1))), 10)
  		# capture the first value of each of the 10 ranges
    	y <- vector()
    	for(n in 1:10) {y <- c(y, as.numeric(strsplit(strsplit(levels(x)[n],",")[[1]][1],"\\(")[[1]][2]))}
    	# capture the second value of each of the 10 ranges
    	z <- vector()
    	for(n in 1:10) {z <- c(z, as.numeric(strsplit((strsplit(levels(x)[n],",")[[1]][2]),"\\]")[[1]][1]))}
    	# chop the internal gRanges data into the 10 bins
    	for(n in 1:10) {bins[[i]] <- c(bins[[i]], avrg(Sample[GenomicRanges::seqnames(Sample) == i & GenomicRanges::start(Sample) >= y[n] & GenomicRanges::end(Sample) < z[n]]))}
    	# calculate std deviation and median for each bin
    	stdDevAll[[i]] <- sd(bins[[i]], na.rm=TRUE)
    	#medianAll[[i]] <- median(bins[[i]], na.rm=TRUE)
    	}
    	
    message("preparing dataframes for chromosome size bias")       
    #mean signal for each chromosome 
    sampleSeq_avrg <- sapply(GenomicRanges::split(Sample, GenomicRanges::seqnames(Sample)),avrg)
    intChrSeq_avrg <- sapply(GenomicRanges::split(intChrAll, GenomicRanges::seqnames(intChrAll)),avrg)
    
    # convert to dataframe
    sampleSeq_avrg <- data.frame(chr=names(sampleSeq_avrg), signal=sampleSeq_avrg, row.names=NULL, stringsAsFactors = F)
    intChrSeq_avrg <- data.frame(chr=names(intChrSeq_avrg), signal=intChrSeq_avrg, row.names=NULL, stringsAsFactors = F)
    
    # remove the small chromosome empty rows
    #sampleSeq_avrg <- na.omit(sampleSeq_avrg)
    #intChrSeq_avrg <- na.omit(intChrSeq_avrg)
    
    # add a column for chrSize   
    for (i in sampleSeq_avrg$chr){
    	sampleSeq_avrg[sampleSeq_avrg[,1]==i,3] <-  chrLen[chrLen[,1]==i,2]
    	sampleSeq_avrg[sampleSeq_avrg[,1]==i,4] <-  stdDevAll[[i]]
    }
    colnames(sampleSeq_avrg)=c("chr", "avrg_signal", "chrSize", "Std Dev")
    
    for (i in intChrSeq_avrg$chr){
    	intChrSeq_avrg[intChrSeq_avrg[,1]==i,3] <-  chrLen[chrLen[,1]==i,2]
    	intChrSeq_avrg[intChrSeq_avrg[,1]==i,4] <-  stdDev[[i]]
    	#intChrSeq_avrg[intChrSeq_avrg[,1]==i,5] <-  variance[[i]]
    }
    colnames(intChrSeq_avrg)=c("chr", "avrg_signal", "chrSize", "Std Dev")
    
    # arrange according to chr size
    sampleSeq_avrg <- dplyr::arrange(sampleSeq_avrg, sampleSeq_avrg$chrSize)
    intChrSeq_avrg <- dplyr::arrange(intChrSeq_avrg, intChrSeq_avrg$chrSize)
    
    message(paste0('\n\nCompleted in ', round((proc.time()[3] - ptm[3]), 2), ' sec.\n'))
    # return data
    return(list(sampleSeq_avrg, intChrSeq_avrg))
}

#######

####### source above script
#source("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/SK1Yue_Bed/noTeloChrSize.R")

####### input bedgraph data
#AH6179Hop1T3_noTelo <- noTeloChrSize(AH6179Hop1T3, genome = "SK1Yue", lengthToCollect = 110000)
#AH6179Hop1T6_noTelo <- noTeloChrSize(AH6179Hop1T6, genome = "SK1Yue", lengthToCollect = 110000)

####### plot ChrSize-EARs
#plot(log2(AH6179Hop1T6_noTelo[[1]][,3]), log2(AH6179Hop1T6_noTelo[[1]][,2]), type="p", ylim=c(-1,1.5), xlim=c(17.7, 20.7), pch = 16, cex = 3, col="#F2AB1E", xlab="Chromosome Size (Kb) log2 scale", ylab="Average binding log2 scale", bty='n')
#points(log2(AH6179Hop1T6_noTelo[[2]][,3]), log2(AH6179Hop1T6_noTelo[[2]][,2]), pch = 16, cex = 3, col="magenta4" )
#abline(lm(log2(AH6179Hop1T6_noTelo[[1]][,2])~log2(AH6179Hop1T6_noTelo[[1]][,3])), col = "#F2AB1E", lwd=3, lty=3)
#abline(lm(log2(AH6179Hop1T6_noTelo[[2]][,2])~log2(AH6179Hop1T6_noTelo[[2]][,3])), col = "magenta4", lwd=3, lty=3)
#arrows( log2(AH6179Hop1T6_noTelo[[1]][,3]), log2(AH6179Hop1T6_noTelo[[1]][,2] - AH6179Hop1T6_noTelo[[1]][,4]), log2(AH6179Hop1T6_noTelo[[1]][,3]), log2(AH6179Hop1T6_noTelo[[1]][,2] + AH6179Hop1T6_noTelo[[1]][,4]), length=0.05, angle=90, code=3, col="#F0AB20")
#arrows( log2(AH6179Hop1T6_noTelo[[2]][,3]), log2(AH6179Hop1T6_noTelo[[2]][,2] - AH6179Hop1T6_noTelo[[2]][,4]), log2(AH6179Hop1T6_noTelo[[2]][,3]), log2(AH6179Hop1T6_noTelo[[2]][,2] + AH6179Hop1T6_noTelo[[2]][,4]), length=0.05, angle=90, code=3, col="magenta4")
#legend(19.2,1.4, c("T=6 ndt80 Hop1 ChIP", "T=6 ndt80 Hop1 ChIP, no EARs"), text.col = c( "#F0AB20", "magenta4"), bty = "n")


