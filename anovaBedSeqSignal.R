##################################################################################################
####                                                                                         #####
####         function: Linear Model testing/ANOVA of signal from chr Size data               #####
####                                                                                         #####
####                        Viji Subramanian: July 8, 2017                                   #####
####                                                                                         #####
##################################################################################################



#2017-07-08

anovaBedSeqSignal = function(Sample1, Sample2) {
    ptm <- proc.time()
    library(hwglabr)
    library(dplyr)
    library(magrittr)
    library(pbapply)
    library(readr)
    source("/Volumes/LabShare/Viji/Scripts/ChIPSeq_R/SK1Yue_Bed/plotChrSize.R")
    
    Sample1_avgChr <- plotChrSize(Sample1)
    Sample2_avgChr <- plotChrSize(Sample2)
    
    a <- mutate(Sample1_avgChr, timepoint='early')
    b <- mutate(Sample2_avgChr, timepoint='late')
    y <- rbind(a,b)
    M <- lm(log2(avrg_signal)~log2(chrSize)*timepoint, y)
	suMMary <- summary.lm(M)
	#anova(M)
    message(suMMary)
    message(paste0('\n\nCompleted in ', round((proc.time()[3] - ptm[3]), 2), ' sec.\n'))

    return(suMMary)
    }
