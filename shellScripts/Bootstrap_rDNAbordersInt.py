"""
Bootstrap rDNA borders: Viji Subramanian 
modified Bootstrap centromeres by Tovah Markowitz: 3/1/17
updated to handle bedgraph files by Viji Subramanian: 6/26/17
excludes 110Kb from chromosome ends(EAR domains) from the sample

inputs: all chromosome unzipped wiggle file, bed/gff file, number of resamples, size of region, root of output file names
calculates: ratio of mean signal of regions versus mean signal of entire mapped genome
"""

###############################
# MODULES

import re
from collections import defaultdict
import random
import optparse
from datetime import datetime
import time, sys

###############################
# FUNCTIONS
                                                                                            
def read_wiggle(allWiggle):
    # read in all chromosome wiggle file
    # output is a dictionary with each chromosome having a separate key
    # and each key having a list of lists where each base has its own list
    # containing position and score
    lineNum = 0
    loop = 0
    chrPattern = re.compile( "chrom=(chr\w+)" )
    wiggleD = defaultdict( list )
    f = open( allWiggle, 'r' )
    for line in iter(f):
        lineNum += 1
        if line.startswith( 'track type' ):
            loop += 1
        elif 'Step' in line: 
            a = chrPattern.search( line )
            chr = a.group(1)
        else:
            wiggleD[chr].append( map( float, (line.strip().split('\t') ) ) )
    f.close()
    return( wiggleD )

def read_bedgraph(allWiggle):
    # Purpose: to read in a bedgraph file and create a sorted dict composed of lists of lists
    # step 1: read in bedgraph
    f = open(allWiggle, 'r')
    bedG = f.readlines()
    f.close()
    # step 2: organize bedgraph into a dict (chromosomes)
    # of lists (each row of bedgraph) of lists (start, end, score)
    for i in range( len(bedG) ):
        bedG[i] = bedG[i].strip().split('\t')        
    bedD = defaultdict(list)
    for i in bedG:
        if len(i) != 4:
            print( "Some rows in this bedgraph are not complete. Cannot create a wiggle file.\n" )
            exit()
        else:
            bedD[i[0]].append( list(map(float, i[1:])) )
    # step 3: for each chromosome do numeric sort by start
    for key in bedD.keys():
        bedD[key].sort()
    # step 4: convert to wiggle
    wiggleD = defaultdict( list )
    for chr in bedD.keys():
        print( chr + ": " + datetime.now().ctime() )
        for row in range( len( bedD[chr] ) ):
            if row != 0:
                    if bedD[chr][row][0] < bedD[chr][row-1][1]:
                            bedD[chr][row][0] = bedD[chr][row-1][1]
            if bedD[chr][row][2] != 0:
                positions = range( int( bedD[chr][row][0] ) + 1, int( bedD[chr][row][1] ) + 1 )
                wiggleD[chr]+= [ [position,bedD[chr][row][2]] for position in positions]
    return(wiggleD)

def readBedGff(inputName):
    # function to read in Gff or Bed files and convert into a universal format
    # read in input file
    f = open(inputName,'r')
    input = f.readlines()
    f.close()
    
    # split into a table
    input = [ row.strip().split('\t') for row in input if not row.startswith("#") ]
        
    # determine type of input file and make four column table
    # column 1: chr#, column 2: start, column 3: end, column 4: name, column 5: orientation
    if len(input[0])==9:          # gff file
        for i in range(len(input)):
            chr = input[i][0]
            input[i] = [ chr, int(input[i][3]), int(input[i][4]), input[i][8] , input[i][6]]
    elif len(input[0])==3:          # simplest bed file
        for i in range(len(input)):
            chr = input[i][0]
            input[i] = [ chr, int(input[i][1]), int(input[i][2]), 'id'+str(i), '+']
    elif len(input[0])==2:          # simplest txt file
        for i in range(len(input)):
            chr = input[i][0]
            input[i] = [ chr, int(0), int(input[i][1]), 'id'+str(i), '+']
    elif len(input[0])==5:          # bed file made for summits
        for i in range(len(input)):
            chr =input[i][0]
            input[i] = [ chr, int(input[i][1]), int(input[i][2]), input[i][3], '+']            
    return (input)

def random_signal_totals( bedFile, nRand, wiggleD, extend):
    # assumes 1 equally sized region per chromosome
    # allows wrapping across telomeres
    # output is two lists: bootstrapped signal for each iteration and the number of bases included in the calculation
    dataRand = [ ]
    randSignal = [ ]
    randLength = [ ]
    wiggleL = [ ]
    extend = int( extend )
    a = [ ]
    x = 0
    b = [ ] 
    bedD = readBedGff(bedFile)
    dataInt = defaultdict( list )
    for ends in bedD:
        chr = ends[0]
        leftStart = ends[1] + extend
        rightEnd = ends[2] - extend
        dataInt[chr] = [ i for i in wiggleD[chr] if (i[0] >= leftStart) & (i[0] <= rightEnd) ]
    # for each iteration of looping/number of randomizations
    for i in dataInt.keys():
        for row in dataInt[i]:
            x += 1
            wiggleL.append( [x, i, row ])
    # repeat the randomized iteration twice for each chromsome to account for 2X ends
    # for each iteration of looping/number of randomizations
    for i in range(0,int(nRand)):
        counter = 0
        while counter < 1:
            # find a random index value in the outer list from the entire genome
            a = random.randint(0,len(wiggleL)-1)
            # define start as the position at the index called
            start = wiggleL[a][0]
            end = start + int(extend)
            # if the chromosome is shorter than end:
            # fix numbering to begin at the beginning of the chromosome
            # get all wiggle data across the two regions
            # save total signal and total number of bases with sequencing data
            if wiggleL[-1][0] < end:
                newEnd = end - wiggleL[-1][0] + 1
                b = [ j[2][1] for j in wiggleL if (j[0] >= start) | (j[0] <= newEnd) ]
            else:
                # otherwise: just get wiggle data across region and save signal total and length
                b = [ j[2][1] for j in wiggleL if (j[0] >= start) & (j[0] <= end) ]
            dataRand.append( [ sum(b), len(b) ] )
            counter += 1
        #if (i % 200) == 0 :
        #    print "\rCompleted " + str(i) + "random samplings. " + datetime.now().ctime()
        # once all chromsomes are complete: total signal and lengths for all 16 "centromeres"
        randSignal = [ i[0] for i in dataRand  ]  
        randLength = [ i[1] for i in dataRand  ] 
    return( randSignal, randLength )

def genome_ave_signal( wiggleD ):
    # calculates the total signal and bases genome-wide
    dataGenome = defaultdict( list )
    for chr in wiggleD.keys():
        a = [ i[1] for i in wiggleD[chr] ]
        dataGenome[chr].append( [ sum(a), len(a) ] )
    genomeSignal = sum( [ i[0][0] for i in dataGenome.values() ] )
    genomeLength = sum( [ i[0][1] for i in dataGenome.values() ] )
    return( genomeSignal, genomeLength )

def real_data_ratio( wiggleD, extend ):
    # defines region "extend" around midpoints of features from bedFile
    # calculates signal and length of each feature and calculates totals
    extend = int( extend )
    a = [ ]
    dataReal = defaultdict( list )
    for chr in wiggleD.keys():
        if chr == 'chrXII':
            start = 447012 - extend
            end = 461699 + extend
            leftarm = [ i[1] for i in wiggleD[chr] if (i[0] >= start) & (i[0] <= 447012) ]
            rightarm = [ i[1] for i in wiggleD[chr] if (i[0] >= 461699) & (i[0] <= end) ]
            a = leftarm + rightarm
            dataReal[chr].append( [ sum(a), len(a) ] )
            endsSignal = sum( [ i[0][0] for i in dataReal.values() ] )
            endsLength = sum( [ i[0][1] for i in dataReal.values() ] )
    return (dataReal, endsSignal, endsLength )

def plot_hist( randRatio, realRatio, outFileRoot ):
    plt.hist(randRatio, bins = 20 )
    plt.axvline( realRatio, col='b', linestyle='dashed', linewidth=2 )
    plt.title( "Bootstrap of ratios" )
    plt.xlabel( "Feature average/ Average genome signal" )
    plt.ylabel( "Count" )
    plt.show()
    plt.savefig( outFileRoot + "_hist.png" )
    plt.close()
    
###############################
# MAIN

def main( allWiggle, extend, nRand, bedFile, outFileRoot ):
    print "Starting"
    print datetime.now().ctime()
    if ("wig" in allWiggle):
        wiggleD = read_wiggle(allWiggle)
        print "Wiggle file read."
    else:
        wiggleD = read_bedgraph(allWiggle)
        print "Bedgraph file read."
    print datetime.now().ctime()
    (genomeSignal, genomeLength) = genome_ave_signal( wiggleD )
    genomeAveSignal = genomeSignal/genomeLength
    print "Genome average calculated."
    print datetime.now().ctime()
    (dataReal, endsSignal, endsLength) = real_data_ratio( wiggleD, extend )
    endsRatio = (endsSignal/endsLength) / genomeAveSignal
    print "Ends ratio calculated."
    print datetime.now().ctime()
    g = open ( outFileRoot + '_realData_ext_' + str(extend) + '.txt','w')
    g.write( 'region\ttotalSignal\tNbases\n' )
    #for chr in chrOrder:
    g.write( 'chrXII' + '\t' + str(endsSignal) + '\t' + str(endsLength) + '\n' )
    g.write( 'genome\t' + str(genomeSignal) + '\t' + str(genomeLength) + '\n' )
    g.write( 'ratio:\t' + str(endsRatio) )
    g.close()
    (randSignal, randLength) = random_signal_totals(bedFile, nRand, wiggleD, extend)
    randRatio = [ ( randSignal[i]/randLength[i] ) / genomeAveSignal for i in range(0,len(randSignal)) ]
    print "Randomized ratios calculated."
    print datetime.now().ctime()
    f = open( outFileRoot + "_randRatio_" + str(nRand) + "_ext_" + str(extend) + ".txt" , 'w' )
    f.write( "\n".join( map( str, randRatio ) ) )
    f.close()


desc="""
A script to a) bootstrap the chromosome ends randomly across the genome
and b) to calculate the ratio of the ratio of binding at the 
chromosome ends versus the entire genome. Bootstrapped ratios are saved
in one file. A second file contains a table of information needed
to calculate the ratio of the real chromosome ends.
"""

# parse object for managing input options
#parser = optparse.OptionParser()
parser = optparse.OptionParser(description=desc)

# essental data, defines commandline options
parser.add_option('-w', dest= "wiggle", default= '', help= "This is the name \
 of the wiggle file with information about all chromosomes.")
parser.add_option('-e', dest= "extend", default= '120000', help= "This is the \
 size of the region to be bootstrapped/used for calculations. Default: 120 kb.")
parser.add_option('-n', dest= "nRand", default= '10000', help= "This is the number \
 of times to run the bootstrapping calculations. Default: 10000.")
parser.add_option('-b', dest= "bedFile", default= '', help= "This is the name \
 of the bed/gff file with information about the lengths of the chromosomes.")
parser.add_option('-o', dest= "outFileRoot", default= '', help= "This is the \
 partial name of the output file. Extension length is automatically included \
 in the file names.")

# load the inputs
(options,args) = parser.parse_args()

# reads the inputs from commandline
wiggle = options.wiggle
extend = options.extend
nRand = options.nRand
bedFile = options.bedFile
outFileRoot = options.outFileRoot

a = main( wiggle, extend, nRand, bedFile, outFileRoot )

# allWiggle="AH6408I-144-183-reps-MACS2/AH6408I-144-183-reps-SacCer3-2mis-PM-M5_MACS2_double_norm/AH6408I-144-183-reps-SacCer3-2mis-PM-M5_Dnorm_all.wig"
# bedFile = "chr_len.txt" """from SK1 Yue data 2017
# outFileRoot = "AH6408I-144-183"

