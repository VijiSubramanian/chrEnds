#!/bin/bash
##SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1                                                                                                          
#SBATCH --cpus-per-task=2
#SBATCH --mem=40GB
#SBATCH --time=40:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=vvs2@nyu.edu
#SBATCH --job-name=boot_chrEnds
#SBATCH --output=/scratch/%u/%x_%j.out

module purge
module load python/intel/2.7.12

# I defined $FOLDER outside of the job, so I could do multiple folders in parallel
# i.e., sbatch --export FOLDER="XX" ~/Bootstrap/Bootstrap_chrEnds110kb.sbatch
cd /scratch/vvs2/endsEnrichment20-110kb/$FOLDER

for infile in echo `ls *.bdg*`
do
     echo "working on: $infile"

     RE="(.+)\.gz"
     if [[ $infile =~ $RE ]]; then
     	  file=${BASH_REMATCH[1]}
       	  gzip -d $infile
     else
	  file=$infile
     fi

     outname=$(echo $file | cut -d "_" -f 1)

	  python Bootstrap_chromosomeEnds.py -w $file -e 110000 -n 5000 -b "/home/vvs2/Library/SK1YueChrLen.txt" -o $outname
     

gzip $file
done

