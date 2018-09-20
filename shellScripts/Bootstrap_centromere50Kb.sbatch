#!/bin/bash
##SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --mem=40GB#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20GB
#SBATCH --time=40:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=vvs2@nyu.edu
#SBATCH --job-name=bootCen50kb
#SBATCH --output=/scratch/%u/%x_%j.out

module purge
module load python/intel/2.7.12

# I defined $FOLDER outside of the job, so I could do multiple folders in parallel
# i.e., sbatch --export FOLDER="XX" ~/Bootstrap_centromere50Kb.sbatch
cd /scratch/vvs2/cenEnrichment50kb/$FOLDER

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

     #for ext in 2000 50000 
     #do
#	  echo $ext
	  python Bootstrap_centromere.py -w $infile -e 50000 -n 5000 -b "/home/vvs2/Library/SK1YueCen.gff" -o $outname
 #    done
#	   python Bootstrap_centromere.py -w $infile -e 50000 -n 5000 -b "/home/vvs2/Library/SK1YueCen.gff" -o $outname
gzip $infile
done




