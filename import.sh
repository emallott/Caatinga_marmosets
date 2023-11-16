#!/bin/bash

#SBATCH -J obi_import
#SBATCH -A b1057
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elizabeth.mallott@northwestern.edu
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=24G
#SBATCH -t 24:00:00
#SBATCH --output=/home/ekm9460/obi_import_invert.out
#SBATCH --error=/home/ekm9460/obi_import_invert.err
#SBATCH -p b1057

module purge all

module load python

source activate /home/ekm9460/obi3-env

PROJECT="/projects/b1057/liz/tiputini_caatinga_diet"

mkdir ${PROJECT}/invert

cd ${PROJECT}/Art11l_ArtR17l

(
for R1 in *_R1.fastq
  do 
    name=`echo $R1 | sed 's/_R1.fastq//g'` 
    R2=${R1%_R1.fastq}_R2.fastq 
    obi import --fastq-input ./$R1 ${PROJECT}/invert/${name}/reads1 
    obi import --fastq-input ./$R2 ${PROJECT}/invert/${name}/reads2 
    obi import --ngsfilter ${PROJECT}/ngs_tiputini_caatinga_invert.txt ${PROJECT}/invert/${name}/ngsfile 
done
)
