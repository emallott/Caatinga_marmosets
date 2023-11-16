#!/bin/bash

#SBATCH -J obi_clean
#SBATCH -A b1057
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elizabeth.mallott@northwestern.edu
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --mem=24G
#SBATCH -t 24:00:00
#SBATCH --output=/home/ekm9460/obi_clean_invert.out
#SBATCH --error=/home/ekm9460/obi_clean_invert.err
#SBATCH -p b1057

module purge all

module load python

source activate /home/ekm9460/obi3-env

PROJECT="/projects/b1057/liz/tiputini_caatinga_diet"

#(
#for identified in ${PROJECT}/invert/*
#  do
#     obi annotate $identified/identified_sequences -S sample:$identified $identified/#tagged_identified_sequences
#done
#)

(
for sample_dms in ${PROJECT}/invert/*.obidms
  do obi clean_dms $sample_dms
done
)

catcmd='obi cat ' 

for sample_dms in ${PROJECT}/invert/*.obidms ; do catcmd="$catcmd-c $sample_dms/tagged_identified_sequences " ; done 

catcmd2="$catcmd ${PROJECT}/combined_invert_samples_dms/combined_samples" 
echo $catcmd2
$catcmd2

obi uniq -m sample ${PROJECT}/combined_invert_samples_dms/combined_samples ${PROJECT}/combined_invert_samples_dms/dereplicated_sequences

obi annotate -k COUNT -k MERGED_sample ${PROJECT}/combined_invert_samples_dms/dereplicated_sequences ${PROJECT}/combined_invert_samples_dms/cleaned_metadata_sequences

obi grep -p "len(sequence)>=50 and sequence['COUNT']>=10" ${PROJECT}/combined_invert_samples_dms/cleaned_metadata_sequences ${PROJECT}/combined_invert_samples_dms/denoised_sequences

obi clean -p 6 -s MERGED_sample -r 0.05 -H ${PROJECT}/combined_invert_samples_dms/denoised_sequences ${PROJECT}/combined_invert_samples_dms/cleaned_sequences
