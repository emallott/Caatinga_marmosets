#!/bin/bash

#SBATCH -J obi_taxa
#SBATCH -A b1057
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elizabeth.mallott@vanderbilt.edu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH -t 720:00:00
#SBATCH --output=/home/ekm9460/obi_taxa_invert_small.out
#SBATCH --error=/home/ekm9460/obi_taxa_invert_small.err
#SBATCH -p b1057

module purge all

module load python

source activate /home/ekm9460/obi3-env

PROJECT="/projects/b1057/liz/tiputini_caatinga_diet"

obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy ${PROJECT}/taxonomy_invert_small/taxonomy/ncbi_tax ${PROJECT}/taxonomy_invert_small/081720_coi_refs ${PROJECT}/taxonomy_invert_small/081720_coi_refs_clean

obi clean_dms ${PROJECT}/taxonomy_invert_small/

obi uniq --taxonomy ${PROJECT}/taxonomy_invert_small/taxonomy/ncbi_tax ${PROJECT}/taxonomy_invert_small/081720_coi_refs_clean ${PROJECT}/taxonomy_invert_small/081720_coi_refs_uniq

obi grep --require-rank=family --taxonomy ${PROJECT}/taxonomy_invert_small/taxonomy/ncbi_tax ${PROJECT}/taxonomy_invert_small/081720_coi_refs_uniq ${PROJECT}/taxonomy_invert_small/081720_coi_refs_uniq_clean

obi build_ref_db -t 0.80 --taxonomy ${PROJECT}/taxonomy_invert_small/taxonomy/ncbi_tax ${PROJECT}/taxonomy_invert_small/081720_coi_refs_uniq_clean ${PROJECT}/taxonomy_invert_small/081720_coi_refs_db_80
