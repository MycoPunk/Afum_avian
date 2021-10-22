#!/usr/bin/bash
#SBATCH -p intel -N 1 -n 6 --mem 16gb  --out logs/IQTREE2.log

module load IQ-TREE/2.1.1

source config.txt
# this will fail but will generate file with fixed sites removed
cd strain_tree
iqtree2 -s $PREFIX.All.SNP.mfa -m GTR+ASC -nt AUTO -B 1000

#iqtree2 -s Baf_Strains3.SNP.combined_selected.mfa.varsites.phy -m GTR+ASC -nt AUTO -B 1000


