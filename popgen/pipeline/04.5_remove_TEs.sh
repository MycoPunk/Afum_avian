#!/usr/bin/bash
#SBATCH -p intel --mem 64gb -N 1 -n 4 --out logs/concat_vcf.log -p short

cd vcf

module load bedtools

#SNPS
bedtools subtract -a Avian_CIFAR_7.All.SNP.combined_selected.vcf.gz -b ../genome/FungiDB-39_AfumigatusAf293_Genome.RM.tab -header > Avian_CIFAR_7.All.SNP.combined_selected.NO.TE.vcf 
#INDELS
bedtools subtract -a Avian_CIFAR_7.All.INDEL.combined_selected.vcf.gz -b ../genome/FungiDB-39_AfumigatusAf293_Genome.RM.tab -header > Avian_CIFAR_7.All.INDEL.combined_selected.NO.TE.vcf 

module load tabix
#compress
bgzip Avian_CIFAR_7.All.SNP.combined_selected.NO.TE.vcf 
bgzip Avian_CIFAR_7.All.INDEL.combined_selected.NO.TE.vcf

#index:
tabix -p vcf Avian_CIFAR_7.All.SNP.combined_selected.NO.TE.vcf.gz
tabix -p vcf Avian_CIFAR_7.All.INDEL.combined_selected.NO.TE.vcf.gz
