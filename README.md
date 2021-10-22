This repo contains the code used for analysis in the A. fumigatus avian project 2021. 

# popgen 
The popgen folder contains code for genotyping and tree building templated on the Stajich Lab Popgen pipeline for A. fumigatus. The code contained in this folder captures the version used for this analysis. The definitive version of the pipeline can be found at https://github.com/stajichlab/PopGenomics_Afumigatus_Global
See the README.md file in the popgen folder for a detailed description of each step

# scripts_r
The scripts_r folder contains the R code used to analyze variant data and generate the figures for the manuscript 

# data
The data folder contains the input data needed to generate the figures used in scripts_r, minus the variant files which are too large to put here and must be generated first using the popgen pipeline
