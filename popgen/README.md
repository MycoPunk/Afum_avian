This is a pipeline templated on the Stajich Lab Popgen pipeline for A.fumigatus. The definitive version can be found at https://github.com/stajichlab/PopGenomics_Afumigatus_Global
The code contained in this repo captures the version used for analysis in the A. fumigatus avian project 2021. 


##

Files that need to be updated
 - samples.csv - this is a list of samples, the first line is a header and will be skipped by the job.
 - population_sets.yaml - all samples in samples.csv with optional population perams. 
 - config.txt  - there are several customization of the species run and data formats (we use CRAM as alignment storage default instead of BAM)
    * There is a PREFIX variable in config.txt which can be used to set the prefix of all output names, change this as you update your datasets so you can generate new results without overwriting old when generating combined VCF files.  Generally the scripts try to not re-run an analysis by checking if an output file already exists. So if you need to re-do an analysis remove files in the vcf, gvcf, or cram/aln folders.

# Steps
These are analysis pipeline steps intended to be run on slurm queueing system. The software is configured on the UCR HPCC system with UNIX modules. 

## Initialization
This step only needs to be run once.

* 00_index.sh - this build index files for alignment. It also will download (example of how to download automatically from FungiDB for those sources). In order for snpEff to work in later step you need to have downloaded the genome annotation in GFF3 format as well

## Alignment and gvcf creation (one per strain)
* 01_align.sh - this is bwa alignment step - this is intended to run one job per strain. i.e. ```sbatch --array=1-10 pipeline/01_align.sh```
* 02_call_gvcf.sh - this is done after alignments are finished - also one job per strain. 

## Genotyping GVCFs
* 03_jointGVCF_call_slice.sh - note: to change the `GVCF_INTERVAL` toa larger number, for example 8 chromosomes, run this as ```sbatch --array=1-8 pipeline/03_joint_GVCF_call_slices.sh```.
This step will run the GATK GenotypeGVCFs step followed by splitting variants into SNP and INDELs, followed by filtering steps to do HardFiltering based on QC cutoffs specified in the scripts 

* 04_combine_vcf.sh - this is a fast running script to gather the results from per-chrom/contig GVCF run into a single VCF file of filtered variants - the output files will be named by the PREFIX variable in the config.txt file.
* 04.5_remove_TEs.sh - this removes Transposable elements from known locations
* 06_make_SNP_tree.sh - this script will generate a FASTA format alignment of variants - using a parallelized approach (which requires the [GNU parallel tool](https://www.gnu.org/software/parallel/)). So a multithreaded slurm job is useful here - it will also launch a FastTree run to generate a phylogenetic tree from the alignment.
* 07_iqtree.sh - This script will generate a IQ-TREE run using the `GTR+ASC` model which is appropriate for SNP-based phylogenetic trees.
* 08_snpEff.sh - This script will generate a custom snpEff database using the GFF3 file in genome - this could be a little fragile so you may need to check that data files are download properly. This script right now has a lot A.fumigatus specific information.
