##This script looks at the distribution of variants between avaian strains and clinical strains of A. fumigatus, and searches for azole resistance variants in avian isoaltes.
#last updated 20.Oct.2021

#load modules
require(data.table)
require(splitstackshape)
library(plyr)
library(dplyr)
library(tidyverse)
library(rlist)
library(gdata)
library(ggtree)
library(ape)
library(phytools)
library(ggplot2)
library(stringr)

#set wd
setwd("~/../data")
options(stringsAsFactors = FALSE)

##set seed for reproducibility
set.seed(666)

##read in data
snpEff_all<-read.delim("Avian_CIFAR_7.snpEff.matrix.tsv", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = FALSE)
resistance_db<-read.delim("Afum_azole_mutations.csv", header = TRUE, sep = ",", fill = TRUE, strip.white = TRUE)

###clean the data for easy processing###
#remove the brackets form ALT
snpEff_all$ALT<-gsub("\\[|\\]", "", snpEff_all$ALT)

#remove the "/" from all the the SNP calls 
snpEff_all[] <- lapply(snpEff_all, gsub, pattern='/', replacement="")

#get set size
dim(snpEff_all[,11:(ncol(snpEff_all) -1)])
#181 strains

#subset to only missense variants 
snpEff_no_intergenic<- snpEff_all[snpEff_all$TYPE == "missense_variant",]
dim(snpEff_no_intergenic)
#76125 positions with missence variants

#in how many genes? 
length(unique(snpEff_no_intergenic$GENE)) #9306

#subset to just Avian strains 
#where AF293 is negative control for Cyp51A mutation and 08-19-02-10 is the positive control for Cyp51A mutation
snpEff_no_intergenic_avian<- snpEff_no_intergenic[grep("USGS|ICMP|AF293", names(snpEff_no_intergenic))]
snpEff_no_intergenic_avian_meta<- cbind(snpEff_no_intergenic[1:10], snpEff_no_intergenic_avian)
names(snpEff_no_intergenic_avian)


#remove rows where all avian isolates match the reference 
#calculate alternative ALT call
unique_calls_per_row<- t(apply(snpEff_no_intergenic_avian, 1, function(x) unique(x)))
new_ALT<- lapply(unique_calls_per_row, function(x) paste(x))
new_ALT2<- lapply(new_ALT, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df <- data.frame(matrix(unlist(new_ALT2), nrow=length(new_ALT2), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df)<- "NEW_ALT_merged"
new_ALT2_sep<- lapply(new_ALT, function(x) paste(x, sep = " ", collapse = ","))
new_ALT2_df_sep <- data.frame(matrix(unlist(new_ALT2_sep), nrow=length(new_ALT2_sep), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_sep)<- "NEW_ALT"

#add back metadata
new_df<- cbind(snpEff_no_intergenic_avian, new_ALT2_df)
new_df_w_ref<- cbind("NEW_ALT" = new_ALT2_df_sep, 
                     "REF" = snpEff_no_intergenic_avian_meta$REF, 
                     "TYPE" = snpEff_no_intergenic_avian_meta$TYPE, 
                     "GENE" = snpEff_no_intergenic_avian_meta$GENE, 
                     "CHANGEDNA" = snpEff_no_intergenic_avian_meta$CHANGEDNA, 
                     "CHANGEPEP" = snpEff_no_intergenic_avian_meta$CHANGEPEP, 
                     "IMPACT" = snpEff_no_intergenic_avian_meta$IMPACT, 
                     "POS" = snpEff_no_intergenic_avian_meta$POS, 
                     "CHROM" = snpEff_no_intergenic_avian_meta$CHROM, 
                     new_df)

#remove those were new alt matches the ref
subset_unique <- function(a,b) FALSE %in%(unlist(strsplit(a,",")) %in% unlist(strsplit(b,",")))
only_in_avian<- new_df_w_ref[apply(new_df_w_ref[,c('NEW_ALT','REF')], 1, function(y) subset_unique(y['NEW_ALT'],y['REF'])),]

#check that this worked
nrow(new_df_w_ref)
nrow(only_in_avian)
all_the_same_so_removed<- new_df_w_ref[!rownames(new_df_w_ref) %in% rownames(only_in_avian),]
nrow(all_the_same_so_removed) #number match up
sum(new_df_w_ref$NEW_ALT == new_df_w_ref$REF) #worked

#across how many genes 
length(unique(only_in_avian$GENE))

##so there are 53187 missence variants (including misscalls/deletions) in the Avian set
#53841 if we count our non-avian 08-19-02-10 Cyp51A positive control

#Are there any variants that were found only in the avian set and NOT in the clinical set?
snpEff_no_intergenic_clinical<- snpEff_no_intergenic[grep("USGS|ICMP|AF293|NEW_ALT|FLANKING|POS|REF|TYPE|GENE|CHANGEDNA|CHANGEPEP|IMPACT|CHROM|ALT|ANN", names(snpEff_no_intergenic), invert = TRUE)]

#any vars unique to avian set?  - do the same for clinical set
#remove rows where all avian isolates match the reference 
#calculate alternative ALT call
unique_calls_per_row_clin<- t(apply(snpEff_no_intergenic_clinical, 1, function(x) unique(x)))
new_ALT_clin<- lapply(unique_calls_per_row_clin, function(x) paste(x))
new_ALT2_clin<- lapply(new_ALT_clin, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df_clin <- data.frame(matrix(unlist(new_ALT2_clin), nrow=length(new_ALT2_clin), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_clin)<- "NEW_ALT_merged_clin"
new_ALT2_sep_clin<- lapply(new_ALT_clin, function(x) paste(x, sep = " ", collapse = ","))
new_ALT2_df_sep_clin <- data.frame(matrix(unlist(new_ALT2_sep_clin), nrow=length(new_ALT2_sep_clin), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_sep_clin)<- "NEW_ALT_CLIN"

#add back metadata
new_df_clin<- cbind(snpEff_no_intergenic_clinical, new_ALT2_df_clin)
new_df_w_ref_clin<- cbind("NEW_ALT_CLIN" = new_ALT2_df_sep_clin, 
                     "REF" = snpEff_no_intergenic_avian_meta$REF, 
                     "TYPE" = snpEff_no_intergenic_avian_meta$TYPE, 
                     "GENE" = snpEff_no_intergenic_avian_meta$GENE, 
                     "CHANGEDNA" = snpEff_no_intergenic_avian_meta$CHANGEDNA, 
                     "CHANGEPEP" = snpEff_no_intergenic_avian_meta$CHANGEPEP, 
                     "IMPACT" = snpEff_no_intergenic_avian_meta$IMPACT, 
                     "POS" = snpEff_no_intergenic_avian_meta$POS, 
                     "CHROM" = snpEff_no_intergenic_avian_meta$CHROM, 
                     new_df_clin)

#remove those were new alt matches the ref
subset_unique <- function(a,b) FALSE %in%(unlist(strsplit(a,",")) %in% unlist(strsplit(b,",")))
only_in_clinical<- new_df_w_ref_clin[apply(new_df_w_ref_clin[,c('NEW_ALT_CLIN','REF')], 1, function(y) subset_unique(y['NEW_ALT_CLIN'],y['REF'])),]

length(unique(only_in_clinical$GENE))
length(unique(only_in_avian$GENE))

length(setdiff(only_in_clinical$GENE, only_in_avian$GENE)) #there are 444 genes with variants in the clinical set that don't have variants in the set 
genes_clin<- setdiff(only_in_clinical$GENE, only_in_avian$GENE)
genes_mut_in_clin_not_avian<- only_in_clinical[only_in_clinical$GENE %in% genes_clin,]
length(unique(genes_mut_in_clin_not_avian$GENE))

length(setdiff(only_in_avian$GENE, only_in_clinical$GENE)) #there are 91 genes in the avian set that don't have variants in the clinical set
genes_avian<- setdiff(only_in_avian$GENE, only_in_clinical$GENE)
genes_mut_in_avian_not_clin<- only_in_avian[only_in_avian$GENE %in% genes_avian,]
length(unique(genes_mut_in_avian_not_clin$GENE))


#for these uniquely mutated genes, are there variants in all of the strains for either set?
counts_clin<- data.frame(names=genes_mut_in_clin_not_avian$NEW_ALT_CLIN,chr=apply(genes_mut_in_clin_not_avian,2,nchar)[,1])
#View(counts_clin) none are in all. 
counts_avian<- data.frame(names=genes_mut_in_avian_not_clin$NEW_ALT,chr=apply(genes_mut_in_avian_not_clin,2,nchar)[,1])
#View(counts_avian)none are in all. 


###sep each gene of interest in the resistance_db, sep snpEff_no_intergenic into their own dfs
#set genes of interest
GOI<- unique(resistance_db$A_fum_gene_name)
length(GOI) #13 genes under consideration. 

#subset df into just genes of interest
snpEff_no_intergenic_GOI <- only_in_avian[only_in_avian$GENE %in% GOI, ]
#how many genes have variants in them?
genes_w_vars<- unique(snpEff_no_intergenic_GOI$GENE)
length(genes_w_vars)
#all 13 of them- interesting. 
#note Afu4g06890 is Cyp51A


#split each gene into it's own df (keep in a list)
listDf_GOI_variants <- split(snpEff_no_intergenic_GOI, f = snpEff_no_intergenic_GOI$GENE)

##split resistance_db by gene
#subset to only genes with variants in the dataset
resistance_db_GOI <- resistance_db[resistance_db$A_fum_gene_name %in% genes_w_vars, ]
#check
unique(resistance_db$A_fum_gene_name)
unique(resistance_db_GOI$A_fum_gene_name) #works - but here, were not actually removing any since there's variants in all 13 GOI

#split each database entry by gene into it's own df (keep in a list)
listDf_resistance_db_GOI <- split(resistance_db_GOI, f = resistance_db_GOI$A_fum_gene_name)

#function to subset all variants in genes of interest to only known resistance mutations
get_resistance_mutations <- function(variants_OI, resistance_db_OI){
  output<- list()
  #for each gene, subset the rows if variants are known resistance mutations 
  for (i in 1:length(variants_OI)){
    cyp51_db<- resistance_db_OI[[i]]
    cyp51_var<- variants_OI[[i]]
    output[[i]] <- cyp51_var[grep(paste(cyp51_db$mutation_code2,collapse="|"), cyp51_var$CHANGEPEP), ]
  }
  names(output)<- names(resistance_db_OI)
  no_empties<- output[sapply(output, function(x) nrow(x)) > 0]
  return(no_empties)
}


#alternate function to find variable mutation at the same position as known resistance variants
get_resistance_mutations_variable <- function(variants_OI, resistance_db_OI){
  output<- list()
  #for each gene, subset the rows if variants are known resistance mutations 
  for (i in 1:length(variants_OI)){
    cyp51_db<- resistance_db_OI[[i]]
    cyp51_var<- variants_OI[[i]]
    output[[i]] <- cyp51_var[grep(paste(cyp51_db$mutation_code2_start_only,collapse="|"), cyp51_var$CHANGEPEP), ]
  }
  names(output)<- names(resistance_db_OI)
  no_empties<- output[sapply(output, function(x) nrow(x)) > 0]
  return(no_empties)
}

#run function
list_of_dfs_resistance_vars_in_dataset<- get_resistance_mutations(variants_OI = listDf_GOI_variants, resistance_db_OI = listDf_resistance_db_GOI)
names(list_of_dfs_resistance_vars_in_dataset)
#Here we only have mutations in Cyp51A
#isolate Cyp51A
Cyp51a<- list_of_dfs_resistance_vars_in_dataset[[1]]
#View(Cyp51a) #canonnical Cyp51 mutation only in our control. - NOTE this will yield an empty if you run w/o the control

#How about if it's a variable change at the same position?
#list_of_dfs_resistance_vars_in_dataset_variable<- get_resistance_mutations_variable(variants_OI = listDf_GOI_variants, resistance_db_OI = listDf_resistance_db_GOI)
#names(list_of_dfs_resistance_vars_in_dataset_variable)
#still only Cyp51A
#isolate this too
#Cyp51a_variable<- list_of_dfs_resistance_vars_in_dataset_variable[[1]]

#check if they're the same
#dim(Cyp51a)
#dim(Cyp51a_variable)
#nope, no new aa changes at positions of known resistance variants.


##make presence/absence dataframe for presence of variants in each gene and map onto tree file. 

#Function to get the number of variants in ea gene for each strain
#over all: for each strain, for each gene, are the calls real variants (they don't match the REF?) If so, TRUE, if no FALSE. 
#n_variants_per_gene_per_strain<- function(gene_df){
#  #shrink the df 
#  REF<- data.frame(gene_df$REF)
#  strains<- gene_df[, 10:(ncol(gene_df) -1)]
#  #create true/false data frame: true if not == reference, and not == "." (possible misscalls)
#  T_F_df<- data.frame(apply(strains, 2, function(x) (x != REF$gene_df.REF) & (x != ".")))
#  if (nrow(strains) >1) {
#    T_F_df<- data.frame(apply(strains, 2, function(x) (x != REF$gene_df.REF) & (x != ".")))
#  } else {
#    T_F_df<- t(data.frame(apply(strains, 2, function(x) (x != REF$gene_df.REF) & (x != "."))))
#  }
#  #get col. sums (the number of TRUE values == the number of variants in that gene for that strain)
#  n_vars_in_gene<- data.frame(colSums(T_F_df))
#  #rename for later
#  df.name<- deparse(substitute(gene_df))
#  colnames(n_vars_in_gene)<- df.name 
#  return(n_vars_in_gene)
#}

#dim(n_variants_per_gene_per_strain(Cyp51a))

#use the above function in it's own function to loop over each mutation to make a df for graphing and get totals across all strains.
#make_by_variant_output <- function(input){
#  output<- data.frame(matrix(NA, nrow = ncol(input[, 10:(ncol(input) -1)]), ncol = 1))
#  #for each row (variant) get presence/ absence for which strains have the variant
#  for (i in 1:nrow(input)){
#    output[[i]]<- n_variants_per_gene_per_strain(input[i,])
#    this_name<- input[i,6]
#    names(output)[[i]]<- this_name
#  }
#  return(output)
#}


#run function
#by_variant_output_df<- make_by_variant_output(input = Cyp51a)
#colSums(by_variant_output_df)

#do the same for all non-syn changes in Cyp51A (not just known resistance aleles)
Cyp51a_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu4g06890",]
Cyp51a_all_but_vir<- Cyp51a_all
#remove the known virulence positions identified above 
#Cyp51a_all_but_vir<- Cyp51a_all[!Cyp51a_all$CHANGEPEP %in% names(by_variant_output_df),]

#make presence/absence df by strain
by_variant_output_df_all_but_vir<- make_by_variant_output(input = Cyp51a_all_but_vir)
colSums(by_variant_output_df_all_but_vir)
colnames(by_variant_output_df_all_but_vir)<- Cyp51a_all_but_vir$CHANGEPEP


#write.table(by_variant_output_df, "variants_in_resistance_aleles.csv", sep=",", row.names = TRUE, col.names = TRUE, quote = FALSE)
#write.table(by_variant_output_df_all_but_vir, "all_nonsyn_variants_in_cyp51a.csv", sep=",", row.names = TRUE, col.names = TRUE, quote = FALSE)


##tree- set colors for groups, fix tip labels etc. 
#rename tips
tree <- read.tree("Avian_iq_tree.tre")
#read in metadata
metadata<-read.delim("Avian_metadata.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

#remove the reference from the tree:
tree<- drop.tip(tree, "Af293-REF", trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0, rooted = is.rooted(tree), collapse.singles = TRUE,
                   interactive = FALSE)

name_map<-read.delim("Avian_name_map.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)


#rename using tree.io
tree_me2<- treeio::rename_taxa(tree, name_map, OG_tree_name, New_tree_name) #%>% write.tree
#clade1_names_fixed <- name_map$name_pop_genome[match(names(clade1_names), name_map$name_Pan_genome)]

#set groups
avian_names<-tree_me2$tip.label[grep(pattern = 'USGS', x = tree_me2$tip.label)]
avian_names_icmp<-tree_me2$tip.label[grep(pattern = 'ICMP', x = tree_me2$tip.label)]
avian_names_all<-tree_me2$tip.label[grep(pattern = 'USGS|ICMP|AF293', x = tree_me2$tip.label)]
outgroup_names<- tree_me2$tip.label[!tree_me2$tip.label %in% avian_names_all]

#subset to only avian strains
tree_me_small <- keep.tip(tree_me2, avian_names_all)

#split data by host ecology
grA_me_small<- split(metadata$STRAIN_NAME, metadata$HOST_ECOLOGY)

#set colors by group
tree_grA_me_small <- ggtree::groupOTU(tree_me_small, grA_me_small)
str(tree_grA_me_small)
levels(attributes(tree_grA_me_small)$group) 
levels(attributes(tree_grA_me_small)$group)[1] <- "raptor"
levels(attributes(tree_grA_me_small)$group) 

attributes(tree_grA_me_small)$group <- factor(x = attributes(tree_grA_me_small)$group, 
                                              levels = c("raptor","sea bird","water bird"))

my_cols_me_small <- c(raptor = "#58788C",
                      sea_bird ="#2A4359",
                      water_bird = "#F2811D")
names(my_cols_me_small) <- levels(attributes(tree_grA_me_small)$group)
scales::show_col(my_cols_me_small); my_cols_me_small

#dark blue: #2A4359 (Kakapoo)
#light blue: #58788C
#yellow: #F2BF27
#orange: #F2811D
#dark orange: #F2490C (Avian)
#almost black: #0D0D0D (Outgroups)

#plot test tree
#plot tree
tree_plot_me_sm <- 
  ggtree(tr = tree_grA_me_small, 
         # color by group attribute, check str(tree_grA_me_small)
         mapping = aes(color = group), 
         layout  = 'rectangular', 
         #branch.length = 'none', 
         #  geom_treescale(x=3, y=NULL, color = "white") +
         # set line thickness
         size = .3, show.legend=FALSE) +
  # adjust coloring of main groups
  scale_color_manual(name = 'Group', values = my_cols_me_small) +
  #geom_text(show.legend = FALSE)
  theme(legend.title=element_text(size=0), # The title of legend 
        legend.text=element_text(size=0))


# plot and add the tip labels
tre_to_plot<- tree_plot_me_sm + geom_tiplab(size = 1, align = TRUE, linesize = .1, linetype = 3, show.legend=FALSE)
tre_to_plot


####


#format resistance allele data for graphing
#colapse df of df because you're an idiot and didn't append these as a list in the first place
#by_variant_output_df_colapse<- bind_rows(lapply(by_variant_output_df,function(i)do.call(cbind,i)))
#row.names(by_variant_output_df_colapse)<- row.names(by_variant_output_df[[1]])

#same for unknown function vars
by_variant_output_df_colapse_unknown<- bind_rows(lapply(by_variant_output_df_all_but_vir,function(i)do.call(cbind,i)))
row.names(by_variant_output_df_colapse_unknown)<- row.names(by_variant_output_df_all_but_vir[[1]])


#order columns by abundance
#colSums(by_variant_output_df_colapse)
#by_abundance<- by_variant_output_df_colapse[,order(colSums(-by_variant_output_df_colapse))]
#colSums(by_abundance) #works
#get total genomes influenced 
#list<- colSums(by_abundance)
#sum(list)

#same for unknowns
colSums(by_variant_output_df_colapse_unknown)
by_abundance_unknown<- by_variant_output_df_colapse_unknown[,order(colSums(-by_variant_output_df_colapse_unknown))]
colSums(by_abundance_unknown) #works

#set presence/absence
#by_variant_output_df_binary_presence<- data.frame(sapply(by_abundance, gsub, pattern = "1", replacement = "present"))
#by_variant_output_df_binary_presence_absence<- data.frame(sapply(by_variant_output_df_binary_presence, gsub, pattern = "0", replacement = "absent"))
#row.names(by_variant_output_df_binary_presence_absence)<- row.names(by_variant_output_df[[1]])

#set presence/absence for unknowns
by_variant_output_df_binary_presence_unknowns<- data.frame(sapply(by_abundance_unknown, gsub, pattern = "1", replacement = "present"))
by_variant_output_df_binary_presence_absence_unknowns<- data.frame(sapply(by_variant_output_df_binary_presence_unknowns, gsub, pattern = "0", replacement = "absent"))
row.names(by_variant_output_df_binary_presence_absence_unknowns)<- row.names(by_variant_output_df_colapse_unknown)


#check that all the names match and you're not missing any data
tree_node_names<- tree_grA_me_small$tip.label
df_names<- rownames(by_variant_output_df_binary_presence_absence_unknowns)
setdiff(tree_node_names, df_names)
#none missing


#plot resistance variants onto tree
#plot

#bind and highlight - can't get width right when adjacent. 
#all_vars<- cbind(by_variant_output_df_binary_presence_absence, by_variant_output_df_binary_presence_absence_unknowns)
cyp51A_tree_plot <-  gheatmap(tre_to_plot, 
                              by_variant_output_df_binary_presence_absence_unknowns,
                              #all_vars,
                              offset=0.14, width=0.2, 
                              low="white", high="black", 
                              colnames = T, 
                              colnames_angle = 45,
                              colnames_position = "top",
                              font.size = 1,
                              colnames_offset_y = -.5,
                              colnames_offset_x = .0)
                              #color="white") +
  #scale_fill_manual(values=c("white", "black")) 
#+ggtitle("Cyp51 variants")

#p<- cyp51A_tree_plot + geom_tiplab(size = .8, align = TRUE, linesize = .25, offset = 1, linetype = 0)
#p<- cyp51A_tree_plot + scale_fill_manual(values=c("white", "black")) 

#export
#ggsave(file="Avian_Cyp51A_variants.pdf",device="pdf", p, width=8, height=8, units="in")


#make supplemental figure of all missence variants from all strains (includng those for which only expression data is availble)
#remove cyp51a results
genes_w_vars



#isolate the four other genes with non-syn changes in resistance genes:


"Afu1g12690" #mdr4
"Afu1g14330" #cdr1B_abcC
"Afu2g03700" #hmg1
"Afu2g14720" #hapE
"Afu3g03500" #mdr3
"Afu4g06890" #Cyp51A
"Afu4g08340" #cox10
"Afu4g10000" #mdr2
"Afu5g06070" #mdr1
"Afu6g04360" #atrF
"Afu6g12400" #fks1
"Afu7g01960" #7G01960
"Afu7g03740" #cyp51B

mdr4_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu1g12690",]
nrow(mdr4_all)# 9 positions

cdr1B_abcC_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu1g14330",]
nrow(cdr1B_abcC_all)# 3 positions

hmg1_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu2g03700",]
nrow(hmg1_all)# 5 positions

hapE_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu2g14720",]
nrow(hapE_all)# 1 position

mdr3_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu3g03500",]
nrow(mdr3_all)# 13 positions

cox10_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu4g08340",]
nrow(cox10_all)# 5 positions

mdr2_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu4g10000",]
nrow(mdr2_all)# 6 positions

mdr1_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu5g06070",]
nrow(mdr1_all)# 4 positions

atrF_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu6g04360",]
nrow(atrF_all)# 8 positions

fks1_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu6g12400",]
nrow(fks1_all)# 4 positions

cyp51B_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu7g03740",]
nrow(cyp51B_all)# 2 positions

AFUA_7G01960_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu7g01960",]
nrow(AFUA_7G01960_all)# 3 positions



#make presence/absence df by strain
mdr4_all_PA<- make_by_variant_output(input = mdr4_all)
colSums(mdr4_all_PA)
cdr1B_abcC_all_PA<- make_by_variant_output(input = cdr1B_abcC_all)
colSums(cdr1B_abcC_all_PA)
hmg1_all_PA<- make_by_variant_output(input = hmg1_all)
colSums(hmg1_all_PA)
hapE_all_PA<- make_by_variant_output(input = hapE_all)
colSums(hapE_all_PA)
mdr3_all_PA<- make_by_variant_output(input = mdr3_all)
colSums(mdr3_all_PA)
cox10_all_PA<- make_by_variant_output(input = cox10_all)
colSums(cox10_all_PA)
mdr2_all_PA<- make_by_variant_output(input = mdr2_all)
colSums(mdr2_all_PA)
mdr1_all_PA<- make_by_variant_output(input = mdr1_all)
colSums(mdr1_all_PA)
atrF_all_PA<- make_by_variant_output(input = atrF_all)
colSums(atrF_all_PA)
fks1_all_PA<- make_by_variant_output(input = fks1_all)
colSums(fks1_all_PA)
cyp51B_all_PA<- make_by_variant_output(input = cyp51B_all)
colSums(cyp51B_all_PA)
AFUA_7G01960_all_PA<- make_by_variant_output(input = AFUA_7G01960_all)
colSums(AFUA_7G01960_all_PA)


#format resistance allele data for graphing
#colapse df of df because you're an idiot and didn't append these as a list in the first place
mdr4_colapse<- bind_rows(lapply(mdr4_all_PA,function(i)do.call(cbind,i)))
row.names(mdr4_colapse)<- row.names(mdr4_all_PA[[1]])

cdr1B_abcC_colapse<- bind_rows(lapply(cdr1B_abcC_all_PA,function(i)do.call(cbind,i)))
row.names(cdr1B_abcC_colapse)<- row.names(cdr1B_abcC_all_PA[[1]])

hmg1_colapse<- bind_rows(lapply(hmg1_all_PA,function(i)do.call(cbind,i)))
row.names(hmg1_colapse)<- row.names(hmg1_all_PA[[1]])

hapE_colapse<- bind_rows(lapply(hapE_all_PA,function(i)do.call(cbind,i)))
row.names(hapE_colapse)<- row.names(hapE_all_PA[[1]])

mdr3_colapse<- bind_rows(lapply(mdr3_all_PA,function(i)do.call(cbind,i)))
row.names(mdr3_colapse)<- row.names(mdr3_all_PA[[1]])

cox10_colapse<- bind_rows(lapply(cox10_all_PA,function(i)do.call(cbind,i)))
row.names(cox10_colapse)<- row.names(cox10_all_PA[[1]])

mdr2_colapse<- bind_rows(lapply(mdr2_all_PA,function(i)do.call(cbind,i)))
row.names(mdr2_colapse)<- row.names(mdr2_all_PA[[1]])

mdr1_colapse<- bind_rows(lapply(mdr1_all_PA,function(i)do.call(cbind,i)))
row.names(mdr1_colapse)<- row.names(mdr1_all_PA[[1]])

atrF_colapse<- bind_rows(lapply(atrF_all_PA,function(i)do.call(cbind,i)))
row.names(atrF_colapse)<- row.names(atrF_all_PA[[1]])

fks1_colapse<- bind_rows(lapply(fks1_all_PA,function(i)do.call(cbind,i)))
row.names(fks1_colapse)<- row.names(fks1_all_PA[[1]])

cyp51B_colapse<- bind_rows(lapply(cyp51B_all_PA,function(i)do.call(cbind,i)))
row.names(cyp51B_colapse)<- row.names(cyp51B_all_PA[[1]])

AFUA_7G01960_colapse<- bind_rows(lapply(AFUA_7G01960_all_PA,function(i)do.call(cbind,i)))
row.names(AFUA_7G01960_colapse)<- row.names(AFUA_7G01960_all_PA[[1]])


#order columns by abundance
mdr4_by_abundance<- mdr4_colapse[,order(colSums(-mdr4_colapse))]
cdr1B_abcC_by_abundance<- cdr1B_abcC_colapse[,order(colSums(-cdr1B_abcC_colapse))]
hmg1_by_abundance<- hmg1_colapse[,order(colSums(-hmg1_colapse))]
hapE_by_abundance<- hapE_colapse[,order(colSums(-hapE_colapse))]
mdr3_by_abundance<- mdr3_colapse[,order(colSums(-mdr3_colapse))]
cox10_by_abundance<- cox10_colapse[,order(colSums(-cox10_colapse))]
mdr2_by_abundance<- mdr2_colapse[,order(colSums(-mdr2_colapse))]
mdr1_by_abundance<- mdr1_colapse[,order(colSums(-mdr1_colapse))]
atrF_by_abundance<- atrF_colapse[,order(colSums(-atrF_colapse))]
fks1_by_abundance<- fks1_colapse[,order(colSums(-fks1_colapse))]
cyp51B_by_abundance<- cyp51B_colapse[,order(colSums(-cyp51B_colapse))]
AFUA_7G01960_by_abundance<- AFUA_7G01960_colapse[,order(colSums(-AFUA_7G01960_colapse))]

#set presence/absence
mdr4_by_abundance_presence<- data.frame(sapply(mdr4_by_abundance, gsub, pattern = "1", replacement = "present"))
mdr4_by_abundance_presence_absence<- data.frame(sapply(mdr4_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(mdr4_by_abundance_presence_absence)<- row.names(mdr4_colapse)

cdr1B_abcC_by_abundance_presence<- data.frame(sapply(cdr1B_abcC_by_abundance, gsub, pattern = "1", replacement = "present"))
cdr1B_abcC_by_abundance_presence_absence<- data.frame(sapply(cdr1B_abcC_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(cdr1B_abcC_by_abundance_presence_absence)<- row.names(cdr1B_abcC_colapse)

hmg1_by_abundance_presence<- data.frame(sapply(hmg1_by_abundance, gsub, pattern = "1", replacement = "present"))
hmg1_by_abundance_presence_absence<- data.frame(sapply(hmg1_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(hmg1_by_abundance_presence_absence)<- row.names(hmg1_colapse)

hapE_by_abundance_presence<- data.frame(sapply(hapE_by_abundance, gsub, pattern = "1", replacement = "present"))
hapE_by_abundance_presence_absence<- data.frame(sapply(hapE_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(hapE_by_abundance_presence_absence)<- row.names(hapE_colapse)

mdr3_by_abundance_presence<- data.frame(sapply(mdr3_by_abundance, gsub, pattern = "1", replacement = "present"))
mdr3_by_abundance_presence_absence<- data.frame(sapply(mdr3_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(mdr3_by_abundance_presence_absence)<- row.names(mdr3_colapse)

cox10_by_abundance_presence<- data.frame(sapply(cox10_by_abundance, gsub, pattern = "1", replacement = "present"))
cox10_by_abundance_presence_absence<- data.frame(sapply(cox10_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(cox10_by_abundance_presence_absence)<- row.names(cox10_colapse)

mdr2_by_abundance_presence<- data.frame(sapply(mdr2_by_abundance, gsub, pattern = "1", replacement = "present"))
mdr2_by_abundance_presence_absence<- data.frame(sapply(mdr2_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(mdr2_by_abundance_presence_absence)<- row.names(mdr2_colapse)

mdr1_by_abundance_presence<- data.frame(sapply(mdr1_by_abundance, gsub, pattern = "1", replacement = "present"))
mdr1_by_abundance_presence_absence<- data.frame(sapply(mdr1_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(mdr1_by_abundance_presence_absence)<- row.names(mdr1_colapse)

atrF_by_abundance_presence<- data.frame(sapply(atrF_by_abundance, gsub, pattern = "1", replacement = "present"))
atrF_by_abundance_presence_absence<- data.frame(sapply(atrF_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(atrF_by_abundance_presence_absence)<- row.names(atrF_colapse)

fks1_by_abundance_presence<- data.frame(sapply(fks1_by_abundance, gsub, pattern = "1", replacement = "present"))
fks1_by_abundance_presence_absence<- data.frame(sapply(fks1_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(fks1_by_abundance_presence_absence)<- row.names(fks1_colapse)

cyp51B_by_abundance_presence<- data.frame(sapply(cyp51B_by_abundance, gsub, pattern = "1", replacement = "present"))
cyp51B_by_abundance_presence_absence<- data.frame(sapply(cyp51B_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(cyp51B_by_abundance_presence_absence)<- row.names(cyp51B_colapse)

AFUA_7G01960_by_abundance_presence<- data.frame(sapply(AFUA_7G01960_by_abundance, gsub, pattern = "1", replacement = "present"))
AFUA_7G01960_by_abundance_presence_absence<- data.frame(sapply(AFUA_7G01960_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(AFUA_7G01960_by_abundance_presence_absence)<- row.names(AFUA_7G01960_colapse)


#mdr4
#cdr1B_abcC
#hmg1
#hapE
#mdr3
#cox10
#mdr2
#mdr1
#atrF
#fks1
#7G01960
#cyp51B
#bind and highlight - can't get width right when adjacent. 
#all_vars<- cbind(by_variant_output_df_binary_presence_absence, by_variant_output_df_binary_presence_absence_unknowns)


#bind and plot together
all_all_vars<- cbind(by_variant_output_df_binary_presence_absence_unknowns,
                     mdr4_by_abundance_presence_absence,
                     cdr1B_abcC_by_abundance_presence_absence,
                     hmg1_by_abundance_presence_absence,
                     hapE_by_abundance_presence_absence,
                     mdr3_by_abundance_presence_absence,
                     cox10_by_abundance_presence_absence,
                     mdr2_by_abundance_presence_absence,
                     mdr1_by_abundance_presence_absence,
                     atrF_by_abundance_presence_absence,
                     fks1_by_abundance_presence_absence,
                     AFUA_7G01960_by_abundance_presence_absence,
                     cyp51B_by_abundance_presence_absence)
  
ncol(cyp51B_by_abundance_presence_absence)
all_tree_plot <-  gheatmap(tre_to_plot, 
                           all_all_vars,
                           offset=0, width=.8, low="white", high="black", 
                           colnames = T, 
                           colnames_angle = 90,
                           colnames_position = "top",
                           font.size = .9,
                           colnames_offset_y = -1,
                           colnames_offset_x = 0,
                           color="white") +
  scale_fill_manual(values=c("white", "black")) 
#+ggtitle("all variants")

all_tree<- all_tree_plot + geom_tiplab(size = .7, align = TRUE, linesize = .25, offset = 1, linetype = 0)
all_tree

#export
#ggsave(file="Avian_all_resistnace_gene_variants.pdf",device="pdf", all_tree, width=7, height=2, units="in")


####
#BASIC STATS/COUNTS ETC. 
metadata<-read.delim("Avian_metadata.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

#how many bird species?
length(unique(metadata$HOST)) #23 + ICMP strain = 24
#how many of each host ecology?
table(metadata$HOST_ECOLOGY)
#how many from multiple death events?
table(metadata$Aspergillosis.Event)

#how many ea event?
how_many_ea_event<- data.frame(table(metadata$Event))
how_many_ea_event_ordered<- how_many_ea_event[order(how_many_ea_event$Freq),]
View(how_many_ea_event_ordered)

#how many events?
as_event_only<- metadata[metadata$Aspergillosis.Event == "TRUE",]
as_event_only_how_many<- length(unique(as_event_only$Event))
as_event_only_how_many

#how many tissue types?
tissue_types<- data.frame(table(metadata$TISSUE))
nrow(tissue_types)
tissue_types
#confirm with authors - table is a little confusing here. 
##YOU ARE HERE - FIX TISSUE CATS
list<- c(140, 6,18,2,15,62,300,140,85,47)
mean(list)

#look at each event individually
#View(table(as_event_only$Event))
table(as_event_only$Event) # seven events with multiple isolates


#how many potential events?
as_event_only_p<- metadata[metadata$Aspergillosis.Event == "UNKNOWN",]
as_event_only_how_many_p<- length(unique(as_event_only_p$Event))
as_event_only_how_many_p
unique(as_event_only_p$Event)

#cause of death primary aspergilliosis?
cause<- sum(metadata$Aspergillosis == "TRUE")
cause


