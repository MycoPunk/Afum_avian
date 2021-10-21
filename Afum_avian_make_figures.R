#make figures for A.fum_avian
#for net tree see:
#https://cran.r-project.org/web/packages/phangorn/vignettes/IntertwiningTreesAndNetworks.html

#set packages 
library(data.table)
library(tidyverse)
library(hrbrthemes)
library(ggtree)
library(ape)
library(phytools)
library(ggplot2)
library(stringr)
library(treeio)
library(dplyr)
library(albersusa)
library(usmap)
library(cowplot)


setwd("~/data")
options(stringsAsFactors = FALSE)

#read in tree
tree_me <- read.tree("Avian_iq_tree.tre")

#remove the reference from the tree:
tree_me<- drop.tip(tree_me, "Af293-REF", trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0, rooted = is.rooted(tree_me), collapse.singles = TRUE,
                   interactive = FALSE)

#rename tips
name_map<-read.delim("Avian_name_map.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

#rename using tree.io
tree_me2<- treeio::rename_taxa(tree_me, name_map, OG_tree_name, New_tree_name) #%>% write.tree
#clade1_names_fixed <- name_map$name_pop_genome[match(names(clade1_names), name_map$name_Pan_genome)]

#set groups
avian_names<-tree_me2$tip.label[grep(pattern = 'USGS', x = tree_me2$tip.label)]
avian_names_icmp<-tree_me2$tip.label[grep(pattern = 'ICMP', x = tree_me2$tip.label)]
avian_names_all<-tree_me2$tip.label[grep(pattern = 'USGS|ICMP', x = tree_me2$tip.label)]
outgroup_names<- tree_me2$tip.label[!tree_me2$tip.label %in% avian_names_all]

avian_names_df<- cbind(strain = avian_names, group = "A")
avian_names_icmp_df<- cbind(strain = avian_names, group = "K")
outgroup_names_df<- cbind(strain = outgroup_names, group = "O")
groups_df<- data.frame(rbind(avian_names_df, avian_names_icmp_df,outgroup_names_df))

grA_me<- split(groups_df$strain, groups_df$group)

#set colors by group
tree_grA_me <- ggtree::groupOTU(tree_me2, grA_me)
str(tree_grA_me)
levels(attributes(tree_grA_me)$group) 
levels(attributes(tree_grA_me)$group)[1] <- "A"

attributes(tree_grA_me)$group <- factor(x = attributes(tree_grA_me)$group, 
                                        levels = c("A", "K", "O"))


#full pallet:
#dark blue: #2A4359 (Kakapoo)
#light blue: #58788C
#yellow: #F2BF27
#orange: #F2811D
#dark orange: #F2490C (Avian)
#almost black: #0D0D0D (Outgroups)


#set colors
my_cols_me <- c(Kakapoo = "#58788C",
                Avian = "#F2490C",
                Other ="#0D0D0D")

#check colors
names(my_cols_me) <- levels(attributes(tree_grA_me)$group)
scales::show_col(my_cols_me); my_cols_me


##make Fig. 2
#plot tree
tree_plot_me <- 
  ggtree(tr = tree_grA_me, 
         # color by group attribute, check str(tree_grA_me)
         mapping = aes(color = group), 
         layout  = 'circular', 
         #branch.length = 'none', 
         #  geom_treescale(x=3, y=NULL, color = "white") +
         # set line thickness
         size = .3) +
  # adjust coloring of main groups
  scale_color_manual(name = 'Group', values = my_cols_me) + 
  theme(legend.title=element_text(size=9), # The title of legend 
        legend.text=element_text(size=7))
#xlim(NA, NA)
#  guides(color = guide_legend(override.aes = list(size = 4))) 

# plot and ddd the tip labels
tree_plot_me + geom_tiplab(size = 1, align = TRUE, linesize = .1, linetype = 3) +
geom_point2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 95), size = 1.5, col= "black") +
geom_point2(aes(label=label, subset = !is.na(as.numeric(label)) & (as.numeric(label) > 80) & (as.numeric(label) < 95)), size = 1.5, col= "dark grey") +
geom_point2(aes(label=label, subset = !is.na(as.numeric(label)) & (as.numeric(label) > 70) & (as.numeric(label) < 85)), size = 1.5, col= "light grey") +
geom_point2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) < 70), size = 1.5, col= "light grey", shape = 1) #make hollow 

#print that shit
#ggsave("Afum_avian_Big_tree.pdf", width=20, height=20, 
#       device = "pdf", units = "cm")


#to add annotations below - note- I ended up just adding these in illustrator though, because it offerend more flexibility for shapes+color mapping.  

#add event annotations 
#event_map<-read.delim("event_map.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
##split byevent
#grA_me_event<- split(event_map$New_tree_name, event_map$event)
####set colors by group for event type 
#tree_grA_me_event <- ggtree::groupOTU(tree_me2, grA_me_event)
#str(tree_grA_me_event)
#tree_grA_me_mat2 <- ggtree::groupOTU(tree_grA_me_mat, grA_me)
#str(tree_grA_me_mat2)
#levels(attributes(tree_grA_me_event)$group) 
#levels(attributes(tree_grA_me_event)$group)[1] <- "A"
#attributes(tree_grA_me_event)$group <- factor(x = attributes(tree_grA_me_event)$group, 
#                                            levels = c("A","B","C","D","E","F","G"))
#my_cols_me_event <- c(A = "#8C824D",
#                      B ="#736A4F",
#                      C = "#705E78", 
#                      D = "#A3A1A8", 
#                      E = "#D95F69",
#                      F = "#F2BE24",
#                      G = "#BACAC0", 
#                      H = "#FFFFFF")

#names(my_cols_me_event) <- levels(attributes(tree_grA_me_event)$group)
#scales::show_col(my_cols_me_event); my_cols_me_event

#tree_plot_me <- 
#  ggtree(tr = tree_grA_me, 
#         # color by group attribute, check str(tree_grA_me)
#         mapping = aes(color = group), 
#         layout  = 'circular', 
#         branch.length = 'none', 
#         #  geom_treescale(x=3, y=NULL, color = "white") +
#         # set line thickness
#         size = .3) +
#  # adjust coloring of main groups
#  scale_color_manual(name = 'Clade', values = my_cols_me) + 
#  theme(legend.title=element_text(size=9), # The title of legend 
#        legend.text=element_text(size=7))
#xlim(NA, NA)
#  guides(color = guide_legend(override.aes = list(size = 4))) 


# plot and ddd the tip labels
#tree_plot_me + geom_tiplab(size = 1, align = TRUE, linesize = .1, linetype = 0)
#ggsave(file="pan_genome_tree_cladogram.png",device="png")


#base plot
#tree_plot_me_event <- 
#  ggtree(tr = tree_grA_me_event, 
#          #color by group attribute, check str(tree_grA_me)
#         mapping = aes(color = group), 
#         layout  = 'circular', 
#         branch.length = 'none') + 
#  geom_treescale(x=3, y=NULL, color = "white") +
  # set line thickness
  #  size = 1)
  # adjust coloring of main groups
#  scale_color_manual(name = 'Clade', values = my_cols_me_event) + 
  #scale_shape_manual("Clade", values = my_cols_me_mat, breaks=c("MAT1", "MAT2", "Unknown"), labels=c("O","R","N"))+
#  theme(legend.title=element_text(size=9), # The title of legend 
#        legend.text=element_text(size=7))
#guides(color = guide_legend(override.aes = list(size = 4))) 


#tree_plot_me_event + geom_tiplab(size = 1, align = TRUE, linesize = 0, linetype = 0) +
#  geom_tippoint(aes(x=x+16, color=group), size=.50)


#combine into one tree
#pi<- tree_plot_me + geom_tiplab(size = .6, align = TRUE, linesize = .20, linetype = NA) 
#p <- pi %<+% event_map + geom_tippoint(aes(x=x+14,color=event), size=.30) + scale_color_manual(values=c("#8C824D",
#                                                                                                         "#736A4F",
#                                                                                                         "#705E78", 
#                                                                                                         "#A3A1A8", 
#                                                                                                         "#D95F69",
#                                                                                                         "#F2BE24",
#                                                                                                         "#BACAC0", 
#                                                                                                         "#FFFFFF",
#                                                                                                         "#FFFFFF",
#                                                                                                         "#FFFFFF"))










#####
##Make Fig. 1B - small tree
#subset
tree_me_small <- keep.tip(tree_me, avian_names_all)

#read in metadata
metadata<-read.delim("Avian_metadata.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)


#split data by host ecology
grA_me_small<- split(metadata$STRAIN_NAME, metadata$HOST_ECOLOGY)

#spllit by year
grA_me_year<- split(metadata$STRAIN_NAME, metadata$year.collected)

#set colors by group
tree_grA_me_small <- ggtree::groupOTU(tree_me_small, grA_me_small)
str(tree_grA_me_small)
levels(attributes(tree_grA_me_small)$group) 
levels(attributes(tree_grA_me_small)$group)[1] <- "raptor"
levels(attributes(tree_grA_me_small)$group) 

attributes(tree_grA_me_small)$group <- factor(x = attributes(tree_grA_me_small)$group, 
                                              levels = c("raptor","sea bird","water bird"))

my_cols_me_small <- c(raptor = "#58788C",
                      sea_bird ="#0D0D0D",
                      water_bird = "#F2811D")
names(my_cols_me_small) <- levels(attributes(tree_grA_me_small)$group)
scales::show_col(my_cols_me_small); my_cols_me_small

#dark blue: #2A4359 (Kakapoo)
#light blue: #58788C
#yellow: #F2BF27
#orange: #F2811D
#dark orange: #F2490C (Avian)
#almost black: #0D0D0D (Outgroups)

####
tree_grA_me_year <- ggtree::groupOTU(tree_me_small, grA_me_year)
str(tree_grA_me_year)
tree_grA_me_year2 <- ggtree::groupOTU(tree_grA_me_small, grA_me)
str(tree_grA_me_year2)
levels(attributes(tree_grA_me_year)$group) 
#levels(attributes(tree_grA_me_mat)$group)[1] <- "1"
attributes(tree_grA_me_year)$group <- factor(x = attributes(tree_grA_me_year)$group, 
                                             levels = c("2014","2015","2016", "2017", "2018", "2019"))
my_cols_me_year <- c("2014" = "#b09e4d",
                     "2015" ="#758e4e",
                     "2016" = "#437a55",
                     "2017" = "#1c6257", 
                     "2018" = "#10494f", 
                     "2019" = "#18313b")

names(my_cols_me_year) <- levels(attributes(tree_grA_me_year)$group)
scales::show_col(my_cols_me_year); my_cols_me_year
####


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


# plot and ddd the tip labels
Fig1B<- tree_plot_me_sm + geom_tiplab(size = 1, align = TRUE, linesize = .1, linetype = 3, show.legend=FALSE)
Fig1B


#plot with year codes
pi<- Fig1B + geom_tiplab(size = .6, align = TRUE, linesize = .20, linetype = NA) 
Fig1B_dot <- pi %<+% metadata + geom_tippoint(aes(x=1.3,color=as.character(year.collected)), size=1.8, shape = 15) + 
  scale_color_manual(values=c("#b09e4d",
                              "#758e4e",
                              "#437a55",
                              "#1c6257", 
                              "#10494f", 
                              "#18313b",
                              "#58788C",
                              "#0D0D0D",
                              "#F2811D"))
plot(Fig1B_dot)

####
##make Fig. 1A - map
#format data
metadata_to_plot<- data.frame(table(metadata$state, metadata$HOST_ECOLOGY))
colnames(metadata_to_plot)<- c("state", "host", "freq")

state_level_df <- data.frame(state = tolower(state.name), 
                             long = state.center$x, 
                             lat = state.center$y,
                             stringsAsFactors = FALSE) %>%
  inner_join(metadata_to_plot, by="state" )

#fix Hawaii as the state.center function does not provide legit location data for HI
hawaii1<- c("hawaii", -157.85809, 21.31560, "raptor", 0)
hawaii2<- c("hawaii", -157.85809, 21.31560, "sea bird", 1)
hawaii3<- c("hawaii", -157.85809, 21.31560, "water bird",1)
state_level_df2<- state_level_df[-c(7,8,9), ]
state_level_df3<- rbind(state_level_df2, hawaii1, hawaii2, hawaii3)

#rescale hawaii points using library(albersusa)
subset<- data.frame(cbind(as.numeric(state_level_df3$long)), as.numeric(state_level_df3$lat), state_level_df3$host, state_level_df3$freq)
colnames(subset)<- c("long", "lat", "host", "freq")
#new.lat <- points_elided(subset)
subset_trans <- usmap_transform(subset)

#remove zero values 
subset_trans_no_zero<- subset_trans[!subset_trans$freq ==0,]

#plot
p<- plot_usmap() + geom_jitter(width = 60000, height = 60000,
  data = subset_trans_no_zero,
  aes(x = long.1, y = lat.1, size = as.numeric(freq), color=host),
  alpha = 0.5)
Fig1A<- p+ scale_color_manual(values=c("#58788C", "#0D0D0D", "#F2811D"))+
  labs(size="n isolates", colour="Host") + scale_size(
    name = waiver(),
    breaks = c(1,4,8,12),
    labels = c(1,4,8,12),
    limits = c(1,12),
    range = c(2, 6),
  ) + theme(legend.direction = "vertical", 
            legend.position = "bottom",
            legend.box = "horizontal"
  )
Fig1A


#my_cols_me_small <- c(raptor = "#58788C",
#                      sea_bird ="#0D0D0D",
#                      water_bird = "#F2811D")

#combine plots
p <- plot_grid(Fig1A,
               Fig1B_dot,
               labels = c("A", "B"), 
               #align = "hv", 
               align = "none",
               axis = "rlbt", 
               nrow = 1, 
               ncol =2, 
               byrow = TRUE,
               rel_widths = c(2.5, 1), 
               rel_heights = 1,2)
p

#print
ggsave("Afum_avian_Fig1.pdf", width=20, height=16, 
       device = "pdf", units = "cm")

