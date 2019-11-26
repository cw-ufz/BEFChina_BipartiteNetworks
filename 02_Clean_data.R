############################### Clean data ####################################

# R script (02/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

###############################

library(vegan)        #vegan_2.5-4        #decostand function
library(phyloseq)     #phyloseq_1.26.1

#sessionInfo()
#R version 3.5.2 (2018-12-20)
###############################

#----- load data

BEF_dom<-readRDS(file="./Cleaned_Data/01_dominant4_phyloseq.RDS")
all_var_sorted<-readRDS(file="./Cleaned_Data/01_dominant4_metadata.RDS")

BEF_dom_fungi <- subset_taxa(BEF_dom, Kingdom=="Fungi") 


#----- check read distribution
#some samples have a high amount of rare OTUs and/or low target fungal sequences

tiff(filename = "./images/02_plot_dom4fungi_samplereads.tiff", 
     width = 4 , height =6 , units = "in", res=300)

  plot(sample_sums(BEF_dom_fungi), ylim=c(0,3500), main="BEF__dom_fungi sample reads",
     ylab="read number", xlab="sample number")

dev.off()

plot(sample_sums(BEF_dom_fungi), ylim=c(0,700)) 
length(which(sample_sums(BEF_dom_fungi)<700)) #21 samples

#----- rarefy my data to even sampling depth 

BEF_rare700<-rarefy_even_depth(BEF_dom_fungi, sample.size = 700,
                               rngseed = 711, replace = TRUE, trimOTUs = TRUE, 
                               verbose = TRUE)
# rngseed=711 is the default value
# function info:21 samples removedbecause they contained fewer reads than `sample.size`.
#               328OTUs were removed

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 5337 taxa and 373 samples ]
#sample_data() Sample Data:       [ 373 samples by 22 sample variables ]
#tax_table()   Taxonomy Table:    [ 5337 taxa by 9 taxonomic ranks ]

tiff(filename = "./images/02_BEF_rare_700_samplereads.tiff", 
     width = 4 , height =6 , units = "in", res=300)

  plot(sample_sums(BEF_rare700), main="BEF_rare_700 sample reads",
     ylab="read number", xlab="sample number")

dev.off()

BEF_10r = prune_taxa(taxa_sums(BEF_rare700) >= 10, BEF_rare700)
#otu_table()   OTU Table:         [ 1926 taxa and 373 samples ]
#sample_data() Sample Data:       [ 373 samples by 22 sample variables ]
#tax_table()   Taxonomy Table:    [ 1926 taxa by 9 taxonomic ranks ]


#----- script output

saveRDS(BEF_10r, file="./Cleaned_Data/02_BEF_10r.RDS")
saveRDS(BEF_rare700, file="./Cleaned_Data/02_BEF_rare700.RDS")

