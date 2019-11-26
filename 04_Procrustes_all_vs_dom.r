####### Procrustes correlation between complete and dominant OTU dataset #######

# R script (04/16) for analyses in Weissbecker et al. 2019
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

######################
library(vegan)            #vegan_2.5-2
library(biomformat)       #biomformat_1.8.0 #read biom file
library(phyloseq)         #phyloseq_1.24.0  
#sessionInfo()
#R version 3.5.1 (2018-07-02)
######################

BEF_all_otu<-readRDS(file="./Cleaned_Data/01_BEF_all_otu_table.RDS")

BEF_10r<-readRDS(file="./Cleaned_Data/02_BEF_10r.RDS") 

Dataset<-BEF_10r
sample_list<-sample_names(Dataset)
BEF_all_otu_subset<-BEF_all_otu[,sample_list]

all_table<-otu_table(BEF_all_otu_subset) 
nrow(all_table) # 21732 taxa
dom_table<-otu_table(Dataset) 
nrow(dom_table)  # 1926 taxa

#----- NMDS calculation for the complete dataset 

otu_all<-data.frame(all_table, check.names = FALSE) 
t_all<-t(otu_all)

otu_hell_all<-decostand(t_all,method="hellinger") 

otu_nmds_all<-metaMDS(otu_hell_all, distance = "bray", k = 2, trymax = 30,
                      autotransform =FALSE)
otu_nmds_all$stress #0.1657157 


tiff(filename = "./images/04_nmds_all.tiff", 
     width = 4 , height =6 , units = "in", res=300)

  ordiplot(otu_nmds_all,type="n",main="BEF_all")
  orditorp(otu_nmds_all,display="sites",cex=1,air=0.01)

dev.off()

#----- NMDS calculation for the rarefied and thresholded OTU dataset 

dom_otu<-data.frame(dom_table, check.names = FALSE)
dom_t<-t(dom_otu) 

otu_hell<-decostand(dom_t,method="hellinger") 

otu_nmds<-metaMDS(otu_hell, distance = "bray", k = 2, trymax = 30, 
                  autotransform =FALSE)
otu_nmds$stress # 0.201795

tiff(filename = "./images/04_nmds_BEF_10r.tiff", 
     width = 4 , height =6 , units = "in", res=300)

  ordiplot(otu_nmds,type="n", main="BEF_10r")
  orditorp(otu_nmds,display="sites",cex=1,air=0.01)

dev.off()

#----- Procrustes correlation analysis of both NMDS ordinations

pro<-protest(otu_nmds_all, otu_nmds)

# Procrustes Sum of Squares (m12 squared):        0.1041 
# Correlation in a symmetric Procrustes rotation: 0.9465 
# Significance:  0.001 
# 
# Permutation: free
# Number of permutations: 999


tiff(filename = "./images/04_Procrustes_all_vs_BEF_10r.tiff", 
     width = 5, height = 6.8, units = "in", res=300) 
  
  plot(pro, cex=2, main="allvs. BEF 10r")

dev.off()





