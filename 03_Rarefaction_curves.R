############################### Rarefaction curves ############################

# R script (03/16) for analyses in Weissbecker et al. 2019
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

###############################

library(extrafont)  # extrafont_0.17    #make available font "Arial"
library(vegan)      # vegan_2.5-4 
library(ggplot2)    # ggplot2_3.1.0
library(dplyr)
library(forcats)

#font_import()      # install fonts if neccessary
#fonts()

#sessionInfo()
#R version 3.5.2 (2018-12-20)
#######################

BEF_10r<-readRDS(file="./output/02_BEF_10r.RDS") 

BEF_10r_high3<-subset_samples(BEF_10r, Div_Plot %in% c(4,8,16)) # 1926 taxa 220 samples
BEF_10r_low2<-subset_samples(BEF_10r, Div_Plot %in% c(2)) #1926 taxa 76 samples
BEF_10r_mono1<-subset_samples(BEF_10r, Div_Plot %in% c(1)) #1926 taxa 77 samples


#-----manual file selection
Dataset<-BEF_10r_mono1 #BEF_10r_high3 or BEF_10r_low2 or BEF_10r_mono1

BEF_otu<-otu_table(Dataset)
BEF_otu_t<-t(BEF_otu)
BEF_otu_df<-data.frame(BEF_otu_t)


#adjust names manually in file and figure!
tiff(filename ="./images/03_Rarecurves_BEF_10r_mono1.tiff",
     width = 6.8, height = 6.8, units = "in", res=300)


  rarecurve(BEF_otu_df, step = 10, sample = 800, label = FALSE, 
          main="BEF_10r_mono1" , ylab="Number of OTUs", xlab="Number of sequence reads")


dev.off()


#sequences per sample
#adjust names manually in file and figure!
tiff(filename = "./images/03_seqs_BEF_10r.tiff", 
     width = 4 , height =6 , units = "in", res=300)

plot(sample_sums(Dataset), ylim=c(0,700), main="BEF_10r",
     ylab="read number", xlab="sample number")

dev.off()

