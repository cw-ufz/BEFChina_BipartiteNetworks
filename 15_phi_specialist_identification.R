################# identify the phi coefficient top 20 specialists ##############

# R script (15/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25


################################################################################

library(phyloseq)

##########################


BEF_data<-readRDS(file="./Cleaned_Data/02_BEF_10r.RDS") 
phi_top20_list<-readRDS(file="./output/14_phi_top20_otulist.RDS")

BEF_reduced<-prune_taxa(phi_top20_list, BEF_data)

phi_top20_specialist_tax<-as.data.frame(BEF_reduced@tax_table@.Data)


#----- output
write.csv2(phi_top20_specialist_tax, file="./output/15_phi_tax_top20spec.csv")







