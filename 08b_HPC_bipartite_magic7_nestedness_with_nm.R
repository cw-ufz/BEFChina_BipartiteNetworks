######## calculate bipartite network values for magic 7 subsampling ############
#
# R script (08b/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

################################################################################
library(bipartite)
library(plyr)
library(vegan)
library(phyloseq)

args<-commandArgs(TRUE)
div<-as.numeric(args[1])
link_mode<-args[2]
input_real<-args[3]
input_nm<-args[4]

#Pfade anpassen
source(file = "/gpfs1/work/lachmann/Nets/scripts/07F_CalcLinkOptions.R")
source(file="/gpfs1/work/lachmann/Nets/scripts/07Fb_NetworkTopologies.R")


####################

#Load datasets: 

print(input_real)
print(input_nm)

combination_samples<-readRDS(file = input_real) 
combination_samples_nm<-readRDS(file = input_nm) #mit oder ohne Pfad, je nach aufruf

net_subsample<-data.frame(matrix(data = NA, nrow = 576, ncol = 4,
                                 dimnames=list(1:576,c("nestedness",
                                                       "nestedness_nm",
                                                       "NODF",
                                                       "NODF_nm"))))
  
read_data<-"pa"

  for(i in 1:576){
#1:576

    #reduce otu table
    pruned_set<-prune_taxa(taxa_sums(combination_samples[[div]][[i]])>0 , combination_samples[[div]][[i]])
    pruned_set_nm<-prune_taxa(taxa_sums(combination_samples_nm[[div]][[i]])>0 , combination_samples_nm[[div]][[i]])
    
     
    #select needed information
    combi_melted<-psmelt(pruned_set)
    combi_melted_nm<-psmelt(pruned_set_nm)
    combi_melted_prep<-combi_melted[,c("OTU","Tree_s","Div_Plot","Abundance")] 
    combi_melted_prep_nm<-combi_melted_nm[,c("OTU","Tree_s","Div_Plot","Abundance")]
    names(combi_melted_prep)<-c("higher", "lower", "webID", "freq")
    names(combi_melted_prep_nm)<-c("higher", "lower", "webID", "freq")
  
    #calculate link options
    combi_melted_prelink<-CalcLinkOptions(combi_melted_prep)
    combi_melted_prelink_nm<-CalcLinkOptions(combi_melted_prep_nm)

    #select the appropriate link type
    selected_link<-combi_melted_prelink[combi_melted_prelink[link_mode]==TRUE,
                           c("higher","lower","webID", read_data,link_mode)]
    selected_link_nm<-combi_melted_prelink_nm[combi_melted_prelink_nm[link_mode]==TRUE,
                                        c("higher","lower","webID", read_data,link_mode)]
    selected_link$webID<-1
    selected_link_nm$webID<-1
    network_data<-unique(selected_link[-5])
    network_data_nm<-unique(selected_link_nm[-5])
    names(network_data)<-c("higher","lower","webID","freq")  
    names(network_data_nm)<-c("higher","lower","webID","freq")
  
    #calculate the network
    network_frame<-frame2webs(network_data,type.out="list") 
    network_frame_nm<-frame2webs(network_data_nm,type.out="list")
  
    #calculate network metrics
    network_topologies<-NetworkTopologies(network_frame)
    network_dataframe<-data.frame(network_topologies)
    
    network_topologies_nm<-NetworkTopologies(network_frame_nm)
    network_dataframe_nm<-data.frame(network_topologies_nm)
    
    net_subsample[i,"nestedness"]<-network_dataframe["nestedness",]
    net_subsample[i,"nestedness_nm"]<-network_dataframe_nm["nestedness",]
    net_subsample[i,"NODF"]<-network_dataframe["NODF",]
	#
    net_subsample[i,"NODF_nm"]<-network_dataframe_nm["NODF",]
	#
    
  }

outfile<- paste0("/gpfs1/work/lachmann/Nets/output/", "net_subsample.nm." ,div, ".", link_mode, ".RDS")

saveRDS(net_subsample, file = outfile) 
















