######## calculate bipartite network values for magic 7 subsampling ############
#
# R script (07b/) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

#this script was only used to test the analysis on the laptop before
#putting the complete analysis onto a high performance computer cluster
#so this script is not fully functional 

#here now added the nestedness calculations


################################################################################
library(bipartite)
library(plyr)
library(vegan)
library(phyloseq)


source(file = "./R_functions/07F_CalcLinkOptions.r") #select in function which 
                                                    #metrics to calculate!
source(file="./R_functions/07Fb_NetworkTopologies.r")
####################


#Load datasets: 
combination_samples<-readRDS(file = "./output/05_tree_combinations_samples.RDS")
combination_samples_nm<-readRDS(file = "./output/05_tree_combinations_samples_nm.RDS")

# [[1]] : Tree monocultures subsample combinations
# [[2]] : Two tree species subsample combinations
# [[3]] : High tree diversity subsample combinations 
#(4,8,16 diversity plots summed)
# [[]][[1]] - [[]][[576]] subsample combinations a tree diversity level

#one loop: div levels (1:3)
#second loop: subsampling combinationa (1:576)
#third loop: percentages valid links (1:5) (20,40,60,80,100)

net_subsample_mono<-matrix(data = NA, nrow = 576, ncol = 4,
                           dimnames=list(1:576, c("nestedness",
                                                  "nestedness_nm",
                                                  "NODF",
                                                  "NODF_nm")))

net_subsample_mono<-data.frame(net_subsample_mono)
net_subsample_di<-net_subsample_mono
net_subsample_high<-net_subsample_mono

net_subsample<-list(net_subsample_mono,net_subsample_di,net_subsample_high)

#head(View(net_subsample_mono))

i=1
j=1
link_mode<-"p100"   
read_data<-"pa"     


for(j in 1:3){
  for(i in 1:10){
#1:576
  
    #reduce otu table
    pruned_set<-prune_taxa(taxa_sums(combination_samples[[j]][[i]])>0 , combination_samples[[j]][[i]])
    pruned_set_nm<-prune_taxa(taxa_sums(combination_samples_nm[[j]][[i]])>0 , combination_samples_nm[[j]][[i]])
    
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
  
    View(network_data_nm)
    #calculate the network
    network_frame<-frame2webs(network_data,type.out="list") 
    network_frame_nm<-frame2webs(network_data_nm,type.out="list")
  
    #calculate network metrics
    network_topologies<-NetworkTopologies(network_frame)
    network_dataframe<-data.frame(network_topologies)
   
    network_topologies_nm<-NetworkTopologies(network_frame_nm)
    network_dataframe_nm<-data.frame(network_topologies_nm)
    
    net_subsample[[j]][i,"nestedness"]<-network_dataframe["nestedness",]
    net_subsample[[j]][i,"nestedness_nm"]<-network_dataframe_nm["nestedness",]
    net_subsample[[j]][i,"NODF"]<-network_dataframe["NODF",]
    net_subsample[[j]][i,"NODF_nm"]<-network_dataframe_nm["NODF",]
    }

}




