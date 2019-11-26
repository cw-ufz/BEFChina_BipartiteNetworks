######## calculate bipartite network values for magic 7 subsampling ############
#
# R script (07/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

#this script was only used to test the analysis on the laptop before
#putting the complete analysis onto a high performance computer cluster
#so this script is not fully functional 

################################################################################
library(bipartite)
library(plyr)
library(vegan)
library(phyloseq)


source(file = "./R_functions/06F_CalcLinkOptions.r")
source(file="./R_functions/06Fb_NetworkTopologies.r")  #select in function which 
                                                      #metrics to calculate!
source(file="./R_functions/06Fc_NetworkTopologiesGroups.r")

####################

#Load datasets: 
combination_samples<-readRDS(file = "./output/06_tree_combinations_samples.RDS")

# [[1]] : Tree monocultures subsample combinations
# [[2]] : Two tree species subsample combinations
# [[3]] : High tree diversity subsample combinations 
#(4,8,16 diversity plots summed)
# [[]][[1]] - [[]][[576]] subsample combinations a tree diversity level

#one loop: div levels (1:3)
#second loop: subsampling combinationa (1:576)
#third loop: percentages valid links (1:5) (20,40,60,80,100)


#select your choice of metrics and adjust the code:

#Calculate network metrics and write to table
net_subsample_mono<-matrix(data = NA, nrow = 576, ncol = 12,
                           dimnames=list(1:576, c("modularity", "connectance",
                           "weighted_connectance",
                           "number_OTUs", "Fungal_generality",
                           "mean_shared_treespecies", "mean_shared_fungalpartners",  "fungal_niche_overlap",
                           "fungal_togetherness", "fungal_Cscore", "richness", "Shannon_diversity")))


net_subsample_mono<-data.frame(net_subsample_mono)
net_subsample_di<-net_subsample_mono
net_subsample_high<-net_subsample_mono

net_subsample<-list(net_subsample_mono,net_subsample_di,net_subsample_high)

i=1
j=1
link_mode<-"p100"
read_data<-"pa"

for(j in 1:3){
  for(i in 1:10){
#1:576
  
    #reduce otu table
    pruned_set<-prune_taxa(taxa_sums(combination_samples[[j]][[i]])>0 , combination_samples[[j]][[i]])
  
    otu_summed<-colSums(otu_table(pruned_set))
    Shannon_summed<-diversity(otu_summed, index = "shannon", MARGIN = 1, base = exp(1))
    richness_summed<-specnumber(otu_summed, MARGIN = 1)
    
    #select needed information
    combi_melted<-psmelt(pruned_set)
    combi_melted_prep<-combi_melted[,c("OTU","Tree_s","Div_Plot","Abundance")] 
    names(combi_melted_prep)<-c("higher", "lower", "webID", "freq")
  
    #calculate link options
    combi_melted_prelink<-CalcLinkOptions(combi_melted_prep)

    #select the appropriate link type
    selected_link<-combi_melted_prelink[combi_melted_prelink[link_mode]==TRUE,
                           c("higher","lower","webID", read_data,link_mode)]
    selected_link$webID<-1
    network_data<-unique(selected_link[-5])
    names(network_data)<-c("higher","lower","webID","freq")  
  
    #calculate the network
    network_frame<-frame2webs(network_data,type.out="list") 
  
    #calculate network metrics
    network_topologies<-NetworkTopologies(network_frame)
    #network_group_topologies <- NetworkTopologiesGroups(network_frame)
 
    #network_group_frame<-data.frame(network_group_topologies)
    network_dataframe<-data.frame(network_topologies)
   
    
    net_subsample[[j]][i,"modularity"]<-sapply(network_frame,function(x) computeModules(x)@likelihood)[1]
    net_subsample[[j]][i,"connectance"]<-network_dataframe["connectance",]
    net_subsample[[j]][i,"weighted_connectance"]<-network_dataframe["weighted connectance",]
    net_subsample[[j]][i,"number_OTUs"]<-network_dataframe["number.of.species.HL",]

    net_subsample[[j]][i,"Fungal_generality"]<-network_group_frame["generality.HL",]
    net_subsample[[j]][i,"mean_shared_treespecies"]<-network_group_frame["mean.number.of.shared.partners.HL",]
    net_subsample[[j]][i,"mean_shared_fungalpartners"]<-network_group_frame["mean.number.of.shared.partners.LL",]
    net_subsample[[j]][i,"fungal_niche_overlap"]<-network_group_frame["niche.overlap.HL",]
    net_subsample[[j]][i,"fungal_togetherness"]<-network_group_frame["togetherness.HL",]
    net_subsample[[j]][i,"fungal_Cscore"]<-network_group_frame["C.score.HL",]
    net_subsample[[j]][i,"richness"]<-richness_summed
    net_subsample[[j]][i,"Shannon_diversity"]<-Shannon_summed

    net_subsample[[j]][i,"no_modules"]<-length(sapply(network_frame,function(x) computeModules(x)@likelihood))
  }

}

















