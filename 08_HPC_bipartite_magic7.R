######## calculate bipartite network values for magic 7 subsampling ############
#
# R script (08/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

# this script was for use on a univa grid engine high performance computer cluster
# (HPC) of the UFZ (named "EVE")

################################################################################
library(bipartite)
library(plyr)
library(vegan)
library(phyloseq)

args<-commandArgs(TRUE)
div<-as.numeric(args[1])
link_mode<-args[2]
input<-args[3]

#Pfade anpassen
source(file = "/gpfs1/work/lachmann/Nets/scripts/07F_CalcLinkOptions.R")
source(file="/gpfs1/work/lachmann/Nets/scripts/07Fb_NetworkTopologies.R")
source(file="/gpfs1/work/lachmann/Nets/scripts/07Fc_NetworkTopologiesGroups.R")





####################

#Load datasets: 
combination_samples<-readRDS(file = input) #mit oder ohne Pfad, je nach aufruf


 net_subsample<-data.frame(matrix(data = NA, nrow = 576, ncol = 12,
                            dimnames=list(1:576, c("modularity", "connectance", 
                            "weighted_connectance",
                            "number_OTUs", "Fungal_generality", 
                            "mean_shared_treespecies", "mean_shared_fungalpartners",  "fungal_niche_overlap",
                            "fungal_togetherness", "fungal_Cscore", "richness", "Shannon_diversity"))))

read_data<-"pa"

  for(i in 1:576){
#1:576

    #reduce otu table
    pruned_set<-prune_taxa(taxa_sums(combination_samples[[div]][[i]])>0 , combination_samples[[div]][[i]])
  
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
    
    net_subsample[i,"modularity"]<-sapply(network_frame,function(x) computeModules(x)@likelihood)[1]
    net_subsample[i,"connectance"]<-network_dataframe["connectance",]
    net_subsample[i,"weighted_connectance"]<-network_dataframe["weighted connectance",]
    net_subsample[i,"number_OTUs"]<-network_dataframe["number.of.species.HL",]

    net_subsample[i,"Fungal_generality"]<-network_group_frame["generality.HL",]
    net_subsample[i,"mean_shared_treespecies"]<-network_group_frame["mean.number.of.shared.partners.HL",]
    net_subsample[i,"mean_shared_fungalpartners"]<-network_group_frame["mean.number.of.shared.partners.LL",]
    net_subsample[i,"fungal_niche_overlap"]<-network_group_frame["niche.overlap.HL",]
    net_subsample[i,"fungal_togetherness"]<-network_group_frame["togetherness.HL",]
    net_subsample[i,"fungal_Cscore"]<-network_group_frame["C.score.HL",]
    net_subsample[i,"richness"]<-richness_summed
    net_subsample[i,"Shannon_diversity"]<-Shannon_summed

 }

outfile<- paste0("/gpfs1/work/lachmann/Nets/output/", "net_subsample." ,div, ".", link_mode, ".RDS")

saveRDS(net_subsample, file = outfile) 
















