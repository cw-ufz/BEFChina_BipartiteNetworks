########### Calculate selected network topological features for organism groups#

# R function for
# R script (07/) for analyses in Weissbecker et al. 2019
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

##########################

#This function calculates selected group level network tolopiological features from 
# frame2web object as input data
library("bipartite") 

NetworkTopologiesGroups<-function(input_data){
  
  BEF_group_topologies <- sapply(input_data,
                                 function(x) grouplevel(x,
                                index = c("mean number of shared partners", "togetherness",
                               "C score","niche overlap",
                               "generality")))
  
  return(BEF_group_topologies)
  
}
