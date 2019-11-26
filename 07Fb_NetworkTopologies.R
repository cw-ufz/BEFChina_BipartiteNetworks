########################## Calculate selected network topological features######

# R function for
# R script (07/) for analyses in Weissbecker et al. 2019
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

##########################
library("bipartite") 
#This function calculates selected network tolopiological features from 
# frame2web object as input data


NetworkTopologies<-function(input_data){

  
  # BEF_topologies<-sapply(input_data, networklevel,
  #                        index=c("connectance",
  #                                "weighted connectance",
  #                                "number of species"))
 
  BEF_topologies<-sapply(input_data, networklevel,
                         index=c("nestedness",
                                 "NODF"))
  
 
  
  return(BEF_topologies)
  
}
