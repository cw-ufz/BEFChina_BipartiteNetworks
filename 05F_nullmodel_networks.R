########################## Nullmodel calculation for networks ##################

# R function for
# R script (05/) for analyses in Weissbecker et al. 2019
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

##########################

#This function calculates nullmodels of input data, count data
#with the implemented method "r2dtable", count sums for each OTU and each
#sample are preserved.


nullmodel_networks<-function(input_data){
  if(class(input_data) == 'phyloseq') input_otu<-data.frame(t(otu_table(input_data)))
  if(class(input_data) == 'data.frame') input_otu<-input_data
  
  (nm_input <- vegan::nullmodel(input_otu, method = "r2dtable"))
  sm_input <- matrix(simulate(nm_input,nsim=1),nrow=nrow(input_otu),ncol=ncol(input_otu))
  
  rownames(sm_input)<-rownames(input_otu)
  colnames(sm_input)<-colnames(input_otu)
  
  if(class(input_data) == 'phyloseq') {
    
  output_phylo<- phyloseq(otu_table(sm_input, taxa_are_rows = FALSE), 
                          tax_table(input_data), sample_data(input_data))
  
  return(output_phylo)
  }
  
  if(class(input_data) == 'data.frame') return(sm_input)

}


