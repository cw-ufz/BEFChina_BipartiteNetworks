###########calculate median phi coefficient for magic7 tree based ##############

# R script (11/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

################################################################################

library(phyloseq)
library(plyr)


#sessionInfo() 
###############################


# get the calculated phi coefficients for each tree species for each fungal OTU
phi_tree_otu<-readRDS(file="./output/10_phi_coeff_trees_otus.RDS")

#get the magic 7 combinations
combination_samples<-readRDS(file = "./output/05_tree_combinations_samples.RDS")

# [[1]] : Tree monocultures subsample combinations
# [[2]] : Two tree species subsample combinations
# [[3]] : High tree diversity subsample combinations 
#(4,8,16 diversity plots summed)
# [[]][[1]] - [[]][[576]] subsample combinations a tree diversity level

phi_median_mono<-data.frame(matrix(data = NA, nrow = 576, ncol = 4,
                           dimnames=list(1:576,
                                  c("median_all", "median_absolut",
                                    "median_positive", "median_negative"))))

phi_median_di<-phi_median_mono
phi_median_high<-phi_median_mono

phi_median<-list(phi_median_mono,phi_median_di,phi_median_high)


j=1
i=1
result<-c(1,2,3)

for(j in 1:3){
  print("result:")
  print(result[j])  
  for(i in 1:576){
    #reduce otu table
    pruned_set<-prune_taxa(taxa_sums(combination_samples[[j]][[i]])>0 , combination_samples[[j]][[i]])
    
    env_data<-data.frame(sample_data(pruned_set))
    otu_data<-data.frame(otu_table(pruned_set))
    
    merge_data<-merge(env_data,otu_data, by="row.names") 
    merge_data_reduce<-merge_data[,c(grep('^Tree_s', colnames(merge_data)),grep('^Otu', colnames(merge_data)))]                       
    merge_data_aggregate<-aggregate(merge_data_reduce[,-1], by=list(merge_data_reduce$Tree_s), FUN=sum)                          
                   
    merge_data_aggregate[merge_data_aggregate == 0] <- NA
    
    View(qhi_reduced)
    
    qhi_reduced<-phi_tree_otu[which(rownames(phi_tree_otu)%in% merge_data_aggregate$Group.1),
                 which(colnames(phi_tree_otu)%in%colnames(merge_data_aggregate))]
    
    put_NA<-which(is.na(merge_data_aggregate)==TRUE)
    NA_list<-apply(merge_data_aggregate[,-1], 1, function(x) which(is.na(x)==TRUE))
   
    for(k in 1:length(NA_list)){
        qhi_reduced[k,names(NA_list[[k]])]<-NA
    }  

    #now I completed the phi coefficient info for the current magic7 subsampling data
    #continue with the calculation of the median values after 3 different approaches:
     
    
    #i) calculate the median of all values, excluding the NAs
    phi_median[[j]][[i,"median_all"]]<-median(qhi_reduced, na.rm = TRUE)
    
    #ii) calculate the median of the absolute values 
    phi_median[[j]][[i,"median_absolut"]]<-median(abs(qhi_reduced), na.rm = TRUE)
    
    #iii)a) calculate the median of the positive values
    phi_positive<-which(as.vector(qhi_reduced)>0)
    phi_median[[j]][[i,"median_positive"]]<-median(as.vector(qhi_reduced)[phi_positive])
    
    #iii)b) calculate the median of the negative values
    phi_negative<-which(as.vector(qhi_reduced)<0)
    as.vector(qhi_reduced)[phi_negative]
    phi_median[[j]][[i,"median_negative"]]<-median(as.vector(qhi_reduced)[phi_negative])
   
    
  }
  
}

saveRDS(phi_median, file="./output/11_phi_median_magic7_trees.RDS")







