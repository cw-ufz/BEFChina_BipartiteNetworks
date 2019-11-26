########################## Generate link options data matrix ##################

# R function for
# R script (07/) for analyses in Weissbecker et al. 2019
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

##########################
library(plyr)

##########################


#This function calculates different options to set links for co-occuring
#count data of two organisms


CalcLinkOptions<-function(input_data){

  # number of replicates
  Dataset_pre<-input_data
  Dataset_pre$count<-rep(1,nrow(input_data))
  agg_countreps <- aggregate(count~higher+lower+webID,data=Dataset_pre,FUN=sum)
  Dataset_pre_treps<-merge(Dataset_pre, agg_countreps, by=c("higher", "lower", "webID"))
  
  # transform to incidence data
  Dataset_pre_treps$pos_reps<-Dataset_pre_treps$freq
  Dataset_pre_treps$pos_reps[Dataset_pre_treps$freq>0]<-1

  # average sequence read numbers
  freqs_mean<-ddply(Dataset_pre_treps, .(higher, lower, webID), summarize,  
                    freq=mean(freq))
  freqs_median<-ddply(Dataset_pre_treps, .(higher, lower, webID), summarize,  
                      freq=median(freq)) 
  Dataset_pre_treps_mean<-merge(Dataset_pre_treps, freqs_mean, by=c("higher", "lower", "webID"))
  Dataset_pre_treps_median<-merge(Dataset_pre_treps_mean, freqs_median, by=c("higher", "lower", "webID"))
 
    # number of replicates with observed co-occurrence fungal_OTU-tree
  agg_posreps <- aggregate(pos_reps~higher+lower+webID,data=Dataset_pre_treps_median,FUN=sum)
  BEF_sel_posreps<-merge(Dataset_pre_treps_median, agg_posreps, by=c("higher", "lower", "webID"))
  names(BEF_sel_posreps)<-c("higher","lower","webID","freq", "number ","reps_tot","pa" ,"mean_read", "median_read","reps_pos")
  
  # select appropriate columns in nice order
  BEF_sel_prelink<-BEF_sel_posreps[,c("higher","lower","webID","freq","mean_read", "median_read","pa", "reps_tot","reps_pos")]
  
  # Number of links:
  # min - minimum number of links - all samples show this co-occurrence
  # mean - half of the samples show this co-occurrence
  # max - maximum number of links: the co-occurrence was observed at least once
  
  #eg 0.3 for 30% of samples need to show co-occurrence to be considered in network analysis
  custom<-0.7
 
  BEF_sel_prelink$min_link<-BEF_sel_prelink$reps_pos==BEF_sel_prelink$reps_tot
  BEF_sel_prelink$mean_link<-BEF_sel_prelink$reps_pos>=0.5*BEF_sel_prelink$reps_tot
  BEF_sel_prelink$max_link<-BEF_sel_prelink$reps_pos>0
  BEF_sel_prelink$custom_link<-BEF_sel_prelink$reps_pos>=custom*BEF_sel_prelink$reps_tot
  BEF_sel_prelink$p20<-BEF_sel_prelink$reps_pos>=0.2*BEF_sel_prelink$reps_tot
  BEF_sel_prelink$p40<-BEF_sel_prelink$reps_pos>=0.4*BEF_sel_prelink$reps_tot
  BEF_sel_prelink$p60<-BEF_sel_prelink$reps_pos>=0.6*BEF_sel_prelink$reps_tot
  BEF_sel_prelink$p80<-BEF_sel_prelink$reps_pos>=0.8*BEF_sel_prelink$reps_tot
  BEF_sel_prelink$p100<-BEF_sel_prelink$reps_pos>=1*BEF_sel_prelink$reps_tot
  
  
  return(BEF_sel_prelink)
  
}
