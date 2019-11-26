############################### Occurrence evaluation ##########################

# R script (06/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

#compare the number and identity of identified frequently occuring fungi
#across the diversity levels

###############################

read_occurmatrices<-readRDS(file="./output/05_BEF_10r_occurMatrix_allsamples.RDS")

fungal_taxonomy<-read.csv2(file="./Data/BEF_vipA_454_ITS_2011_curated_taxonomy_guild.csv",
                           row.names = 1)

frequents_magic7<-list()

for (j in 1:length(read_occurmatrices)){
    
  max_results<-apply(read_occurmatrices[[j]], MARGIN = 2, max)
  frequent_list<-which(max_results == 7)
  frequents_magic7[[names(read_occurmatrices)[j]]]<-names(max_results[frequent_list]) 
     
}


length(frequents_magic7$occmatMono) #63 (26 unique)
length(frequents_magic7$occmatLow)  #47 (11 unique)
length(frequents_magic7$occmatHigh) #29 (3 unique)

# length(common_all)                  #15 common frequent fungal OTUs across
#                                     all three diversity levels


#----- output 1/3
saveRDS(frequents_magic7, file="./output/06_BEF10r_frequentsMagic7_allsamples.RDS")





uniq_high<-setdiff(frequents_magic7$occmatHigh,
                   unique(c(frequents_magic7$occmatMono, frequents_magic7$occmatLow)))
#3

uniq_mono<-setdiff(frequents_magic7$occmatMono,
                   unique(c(frequents_magic7$occmatHigh, frequents_magic7$occmatLow)))
#26
uniq_low<-setdiff(frequents_magic7$occmatLow,
                  unique(c(frequents_magic7$occmatHigh, frequents_magic7$occmatMono)))
#11

common_low<-intersect(frequents_magic7$occmatMono, frequents_magic7$occmatLow)
common_all<-intersect(common_low, frequents_magic7$occmatHigh)

length(common_all)

all_OTUs<-length(unique(c(frequents_magic7$occmatMono, frequents_magic7$occmatLow,
                          frequents_magic7$occmatHigh)))
#82 frequently occurring OTUs (appearing in at least one diversity level)
#

frequent_table<-matrix(data=c(length(frequents_magic7$occmatMono),
                                length(uniq_mono),
                                length(frequents_magic7$occmatLow),
                                length(uniq_low),
                                length(frequents_magic7$occmatHigh),
                                length(uniq_high),
                                all_OTUs, 
                                NA),
                         ncol=4, 
                         nrow = 2,
                         dimnames = list(c("Total number frequents", "Number unique frequents"),
                                         c("Monocultures", "Di_div", "High_div", "All")))

#----- output 2/3
write.csv(frequent_table, file="./output/06_frequent7_counts_allsamples.csv")





frequent_sets<-list(uniq_mono, uniq_low, uniq_high, common_low, common_all)
frequent_titles<-c("Unique_Mono", "Uniq_Low", "Uniq_High", "Common_Low", "Common_all")
  
frequent_taxonomy<-list()

for(i in 1:length(frequent_sets)){

  taxonomy_pos<-which(rownames(fungal_taxonomy) %in% frequent_sets[[i]])
  frequent_taxonomy[[frequent_titles[i]]]<-fungal_taxonomy[taxonomy_pos, 2:7]


}
View(frequent_taxonomy$Common_all)

#----- output 3/3
saveRDS(frequent_taxonomy, file="./output/06_taxonomy_frequent7_levels_allsamples.RDS")
write.csv2(frequent_taxonomy$Common_all, file="./output/06_taxonomy_frequent7_levels_allsamples.csv")



#which counts do unique Mono/Low/High frequents have in the other tree
#diversity treamtents?

#need to select the unique dataset and the dataset to compare to manually!!!

  gen_cross_selection<-which(colnames(read_occurmatrices$occmatMono) %in% uniq_high)
  
  occ_mat<-read_occurmatrices$occmatMono[,gen_cross_selection]
  
  no_count_OTU<-which(apply(occ_mat, MARGIN=2, max)==0)
  no_count_OTU
  one_count_OTU<-which(apply(occ_mat, MARGIN=2, max)==1)
  one_count_OTU
  two_count_OTU<-which(apply(occ_mat, MARGIN=2, max)==2)
  two_count_OTU
  three_count_OTU<-which(apply(occ_mat, MARGIN=2, max)==3)
  three_count_OTU
  
  View(occ_mat)
  
#Low   
#unique frequents mono: OTU418, OTU949 appear only 2 times in the two species mixtures

#High
#unique frequent mono: OTU418 appear only 2 times in the high mixtures
#unique frequent low: OTU616 appear only 2 times in the high mixtures...

  
#OTU418:Saprotroph -  Fungi -  Ascomycota -   Sordariomycetes -  Chaetosphaeriales
#  Chaetosphaeriaceae -  Chaetosphaeria

#OTU616: Saprotroph -  Fungi -  Ascomycota -  Archaeorhizomycetes -  Archaeorhizomycetales
#  Archaeorhizomycetaceae -  Archaeorhizomyces  
   
#OTU949: Saprotroph -  Fungi -  Ascomycota - Eurotiomycetes - Eurotiales
#  Aspergillaceae - Penicillium

   
  
tiff(filename ="./images/06_BEF_10r_hist_highunique_in_Mono.tiff",
       width = 6.8, height = 6.8, units = "in", res=300)
  
   hist_result<-hist(occ_mat, plot=TRUE)
 
dev.off()  

