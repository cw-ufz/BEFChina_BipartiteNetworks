############## summarize bipartite tables and do anova statistics ##############
#
# R script (09/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

#real data NODF always lower (less nested) than nullmodel NODF
################################################################################

library(pgirmess)

########################


#manually select result files of either real or simulated (nm) data:
Linkages_p20<-lapply(1:3,
                    # function(x) readRDS(file = paste0("./output/net_subsample.", x, ".p20_all.RDS")))
                    function(x) readRDS(file = paste0("./output/net_subsample.nm.", x, ".p20.RDS")))
names(Linkages_p20)<-c("Mono", "Low", "High")
Linkages_p40<-lapply(1:3,
                     #function(x) readRDS(file = paste0("./output/net_subsample.", x, ".p40_all.RDS")))
                    function(x) readRDS(file = paste0("./output/net_subsample.nm.", x, ".p40.RDS")))
names(Linkages_p40)<-c("Mono", "Low", "High")
Linkages_p60<-lapply(1:3,
                    # function(x) readRDS(file = paste0("./output/net_subsample.", x, ".p60_all.RDS")))
                    function(x) readRDS(file = paste0("./output/net_subsample.nm.", x, ".p60.RDS")))
names(Linkages_p60)<-c("Mono", "Low", "High")
Linkages_p80<-lapply(1:3,
                    # function(x) readRDS(file = paste0("./output/net_subsample.", x, ".p80_all.RDS")))
                    function(x) readRDS(file = paste0("./output/net_subsample.nm.", x, ".p80.RDS")))
names(Linkages_p80)<-c("Mono", "Low", "High")
Linkages_p100<-lapply(1:3,
                     # function(x) readRDS(file = paste0("./output/net_subsample.", x, ".p100_all.RDS")))
                      function(x) readRDS(file = paste0("./output/net_subsample.nm.", x, ".p100.RDS")))
names(Linkages_p100)<-c("Mono", "Low", "High")


head(View(Linkages_p60[[1]]))
#short boxplot plotting of fungal richness and fungal Shannon diversity

tiff(filename ="./images/09_boxplot_richness.tiff",
     width = 6.8, height = 6.8, units = "in", res=300)

boxplot(Linkages_p20[[1]]$richness, 
        Linkages_p20[[2]]$richness,
        Linkages_p20[[3]]$richness,
        col=c("green", "blue", "orange"),
        names=c("Monoculture", "Two tree species", "High tree diversity"),
        ylab="Fungal richness")
stripchart(Linkages_p20[[1]]$richness, vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=1)
stripchart(Linkages_p20[[2]]$richness, vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=2)
stripchart(Linkages_p20[[3]]$richness, vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=3)

dev.off()




tiff(filename ="./images/09_boxplot_Shannon.tiff",
     width = 6.8, height = 6.8, units = "in", res=300)

boxplot(Linkages_p20[[1]]$Shannon_diversity, 
        Linkages_p20[[2]]$Shannon_diversity,
        Linkages_p20[[3]]$Shannon_diversity,
        col=c("green", "blue", "orange"),
        names=c("Monoculture", "Two tree species", "High tree diversity"),
        ylab="Fungal Shannon diversity")
stripchart(Linkages_p20[[1]]$Shannon_diversity, vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=1)
stripchart(Linkages_p20[[2]]$Shannon_diversity, vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=2)
stripchart(Linkages_p20[[3]]$Shannon_diversity, vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=3)


dev.off()



datasets_to_merge<-list(Linkages_p20,Linkages_p40,Linkages_p60,Linkages_p80,Linkages_p100)
names(datasets_to_merge)<-c("p20", "p40", "p60", "p80","p100")


for (i in 1:5){
  for(j in 1:3){
    compare_data[i,j]<-max(datasets_to_merge[[i]][[j]]$number_OTUs)
   }
}


# Statsitic evaluation between diversity levels: 


for(i in 1:length(datasets_to_merge)){

  datasets_to_merge[[i]][[1]]$div<-1
  datasets_to_merge[[i]][[2]]$div<-2
  datasets_to_merge[[i]][[3]]$div<-3
}

data_merged<-lapply(1:5, function(x) rbind(datasets_to_merge[[x]][[1]],
                                           datasets_to_merge[[x]][[2]],
                                           datasets_to_merge[[x]][[3]]))
names(data_merged)<-c("p20", "p40", "p60", "p80","p100")


#visualize the data  
plot(data_merged$p20$Shannon_diversity)
plot(data_merged$p20$richness)
boxplot(data_merged$p20$Shannon_diversity~data_merged$p20$div)
# data ist not normally distributed and also later residuals are not!...


#test for statistical differences --> non-paramteric tests

#create a result table
net_stats<-matrix(data = NA, nrow=5*3, ncol=12, 
                  dimnames = list(c("p20", "p20.median:mono_low_high", "p20.groups:1-2.1-3.2-3",
                                  "p40", "p40.median:mono_low_high", "p40.groups:1-2.1-3.2-3",
                                  "p60", "p60.median:mono_low_high", "p60.groups:1-2.1-3.2-3",
                                  "p80", "p80.median:mono_low_high", "p80.groups:1-2.1-3.2-3",
                                  "p100", "p100.median:mono_low_high", "p100.groups:1-2.1-3.2-3"),
                                  c(colnames(data_merged$p20)[-13])))
 #                                 c(colnames(data_merged$p20)[-5])))

i=1
j=1

for(i in 1:5){
  for(j in 1:12){
#  for(j in 1:4){
    print(j)
    col_sel<-colnames(data_merged[[i]])[j]
    lm_model<-lm(data_merged[[i]][[col_sel]]~data_merged[[i]]$div)
    kruskal_result<-kruskal.test(data_merged[[i]][[col_sel]]~div, data=data_merged[[i]])

    
    if(is.na(kruskal_result$p.value)==TRUE){
      p_kruskal<-"NA"
    }else if(kruskal_result$p.value<0.001){
      p_kruskal<-"<0.001"
    }else{
      p_kruskal<-round(kruskal_result$p.value, digits=4)
      }

    row_position<-c(1,4,7,10,13)
    net_stats[row_position[i],j]<-p_kruskal

    median_results<-aggregate(data_merged[[i]][[col_sel]]~data_merged[[i]]$div, FUN=median)
    round(median_results$`data_merged[[i]][[col_sel]]`,2)
    median_results_collapsed<-paste(c(round(median_results$`data_merged[[i]][[col_sel]]`,2)), collapse="_")
  
    net_stats[(row_position[i]+1),j]<-median_results_collapsed
  
  
    if(is.na(kruskal_result$p.value)==TRUE){
      net_stats[(row_position[i]+2),j]<-"no difference"
    }else if(kruskal_result$p.value<0.05){

      kruskalmc_results<-kruskalmc(data_merged[[i]][[col_sel]]~div, data=data_merged[[i]],
                                 probs=0.05/(12*5))
   
      kruskalmc_collapsed<-paste(kruskalmc_results$dif.com$difference, collapse="_")
      net_stats[(row_position[i]+2),j]<-kruskalmc_collapsed
    
      #Multiple comparison test after Kruskal-Wallis 
      #p.value: 0.05 
      #Comparisons
      #obs.dif critical.dif difference
      #1-2 322.1736     70.38871       TRUE
      #1-3 594.2153     70.38871       TRUE
      #2-3 916.3889     70.38871       TRUE
    }else{
      net_stats[(row_position[i]+2),j]<-"no difference"
    } 
  }
}


#----- output
saveRDS(net_stats, file = "./output/09_netstats.nestedness.RDS")
write.csv2(net_stats, file = "./output/09_netstats.nestedness.csv")
#-----



#----- Statistic evaluation nestedness: observed vs. null model


#the comparisons need to be made for each linkage threshold and each diversity level 
#separately!

datasets_to_merge<-list(Linkages_p20,Linkages_p40,Linkages_p60,Linkages_p80,Linkages_p100)
names(datasets_to_merge)<-c("p20", "p40", "p60", "p80","p100")

net_stats.nm<-matrix(data = NA, nrow=5*3, ncol=2, 
                  dimnames = list(c("p20.null-mono,null-low.null-high", "p20.median:mono_low_high", "p20.nm.median:mono_low_high",
                                    "p40.null-mono,null-low.null-high", "p40.median:mono_low_high", "p40.nm.median:mono_low_high",
                                    "p60.null-mono,null-low.null-high", "p60.median:mono_low_high", "p60.nm.median:mono_low_high",
                                    "p80.null-mono,null-low.null-high", "p80.median:mono_low_high", "p80.nm.median:mono_low_high",
                                    "p100.null-mono,null-low.null-high", "p100.median:mono_low_high", "p100.nm.median:mono_low_high"),
                                 c("nestedness", "NODF")))

#we perform a two sided rank test : Mann Whitney Wilcoxon
# we adjust the p value (confidence level) to 
# the total of 30 individual tests: p=0.05/30 = 0.001666667
#--> confidence level = 0.9983333


p_wilkox_list<-list()
median_results_metric<-list()
median_results_metric_nm<-list()


for(metric in c("nestedness", "NODF")){
  
  for(i in 1:5){    #link thresholds
    
    for (div in c("Mono","Low","High")){
    
      null_metric<-paste0(metric,"_nm")
      lm_model<-lm(unlist(datasets_to_merge[[i]][[div]][metric])~
                   unlist(datasets_to_merge[[i]][[div]][null_metric]))
    
      p<-0.05/(3*5*2) #p correction for multiple comparisons... 0.001666667
    
      wilkox_result<-wilcox.test(unlist(datasets_to_merge[[i]][[div]][metric]),
                unlist(datasets_to_merge[[i]][[div]][null_metric]),
                alternative = "two.sided",
                conf.level=1-p)
  
    
      wilkox_result$p.value
   
      if(is.na(wilkox_result$p.value)==TRUE){
       p_wilkox<-"NA"
      }else if(wilkox_result$p.value<0.001){
       p_wilkox<-"<0.001"
      }else{
        p_wilkox<-"n.s"
      }
    
      p_wilkox_list[div]<-p_wilkox
    
    
      median_results_metric[div]<-round(median(unlist(datasets_to_merge[[i]][[div]][metric])),2)
      median_results_metric_nm[div]<-round(median(unlist(datasets_to_merge[[i]][[div]][null_metric])),2)
    
    }
  
  
  counter<-c(1,4,7,10,13)  
  line=counter[i]
    
  net_stats.nm[line,metric]<-paste(c(p_wilkox_list["Mono"], p_wilkox_list["Low"], p_wilkox_list["High"]), 
        collapse="_")
  
  net_stats.nm[line+1,metric]<-paste(c(median_results_metric["Mono"], 
                                       median_results_metric["Low"], 
                                       median_results_metric["High"]), 
                                     collapse="_")
  
  net_stats.nm[line+2,metric]<-paste(c(median_results_metric_nm["Mono"], 
                                       median_results_metric_nm["Low"], 
                                       median_results_metric_nm["High"]), 
                                     collapse="_")
  
  
  }

}


#----- output
saveRDS(net_stats.nm, file = "./output/09_net_stats.nestedness.nm.RDS")
write.csv2(net_stats.nm, file = "./output/09_net_stats.nestedness.nm.csv")

################################################################################

boxplot(unlist(datasets_to_merge[[1]][["Mono"]]["NODF"]), 
        unlist(datasets_to_merge[[1]][["Mono"]]["NODF_nm"]),
        names=c("NODF","NODF_nm"),
        ylab="NODF (nestedness)")


# boxplot(unlist(datasets_to_merge[[3]][["Mono"]]["nestedness"]), 
#         unlist(datasets_to_merge[[3]][["Mono"]]["nestedness_nm"]),
#         names=c("nestedness","nestedness_nm"),
#         ylab="NODF (nestedness)")







