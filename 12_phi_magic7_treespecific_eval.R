### evaluate the tree specific phi coeffiecient of the magic 7 subsamplings ####

# R script (12/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

#select the positive phi correlation coefficient values

################################################################################

library(dunn.test)
library(pgirmess) # #kruskal wallis test multiple group comparisons, post hoc test

############################


#load the phi coefficient data
phi_magic7_trees<-readRDS(file="./output/11_phi_median_magic7_trees.RDS")
names(phi_magic7_trees)<-c("mono", "low", "high")


tiff(filename ="./images/12_phi_median_positive_boxplot_text.tiff",
     width = 6.8, height = 6.8, units = "in", res=300)

boxplot(phi_magic7_trees$mono$median_positive, 
        phi_magic7_trees$low$median_positive,
        phi_magic7_trees$high$median_positive,
        col=c("green", "blue", "orange"),
        cex.lab=1.5,
        cex.axis=1.5,
        names=c("", "", ""),
        ylab="Median phi coefficient")

axis(side = 1, line = 0, at = c(1, 2, 3), cex.axis=1.5,
     labels = c("Monocultures", "Two tree","High tree"),
     tick = F)
axis(side = 1, line = 1.5, at = c(1, 2, 3), cex.axis=1.5,
     labels = c("", "species", "diversity"), tick = F)

stripchart(phi_magic7_trees$mono$median_positive,vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=1)
stripchart(phi_magic7_trees$low$median_positive,
           vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=2)
stripchart(phi_magic7_trees$low$median_positive,
           vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=3)

text(0.7, 0.06, labels="a", cex=1.5)
text(1.7, 0.057, labels="b", cex=1.5)
text(2.7, 0.059, labels="c", cex=1.5)
text(0.7, 0.067, labels="p<0.05", cex=1.5)
dev.off()




Dataset_phi<-phi_magic7_trees    

result_table<-data.frame(matrix(data=NA, ncol = 5 ,nrow = 3*576,
                     dimnames=list(as.character(1:(3*576)), 
                      c("median_all", "median_absolute",
                        "median_positive", "median_negative",
                        "level"))))

result_table[,5]<-c(rep("monoculture", 576), rep("Low", 576), rep("High", 576))
result_table[,1]<-c(Dataset_phi[[1]]$median_all, Dataset_phi[[2]]$median_all,
                    Dataset_phi[[3]]$median_all)
result_table[,2]<-c(Dataset_phi[[1]]$median_absolut, Dataset_phi[[2]]$median_absolut,
                    Dataset_phi[[3]]$median_absolut)
result_table[,3]<-c(Dataset_phi[[1]]$median_positive, Dataset_phi[[2]]$median_positive,
                    Dataset_phi[[3]]$median_positive)
result_table[,4]<-c(Dataset_phi[[1]]$median_negative, Dataset_phi[[2]]$median_negative,
                    Dataset_phi[[3]]$median_negative)


#median_positive
kruskal_result_median_positive<-kruskal.test(result_table$median_positive ~ as.factor(result_table$level), data = result_table)
#data:  result_table$median_positive by as.factor(result_table$level)
#Kruskal-Wallis chi-squared = 75.082, df = 2, p-value < 2.2e-16

post_hoc_result_median_positive<- dunn.test(as.numeric(as.character(result_table$median_positive)), 
                           result_table$level, method="bh")    
#Bonferroni corrected, but also correct for the 4 different kinds of phi correlation
#coefficient analysis: set p=0.05/4 .... significant p= 0.0125

post_hoc_result_median_positive$P.adjusted
post_hoc_result_median_positive$comparisons
post_hoc_result_median_positive$Z

#"High - Low" : 1.864776e-05        # 4.123625 more than low
#"High - monoculture": 4.258036e-06 #-4.538061 less than mono
#"Low - monoculture" : 6.972525e-18 #-8.661687 less than monoculture


kruskalmc_results<-kruskalmc(result_table$median_positive~result_table$level,
                             data=result_table,
                             probs=0.05/(3))

# Multiple comparison test after Kruskal-Wallis 
# p.value: 0.01666667 
# Comparisons
# obs.dif critical.dif difference
# High-Low         121.2135     81.53049       TRUE
# High-monoculture 133.3958     81.53049       TRUE
# Low-monoculture  254.6094     81.53049       TRUE

