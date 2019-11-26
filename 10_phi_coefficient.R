####Calculate phi coefficient - degree of specialization to a tree species #####
#
# R script (10/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25


################################################################################

library(phyloseq)
library(vegan)

###################################

source(file="./R_functions/05F_nullmodel_networks.r")

# Define the phi coeffient based on vectors with abundances
# The phi coefficient is based on presence/absence data, 

phi.coeff <- function(x,y){
  # x is the vector of absolute frequencies of the species in a subset, for which 
  #    the degree of specialisation is calculated
  # y is the vector of absolute frequencies of species in total
  sx <- sum(x)
  sy <- sum(y)
  if (sx!=0 & sy!=0){
    phi <- unlist(sapply(1:length(x),function(i){
      a <- as.double(x[i])
      b <- as.double(sx-x[i])
      c <- as.double(y[i]-x[i])
      d <- as.double(sy-y[i]-(sx-x[i]))
      (a*d-b*c)/sqrt((a+c)*(b+d)*(a+b)*(c+d))
    }))
    return(phi)
  } else {
    return(rep(0,length(x)))
  }
}


#load the data
BEF_sel<-readRDS(file="./output/02_BEF_10r.RDS") # 1926 taxa and 373 samples

BEF_otu<-data.frame(t(otu_table(BEF_sel)))
BEF_otu_pa<-decostand(BEF_otu, method="pa")

#get a column with the tree species names
sample_info<-data.frame(sample_data(BEF_sel))
BEF_otu$Tree_s<-as.factor(sample_info$Tree_s)

BEF_otu_pa_nm_r2dtable<-nullmodel_networks(BEF_otu_pa)
#currently method "r2dtable"




#manually select dataset
Dataset<- BEF_otu_pa #or BEF_otu_pa or BEF_otu_pa_nm_r2dtable


tree.classes <- as.character(sample_info$Tree_s)
otu_x1 <- aggregate(Dataset, by=list(tree.classes), FUN=sum)


otus_onetree<-colnames(otu_x1[which(colSums(otu_x1[,-1])==1)])

one_tree_list<-which(colnames(BEF_otu) %in% otus_onetree)
one_sample_otus<-names(which(colSums(BEF_otu_pa[,one_tree_list])==1))
#[1] "Otu00305" "Otu00403" "Otu00595" "Otu00927" "Otu00932" "Otu02176"


dimnames(otu_x1)[[1]] <- otu_x1[,1]
otu_x1 <- otu_x1[,-1]
otu_y1 <- table(tree.classes)


result_otu <- matrix(NA,nrow=length(unique(tree.classes)), ncol=ncol(otu_x1))
dimnames(result_otu)[[1]] <- unique(tree.classes)
#head(View(result_otu))
result_otu <- apply(otu_x1,2,function(x) phi.coeff(x,otu_y1))
dimnames(result_otu)[[1]] <- unique(tree.classes)
View(result_otu)



#----- output
saveRDS(result_otu, file="./output/10_phi_coeff_trees_otus.RDS")


