################### phi coefficient heatmap calculation ########################
#
# R script (14/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

################################################################################

library(gplots)         #gplots_3.0.1.1
library(colorspace)     #colorspace_1.4-0 
library(dendextend)     #dendextend_1.9.0

#sessionInfo()

#########################

phi_tree_otu<-readRDS(file="./output/10_phi_coeff_trees_otus.RDS")

#try heatmap with top 200max phi value OTUs:
phi_maxes<-apply(phi_tree_otu, MARGIN=2, max)
phi_maxes_ordered<-sort(phi_maxes, decreasing=TRUE)
phi_top200_list<-rownames(data.frame(phi_maxes_ordered[1:200]))

phi_top200<-phi_maxes_ordered[1:200]
#range(phi_top200) #0.2335616 to 0.5694651

phi_top20_list<-rownames(data.frame(phi_maxes_ordered[1:20]))

#-----output

saveRDS(phi_top20_list, file="./output/14_phi_top20_otulist.RDS")

###############


phi_top200<-phi_tree_otu[,phi_top200_list]

Dataset_phi_heatmap<-phi_top200      

phi_pre_heatmap<-t(Dataset_phi_heatmap)

d_occurr_phi<-dist(phi_pre_heatmap, method="euclidean", diag=FALSE)
hc_occurr_phi<-hclust(d_occurr_phi, method = "complete")

dend_occurr_phi <- as.dendrogram(hc_occurr_phi)


# Color the branches 
dend_occurr_phi <- color_branches(dend_occurr_phi, k=1, col="darkgrey")

some_col_func <- function(n) rev(colorspace::diverging_hsv(n, h = c(0, 260), s = 2, v = 1.3, power = 1,
                                                           gamma = NULL, fixup = TRUE, alpha = 1))

positions <- c("Cyclobalanopsis glauca (EcM)","Cyclobalanopsis myrsinifolia (EcM)",
               "Castanopsis sclerophylla (EcM)",
               "Castanopsis carlesii (EcM)",
               "Quercus serrata (EcM)","Quercus fabri (EcM)",
               "Castanea henryi (EcM)","Lithocarpus glaber (EcM)",
               "Choerospondias axillaris (AM)","Rhus chinensis (AM)",
               "Koelreuteria bipinnata (AM)","Sapindus saponaria (AM)",
               "Nyssa sinensis (AM)","Liquidambar formosana (AM)",
               "Triadica sebifera (AM)","Schima superba (AM)")


phi_pre_heatmap_ordered<-phi_pre_heatmap[,positions]

tiff(filename ="./images/25_phi_heatmap_treespecies_top200_2.tiff",
#     width = 10, height = 6.8, units = "in", res=300)
      width = 10, height = 6.8, units = "in", res=300)


#par(mar=c(5.1, 4.1, 4.1, 2.1))
#sets the margin sizes in the following order: bottom, left, top, and right

gplots::heatmap.2(as.matrix(phi_pre_heatmap_ordered), 
                  main = "",
                  srtCol = 16.5,
                  dendrogram = "row",
                  Rowv = dend_occurr_phi,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(5,0.1),      
                  key.xlab = "Phi coefficient",
                  denscol = "black",
                  density.info = "density",
                  #RowSideColors = rev(labels_colors(dend_occurr[[1]])), # to add nice colored strips        
                  col = some_col_func
                
)


dev.off()


















