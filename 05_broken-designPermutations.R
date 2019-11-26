############################### Broken stick design permutations ################

# R script (05/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

###############################

library(dplyr)          # dplyr_0.8.0.1   #%>% 
library(igraph)         # igraph_1.2.2    #graph_from_adjacency_matrix()
library(scales)         # scales_1.0.0
library(phyloseq)       # phyloseq_1.26.1

library(vegan)          #vegan_2.5-4      #nullmodel()
library(TeachingDemos)  #char2seed()      #TeachingDemos_2.10
#library(qgraph)

#sessionInfo()
#R version 3.5.2 (2018-12-20)
###############################

source(file="./R_functions/05F_nullmodel_networks.r")

char2seed('Weissbecker', set = TRUE)
#check what actually is the seed number:
#tmp <- char2seed('Weissbecker',set=FALSE)

# Input data

BEF_sel<-readRDS(file="./output/02_BEF_10r.RDS") # 1926 taxa and 373 samples

BEF_sel_otu<-t(otu_table(BEF_sel))

#reduce to fewer variables to reduce later melting time...
BEF_sel_t<- phyloseq(otu_table(BEF_sel_otu, taxa_are_rows = FALSE), 
                             tax_table(BEF_sel), 
                             sample_data(BEF_sel)[,c("Plot","Tree_s","Div_Plot",
                                                     "Div_local")])

#-----overall null model
BEF_sel_nm<-nullmodel_networks(BEF_sel)

#----select overall dataset manually and also final name in saving the
#    results at the end

overall_data<-BEF_sel_nm    #BEF_sel_t or BEF_sel_nm

pos_new_var<-length(names(sample_data(overall_data)))+1
# Assign running numbers to tree species according to plantation scheme

sample_data(overall_data)[[pos_new_var]]<-rep(0,dim(sample_data(overall_data))[1])
names(sample_data(overall_data))[pos_new_var]<-"Tree_int"

list_1<-which(sample_data(overall_data)[["Tree_s"]]=="Castanea henryi (EcM)")
list_2<-which(sample_data(overall_data)[["Tree_s"]]=="Nyssa sinensis (AM)")
list_3<-which(sample_data(overall_data)[["Tree_s"]]=="Liquidambar formosana (AM)")
list_4<-which(sample_data(overall_data)[["Tree_s"]]=="Sapindus saponaria (AM)")
list_5<-which(sample_data(overall_data)[["Tree_s"]]=="Choerospondias axillaris (AM)")
list_6<-which(sample_data(overall_data)[["Tree_s"]]=="Triadica sebifera (AM)")
list_7<-which(sample_data(overall_data)[["Tree_s"]]=="Quercus serrata (EcM)")
list_8<-which(sample_data(overall_data)[["Tree_s"]]=="Castanopsis sclerophylla (EcM)")
list_9<-which(sample_data(overall_data)[["Tree_s"]]=="Cyclobalanopsis glauca (EcM)")
list_10<-which(sample_data(overall_data)[["Tree_s"]]=="Quercus fabri (EcM)")
list_11<-which(sample_data(overall_data)[["Tree_s"]]=="Rhus chinensis (AM)")
list_12<-which(sample_data(overall_data)[["Tree_s"]]=="Schima superba (AM)")
list_13<-which(sample_data(overall_data)[["Tree_s"]] %in% c("Castanopsis eyrei (EcM)", 
                                                       "Castanopsis carlesii (EcM)"))
list_14<-which(sample_data(overall_data)[["Tree_s"]]=="Cyclobalanopsis myrsinifolia (EcM)")
list_15<-which(sample_data(overall_data)[["Tree_s"]]=="Lithocarpus glaber (EcM)")
list_16<-which(sample_data(overall_data)[["Tree_s"]]=="Koelreuteria bipinnata (AM)")

list_collection<-list(list_1,list_2, list_3, list_4, list_5, list_6, list_7,list_8,
                      list_9, list_10, list_11, list_12,list_13, list_14, list_15, list_16)

j=1
for (i in list_collection){
  sample_data(overall_data)[["Tree_int"]][i]<-j
  j=j+1
}


# Express the experimental design in numbers and tables and derive all possible
# tree species combinations of the lower tree diversity level and for each of these
# one random assignment of tree species data drawn from the higher tree species plots

#----- a) get all 2 species mixture tree combinations of 7 trees, selecting only
#         one tree species of each plot

# species composition on the eight two-species-mixture plots
# one row is a plot and the two columns are the two species mixtures,
# currently represented by ongoing numbers

specs <- cbind((1:8*2)-1,1:8*2) 

# list all possible combinations of chosing 7 tree species from the eight plots
# 2^7 = 128 possible combinations for one subset of 7 plots
# 2^7 * 8 = 1024 possible combinations for all 8 plot combinations of 7 chosen plots
# (2^7 * 2^3= 2^10 = 1024)

# thus creating an empty matrix with 7 cols for the seven target species and 
# and 1024 possible species combinations

#intToBits transforms an integer number to 32 length bit coding(base 2)
# we are only interested in the first 7 bits as equivalent of our 7 plots
#eg. intToBits(1)[1:7]  # 01 00 00 00 00 00 00  # 2^n start with 0--> 2^0 = 1
# intToBits(2)[1:7]     # 00 01 00 00 00 00 00  # 2= 2^1

#as.integer() translates the 00 and 01 into 0 and 1 ...
#ifelse evaluates 0 to FALSE and 1 to TRUE

# if 0, then the first colum of the species coding matrix is chosen
# (so to say 0 as the default selection of the first numbered species on this plot
# and 1 as the second present species on the plot)

#so starting with j=1 and i=0 selects resulting in intoBits= 00 00 00 00 00 00 00
# selects all the first listed species from all plots but the first plot, 
# which is excluded


speciesSelection <- matrix(0,ncol=7,nrow=128*8)
for(j in 1:8){
  for(i in 0:127){
    speciesSelection[(i+1)+(j-1)*128,] <- ifelse(as.integer(intToBits(i)[1:7]),
                                                 specs[-j,][,2], #TRUE: 1
                                                 specs[-j,][,1]) #FALSE:0
   }
}

#----- erase all speciesSelection combinations with tree species number 13!
#as data is insufficient for this tree species!

speciesSelection<-speciesSelection[-(which(speciesSelection == 13, arr.ind=TRUE)[,1]),]
#str(speciesSelection)#576

# access the fungal data for the above calculated data combinations of 2 species mixtures


BEF_2specs<-subset_samples(overall_data, Div_Plot == 2)
selection_data_2<-list()
selection_2_samples<-list()

for( i in 1:nrow(speciesSelection)){
    selection_data_2[[i]]<-otu_table(subset_samples(BEF_2specs, Tree_int %in% speciesSelection[i,]))
    selection_2_samples[[i]]<-subset_samples(BEF_2specs, Tree_int %in% speciesSelection[i,])
    }

# access the fungal data for the above calculated data combinations of monocultures
BEF_mono<-subset_samples(overall_data, Div_Plot == 1)


#--- Start "magic seven" permutations: broken stick data selection

selection_data_1<-list()
selection_1_samples<-list()

for( i in 1:nrow(speciesSelection)){
  selection_data_1[[i]]<-otu_table(subset_samples(BEF_mono, Tree_int %in% speciesSelection[i,]))
  selection_1_samples[[i]]<-subset_samples(BEF_mono, Tree_int %in% speciesSelection[i,])
  }


#----- b) chosing one associated high tree diversity tree species-plot combination
#         for each of the above combinations


#select data

BEF_vierer<-subset_samples(overall_data, Div_Plot == 4)
BEF_achter<-subset_samples(overall_data, Div_Plot == 8)
BEF_16s<-subset_samples(overall_data, Div_Plot == 16)

vierer <- list(c(1:4),c(5:8),c(9:12),c(13:16))
achter <- list(c(1:8),c(9:16))

selection_data_high<-list()
selection_high_samples<-list()

# for all two species mixtures tree combinations deduce the corresponding higher tree 
# species-plot combinations, takes about 3 min


for(i in 1:nrow(speciesSelection)){  
  wievieleVierer <- sapply(vierer, function(x) length(which(x %in% speciesSelection[i,])))
  doppelvierer <- sapply(vierer, function(x) if(length(which(x %in% speciesSelection[i,]))>1) x)
  einzelviererAuswahl <- intersect(speciesSelection[i,],vierer[[which(wievieleVierer==1)]])
  viererAuswahl <- c(einzelviererAuswahl,
                   unlist(sapply(doppelvierer,function(x) if(!is.null(x)) sample(speciesSelection[i,speciesSelection[i,] %in% x],1))))

  speciesSelectionRest <- setdiff(speciesSelection[i,],viererAuswahl)
  wievieleAchter <- sapply(achter, function(x) length(which(x %in% speciesSelectionRest)))
  doppelachter <- sapply(achter, function(x) if(length(which(x %in% speciesSelectionRest))>1) x)
  einzelachterAuswahl <- intersect(speciesSelectionRest,achter[[which(wievieleAchter==1)]])
  achterAuswahl <- c(einzelachterAuswahl,
                   unlist(sapply(doppelachter,function(x) if(!is.null(x)) sample(speciesSelectionRest[speciesSelectionRest %in% x],1))))

  speciesSelectionFinalist <- setdiff(speciesSelectionRest,achterAuswahl)

  # select from OTU table based on this

  vierer_table<-otu_table(subset_samples(BEF_vierer, Tree_int %in% viererAuswahl))
  vierer_samples<-subset_samples(BEF_vierer, Tree_int %in% viererAuswahl)
  
  achter_table<-otu_table(subset_samples(BEF_achter, Tree_int %in% achterAuswahl))
  achter_samples<-subset_samples(BEF_achter, Tree_int %in% achterAuswahl)
  
  max_table<-otu_table(subset_samples(BEF_16s, Tree_int %in% speciesSelectionFinalist))
  max_samples<-subset_samples(BEF_16s, Tree_int %in% speciesSelectionFinalist)
  
  table_4a8<-rbind(vierer_table, achter_table[, colnames(vierer_table)])
  selection_data_high[[i]]<-rbind(table_4a8, max_table[,colnames(table_4a8)])
  
  selection_high_samples[[i]]<-merge_phyloseq(vierer_samples,achter_samples,max_samples)
  
  
}

#----- selection of otu_tables for all combinations

#only OTUs
input_treediv_otumatrices<-list( selection_data_1,selection_data_2,selection_data_high )

#whole phyloseq objects
treediv_samples<-list(selection_1_samples, selection_2_samples,selection_high_samples  )


#--- output 1/2
saveRDS(input_treediv_otumatrices, file = "./output/05_tree_combinations_otu.RDS")
saveRDS(treediv_samples, file = "./output/05_tree_combinations_samples_nm.RDS")
###########




#----- now analyses for frequently occurring fungi: occurrence matrix
# prepare adjacency matrix for network calculations and visualization

occmatMono <- matrix(0,nrow=nrow(speciesSelection),ncol=ntaxa(BEF_sel),
                     dimnames=list(c(1:nrow(speciesSelection)),taxa_names(BEF_sel))) 
occmatLow <- matrix(0,nrow=nrow(speciesSelection),ncol=ntaxa(BEF_sel),
                    dimnames=list(c(1:nrow(speciesSelection)),taxa_names(BEF_sel)))
occmatHigh <- matrix(0,nrow=nrow(speciesSelection),ncol=ntaxa(BEF_sel),
                     dimnames=list(c(1:nrow(speciesSelection)),taxa_names(BEF_sel)))

result_occurrence_matrices<-list(occmatMono,occmatLow,occmatHigh)
names(result_occurrence_matrices)<-c("occmatMono", "occmatLow", "occmatHigh")

result_index=1

for(otu_matrices_collection in input_treediv_otumatrices){  
  print("result_index")
  print(result_index)
  
  for(i in 1:length(otu_matrices_collection)){
     
      single_permutation<-otu_matrices_collection[[i]]
      Plot_info<-sample_data(overall_data)[["Plot"]]
      Plot_info<-data.frame(Plot_info)
      rownames(Plot_info)<-rownames(sample_data(overall_data))

      otu_plot<-merge(single_permutation, Plot_info,by="row.names",all.x=TRUE)

      otu_plot = otu_plot %>% 
        select(Plot_info, everything())

      otu_plot<-otu_plot[,-which(colnames(otu_plot)=="Row.names")]
      
      
      #determine presence of each OTU on each plot 
      #(for a preselected single tree species as determined above)
      otu_plot<-aggregate(data.frame(otu_plot[,-1]), 
                          by=list(as.numeric(otu_plot$Plot_info)), FUN=sum)
      
      #transform to presence/absence
      otu_plot[otu_plot>0]<-1
      otu_occ<-colSums(otu_plot)[-1]
      
      
     result_occurrence_matrices[[result_index]][i,names(otu_occ)]<- otu_occ
    
  }
  
  result_index=result_index+1
  
}

#----- output 2/2
saveRDS(result_occurrence_matrices, file="./output/05_BEF_10r_occurMatrix_allsamples.RDS")



