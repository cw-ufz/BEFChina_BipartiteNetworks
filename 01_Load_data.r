############################### Load data #####################################

# R script (01/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

#########################

library(phyloseq)         #phyloseq_1.26.0        #merge otu,tax,metadata info
library(vegan)            #vegan_2.5-4            #diversity calculations
library(data.table)       #data.table_1.12.0      #dcast function
library(biomformat)       #biomformat_1.10.0       #read biom file
library(BiocManager)      #BiocManager_1.30.4
#sessionInfo()
#R version 3.5.2 (2018-12-20)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)

######################### install biomformat package
#BiocManager::install(c("phyloseq", "biomformat", "rhdf5"))
######################### 

#----- Load otu_table, prune to dominant OTU dataset (reads>=4)

x = read_biom("./Data/BEF_vipA_454_ITS_1283866816.json_Nov2017.biom") 
otu_matrix<-as(biom_data(x), "matrix")                               
BEF_all_otu<- phyloseq(otu_table(otu_matrix, taxa_are_rows = TRUE))
BEF_dom_otu = prune_taxa(taxa_sums(BEF_all_otu) > 3, BEF_all_otu)

#----- Load taxonomic information and fungal group assignment, 
#      create phyloseq object

fungal_taxonomy<-read.csv2(file="./Data/BEF_vipA_454_ITS_2011_curated_taxonomy_guild.csv",
                           row.names = 1)
taxonomy_fungi_pre<-as.matrix(fungal_taxonomy)
row.names(taxonomy_fungi_pre)<-rownames(fungal_taxonomy)
taxonomy_fungi<-taxonomy_fungi_pre 

BEF_dom<- phyloseq(otu_table(BEF_dom_otu, taxa_are_rows = TRUE),
                   tax_table(taxonomy_fungi))

#----- Load metadata

Env_var<-read.csv("./Data/BEF_vipA_2011_soil_typo_metadata.csv",row.names=1, 
                  dec=",", sep=";") 
AddEnv_var<-read.csv("./Data/AddEnvData.csv",row.names=1, 
                  dec=",", sep=";") 

#correct typos in the Tree species names !
# Add the myco tree information and check the correct spelling of the trees!
AddEnv_var$Tree_s<-gsub("Castanea_henryi", "Castanea henryi (EcM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Quercus_fabri", "Quercus fabri (EcM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Cyclobalanopsis_myrsinaefolia", "Cyclobalanopsis myrsinifolia (EcM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Lithocarpus_glaber", "Lithocarpus glaber (EcM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Castanopsis_sclerophylla", "Castanopsis sclerophylla (EcM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Cyclobalanopsis_glauca", "Cyclobalanopsis glauca (EcM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Quercus_serrata", "Quercus serrata (EcM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Castanopsis_fargesii", "Castanopsis carlesii (EcM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Castanopsis_eyrei", "Castanopsis eyrei (EcM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Nyssa_sinensis", "Nyssa sinensis (AM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Liquidambar_formosana", "Liquidambar formosana (AM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Koelreuteria_bipinnata","Koelreuteria bipinnata (AM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Schima_superba","Schima superba (AM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Sapindus_mukorossi","Sapindus saponaria (AM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Sapium_sebiferum", "Triadica sebifera (AM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Choerospondias_axillaris", "Choerospondias axillaris (AM)", AddEnv_var$Tree_s)
AddEnv_var$Tree_s<-gsub("Rhus_chinensis", "Rhus chinensis (AM)", AddEnv_var$Tree_s)
unique(AddEnv_var$Tree_s)


#----- Assure correct data order
reference<-t(otu_matrix)
var_list<-list(Env_var,AddEnv_var)
names(var_list)<-list("Env_var","AddEnv_var")
all_var_sorted<-lapply(var_list, function(x) x[match(rownames(reference),
                                                     rownames(x)),]) 

#----- Final phyloseq object

BEF_dom_addmeta<- phyloseq(otu_table(BEF_dom_otu, taxa_are_rows = TRUE), 
                           tax_table(taxonomy_fungi), 
                           sample_data(cbind(all_var_sorted$Env_var,
                                             all_var_sorted$AddEnv_var[-2])))
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6693 taxa and 394 samples ]
#sample_data() Sample Data:       [ 394 samples by 22 sample variables ]
#tax_table()   Taxonomy Table:    [ 6693 taxa by 9 taxonomic ranks ]


#----- script output

saveRDS(BEF_dom_addmeta,
       file="./Cleaned_Data/01_dominant4_phyloseq.RDS")
saveRDS(all_var_sorted,
       file="./Cleaned_Data/01_dominant4_metadata.RDS")
saveRDS(BEF_all_otu,
       file="./Cleaned_Data/01_BEF_all_otu_table.RDS")