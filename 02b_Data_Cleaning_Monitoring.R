############################### Monitor data along analysis#####################

# R script (02b/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

#########################
library(biomformat)       #biomformat_1.10.1       #read biom file
library(BiocManager)      #BiocManager_1.30.4
#sessionInfo()
#R version 3.5.2 (2018-12-20)
#########################


#----- load datasets (raw+subsets)

# 1) raw data
x = read_biom("./Data/BEF_vipA_454_ITS_1283866816.json_Nov2017.biom") 
otu_matrix<-as(biom_data(x), "matrix")    

Env_var<-read.csv("./Data/BEF_vipA_2011_soil_typo_metadata.csv",row.names=1, 
                  dec=",", sep=";") 

BEF_raw<- phyloseq(otu_table(otu_matrix, taxa_are_rows = TRUE),
                   sample_data(Env_var)) #21732 taxa 

# 2) dominant OTU data (OTU with at least 4 reads ) #6693 taxa
BEF_dom<-readRDS(file="./Cleaned_Data/01_dominant4_phyloseq.RDS")

# 3) normalization to 700 sequences per sample 

BEF_rare700<-readRDS(file="./Cleaned_Data/02_BEF_rare700.RDS")


# 4) exclude OTUs with less than 10 reads
BEF_10r<-readRDS(file="./Cleaned_Data/02_BEF_10r.RDS")

Datacheck<-list(BEF_raw, BEF_dom, BEF_rare700, BEF_10r)

#----- calculation data characteristics

#input for calculations: phyloseq object

#for each dataset calculate 
  #total number of OTUs
  #total number of sequences
  #total number of OTUs for each of the 4 main fungal functional groups
  #% of OTUs of the main 4 functional groups
  #% of sequence reads of the main 4 functional groups

data_monitor<-matrix(data=NA, ncol = 14, nrow = 4,
                     dimnames = list(c("BEF_raw", "BEF_dominant","BEF_rare700" ,"BEF_10r"),
                                     c("No. of OTUs", "No. of reads",
                                        "Sapro OTUs", "Sapro OTU %", "Sapro sequence %",
                                       "EcM OTUs", "EcM OTU %", "EcM sequence %",
                                       "AM OTUs", "AM OTU %", "AM sequence %",
                                       "Patho OTUs", "Patho OTU %", "Patho sequence %"
                                       )))
              
fungalguilds<-list("Saprotroph", "Ectomycorrhizal", "Arbuscular Mycorrhizal",
                   "Plant Pathogen")



for(i in 1:length(Datacheck)){
  data_monitor[i,1]<-ntaxa(Datacheck[[i]])
  data_monitor[i,2]<-sum(taxa_sums(Datacheck[[i]]))
  if(i != 1){
    for(j in 1:length(fungalguilds)){
    Guild_phylo = subset_taxa(Datacheck[[i]], Fungal_guild==fungalguilds[[j]])#185  #121 #192
    data_monitor[i,3*j]<-ntaxa(Guild_phylo)
    data_monitor[i,3*j+1]<-round(ntaxa(Guild_phylo)/ntaxa(Datacheck[[i]]), digits = 2)
    data_monitor[i,3*j+2]<-round(sum(taxa_sums(Guild_phylo))/sum(taxa_sums(Datacheck[[i]])),digits = 2) 
    }
  }  
}  

#----- output

write.csv2(data_monitor, file="./output/02b_data_pruning_development.csv")



