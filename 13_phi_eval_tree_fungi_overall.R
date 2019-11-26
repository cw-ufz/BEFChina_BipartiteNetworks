### evaluating phi coefficient differences tree and fungal species/groups  #####

# R script (13/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25



################################################################################

library(phyloseq)
library(pgirmess)   #kruskal wallis test multiple group comparisons, post hoc test

############################

phi_tree_otu<-readRDS(file="./output/10_phi_coeff_trees_otus.RDS")
BEF_data<-readRDS(file="./Cleaned_Data/02_BEF_10r.RDS") 


#Sort for the different fungal functional groups

BEF_FungalFunction_sel<-list()
fungalfunctions<-c("Saprotroph","Ectomycorrhizal",
                   "Arbuscular Mycorrhizal","Plant Pathogen")

for(i in 1:length(fungalfunctions)){

  Group_sel<-subset_taxa(BEF_data, Fungal_guild==fungalfunctions[i])
  Group_sel<-data.frame(t(otu_table(Group_sel)))
  BEF_FungalFunction_sel[[i]] = colnames(Group_sel)
}
BEF_FungalFunction_sel[[1]]
names(BEF_FungalFunction_sel)<-fungalfunctions

#----- output 1/2
#saveRDS(BEF_FungalFunction_sel, file="./output/13_fungal_OTUnames.RDS")


# Get the respective phi coefficients

phi_fungal_spec<-list()

for( i in 1:length(fungalfunctions)){

  col_list<-which(colnames(phi_tree_otu) %in% BEF_FungalFunction_sel[[i]])
  phi_fungal_selection<-phi_tree_otu[,col_list]
  phi_fungal_spec[[i]]<-phi_fungal_selection[phi_fungal_selection>0]
}

names(phi_fungal_spec)<-fungalfunctions

tiff(filename ="./images/13_phi_fungal_groups_boxplot_text.tiff",
     width = 11, height = 8, units = "in", res=300)

par(cex.lab=1.2) 

boxplot(phi_fungal_spec$Ectomycorrhizal,
        phi_fungal_spec$`Arbuscular Mycorrhizal`,
        phi_fungal_spec$`Plant Pathogen`,
        phi_fungal_spec$Saprotroph, 
        rep(10,576),
        ylim = c(0, 0.5),
        col=c("green", "blue","red","brown","white"),
        ylab="positive phi coefficients",
        cex.lab=1.5,
        cex.axis=1.5,
        names=c("","","",
                "",""))
axis(side = 1, line = 0, at = c(1, 2, 3,4), cex.axis=1.5,
     labels = c("EcM fungi", "AM fungi","Plant pathogen", "Saprotrophic"),
     tick = F)
axis(side = 1, line = 1.5, at = c(1, 2, 3,4), cex.axis=1.5,
     labels = c("", "", "fungi", "fungi"), tick = F)

stripchart(phi_fungal_spec$Ectomycorrhizal,vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=1)
stripchart(phi_fungal_spec$`Arbuscular Mycorrhizal`,vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=2)
stripchart(phi_fungal_spec$`Plant Pathogen`,vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=3)
stripchart(phi_fungal_spec$Saprotroph,vertical = TRUE, 
           method = "stack", add = TRUE, pch = 20, cex=0.5, col = 'black', at=4)
text(0.8, 0.15, labels="a", cex=1.5)
text(1.8, 0.13, labels="b", cex=1.5)
text(2.8, 0.13, labels="b", cex=1.5)
text(3.8, 0.13, labels="b", cex=1.5)
text(5.3, 0.5, labels="p<0.05", cex=1.5)

dev.off()



phi_fungal_spec_vet<-unlist(phi_fungal_spec)

phi_fungal_spec_table<-data.frame(matrix(data=NA, ncol=2, nrow=length(phi_fungal_spec_vet),
       dimnames = list(as.character(1:length(phi_fungal_spec_vet)), 
                       c("phi_coefficient", "Fungal_guild"))))
phi_fungal_spec_table[,1]<-phi_fungal_spec_vet

     phi_fungal_spec_table[grep('^Saprotroph', names(phi_fungal_spec_vet)), 2]<-"Saprotroph"     
     phi_fungal_spec_table[grep('^Ectomycorrhizal', names(phi_fungal_spec_vet)), 2]<-"Ectomycorrhizal"
     phi_fungal_spec_table[grep('^Plant', names(phi_fungal_spec_vet)), 2]<-"Plant Pathogen"
     phi_fungal_spec_table[grep('^Arbuscular', names(phi_fungal_spec_vet)), 2]<-"Arbuscular Mycorrhizal"

 
krukal_phi_fungi<-kruskal.test(phi_fungal_spec_table$phi_coefficient  ~ as.factor(phi_fungal_spec_table$Fungal_guild),
                               data = phi_fungal_spec_table)

#data:  phi_fungal_spec_table$phi_coefficient by as.factor(phi_fungal_spec_table$Fungal_guild)
#Kruskal-Wallis chi-squared = 89.653, df = 3, p-value < 2.2e-16


mc_phi_fungi<-kruskalmc(phi_fungal_spec_table$phi_coefficient~phi_fungal_spec_table$Fungal_guild,
                             data=phi_fungal_spec_table,
                             probs=0.05/(4))
                             
# Multiple comparison test after Kruskal-Wallis 
# p.value: 0.0125 
# Comparisons
# obs.dif critical.dif difference
# Arbuscular Mycorrhizal-Ectomycorrhizal 574.15485     309.7858       TRUE
# Arbuscular Mycorrhizal-Plant Pathogen  121.64223     311.3567      FALSE
# Arbuscular Mycorrhizal-Saprotroph       40.82844     257.3146      FALSE
# Ectomycorrhizal-Plant Pathogen         695.79709     272.2300       TRUE
# Ectomycorrhizal-Saprotroph             614.98329     208.2713       TRUE
# Plant Pathogen-Saprotroph               80.81380     210.6009      FALSE                             
                             
#---> specialization of EcM is significantly higher than for other fungal 
#functional groups that do not further differ
                             
#ok, there is no general difference among the degree of specialization
#among the fungal functional groups of Saprotroph, AM , plant pathogens
