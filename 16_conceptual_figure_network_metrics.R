#################### Conceptual figure of network metrics ######################

# R script (16/16) for analyses in Weissbecker et al. 2019 
# Bipartite Network Analysis Soil Fungi in Chinese Experimental Forests

# version: 2019-11-25

################################################################################

pdf("~/Desktop/test.pdf",width=17/2.54,height=17/2.54,pointsize=6)
par(mar=rep(0.5,4),mfcol=c(4,2))
conmat <- matrix(sample(c(0,1),700,T,c(3/7,4/7)),nrow=7,ncol=100)
no_conmat <- matrix(sample(c(0,1),700,T,c(5/7,2/7)),nrow=7,ncol=100)
plotweb(conmat,bor.col.interaction = NA,bor.col.high = NA,bor.col.low = NA,labsize = 0.01,col.interaction = alpha("grey50",0.5))
plotweb(no_conmat,bor.col.interaction = NA,bor.col.high = NA,bor.col.low = NA,labsize = 0.01,col.interaction = alpha("grey50",0.5))

networklevel(conmat,"generality")
networklevel(no_conmat,"generality")
networklevel(conmat,"connectance")
networklevel(no_conmat,"connectance")


nestmat <- matrix(0,nrow=7,ncol=100)
nestmat[1,] <- sample(c(0,1),100,replace=T,prob=c(0.1,0.9))
nestmat[2,nestmat[1,]==1] <- sample(c(0,1),length(which(nestmat[1,]==1)),replace=T,prob=c(0.2,0.8))
nestmat[2,nestmat[1,]!=1] <- sample(c(0,1),length(which(nestmat[1,]!=1)),replace=T,prob=c(0.9,0.1))
nestmat[3,nestmat[2,]==1] <- sample(c(0,1),length(which(nestmat[2,]==1)),replace=T,prob=c(0.3,0.7))
nestmat[3,nestmat[2,]!=1] <- sample(c(0,1),length(which(nestmat[2,]!=1)),replace=T,prob=c(0.95,0.05))
nestmat[4,nestmat[3,]==1] <- sample(c(0,1),length(which(nestmat[3,]==1)),replace=T,prob=c(0.3,0.7))
nestmat[4,nestmat[3,]!=1] <- sample(c(0,1),length(which(nestmat[3,]!=1)),replace=T,prob=c(0.95,0.05))
nestmat[5,nestmat[4,]==1] <- sample(c(0,1),length(which(nestmat[4,]==1)),replace=T,prob=c(0.3,0.7))
nestmat[5,nestmat[4,]!=1] <- sample(c(0,1),length(which(nestmat[4,]!=1)),replace=T,prob=c(0.97,0.03))
nestmat[6,nestmat[5,]==1] <- sample(c(0,1),length(which(nestmat[5,]==1)),replace=T,prob=c(0.4,0.6))
nestmat[6,nestmat[5,]!=1] <- sample(c(0,1),length(which(nestmat[5,]!=1)),replace=T,prob=c(0.98,0.02))
nestmat[7,nestmat[6,]==1] <- sample(c(0,1),length(which(nestmat[6,]==1)),replace=T,prob=c(0.6,0.4))
nestmat[7,nestmat[6,]!=1] <- sample(c(0,1),length(which(nestmat[6,]!=1)),replace=T,prob=c(0.99,0.01))
nestmat <- nestmat[order(rowSums(nestmat)),order(colSums(nestmat))]
nu_nestmat <- matrix(sample(c(0,1),700,T,prob=c(700-sum(nestmat),sum(nestmat))),nrow=7,ncol=100)

networklevel(nestmat,"nestedness")
networklevel(nu_nestmat,"nestedness")

plotweb(nestmat,bor.col.interaction = NA,bor.col.high = NA,bor.col.low = NA,labsize = 0.01,col.interaction = alpha("grey50",0.5))
plotweb(nu_nestmat,bor.col.interaction = NA,bor.col.high = NA,bor.col.low = NA,labsize = 0.01,col.interaction = alpha("grey50",0.5))

networklevel(nestmat,"generality")
networklevel(nu_nestmat,"generality")
networklevel(nestmat,"connectance")
networklevel(nu_nestmat,"connectance")

genmat <- sapply(1:100, function(x) sample(c(0,1),7,replace=T,c(3/7,4/7)))
sum(genmat)
no_genmat <- matrix(sample(c(0,1),700,T,c(6/7,1/7)),nrow=7,ncol=100)
no_genmat[1,colSums(no_genmat)==0] <- 1
no_genmat[1,no_genmat[1,]==0] <- 1
sum(no_genmat)
no_genmat[2,no_genmat[2,]==0] <- 1
sum(no_genmat)
no_genmat[3,no_genmat[3,]==0] <- 1
sum(no_genmat)
no_genmat[4,no_genmat[4,]==0][1:(398-343)] <- 1
sum(no_genmat)

networklevel(genmat,"generality")
networklevel(no_genmat,"generality")
networklevel(genmat,"connectance")
networklevel(no_genmat,"connectance")

plotweb(genmat,bor.col.interaction = NA,bor.col.high = NA,bor.col.low = NA,labsize = 0.01,col.interaction = alpha("grey50",0.5))
plotweb(no_genmat,bor.col.interaction = NA,bor.col.high = NA,bor.col.low = NA,labsize = 0.01,col.interaction = alpha("grey50",0.5))


modmat <- matrix(0,nrow=7,ncol=100)
modmat[1,] <- sample(c(0,1),100,replace=T,c(5/7,2/7))
modmat[2,modmat[1,]==1] <- sample(c(0,1),length(which(modmat[1,]==1)),replace=T,c(0.05,0.95))
modmat[2,modmat[1,]!=1] <- sample(c(0,1),length(which(modmat[1,]!=1)),replace=T,c(0.95,0.05))
modmat[3,modmat[1,]==1] <- sample(c(0,1),length(which(modmat[1,]==1)),replace=T,c(0.05,0.95))
modmat[3,modmat[1,]!=1] <- sample(c(0,1),length(which(modmat[1,]!=1)),replace=T,c(0.95,0.05))
modmat[4,colSums(modmat[1:3,])==0] <- sample(c(0,1),length(which(colSums(modmat[1:3,])==0)),replace=T,c(0.35,0.65))
modmat[4,colSums(modmat[1:3,])!=0] <- sample(c(0,1),length(which(colSums(modmat[1:3,])!=0)),replace=T,c(0.95,0.05))
modmat[5,modmat[4,]!=1] <- sample(c(0,1),length(which(modmat[4,]!=1)),replace=T,c(0.95,0.05))
modmat[5,modmat[4,]==1] <- sample(c(0,1),length(which(modmat[4,]==1)),replace=T,c(0.05,0.95))
modmat[6,colSums(modmat[1:5,])==0] <- 1
modmat[7,modmat[6,]==1] <- sample(c(0,1),length(which(modmat[6,]==1)),replace=T,c(0.05,0.95))
modmat[7,modmat[6,]!=1] <- sample(c(0,1),length(which(modmat[6,]!=1)),replace=T,c(0.95,0.05))

nu_modmat <- matrix(sample(c(0,1),700,T,prob=c(700-sum(modmat),sum(modmat))),nrow=7,ncol=100)

plotweb(modmat,bor.col.interaction = NA,bor.col.high = NA,bor.col.low = NA,labsize = 0.01,col.interaction = alpha("grey50",0.5))
plotweb(nu_modmat,bor.col.interaction = NA,bor.col.high = NA,bor.col.low = NA,labsize = 0.01,col.interaction = alpha("grey50",0.5))

computeModules(modmat)
computeModules(nu_modmat)
dev.off()
