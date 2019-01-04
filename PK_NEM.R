setwd("/home/bit/praveen/LunCanceRupi/Res/NEM")
load("rupiInputNEM.RData")
D= denMat.Rupi
control = set.default.parameters(unique(colnames(D)), type="CONTmLLMAP",para=c(0.13, 0.05))
### without prior
res5=nem.bootstrap(D, nboot=10000,inference="ModuleNetwork", control=control)
plot(res5$graph)

### with prior

control$Pe=EgenePrior
control$Pm=prior.Sgene
res.prior2=nem.bootstrap(D, nboot=10000,inference="ModuleNetwork", control=control, verbose=F)
save(res.Eprior, file='rupiResWithPrior2.RData')






plot(res5, plot.probs=TRUE, transitiveReduction=TRUE,  PDF=F, thresh=0.5)
plot(newres, plot.probs=TRUE, transitiveReduction=TRUE,  PDF=F)

plot(Sgra, plot.probs=TRUE, transitiveReduction=TRUE, PDF=F)

library(org.Hs.eg.db)
library(annotate)
no=lookUp(nodes(res.prior$graph), 'org.Hs.eg', 'SYMBOL')
namesN=as.character(unlist(no))
nodes(res.prior$graph)=namesN

geneSet1=c(testGenes, Egenes[1:80])
geneSet2=c(testGenes, Egenes[81:160])
geneSet3=c(testGenes, Egenes[161:240])
geneSet4=c(testGenes, Egenes[241:320])
geneSet5=c(testGenes, Egenes[321:344])


for(k in 1:5){
genesIntersect=intersect(geneSets[[k]], allNodesKegg)
PATHWAYs[[k]] = dkSP(keggAll, genesIntersect)
colnames(PATHWAYs[[k]])=genesIntersect
rownames(PATHWAYs[[k]])=genesIntersect
}

PKNOMs=list()
for(k in 1:length(geneSets)){
#print(dim(PATHWAYs[[k]]))
#GOs[[k]] = mergeMatrix(GO[[k]],geneSets[[k]])
#DOMAINs[[k]] = mergeMatrix(DOMAINs[[k]],geneSets[[k]])
#DOMAIN.INTs[[k]] = mergeMatrix(DOMAIN.INTs[[k]],geneSets[[k]])
#PATHWAYs[[k]] =mergeMatrix(PATHWAYs[[k]] , geneSets[[k]])
#PCOMs[[k]] = mergeMatrix(PCOMs[[k]],geneSets[[k]])
#DBASEs[[k]] = mergeMatrix(DBASEs[[k]],geneSets[[k]])
diag(PATHWAYs[[k]]) = 0
diag(GOs[[k]]) = 0
diag(DBASEs[[k]]) = 0
diag(PCOMs[[k]]) = 0
diag(DOMAINs[[k]]) = 0
diag(DOMAIN.INTs[[k]]) = 0
}

for(k in 1:length(geneSets)){
print(dim(GOs[[k]]))
print(dim(DOMAINs[[k]]))
print(dim(DOMAIN.INTs[[k]]))
print(dim(PATHWAYs[[k]]))
print(dim(PCOMs[[k]]))
print(dim(DBASEs[[k]]))
}
save(GOs, DOMAINs, DOMAIN.INTs, PATHWAYs, PCOMs, DBASEs, file='EgeneInfo.RData')

SourceS=list()
OrigData=list()
for(k in 1:5){

OrigData[[k]]=data.frame(GO=c(GOs[[k]]), DOMAIN=c(DOMAINs[[k]]), DOMAIN.INT=c(DOMAIN.INTs[[k]]), PATHWAY=c(PATHWAYs[[k]]), PCOM=c(PCOMs[[k]]), DBASE=c(DBASEs[[k]]))
testGenes=colnames(GOs[[k]])
correlNet=cor.net(origData)
correlNet[is.na(correlNet)]=0
####################
diag(correlNet)=1
dissimilarity = 1 - correlNet
distance = as.dist(dissimilarity)
clust = hclust(distance)
labels=cutree(clust, k=2)
acceptedNod=clust$label[labels==1]
acceptedNod = unique(acceptedNod)
print(acceptedNod)
PKdata=origData[acceptedNod]
data=PKdata
data[which(is.nan(data))] <- 0 
data[is.nan(data)] <- 0 
data[is.na(data)] <- 0
data[is.null(data)] <- 0
SourceS[[k]]=list()
for(i in 1:ncol(data)){
		SourceS[[k]][[i]]=matrix(data[,acceptedNod[i]], length(testGenes), length		(testGenes),byrow=T)
		rownames(SourceS[[k]][[i]])=testGenes
		colnames(SourceS[[k]][[i]])=testGenes
	}
}


 for(i in 1:5){
 print(length(SourceS[[i]][[1]]))
  print(length(SourceS[[i]][[2]]))
   print(length(SourceS[[i]][[3]]))
    print(length(SourceS[[i]][[4]]))
     print(length(SourceS[[i]][[5]]))
     cat("#####")
 }


origDataS=

source=list(GOs[[i]], DOMAINs[[i]], DOMAIN.INTs[[i]], PATHWAYs[[i]], PCOMs[[i]], DBASEs[[i]] )
data=matrix(0, nrow=ncol(source[[1]])*ncol(source[[1]]), ncol=length(source))

PKNOMs=list()
for(k in 1:5){
sorce=SourceS[[k]]
testGenes=colnames(SourceS[[k]][[1]])
datas=data.frame(matrix(0, nrow=length(testGenes)*length(testGenes), ncol=length(source)))
for(i in 1:length(source)){
	datas[,i]=c(source[[i]])
}
testGenes=colnames(SourceS[[k]][[1]])
PKNOMs[[k]]=calcNO(datas, testGenes)
diag(PKNOMs[[k]])=0
colnames(PKNOMs[[k]])=testGenes
rownames(PKNOMs[[k]])=testGenes
print(dim(PKNOMs[[k]]))
}
PKNOM1=PKNOMs[[1]][-c(1:20), -c(21:100)]
PKNOM2=PKNOMs[[2]][-c(1:20), -c(21:100)]
PKNOM3=PKNOMs[[3]][-c(1:20), -c(21:100)]
PKNOM4=PKNOMs[[4]][-c(1:20), -c(21:100)]
PKNOM5=PKNOMs[[5]][-c(1:20), -c(21:44)]
PKNOMEgene=rbind(PKNOM1, PKNOM2,PKNOM3,PKNOM4, PKNOM5)


EG.prior=t(apply(PKNOMEgene, 1, function(x) x/sum(x)))
illD=rownames(denMat.Rupi)
DM=denMat.Rupi
newrow=Sel.Genes[Sel.Genes$Probe_Id %in% illD,]$Entrez_Gene_ID
length(newrow)
rownames(DM)=newrow
newP=data.frame(matrix(0, 585,20))
colnames(newP)=colnames(DM)
rownames(newP)=rownames(DM)
for(i in 1:length(rownames(DM))){
x=rownames(DM)[i]
newP[i,]=EG.prior[x,]
}



dimnames(sGeneP.NOM)=dimnames(prior.Sgene)
D= denMat.Rupi
control = set.default.parameters(unique(colnames(D)), type="CONTmLLMAP",para=c(0.13, 0.05))
control1=control
control2=control
control3=control
control4=control
control5=control
control1$Pm=sGeneP.NOM
control2$Pm=prior.Sgene
control3$Pe=EgenePrior
control4$Pe=EgenePrior
control4$Pm=sGeneP.NOM
control5$Pe=EgenePrior
control5$Pm=prior.Sgene


resNOM.SG=nem.bootstrap(D, nboot=5000,inference="ModuleNetwork", control=control1, verbose=F)
resLFM.SG=nem.bootstrap(D, nboot=5000,inference="ModuleNetwork", control=control2, verbose=F)
res.EG=nem.bootstrap(D, nboot=5000,inference="ModuleNetwork", control=control3, verbose=F)
resNOM.SG_EG=nem.bootstrap(D, nboot=5000,inference="ModuleNetwork", control=control4, verbose=F)
resLFM.SG_EG=nem.bootstrap(D, nboot=5000,inference="ModuleNetwork", control=control5, verbose=F)


nemModelSelection(lambdas,D,inference="nem.greedy",models=NULL,control=control,verbose=FALSE)
plot(res.EG, plot.probs=TRUE, transitiveReduction=TRUE, filename='resEG', PDF=TRUE)

plot(netnewSP, plot.probs=TRUE, transitiveReduction=TRUE, filename='resEG', PDF=FALSE)

("BCL10",  "EGFR" ,  "ESPL1" , "GSK3A" , "GSK3B" , "ITGB4" , "LEPR"  , "LKB1"  , "mTOR"  , "p70S6K", "P90S6K", "PIK3C3", "PRKAA1", "PRKAB1", "RAF1" ,  "SRC" ,  "TRUB2" , "TSC1"  , "TSC2"  , "WDR3" ) 

