load("newdynoNEMRupiall.rda")
allProbes=data.frame(s$ProbeID, s$Symbol, s$Probe_Id,s$Entrez_Gene_ID)
rownames(allProbes)=s$ProbeID
mapposSym=list()
mapposEnt=list()
for(i in 1:length(mynemlfmBS$mappos)){
	mapposSym[[i]]=allProbes[as.character(mynemlfmBS$mappos[[i]]),2]
	mapposEnt[[i]]=allProbes[as.character(mynemlfmBS$mappos[[i]]),3]	
}
rownames(DataF)
allProbes[allProbes$s.Probe_Id==rownames(DataF),]
y=data.frame()
for(i in 1:length(rownames(denMat.Rupi))){
x=subset(allProbes, Probe_Id==rownames(denMat.Rupi)[i], select = c(ProbeID, Symbol, Probe_Id))
y=rbind(y,x)
}
names(mapposSym)=names(mapposEnt)=names(mynemlfmBS$mappos)

load("FITERPPA.RData")

KOs = c("EGFR", "MTOR", "RPS6KB1", "PRKAA1", "PRKAB1")
T = topTable(fitE, coef=7, number=Inf)
Pr = as.matrix(T[,"P.Value"])
Padjr = as.matrix(T[,"adj.P.Val"])
rppafc=as.matrix(T[,"logFC"])
rownames(Pr) = T$ID
rownames(Padjr) = T$ID
rownames(rppafc) = T$ID
for(con in 8:11){
	T = topTable(fitE, coef=con, number=Inf)
	Pr = cbind(Pr, as.matrix(T[match(rownames(Pr), T$ID),"P.Value"]))
	Padjr = cbind(Padjr, as.matrix(T[match(rownames(Padjr), T$ID),"adj.P.Val"]))
	rppafc = cbind(rppafc, as.matrix(T[match(rownames(rppafc), T$ID),"logFC"]))
	
}
colnames(Pr) = colnames(rppafc)=colnames(Padjr) = KOs


rowsN=sapply(strsplit(rownames(rppafc), "_"), "[", 1)
rowsN=strsplit(rowsN, ", ")

eid=c()
for(i in 1:length(mapposEnt)){
eid=c(eid, as.character(mapposEnt[[i]]))
}
eid=unique(eid)

genes=unique(genes)
intersect(genes, rowsN)
rowEID=unlist(mget(rowsN, org.Hs.SYMBOL2EG, ifnotfound=NA))

V(igraph)$name[graph.bfs(igraph, RPnodes[1], neimode="out", unreachable=FALSE)$order]


dsGenes=list()
for(i in 1:5){
	x=V(igraph)$name[graph.bfs(igraph, RPnodes[i], neimode="out", unreachable=FALSE)$order]
	y= c(as.character(mapposSym[[RPnodes[i]]]))
	z=c(x,y)
	dsGenes[[i]]= z[!is.na(z)]
	rm(x,y,z)
} 
names(dsGenes)= RPnodes

xeid=c()
for(i in 1:107){
xeid=c(xeid, mget(as.character(rowsN)[i],  org.Hs.egALIAS2EG,ifnotfound=NA ))
}


#### Fix the RPPA data
rowsN[52]="EGFR"
rowsN[102]="EGFR"
newrppafc=matrix(0, 5, 121)
rppafc = cbind(rppafc[,1:31], rppafc[,32], rppafc[,32], rppafc[,33:35], rppafc[,36],rppafc[,36],rppafc[,36], rppafc[,37:47], rppafc[,48], rppafc[,48], rppafc[,49:52], rppafc[,53], rppafc[,53], rppafc[,53], rppafc[,54:65], rppafc[,66], rppafc[,66], rppafc[,66], rppafc[,67:76], rppafc[,77], rppafc[,77], rppafc[,77], rppafc[,78:79], rppafc[,80], rppafc[,80], rppafc[,80], rppafc[,81:93], rppafc[,94], rppafc[,94], rppafc[,95:104], rppafc[,105], rppafc[,105], rppafc[,106:107])  
rppafc=t(rppafc)
colnames(newrppafc)=c(colnames(rppafc)[1:31], "MAPK1", "MAPK3", colnames(rppafc)[33:35], "AKT1","AKT2","AKT3", colnames(rppafc)[37:47], "SMAD2", "SMAD3", colnames(rppafc)[49:52], "AKT1","AKT2","AKT3", colnames(rppafc)[54:65], "AKT1","AKT2","AKT3", colnames(rppafc)[67:76], "AKT1","AKT2","AKT3", colnames(rppafc)[78:79], "RPS6KA1", "RPS6KA2", "RPS6KA3", colnames(rppafc)[81:93], "MAPK1", "MAPK3", colnames(rppafc)[95:104], "STAT5A", "STAT5B", colnames(rppafc)[106:107]) 

newrppafc[, 1:121] = cbind(rppafc[,1:31], rppafc[,32], rppafc[,32], rppafc[,33:35], rppafc[,36],rppafc[,36],rppafc[,36], rppafc[,37:47], rppafc[,48], rppafc[,48], rppafc[,49:52], rppafc[,53], rppafc[,53], rppafc[,53], rppafc[,54:65], rppafc[,66], rppafc[,66], rppafc[,66], rppafc[,67:76], rppafc[,77], rppafc[,77], rppafc[,77], rppafc[,78:79], rppafc[,80], rppafc[,80], rppafc[,80], rppafc[,81:93], rppafc[,94], rppafc[,94], rppafc[,95:104], rppafc[,105], rppafc[,105], rppafc[,106:107])  


###### For every perturbation extract relevant row and columns
setwd("OLD")
for(a in 1:length(allPriors)){
control = set.default.parameters(unique(colnames(denMat.Rupi)), para=c(0.13, 0.05), type="CONTmLLMAP")
prior=allPriors[[a]]
colnames(prior)=rownames(prior)=priorgenes
control$Pm=prior
mynem=nem(denMat.Rupi, inference="nem.greedy", control=control)
#mynem= mynemlfmBS
igraph=igraph.from.graphNEL(mynem$graph)
M = transitive.reduction((as(mynem$graph, "matrix") > 0.5)*1)
igraph2=as(M, "graphNEL")
#igraph2=igraph.from.graphNEL(igraph2)
#par(mfrow = c(2, 2))
#igraph=igraph2

#allProbes=data.frame(s$ProbeID, s$Symbol, s$Probe_Id,s$Entrez_Gene_ID)
#rownames(allProbes)=s$ProbeID
mapposSym=list()
mapposEnt=list()
for(i in 1:length(mynem$mappos)){
	mapposSym[[i]]=allProbes[as.character(mynem$mappos[[i]]),2]
	mapposEnt[[i]]=allProbes[as.character(mynem$mappos[[i]]),3]
}

names(mapposSym)=names(mapposEnt)=names(mynemlfmBS$mappos)


dsGenes=list()
for(i in 1:5){
	x=V(igraph)$name[graph.bfs(igraph, RPnode[i], neimode="out", unreachable=FALSE)$order]
	y=c()
	for(j in 1:)
	y= c(as.character(mapposSym[[i]]))
	
	z=c(x,y)
	dsGenes[[i]]= z[!is.na(z)]
	rm(x,y,z)
} 
names(dsGenes)= RPnode

newrppafc =validationData$proteiomic
fileV=paste(names(allPriors)[[a]],"_validationBarO.pdf", sep="" )
pdf(fileV)
for (i in 1:5){
	x=names(dsGenes)[i]
	y=dsGenes[[i]]
	y=intersect(colnames(newrppafc), y)
	barplot(newrppafc[as.character(x),as.character(y)], ylim=c(-3,3),xlab="Downstream proteins", ylab="log fold change", las = 2, main=paste("Knock-down for ", x, sep=" "))
}
dev.off()
fileG=paste(names(allPriors)[[a]],"_graphO.pdf", sep="" )
pdf(fileG)
plot(igraph2)
dev.off()
}

> colnames(allPriors$LFM)
 [1] "6794"  "5289"  "26995" "3953"  "7249"  "7248"  "2932"  "2931"  "5894"  "3691" 
[11] "5564"  "8915"  "1956"  "9700"  "6714"  "10885" "5562"  "6195"  "2475"  "6198" 
> priorgenes=c("STK11", "PIK3C3","TRUB2", "LEPR", "TSC2", "TSC1", "GSK3B", "GSK3A", "RAF1", "ITGB4", "PRKAB1", "BCL10", "EGFR", "ESPL1", "SRC", "WDR3", "PRKAA1", "RPS6KA1", "MTOR", "RPSKB1")

e.type=infer.edge.type(mynem, fcData)
MAPK1, MAPK3
AKT1, AKT2, AKT3
SMAD2, SMAD3
AKT1, AKT2, AKT3
RPS6KA1, RPS6KA2, RPS6KA3
STAT5A, STAT5B 



KOs = c("EGFR", "MTOR", "RPS6KB1", "PRKAA1", "PRKAB1")
T = topTable(fitE, coef=7, number=Inf)
Pr = as.matrix(T[,"P.Value"])
Padjr = as.matrix(T[,"adj.P.Val"])
rppafc=as.matrix(T[,"logFC"])
rownames(Pr) = rownames(T)
rownames(Padjr) = rownames(T)
rownames(rppafc) = rownames(T)
for(con in 8:11){
	T = topTable(fitE, coef=con, number=Inf)
	Pr = cbind(Pr, as.matrix(T[match(rownames(Pr), rownames(T)),"P.Value"]))
	Padjr = cbind(Padjr, as.matrix(T[match(rownames(Padjr), rownames(T)),"adj.P.Val"]))
	rppafc = cbind(rppafc, as.matrix(T[match(rownames(rppafc), rownames(T)),"logFC"]))
	
}
colnames(Pr) = colnames(rppafc)=colnames(Padjr) = KOs


######################################################
######################################################
KOs = c("EGFR", "MTOR", "RPS6KB1", "PRKAA1", "PRKAB1")
T = topTable(fitE, coef=7, number=Inf)
P = as.matrix(T[,"P.Value"])
Padj = as.matrix(T[,"adj.P.Val"])
lfc= as.matrix(T[,"logFC"])
rownames(P) = rownames(T)
rownames(Padj) =rownames(T)
rownames(lfc) = rownames(T)#
for(con in 8:11){
   T = topTable(fitE, coef=con, number=Inf)
   P = cbind(P, as.matrix(T[match(rownames(P), rownames(T)),"P.Value"]))
   Padj = cbind(Padj, as.matrix(T[match(rownames(Padj), rownames(T)),"adj.P.Val"]))
   lfc = cbind(lfc, as.matrix(T[match(rownames(lfc), rownames(T)),"logFC"])) #
}
colnames(P) = colnames(Padj) =colnames(lfc) = KOs

#> dat=dat[-c(108:113),]
#> dat.m=dat.m[-c(108:113),]
#> anno=anno[-c(108:113),]


# extract phopho proteins
# pattern __ in row names
lfc.m=lfc[grep("__",rownames(lfc)),]
Padj.m=Padj[grep("__",rownames(Padj)),]

pdf("logfold_changes.pdf", height=8, width=18)
g=list()
SG=names(mynemlfmBS$mappos)[-21]
EG=c()
for(i in 1:20){
	EG= c(EG, as.character(mapposSym[[i]]))
}
EG=unique(EG)
## Check for no names i.e "" was at 6th position so EG=EG[-6]
for(i in 1:5){
plotThis=data.frame(PV=Padj[,i], FC=lfc[,i])
Phospho=rep("FALSE",nrow(plotThis))
Phospho[grep("__",rownames(plotThis))]="TRUE"
Signi=rep(" ",nrow(plotThis))
Signi[which(plotThis$PV<0.05)]="*"
Category1=c()

for(j in 1:length(SG)){
	Category1=c(Category1, grep(SG[j], rownames(plotThis)))
}
Category2=c()
for(k in 1:length(SG)){
	Category2=c(Category2, grep(EG[j], rownames(plotThis)))
}

Category=rep("Others", nrow(plotThis))
Category[Category1]="S-genes"
Category[Category2]="E-genes"
Category[intersect(Category1, Category2)]="Both"

plotThis=data.frame(plotThis, Phospho=Phospho, Significant=Signi, Category=Category)
tit=paste("Knockdown for ", colnames(lfc)[i])
## rownames(plotThis)=sapply(strsplit(rownames(plotThis), "__"), "[",1)
g[[i]]=ggplot(data=plotThis, aes(x=factor(rownames(plotThis)) ,y=FC, fill=Phospho)) + geom_bar(position = 'dodge') + geom_text(aes(label=Significant))+opts(axis.text.x=theme_text(angle=90)) + xlab("RPPA proteins") +ylab("log Fold change") + ggtitle(tit)+facet_wrap(~ Category)
}
dev.off()


getSEgenes=function(nemres, mapposSym,genes){
	igraph=igraph.from.graphNEL(nemres$graph)
	sgenes=list()
	egenes=list()
	for(i in 1:length(genes)){
		x=V(igraph)$name[graph.bfs(igraph, genes[i], neimode="out", unreachable=FALSE)$order]
		x=x[-which(is.na(x))]
		sgenes[[i]]=x
		eg=c()
		for(j in 1:length(x)){
			y=as.character(mapposSym[x[j]][[1]])
			eg=c(eg,y)	
		}
		eg=unique(eg)
		eg=eg[-which(is.na(eg))]
		egenes[[i]]=eg
	}
	names(sgenes)=names(egenes)=genes
	efects=list(SGENE=sgenes, EGENE=egenes)
	return(efects)
}

computeData=function(pv, eft){
	sgrow=c()
	egrow=c()
	sg=eft[[1]]
	eg=eft[[2]]
	for (i in 1:length(sg)){
		sgrow=c(sgrow, grep(sg[i],rownames(pv)))
	}
	for(j in 1:length(eg)){
		egrow=c(egrow, grep(eg[j],rownames(pv)))
	}
	sgrow=unique(sgrow)
	egrow=unique(egrow)
	allrows=list(SGrow=sgrow, EGrow=egrow)
	return(allrows)
}



Padj=Padj[,7:11]
lfc=lfc[,7:11]

g=list()
for(i in 1:5){
plotThis=data.frame(PV=Padj[,i], FC=lfc[,i])
Phospho=rep("FALSE",nrow(plotThis))
Phospho[grep("__",rownames(plotThis))]="TRUE"
Signi=rep(" ",nrow(plotThis))
Signi[which(plotThis$PV<0.05)]="*"
category=rep("Others", nrow(plotThis))
rowval=computeData(plotThis, list(eft[[1]][[i]], eft[[2]][[i]]))
category[rowval[[1]]]="S-gene"
category[rowval[[2]]]="E-gene"
ro=intersect(rowval[[1]], rowval[[2]])
category[ro]="Both"
plotThis=data.frame(plotThis, Phospho=Phospho, Significant=Signi, Category=category)
tit=paste("Knockdown for ", KN[i])
#tit=tits[i]
# rownames(plotThis)=sapply(strsplit(rownames(plotThis), "__"), "[",1)
g[[i]]=ggplot(data=plotThis, aes(x=factor(rownames(plotThis)) ,y=FC, fill=Phospho)) + geom_bar(position = 'dodge') + geom_text(aes(label=Significant))+opts(axis.text.x=theme_text(angle=90)) + xlab("RPPA proteins") +ylab("log Fold change") + ggtitle(tit)+facet_grid(. ~ Category)

}

pdf("logfold_changes.pdf", height=8, width=18)
dev.off()

 eg=c()
for(i in 1:20){
 y=mapposSym[[i]]
 y=as.character(y)
 eg=c(eg,y)
 }
 
 egrow=c()
for(j in 1:length(eg)){
		egrow=c(egrow, grep(eg[j],rownames(lfc)))
	}

 5 19 22 44 57 58 61 69 74 25 27 25 27 25 27 72 85 49 50 56 93
 
hmd=lfc[egrow,1:6])

library(nem)
dens = getDensityMatrix(P, dirname=file.path(DIR, "QualityControlRPPA")) # mäßiger fit
library(EMA)
cl.sample = clustering(dens, method="ward", metric="pearson")
cl.gene = clustering(t(dens), method="ward", metric="pearson")
pdf(file=file.path(DIR, "HeatmapRPPAEffects.pdf"))
clustering.plot(tree=cl.sample, tree.sup=cl.gene, data=dens, names.sup=FALSE, lab=data.frame(KOs), title="log-densities")
dev.off()




hmd=Padj[,1:6]
cl.sample = clustering(hmd, method="ward", metric="pearson")
cl.gene = clustering(t(hmd), method="ward", metric="pearson")
clustering.plot(tree=cl.sample, tree.sup=cl.gene, data=hmd, names.sup=FALSE, lab=data.frame(colnames(hmd)), title="log-densities")

phosp=grep("__", rownames(hmd))
hmdAct=lfc[,1:3]
hmdInh=lfc[,4:6]
hmpAct=Padj[,1:3]
hmpInh=Padj[,4:6]
colnames(hmdAct)=colnames(hmpAct)=c("Activation_low", "Activation_medium", "Activation_high")
colnames(hmdInh)=colnames(hmpInh)=c("Inhibition_low", "Inhibition_medium", "Inhibition_high")
phosp=grep("__", rownames(hmd))

hmdAct=lfc[,1:3]
hmdInh=lfc[,4:6]
hmpAct=Padj[,1:3]
hmpInh=Padj[,4:6]

hmdAct.net=lfc[egrow,1:3]
hmdInh.net=lfc[egrow,4:6]
hmpAct.net=Padj[egrow,1:3]
hmpInh.net=Padj[egrow,4:6]

colnames(hmdAct.net)=colnames(hmpAct.net)=c("Activation_low", "Activation_medium", "Activation_high")
colnames(hmdInh.net)=colnames(hmpInh.net)=c("Inhibition_low", "Inhibition_medium", "Inhibition_high")
phosp=grep("__", rownames(hmd))


phosp=grep("__", rownames(lfc))
hmdAct.ph=lfc[phosp,1:3]
hmdInh.ph=lfc[phosp,4:6]
hmpAct.ph=Padj[phosp,1:3]
hmpInh.ph=Padj[phosp,4:6]

colnames(hmdAct.ph)=colnames(hmpAct.ph)=c("Activation_low", "Activation_medium", "Activation_high")
colnames(hmdInh.ph)=colnames(hmpInh.ph)=c("Inhibition_low", "Inhibition_medium", "Inhibition_high")
phosp=grep("__", rownames(hmd))

phosp=grep("__", rownames(hmdAct.net))
hmdAct.ph.net= hmdAct.net[phosp,]
hmdInh.ph.net= hmdInh.net[phosp,]
hmpAct.ph.net=hmpAct.net[phosp,]
hmpInh.ph.net=hmpInh.net[phosp,]

colnames(hmdAct.ph.all)=colnames(hmpAct.ph.all)=c("Activation_low", "Activation_medium", "Activation_high")
colnames(hmdInh.ph.all)=colnames(hmpInh.ph.all)=c("Inhibition_low", "Inhibition_medium", "Inhibition_high")



pdf("LFheatmapActAll.pdf", height=28, width=14)
heatmap(hmdAct)
title("(All proteins) log fold change- Activation")
dev.off()

pdf("LFheatmapActNet.pdf", height=28, width=14)
heatmap(hmdAct.net)
title("(Network proteins) log fold change- Activation")
dev.off()

pdf("LFheatmapActNetPH.pdf", height=28, width=14)
heatmap(hmdAct.ph.net)
title("(Network proteins - phospho) log fold change- Activation")
dev.off()

pdf("LFheatmapActAllPH.pdf", height=28, width=14)
heatmap(hmdAct.ph.net)
title("(All proteins - phospho) log fold change- Activation")
dev.off()

####

pdf("LFheatmapActAll.pdf", height=28, width=14)
heatmap(hmdAct)
title("(All proteins) log fold change- Activation")
dev.off()

pdf("LFheatmapActNet.pdf", height=28, width=14)
heatmap(hmdAct.net)
title("(Network proteins) log fold change- Activation")
dev.off()

pdf("LFheatmapActNetPH.pdf", height=28, width=14)
heatmap(hmdAct.ph.net)
title("(Network proteins - phospho) log fold change- Activation")
dev.off()

pdf("LFheatmapActAllPH.pdf", height=28, width=14)
heatmap(hmdAct.ph.net)
title("(All proteins - phospho) log fold change- Activation")
dev.off()

#######

pdf("LFheatmapInhAll.pdf", height=28, width=14)
heatmap(hmdInh)
title("(All proteins) log fold change- Inhibition")
dev.off()

pdf("LFheatmapInhNet.pdf", height=28, width=14)
heatmap(hmdInh.net)
title("(Network proteins) log fold change- Inhibition")
dev.off()

pdf("LFheatmapInhNetPH.pdf", height=28, width=14)
heatmap(hmdInh.ph.net)
title("(Network proteins - phospho) log fold change- Inhibition")
dev.off()

pdf("LFheatmapInhAllPH.pdf", height=28, width=14)
heatmap(hmdInh.ph.net)
title("(All proteins - phospho) log fold change- Inhibition")
dev.off()




myMap=lfc[,7:11]

cl.sample = clustering(hmdAct.net, method="ward", metric="pearson")
cl.gene = clustering(t(hmdAct.net), method="ward", metric="pearson")
clustering.plot(tree=cl.sample, tree.sup=cl.gene, data=hmdAct.net, names.sup=FALSE, lab=data.frame(colnames(hmdAct.net)), scale="none", title="log fold change")

pdf("inhNetEMA.pdf")
cl.sample = clustering(hmdInh.net, method="ward", metric="pearson")
cl.gene = clustering(t(hmdInh.net), method="ward", metric="pearson")
clustering.plot(tree=cl.sample, tree.sup=cl.gene, data=hmdInh.net, names.sup=TRUE, lab=data.frame(colnames(hmdInh.net)),scale="none",  title="log fold change (Inhibition)")
dev.off()
null device 
          1 
pdf("InhAllEMA.pdf")
cl.sample = clustering(hmdInh, method="ward", metric="pearson")
cl.gene = clustering(t(hmdInh), method="ward", metric="pearson")
clustering.plot(tree=cl.sample, tree.sup=cl.gene, data= hmdInh, names.sup=TRUE, lab=data.frame(colnames(hmdInh)), scale="none", title="log fold change (Activation)")
dev.off()


for(i in 1:3){
plotThis=data.frame(Proteins=rownames(hmpAct.net),PV=hmpAct.net[,i], FC=hmdAct.net[,i])
Phospho=rep("FALSE",nrow(plotThis))
Phospho[grep("__",plotThis$Proteins)]="TRUE"
Signi=rep(" ",nrow(plotThis))
Signi[which(plotThis$PV<0.05)]="*"
plotThis=data.frame(plotThis, Phospho=Phospho, Significant=Signi)
tit=paste("colnames(hmpAct.net)[i]")
## rownames(plotThis)=sapply(strsplit(rownames(plotThis), "__"), "[",1)
g[[i]]=ggplot(data=plotThis, aes(x=Proteins ,y=FC, fill=Phospho)) + geom_bar(position = 'dodge') + geom_text(aes(label=Significant))+opts(axis.text.x=theme_text(angle=90)) + xlab("RPPA proteins") +ylab("log Fold change") + ggtitle(tit)
}
dev.off()

pdf("knockdownplot.pdf")
for(i in 7:11){
	plotThis=data.frame(Proteins=rownames(hmpAct.net),PV=Padj[rownames(hmpInh.net),i], FC=lfc[rownames(hmpInh.net),i])
Phospho=rep("FALSE",nrow(plotThis))
Phospho[grep("__",plotThis$Proteins)]="TRUE"
Signi=rep(" ",nrow(plotThis))
Signi[which(plotThis$PV<0.05)]="*"
plotThis=data.frame(plotThis, Phospho=Phospho, Significant=Signi)
j=i-6
tit=paste("Knockdown ", KOs[j], sep="")
g[[i]]=ggplot(data=plotThis, aes(x=Proteins ,y=FC, fill=Phospho)) + geom_bar(position = 'dodge') + geom_text(aes(label=Significant))+opts(axis.text.x=theme_text(angle=90)) + xlab("RPPA proteins") +ylab("log Fold change") + ggtitle(tit)

}
dev.off()
