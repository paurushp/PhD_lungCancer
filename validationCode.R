load("newdynoNEMRupiall.rda")
allProbes=data.frame(s$ProbeID, s$Symbol, s$Entrez_Gene_ID)
rownames(allProbes)=s$ProbeID
mapposSym=list()
mapposEnt=list()
for(i in 1:length(mynemlfmBS$mappos)){
	mapposSym[[i]]=allProbes[as.character(mynemlfmBS$mappos[[i]]),2]
	mapposEnt[[i]]=allProbes[as.character(mynemlfmBS$mappos[[i]]),3]
}

names(mapposSym)=names(mapposEnt)=names(mynemlfmBS$mappos)

load("FITERPPA.RData")

KOs = c("EGFR", "MTOR", "RPS6KB1", "PRKAA1", "PRKAB1")
T = topTable(fitE, coef=7, number=Inf)
P = as.matrix(T[,"P.Value"])
Padj = as.matrix(T[,"adj.P.Val"])
rppafc=as.matrix(T[,"logFC"])
rownames(P) = T$ID
rownames(Padj) = T$ID
rownames(rppafc) = T$ID
for(con in 8:11){
	T = topTable(fitE, coef=con, number=Inf)
	P = cbind(P, as.matrix(T[match(rownames(P), T$ID),"P.Value"]))
	Padj = cbind(Padj, as.matrix(T[match(rownames(Padj), T$ID),"adj.P.Val"]))
	rppafc = cbind(rppafc, as.matrix(T[match(rownames(rppafc), T$ID),"logFC"]))
	
}
colnames(P) = colnames(rppafc)=colnames(Padj) = KOs


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
pdf("ValidationBarplots.pdf")
for (i in 1:5){
	x=names(dsGenes)[i]
	y=dsGenes[[i]]
	y=intersect(colnames(newrppafc), y)
	barplot(newrppafc[as.character(x),as.character(y)], ylim=c(-3,3),xlab="Downstream proteins", ylab="log fold change", las = 2, main=paste("Knock-down for ", x, sep=" "))
}
dev.off()




MAPK1, MAPK3
AKT1, AKT2, AKT3
SMAD2, SMAD3
AKT1, AKT2, AKT3
RPS6KA1, RPS6KA2, RPS6KA3
STAT5A, STAT5B 