known = read.csv("metacoreNetSgene.csv", sep="\t")
GMetaRupi = new("graphNEL", edgemode="directed", nodes=union(as.character(known[,2]), as.character(known[,4])))
GMetaRupi = addEdge(as.character(known[,2]), as.character(known[,4]), GMetaRupi)
sgl=read.table("/home/bit/praveen/LunCanceRupi/sgenelistName.txt")
sgl2=c("Bcl-10", "EGFR", "Separase", "GSK3 alpha", "GSK3 beta", "ITGB4", "Leptin receptor", "mTOR", "PI3K cat class III (Vps34)","AMPK alpha 1 subunit", "AMPK beta subunit", "c-Raf-1", "p90RSK1", "p70 S6 kinases", "c-Src", "LKB1", "LKB1", "TRUB2","Hamartin", "Tuberin", "WDR3")


##### Calculate significance
load("/home/bit/praveen/LunCanceRupi/Res/NEM/rupiInputNEM.RData")
pall=denMat.Rupi
library(doMC)
p1 = nem.calcSignificance(pall, res5) # [1] 0.0005809816
print(p1)
p1 = nem.calcSignificance(pall, res.prior2) #[1] 0.0005809816
print(p1)
save(net, p, file=file.path(DIR, "../AllGenes.rda"))

###### change node names corresponding to metacore

sgl2=c("Bcl-10", "EGFR", "Separase", "GSK3 alpha", "GSK3 beta", "ITGB4", "Leptin receptor", "mTOR", "PI3K cat class III (Vps34)","AMPK alpha 1 subunit", "AMPK beta subunit", "c-Raf-1", "p90Rsk", "p70 S6 kinase2", "c-Src", "LKB1", "TRUB2","Hamartin", "WDR3", "Tuberin")

nodes(res5$graph)=sgl2
nodes(res.prior2$graph)=sgl2

###### Calculate the comparison with metacore

intersect(nodes(GMetaRupi), sgl2)
nodes(GMetaRupi)
mm=dkSP(GMetaRupi, sgl2)
paths=sp.path.compute(sgl2, G=GMetaRupi)

# a) How many edges in inferred network can be explained by the literature?
em = edgeMatrix(NemNet)
eval = data.frame()
for(i in 1:NCOL(em)){
	cand.path = paste(nodes(NemNet)[em["from", i]], "->",nodes(NemNet)[em["to", i]], sep="")
	path = sp.between(GMetaRupi, nodes(NemNet)[em["from", i]], nodes(NemNet)[em["to", i]])
	if(!is.na(path[[1]]$length)){
		eval = rbind(eval, data.frame(edge=cand.path, explained.by=paste(path[[1]]$path_detail, collapse="->")))
	}
	else
		eval = rbind(eval, data.frame(edge=cand.path, explained.by="none"))
}
library(xlsx)
write.xlsx(eval, file="EvaluationAllFPR2.xlsx")


# b) Which part of the literature is explained by the inferred network?
paths2 = paths[sapply(paths, function(p) !is.na(p$length))]
eval2 = data.frame()
for(p in 1:length(paths2)){
	mypath = names(paths2)[p]
	mynodes = strsplit(mypath, ":")[[1]]
	path = sp.between(NemNet, mynodes[1],  mynodes[2])
	if(!is.na(path[[1]]$length)){
		eval2 = rbind(eval2, data.frame(literature.path=paste(mynodes[1], "->", mynodes[2], sep=""), explained.by=paste(path[[1]]$path_detail, collapse="->")))
	}
	else{
		eval2 = rbind(eval2, data.frame(literature.path=paste(mynodes[1], "->", mynodes[2], sep=""), explained.by="none"))
	}
}
write.xlsx(eval2, file="EvaluationAllTPR2.xlsx")


a	b
thresholding=function(matr, thresh){
M = matr
M[matr < thresh] <- 0.0
return(M)
}


