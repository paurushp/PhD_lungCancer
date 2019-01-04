library(graphite)
library(igraph)
library(RBGL)
nsclc=convertIdentifiers(kegg$"Non-small cell lung cancer", type="symbol")
stringnet=read.csv("/home/praveen/Paurush/NSLC/stringnetAll.txt", head=TRUE, sep="\t") 
gdata=data.frame(from=stringnet[,1], to=stringnet[,2],weight=stringnet[,15])
g = graph.data.frame(gdata, directed=FALSE)
g2=igraph.to.graphNEL(g)


cancer_cen=c("AKT1", "ALK", "BAP1", "BRAF", "CCDC6", "CD74", "EGFR", "EML4", "ERBB2", "EZR", "FGFR2", "HIP1", "KDR", "KIF5B", "LRIG3", "MAP2K1", "MAP2K2", "NFE2L2", "NKX2-1", "NRG1", "RET", "ROS1", "SDC4", "SLC34A2", "SMARCA4", "SOX2", "STK11", "TFG", "TPM3", "TPR")


sp.path.compute = function(genelist, genelist2=NULL, G,verbose=TRUE){
	require(RBGL)
	if(is.null(genelist2)){
		paths = c()
		pathL = matrix(0, nrow=length(genelist), ncol=length(genelist))
		pathLweight = matrix(0, nrow=length(genelist), ncol=length(genelist))	
		for(i in 1:length(genelist)){			
			for(j in 1:length(genelist)){
				if(i != j){
					path = sp.between(G, genelist[i], genelist[j])					
					paths = c(paths, path)
					if(!is.na(path[[1]]$length)){
						pathL[i,j] = length(path[[1]][[3]][[1]])
						pathLweight[i,j] = path[[1]][[1]]				
					}	
				          
				}
				if(verbose==TRUE){				
					cat(".")	
				}		
			}	
			if(verbose==TRUE){			
				cat(i)	
			}
		}
		res=list(paths, pathL, pathLweight)
		return(res)				
	}	
	else{		
		paths = sp.between(G, genelist, genelist2)
		return(paths)
	}
}

pp=sp.path.compute(genelist=nodes(Gs7), G=Gs7)
pp_sgene=sp.path.compute(genelist=Sgene, G=Gs7)
pp_kegg=sp.path.compute(genelist=nsclcgene, G=Gs7)
pp_cancerCen=sp.path.compute(genelist=cancercen, G=Gs7)

hippie=read.csv("/Volumes/MacintoshMyApp/Downloads/hippie_current.txt", sep="\t", head=FALSE)
hippie=hippie[,c(1:5)]
x=which(hippie[,3]=="")
y=which(hippie[,1]=="")
hippie_ed=hippie[-c(x,y),]
nodes=unique(c(as.character(hippie_ed[,1]), as.character(hippie_ed[,3])))
G = new("graphNEL", nodes=nodes, edgemode="directed")
G = addEdge(as.character(hippie_ed[,1]), as.character(hippie_ed[,3]), G)
nodes(G)=sapply(strsplit(nodes(G), "_"), "[[", 1)
pp_sgene=sp.path.compute(genelist=intersect(Sgene, nodes(G)), G=G)




pdetail=data.frame(Path=factor(),Path.Literature=factor())

for(i in 1:length(pp_sgene[[1]])){
	#print(paste(pp_sgene[[1]][[i]][[2]],sep="->"))
	pdetail=rbind(pdetail, data.frame(Path=factor(names(pp_sgene[[1]])[[i]]),Path.Literature=paste(pp_sgene[[1]][[i]][[2]],collapse="->")))
}

dimnames(pp_sgene[[2]])list(Sgene,Sgene)

new_permutation_test_network=function(net, sgenes, permut=10000, refSP){
	mat=matrix(0,nrow=length(nodes(net)),ncol=length(nodes(net)))
	rownames(mat)=colnames(mat)=nodes(net)
  #     print(mat)	
	for(i in 1:permut){
		pbar=100*i/permut
		if(pbar%%1==0|pbar%%0.5==0){	
			cat(pbar)
		}
		net_here=permute.net(net)
		net_here = as(net_here, "graphNEL")
		net=net_here
		sp_mat_here=sp.path.compute(genelist=nodes(net_here), genelist2=NULL, net_here,verbose=FALSE)[[2]]
		#sp_mat_here[sp_mat_here==0]<-1/0
		dimnames(sp_mat_here)=dimnames(mat)
cat(".")
		#print(sp_mat_here[sgenes,sgenes] - refSP[sgenes,sgenes])
		for(j in 1:nrow(sp_mat_here)){
			for(k in 1:nrow(sp_mat_here)){
				if(sp_mat_here[sgenes[j],sgenes[k]] < refSP[sgenes[j],sgenes[k]]){
					mat[sgenes[j],sgenes[k]]=mat[sgenes[j],sgenes[k]]+1				
				}	
			}
		}
		
	}
dimnames(mat)=list(sgenes,sgenes)
return(mat)

}

load("finalNet.RData")
ptest=new_permutation_test_network(net=myres$graph, sgenes=nodes(myres$graph), permut=10000, refSP)

ptest=new_permutation_test_network(net=myres$graph, sgenes=nodes(myres$graph), permut=10000, pp_sgene[[2]])
ptest3=new_permutation_test_network(net=gg2, sgenes=nodes(myres$graph), permut=5000, pp_sgene[[2]])

ptest3[Sgene,Sgene]


sp_mat_here1=sp.path.compute(genelist=nodes(myres$graph), genelist2=NULL, myres$graph,verbose=FALSE)[[2]]
permute.net = function(net){ # Permute the graph
	genes=nodes(net)
	LV=as(net,"matrix") 
	LV[LV >= 0.5]<-1       
	operation=c(0,0,0)        # operations possible for the node
        newLV=LV                # updated matrix
        i=sample(1:length(genes), 1)
        j=sample(1:length(genes), 1) # randomly select indices
        if(i!=j){
                if(LV[i,j] == 0)
                        operation[1] = 1
                if(LV[i,j] == 1 )
                        operation[2] = 1
                if(LV[i,j]!= newLV[j,i])
                        operation[3] = 1
                op = sample(1:3, 1)        # Random selection of operation
                if(operation[op]==1 && op==1 )
                        newLV [i,j]= 1 # add edge
                if(operation[op]==1 && op==2)
                        newLV [i,j]= 0 # remove edge
                if(operation[op]==1 && op==3){
                        tempLV= newLV[i,j]
                        newLV[i,j] = newLV[j,i] # swap direction
                        newLV[j,i] = tempLV
                }
        }
#print(newLV-LV)
        return (newLV)

}

ptest=new_permutation_test_network(net=myres$graph, sgenes=nodes(myres$graph), permut=100, refSP)

permute.net = function(net){
	genes=nodes(net)
	LV=as(net,"matrix") 
	LV[LV >= 0.5]<-1 
	newLV=LV      
        i=sample(1:length(genes), 1)
        j=sample(1:length(genes), 1) # randomly select indices
	if(i!=j){
		if(LV[i,j] == 0){
			newLV[i,j]=1	
		}
		if(LV[i,j]==1){
			newLV[i,j]=0
		}
	}
	return (newLV)
}

pdf("esdfKEGG.pdf")
plot(ecdf(c(pp_kegg[[2]])))
dev.off()
pdf("esdfsgenes.pdf")
plot(ecdf(c(pp_sgene[[2]])))
dev.off()
pdf("esdfcancerCen.pdf")
plot(ecdf(c(pp_cancerCen[[2]])))
dev.off()

Sgene=c("PRKAA1","PRKAB1","ESPL1","WDR3","RAF1","ITGB4","SRC","MTOR","PIK3C3","TSC1","TSC2","BCL10","EGFR","STK11","GSK3A","GSK3B","TRUB2","LEPR","RPS6KA1","RPS6KB1")
nsclcgene=c(
"AKT3","CDK4","CDK6","CDKN2A","RASSF1","E2F1","E2F2","EGF","EGFR","ERBB2","AKT1","AKT2","PIK3R5","GRB2","HRAS","ARAF","KRAS","NRAS","PDPK1","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PLCG1","PLCG2","PRKCA","PRKCB","PRKCG","MAPK1","MAPK3","MAP2K1","MAP2K2","RAF1","RARB","CCND1","RXRA","RXRB","SOS1","SOS2","BRAF","TGFA","RASSF5","PIK3R3","FOXO3","BAD","CASP9","RB1","STK4","RXRG")
cancercen=c("AKT1","ALK","BAP1","BRAF","CCDC6","CD74","EGFR","EML4","ERBB2","EZR","FGFR2","HIP1","KDR","KIF5B","LRIG3","MAP2K1","MAP2K2","NFE2L2","NKX2-1","NRG1","RET","ROS1","SDC4","SLC34A2","SMARCA4","SOX2","STK11","TFG","TPM3","TPR")



enrich=function(gr, annot, geneset){
x=list()
y=list()
for(i in 1:length(nodes(gr$graph))){
	redGraph=as(transitive.reduction(as(gr$graph, "matrix")),"graphNEL")
	downstream=dfs(redGraph, nodes(gr$graph)[i], checkConn=TRUE)$finish
	egenes=unlist(gr$mappos[downstream])
	egenes=annot[annot$Probe_Id %in% as.character(egenes),]$Symbol
	allg=union(downstream,egenes)
	#print(allg)
	en=length(intersect(allg,geneset))
	cat(nodes(gr$graph)[i])	
	x[[i]]=length(allg)
	y[[i]]=en
		
	print(en)
}
names(x)=names(y)=nodes(gr$graph)
return(list(x,y))
}

#enrich(myres, annot, nsclcgene)
#enrich(myres, annot, cancercen)
keggenrich=enrich(myres, annot, nsclcgene)
cancerenrich=enrich(myres, annot, cancercen)

enrichtcga=function(mat,gene,geneset){
p=list()
for(i in 1:ncol(dematrix)){
degenes=sapply(strsplit(DEgene[[i]], "[|]"), "[[",1)
p[[i]]=length(intersect(degenes,geneset))
}
names(p)=colnames(dematrix)
return(p)
}


keggtcga=enrichtcga(dematrix,DEgene,nsclcgene)
cancercentcga=enrichtcga(dematrix,DEgene,cancercen)


barplot(unlist(keggtcga), xlab = "somatic mutation genes",  ylab= "#DE genes shared with KEGG NSCLC pathway")

barplot(unlist(cancercentcga), xlab = "somatic mutation genes",  ylab= "#DE genes shared with Cancer Census data")

for(i in 1:length(keggenrich[[1]])){
x=c(x, length(keggenrich[[1]][[1]])
y=
)
}

cl.sample = clustering(attachprob, method="ward", metric="pearson")
cl.gene = clustering(t(attachprob), method="ward", metric="pearson")
pdf("heatmapAttachProb.pdf")
clustering.plot(tree=cl.sample, tree.sup=cl.gene, data=attachprob, names.sup=FALSE, scale="none", lab=data.frame(colnames(attachprob)), title="P-values densities")


plotEffects(DataF,myres,border=TRUE,legend=TRUE,order=NULL,orderSCC=FALSE,palette="BlueRed")



cl.sample = clustering(attachprob, method="ward", metric="pearson")
cl.gene = clustering(t(attachprob), method="ward", metric="pearson")
pdf("attachprob.pdf")
clustering.plot(tree=cl.sample, tree.sup=cl.gene, data=attachprob, names.sup=FALSE, scale="none", lab=data.frame(colnames(attachprob)), title="Attachment probability")



