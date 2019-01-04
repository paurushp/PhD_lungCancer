edgeLab=matrix(0, nrow(NemMat), ncol(NemMat))
for(i in 1:nrow(NemMat)){
	for(j in 1:ncol(NemMat)){
		r=rownames(NemMat)[i]
		x=colnames(NemMat)[j]
		Ddat=data.frame(cbind(Da=as.character(DFC2[,1]), Db=DFC2[,r]))
		COL=subset(Ddat, Ddat$Da!=r)
		COL2=subset(COL, DFC2$Da==x)
		S=mean(COL2)
		if(S<0){
			edgeLab[i,j]= -1
			}
		if(S>0){
			edgeLab[i,j]= 1
			}
	}
}

		
toList=(unique(as.character(DFC2$Gsm)))
EdgeLab=matrix(0, 20, 16)				
for(i in 1:20){
	for(j in 1:16){
		r=c()
		DFC3=DFC2[,c(1,i)]
		s=c(r, DFC3[DFC3$Gsm %in% toList[j],][,2]))
		EdgeLab[i,j]=mean(s)
	}
}
edgeLa=function(mat){
	for(i in 1:20){
		for(j in 1:16){
		if(mat[i,j]>0)
			mat[i,j]=1
		if(mat[i,j]<0)
			mat[i,j]=-1
		return(mat)
}
}

thresholding4=function(matr, thresh){
M = matr
M [matr > thresh] <- 1
M[matr < thresh] <- -1
return(M)
# > M2 <- M
# > M2[M < thresh] <- 0
# > M2[M >= thresh] <- 1
}



EdgeLab=data.frame(EdgeLab)
colnames(EdgeLab)=toList
rownames(EdgeLab)=colnames(DFC2)[2:21]
for(i in 1:length(toList))
print(toList[i])


egdesL=names(edgeData(NEMG))
for(i in 55:59){
s=strsplit(egdesL[i], "\\|")
ED[[i]]$label=et[s[[1]][1],s[[1]][2]]
}

length(egdesL)

