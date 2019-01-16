mygraph=trialNEt[[1]]$graph
nodes(mygraph)=nodes(mynemNOMBS.388)
tau=seq(0.0,1, 0.1)
DIR="/Volumes/MacintoshHD1/FinalResults/newGraphTCNOM/"
setwd(DIR)
exp2=data.frame()
for(i in 1:length(tau)){
	M=transitive.reduction((as(mygraph, "matrix")>tau[i])*1)
	#M2=as(trialNEt[[1]]$graph, "matrix")
	#M2[M2>tau[i]]<-1
	pdf(paste("Histogram for Threshold =",tau[i], ".pdf",sep="" ))
	hist(M)
	dev.off()
	G_here=as(M, "graphNEL")
	#G_here2=as(M2, "graphNEL")
	pdf(paste("_for Threshold =",tau[i], ".pdf",sep="" ))
	plot(G_here)
	dev.off()
	fileL=paste("Boot-Network-by-the-Literature for threshold=", tau[i],".csv", sep="" )
	fileM=paste("Boot-Literature-by-the-Network for threshold=", tau[i],".csv", sep="" )
	cat("Validating Model against Literatur \n")
	Valid_Lit=validateLit(G_here)
	allNone=sum(Valid_Lit$explained.by=="none")
	allNF=sum(Valid_Lit$explained.by=="NF")
	allEnt=nrow(Valid_Lit)-allNF
	allExp=nrow(Valid_Lit)-(allNone+allNF)
	Frac1=(allExp/allEnt)*100
	cat("Validating Literature against Model \n")
	Valid_Mode=validateMod(G_here)
	allNone=sum(Valid_Mode$explained.by=="none")
	allNF=sum(Valid_Mode$explained.by=="NF")
	allEnt=nrow(Valid_Mode)-allNF
	allExp=nrow(Valid_Mode)-(allNone+allNF)
	Frac2=(allExp/allEnt)*100
	exp2=rbind(exp2, data.frame(Threshold=tau[i], Model.View=Frac1,Know.View=Frac2))
	write.csv(Valid_Lit, file=fileL)
	write.csv(Valid_Mode, file=fileM)
}


DIR="/Volumes/MacintoshHD1/FinalResults/newGraphNTNOM/"
tau=seq(0.0,1, 0.1)
sp.path.compute.R
setwd(DIR)
exp3=data.frame()
for(i in 1:length(tau)){
	#M=transitive.reduction((as(mygraph, "matrix")>tau[i])*1)
	M=as(mygraph, "matrix")
	M[M<tau[i]]<-0
	M[M>tau[i]]<-1
	pdf(paste("Histogram for Threshold =",tau[i], ".pdf",sep="" ))
	hist(M)
	dev.off()
	G_here=as(M, "graphNEL")
	#G_here2=as(M2, "graphNEL")
	pdf(paste( names(allFinalBS)[k], "_for Threshold =",tau[i], ".pdf",sep="" ))
	plot(G_here)
	dev.off()
	fileL=paste("Boot-Network-by-the-Literature for threshold=", tau[i],".csv", sep="" )
	fileM=paste("Boot-Literature-by-the-Network for threshold=", tau[i],".csv", sep="" )
	cat("Validating Model against Literatur \n")
	Valid_Lit=validateLit(G_here)
	allNone=sum(Valid_Lit$explained.by=="none")
	allNF=sum(Valid_Lit$explained.by=="NF")
	allEnt=nrow(Valid_Lit)-allNF
	allExp=nrow(Valid_Lit)-(allNone+allNF)
	Frac1=(allExp/allEnt)*100
	cat("Validating Literature against Model \n")
	Valid_Mode=validateMod(G_here)
	allNone=sum(Valid_Mode$explained.by=="none")
	allNF=sum(Valid_Mode$explained.by=="NF")
	allEnt=nrow(Valid_Mode)-allNF
	allExp=nrow(Valid_Mode)-(allNone+allNF)
	Frac2=(allExp/allEnt)*100
	exp3=rbind(exp3, data.frame(Threshold=tau[i], Model.View=Frac1,Know.View=Frac2))
	write.csv(Valid_Lit, file=fileL)
	write.csv(Valid_Mode, file=fileM)
}






