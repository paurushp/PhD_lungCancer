##########################################################
# Script to analyze Illumina Lung Cancer data from DKFZ  #
# 1. Preprocessing of data				 #
# 2. Limma analysis of data				 #
# 3. Prepare to rum NEM					 #
#							 #
#							 # 
##########################################################

library(lumi)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(limma)
DIR="/home/bit/praveen/BooleanNET/RupiData"
library(arrayQualityMetrics)
library(sva)
#pdat = read.csv(file.path(DIR, "phenoData", "sample_sheet_final.txt"),sep="\t")
#D=read.txt(file="")
#ann=colnames(D)[147:177] # get annotation column names
#dat.all = lumiR(file.path(DIR, "combinednew.txt"), columnNameGrepPattern = list(exprs='mean', se.exprs='sd',detection='detection', beadNum=NA), annotationColumn=ann, lib.mapping='lumiHumanIDMapping', dec=".", sep="\t")
#arrayQualityMetrics(dat.all2, outdir=file.path(DIR, "QualityControlRaw/"), do.logtransform=TRUE) # one chip as an outlier detected


file1=read.table("/home/bit/praveen/BooleanNET/RupiData/Illumina_Array_Knockdown_H1650_Serial1_all.qnorm_Kuner.txt", head=T, sep="\t") #Read file to look at data and extract column annotation
ann=colnames(s)[147:177]
dat.all = lumiR(file.path(DIR, "combinednew.txt"), columnNameGrepPattern = list(exprs='mean', se.exprs='sd',detection='detection', beadNum=NA),annotationColumn=ann, lib.mapping='lumiHumanIDMapping', dec=".", sep="\t") # Read the combined batch data
arrayQualityMetrics(dat.all, outdir=file.path(DIR, "QualityControlRaws/"), do.logtransform=TRUE) #Check the array quality
dat.norm = lumiExpresso(dat.all, varianceStabilize.param=list(method="vst")) # Do vst normalization
dat.normQ = lumiExpresso(dat.all, varianceStabilize.param=list(method="log2")) # Do log transformation
arrayQualityMetrics(dat.normQu, outdir=file.path(DIR, "QualityControlNorm/"), force=TRUE) #Checkk array quality

dat.normQu=lumiExpresso(dat.normQ, variance.stabilize=FALSE, normalize.param = list(method='quantile')) # Perform quantile normalization

ex=read.csv("/home/bit/praveen/BooleanNET/RupiData/Experiment_description_Knockdown_H1650_2012_Kuner.csv", head=T,sep="\t") #read pheno data
pData(dat.normQu)=ex # assign pheno data
mod = model.matrix(~Gene.Symbol, data=pData(dat.normQu)) # model the matrix to get rid of batch effect
mod = model.matrix(~factor(Gene.Symbol), data=pData(dat.normQu)) # model the matrix to get rid of batch effect
factor(pData(dat.normQu)$Cell.culture.experimental.batch) # check batches
ex2 = ComBat(exprs(dat.normQu), factor(pData(dat.normQu)$Cell.culture.experimental.batch), mod) # Run to combine batches removin batch effect
exprs(dat.normQu)=ex2 # Assign new data as expression

dat.nAP=detectionCall(dat.normQu, type="matrix") # Probe selection
take = c()
for(s in unique(pData(dat.normQu)$Hybridization.Name)){
     take = union(take, which(rowMeans(dat.nAP[, pData(dat.normQu)$Hybridization.Name==s, drop=F] == "P") == 1)) # take only probes, which are present in each patient for one group
}
dat.normal = dat.normQu[take,] # After selecting probes

########################################
# dim(dat.normal) #after filtering     #
# Features  Samples		       # 
# 47849       72 		       #
# dim(dat.normQu) 		       #
# Features  Samples   		       # 
# 48097       72 		       #
########################################


mappingInfo = nuID2EntrezID(featureNames(dat.normal), lib.mapping='lumiHumanIDMapping', returnAllInfo=TRUE) # Get mapping info for the data

dat.normal = dat.normal[mappingInfo[,"EntrezID"] != "",] # Get the data for Entrez ID
dat.normal1 = dat.normal[mappingInfo[,"EntrezID"] != NA,]

########################################
# dim(dat.normal)   #after selection   #
# Features  Samples 		       #
#   44502       72 		       #
########################################

mp = mappingInfo[mappingInfo[,"EntrezID"] != "",] # Filter-in data which have Entrez IDs

########################################
# Begin lima analysis
########################################
setwd("/home/bit/praveen/LunCanceRupi/Res")
sampleNames = paste("T", pData(dat.normal)$Hybridization.Name, sep=".")


###########################################################
#  Sample names for contrast matrix and design matrix     #
###########################################################
sn=as.character(sampleNames)
nam=c()
 for(i in 1:72){
 name=unlist(strsplit(sn[i],"_"))
 nam[i]=name[2]
 }


sn=as.character(sampleNames)
Cnam=c()
 for(i in 1:23){
 name=unlist(strsplit(colnames(design)[i],")"))
 Cnam[i]=name[2]
 }
cont=paste("Mock", Cnam, sep="-")

############################################################################
design = model.matrix(~ -1 + factor(nam),levels=unique(nam))
colnames(design)=Cnam

contrast=makeContrasts(cont[2], cont[3],cont[4],cont[5],cont[6],cont[7],cont[8],cont[9],cont[11],cont[13],cont[14],cont[15],cont[16],cont[17],cont[18],cont[19],cont[20],cont[21],cont[22],cont[23],levels=design)




############################################################################
colnames(design) = unique(nam)

fit = lmFit(dat.normal, design=design)
fit2 = contrasts.fit(fit, contrast)
fitE = eBayes(fit2)
Sel.Genes = topTable(fitE, number="Inf")
Sel.Genes$P.Value
length(unique(Sel.Genes$Entrez_Gene_ID))
fit = lmFit(dat.normal,design=design)
fit = eBayes(fit)


data=list()
ListGenes=list()
for(i in 1: ncol(contrast)){
Sel.Genes = topTable(fitE, number="Inf", coef=i)

num.5p=sum(Sel.Genes$adj.P.Val < 0.05, na.rm=T)
num.1p=sum(Sel.Genes$adj.P.Val < 0.01, na.rm=T)
num.5pFC=sum(Sel.Genes$adj.P.Val < 0.05 & abs(Sel.Genes$logFC) > 1, na.rm=T)
num.1pFC=sum(Sel.Genes$adj.P.Val < 0.01 & abs(Sel.Genes$logFC) > 1, na.rm=T)
genes.5p=Sel.Genes$Entrez_Gene_ID[1:num.5p]
genes.1p=Sel.Genes$Entrez_Gene_ID[1:num.1p]
genes.5pFC=Sel.Genes$Entrez_Gene_ID[1:num.5pFC]
genes.1pFC=Sel.Genes$Entrez_Gene_ID[1:num.1pFC]
PV.5p=Sel.Genes$adj.P.Val[1:num.5p]
PV.1p=Sel.Genes$adj.P.Val[1:num.1p]
PV.5pFC=Sel.Genes$adj.P.Val[1:num.5pFC]
PV.1pFC=Sel.Genes$adj.P.Val[1:num.1pFC]
Data.5p=data.frame(E.ID=genes.5p, P.Val=PV.5p)
Data.1p=data.frame(E.ID=genes.1p, P.Val=PV.1p)
Data.5pFC=data.frame(E.ID=genes.5pFC, P.Val=PV.5pFC)
Data.1pFC=data.frame(E.ID=genes.1pFC, P.Val=PV.1pFC)
data[[i]]=list(Data5p=Data.5p, Data1p=Data.1p, Data5pFC=Data.5pFC, Data1pFC=Data.1pFC)
ListGenes[[i]]=list(genes5p=genes.5p, genes1p=genes.1p, genes5pFC=genes.5pFC, genes1pFC=genes.1pFC)
}


for(i in 1:20){
	print(c(length(ListGenes[[i]]$genes5p),length(ListGenes[[i]]$genes1p),length(ListGenes[[i]]$genes5pFC),length(ListGenes[[i]]$genes1pFC)))
}
 

############################################################################
# get uniqe differentially expressed genes for each cut-off
############################################################################
genelist5p=c() 
genelist1p=c()
genelist1pFC=c()
genelist5pFC=c()
for(i in 1:20){
	genelist5p=c(genelist5p, ListGenes[[i]]$genes5p)
	genelist1p=c(genelist1p, ListGenes[[i]]$genes1p)
	genelist5pFC=c(genelist5pFC, ListGenes[[i]]$genes5pFC)
	genelist1pFC=c(genelist1pFC, ListGenes[[i]]$genes1pFC)
}
length(unique(genelist5p))
length(unique(genelist1p))
length(unique(genelist5pFC))
length(unique(genelist1pFC))

genes=unique(genelist1pFC)


Mat=matrix(0, nrow= length(genes), ncol=ncol(contrast))
for(i in 1:20){
Sel.Genes = topTable(fitE, number="Inf", coef=i)
Mat[,i]=Sel.Genes[Sel.Genes$Entrez_Gene_ID %in% genes,]$P.Value
##################################
# 43723 genes after removing NAs #
##################################

############################################################################
ProbeID=Sel.Genes$Probe_Id[1:43723]
dataFull2=matrix(0, nrow=length(ProbeID), ncol=20)
dataFull2=data.frame(dataFull2)
colnames(dataFull)=SgeneName
rownames(dataFull)=ProbeID
for(i in 1:20){
Sel.Genes = topTable(fitE, number="Inf", coef=i)
dataFull[,i]=Sel.Genes[Sel.Genes$Probe_Id %in% ProbeID,]$P.Value
}

DataF=getDensityMatrix(dataFull, dirname="/home/bit/praveen/LunCanceRupi/Res/DensityMat", startab=c(0.3,10), startlam=c(0.6,0.1,0.3), tol=1e-4)

for(i in 1:20){
Sel.Genes = topTable(fitE, number="Inf", coef=i)
dataFull[,i]=Sel.Genes[Sel.Genes$Probe_Id %in% ProbeID,]$P.Value
}

data.mRNA=cbind(Symbol=Sel.Genes[Sel.Genes$Probe_Id %in% ProbeID,]$Symbol, dataFull)
pertGene=(colnames(data.mRNA)[2:ncol(data.mRNA)])
genes.inter=intersect(data.mRNA$Symbol, pertGene)
selected.data=data.mRNA[Sel.Genes$Symbol %in% genes.inter, ]


dataFull2=cbind(Symbol=glist,dataFull1)
Sel.Genes[Sel.Genes$Probe_Id %in% rownames(dataFull),]$Symbol

dataFull2.rel=dataFull2[dataFull2$Symbol %in% pertGene, ]


ProbeID=Sel.Genes$Probe_Id[1:43723]
dataLFC=matrix(0, nrow=length(ProbeID), ncol=20)
dataLFC=data.frame(dataLFC)
colnames(dataLFC)=SgeneName
rownames(dataLFC)=ProbeID
for(i in 1:20){
Sel.Genes = topTable(fitE, number="Inf", coef=i)
dataLFC[,i]=Sel.Genes[Sel.Genes$Probe_Id %in% ProbeID,]$logFC
}
Gsm=Sel.Genes[Sel.Genes$Probe_Id %in% ilID,]$Symbol
ilID=Sel.Genes[Sel.Genes$Symbol %in% colnames(dataLFC),]$Probe_Id

dataLFC2=matrix(0, nrow=15, ncol=20)
for(i in 1:20){
Sel.Genes = topTable(fitE, number="Inf", coef=i)
dataLFC2[,i]=Sel.Genes[Sel.Genes$Entrez_Gene_ID %in% newGe,]$logFC
}

Gsm=Sel.Genes[Sel.Genes$Probe_Id %in% ilID2,]$Symbol
ilID2=Sel.Genes[Sel.Genes$Entrez_Gene_ID %in% newGe,]$Probe_Id
ilID2=Sel.Genes[Sel.Genes$Entrez_Gene_ID %in% newGe,]$Probe_Id

dataforLM=data.frame()#(0, nrow=15, ncol=20)
for(i in 1:20){
Sel.Genes = topTable(fitE, number="Inf", coef=i)
dataforLM=rbind(dataforLM, Sel.Genes[Sel.Genes$Symbol %in% SgeneName, ]$adj.P.Val)
}

#############################################
#
#PLOT Heat Maps
#############################################
DataF.frame=data.frame(DataF)
rownames(DataF.frame)=ProbeID
DD=DataF.frame[rownames(DataF.frame) %in% Probe_ID,]
DD=data.matrix(DD)
rownames(DD)=Probe_ID
heatmap(DD)
ag=clustering(DD)
clustering.plot(ag, data=DD)
pdf("clustering.pdf", height=6, width=12)
clustering.plot(ag, data=DD)
dev.off()

DFC=dataLFC[rownames(dataLFC) %in% Probe_ID,]
DFC=data.matrix(DFC)
#############################################

for(i in 1:20){
Sel.Genes = topTable(fitE, number="Inf", coef=i)
EntrezID=Sel.Genes[Sel.Genes$Entrez_Gene_ID %in% genes,]$Entrez_Gene_ID
ProbeID=Sel.Genes[Sel.Genes$Entrez_Gene_ID %in% genes,]$Probe_Id
P.Value=Sel.Genes[Sel.Genes$Entrez_Gene_ID %in% genes,]$P.Value
}
#design = model.matrix(~ -1 + factor(c(rep(rep(1,3),72))))
#colnames(design)=c("NTC", "MOCK", "tumor")

#fit=lmFit(dat.normal, design=design)
#fit2=contrasts.fit(fit, contrast)
#fitE=eBayes(fit2)
#Sel.Genes=topTable(fitE, number="1000")
Sel.Genes=topTable(fitE, number=1000, adjust="BH",p.value=0.01)
# 267 genes at 0.01
# topTable(fit, number="Inf", coef=Cnam, adjust="BH", p.value=0.05,sort.by ="F")

###########################################################################
###########################################################################
# Convert Symbols to EIDs
#
###########################################################################
###########################################################################
IDs=mget(newSet,org.Hs.egALIAS2EG, ifnotfound=NA)
save(testGenes, IDs, geneSet, file='sgeneData.RData')


###################
dat.norm = lumiExpresso(dat.all[,-1],varianceStabilize.param=list(method="vst")) # uses log2 transformation here, because BEAD_STDERR is missing!
arrayQualityMetrics(dat.norm, outdir=file.path(DIR, "QualityControlNorm/"), force=TRUE)

pData(dat.norm)=pdat[-1,]
dat.nAP = detectionCall(dat.norm, type="matrix")
take = c()
for(s in unique(pData(dat.norm)$type)){
     take = union(take, which(rowMeans(dat.nAP[, pData(dat.norm)$type ==s, drop=F] == "P") == 1)) # take only probes, which are present in each patient for one group
}
dat.norm = dat.norm[take,]

mappingInfo = nuID2EntrezID(featureNames(dat.norm), lib.mapping='lumiHumanIDMapping', returnAllInfo=TRUE)
dat.norm = dat.norm[mappingInfo[,"EntrezID"] != "",] # only consider probes mapping to Entrez gene IDs
mp = mappingInfo[mappingInfo[,"EntrezID"] != "",]

# limma-analysis
library(limma)
sampleNames = paste("T", pData(dat.norm)$type, sep=".")
design = model.matrix(~ -1 + factor(sampleNames),levels=unique(sampleNames))
colnames(design) = unique(sampleNames)
fit = lmFit(dat.norm,design=design)
contrasts = makeContrasts("T.500 - T.200", "T.500 - T.hBC", "T.200 -T.hBC", levels=unique(sampleNames))
fit2 = contrasts.fit(fit, contrasts)
fitE = eBayes(fit2)


#######################################


Ser1=D[,1:177]
Ser2=D3[,180:323]
D4=cbind(Ser1,Ser2)


DIR="/home/bit/praveen/LunCanceRupi/Data"

dat.allCom=lumiR(file.path(DIR, "combined.txt"), columnNameGrepPattern = list(exprs='mean', se.exprs='sd',detection='p', beadNum=NA), annotationColumn=ann, lib.mapping='lumiHumanIDMapping', dec=".", sep="\t")
DIR="/home/bit/praveen/LunCanceRupi/Res"
arrayQualityMetrics(dat.allCom, outdir=file.path(DIR, "QualityControlRaw/"), do.logtransform=TRUE)




