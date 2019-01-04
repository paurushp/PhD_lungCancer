D2=read.table(file="/home/abidata/datasets/ModelSystems/NSLC/Kuner_Data_June2012/Illumina_Array_Knockdown_H1650_Serial2_all.qnorm_Kuner.txt", head=T, sep="\t")
Sg1=colnames(D)[3:146]
name1=list()
for(i in 1:144){
nam=unlist(strsplit(Sg1[i], "_"))
name1[[i]]=nam[2]
name1=unlist(name1)
}

Sg2=colnames(D2)[3:146]
name2=list()
for(i in 1:144){
nam=unlist(strsplit(Sg2[i], "_"))
name2[[i]]=nam[2]
name2=unlist(name2)
}

# ######################################
# S-gene list from q-normalized data files

SGs1=unique(name1)
SGs2=unique(name2)
S.G=union(SGs1,SGs2)

# S.G=S.G[3:22] # first two names removed

# ######################################
# S gene list from the Experiment description txt file

D3=read.table(file="/home/abidata/datasets/ModelSystems/NSLC/Kuner_Data_June2012/Experiment_description_Knockdown_H1650_2012_Kuner.txt", head=T, sep="\t")
colnames(D3)
S.Genes=unique(D3[,4])[3:22]

# ######################################
# Get Entrez ID for S genes
#
library(org.Hs.eg.db)
get(S.G[i], org.Hs.egENTREZID)
