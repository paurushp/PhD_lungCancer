load("/home/abidata/datasets/ModelSystems/NSLC/Kuner_Data_June2012/EgenesRPPA.rda")
egenes[50]="EGFR" # Originally "EGFR     "
egenes[56]="PTEN" # Originally "PTEN "
# "RAS" NO ID FOUND so NA
library(org.Hs.eg.db)
EgeneIDs=mget(egenes, org.Hs.egALIAS2EG, ifnotfound=NA)
for(i in 1:length(EgeneIDs)){
EgeneIDs[[i]]=EgeneIDs[[i]][1]
}
EG=as.numeric(EgeneIDs)
load("/home/bit/praveen/prior.RData")
SG=rownames(s.prior)
testGenes=c(SG, EG)
