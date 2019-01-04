
cd /home/praveen/Paurush/COSBI/Clustering/DATA/TGCA2
~/MyR

cli_aliquot_gbm<-read.csv(file='/home/praveen/Paurush/COSBI/Clustering/DATA/TGCA2/Clinical/Biotab/nationwidechildrens.org_biospecimen_aliquot_gbm.txt', head=TRUE, sep="\t")
dim(cli_aliquot_gbm)
colnames(cli_aliquot_gbm)
DIR="/home/praveen/Paurush/COSBI/Clustering/DATA/TGCA2/Clinical/Biotab/"
Clinical<-list.files(DIR)
clinicalFeatures=list()
for(i in 1:length(Clinical)){
	clinicalFeatures[[i]]=read.csv(file=paste(DIR,Clinical[[i]],sep=""),head=TRUE, sep="\t")
}
intersect(clinicalFeatures[[1]][,1], clinicalFeatures[[2]][,1])

patients=intersect(clinicalFeatures[[1]][,1], clinicalFeatures[[2]][,1])
clinF=clinicalFeatures[-c(3,4,10:15)]
cf=merge(clinF[[1]],clinF[[2]], by = 1, incomparables = NA)
for(i in 3:length(clinF)){
cf=merge(cf,clinF[[3]], by = "SampleID", incomparables = NA)
}

cf[,-c(5,9,14, 15,21, 28:32, 36:40, 43:47, 50:54)]
cf2=cf[, c(4,7,17,20,25, 26, 36)]
