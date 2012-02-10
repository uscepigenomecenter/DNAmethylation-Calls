## Author: Houtan Noushmehr, PhD
## houtana@gmail.com
## 310.570.2362
## date: Jan 12, 2012
## first of 3 parts to make necessary calls.
## Probe classification by gene level

## first, categorize each gene and probe to three classes
## SNC (strong negative correlation), WNC (weak negative correlation), NNC (no negative correlation).

load(file="CpGProbes_collapsed_to_one_gene.rd")
genenames$correlation_call <- c("UNK")
for(i in 1:dim(genenames)[1]){
	if(genenames[i,"corvalue_round"] < c(-0.4)){

	genenames[i,"correlation_call"]<-c("SNC")
} else {
if(genenames[i,"corvalue_round"] >= c(-0.4) && genenames[i,"corvalue_round"] <= c(-0.2) ){
	
	genenames[i,"correlation_call"]<-c("WNC")
} else {
genenames[i,"correlation_call"]<-c("NNC")
}


}
}

genenames[,"correlation_call"] <- as.factor(genenames[,"correlation_call"])

dimnames(genenames)[[1]] <- genenames[,"probeID"]
### categorize each gene to Tumor Probe to 10th, 50th, and 90th percentile and the same for the normals
## tumors
load(file="dataFreeze_12_15_2011_dnamethylation_geneexpression.rda")

dimnames(y.drs)[[1]] <- y.drs[,2]
dat.t <- (y.drs[as.character(genenames[,"probeID"]),])
tum.per <- apply(dat.t[,3:281], 1, quantile, c(.1,.5,.9), na.rm=T)
tum.per <- t(tum.per)
dimnames(tum.per)[[2]] <- c("T10","T50", "T90")
## normals
gbm.normal <- read.delim(file="non_tumor_brain_noushmehr_et_al_2010.txt")
dimnames(gbm.normal)[[1]] <- gbm.normal[,1]
gbm.normal <- gbm.normal[,-1]
dat.n <- (gbm.normal[as.character(genenames[,"probeID"]),])
nor.per <- apply(dat.n, 1, quantile, c(.1,.5,.9), na.rm=T)
nor.per <- t(nor.per)
dimnames(nor.per)[[2]] <- c("N10","N50", "N90")


genenames <- cbind(genenames, tum.per, nor.per)


### now, make conditions in the following manner

genenames$N.METH <- c("UNK")
genenames$T.METH <- c("UNK")
## for normals:


for(n in 1:dim(genenames)[1]){
if(is.na(genenames[n, "N10"] < 0.25 && genenames[n, "N50"] < 0.25 && genenames[n, "N90"] < 0.25)){
	genenames[n, "N.METH"] <- c("UNK")
} else{

	if(genenames[n, "N90"] < 0.25){

	genenames[n, "N.METH"] <- c("CUN")


} else {
	if(genenames[n, "N10"] > 0.75){

	genenames[n, "N.METH"] <- c("CMN")

} else {

if(genenames[n, "N10"] > 0.25 && genenames[n, "N90"] < 0.75){
	
	genenames[n, "N.METH"] <- c("IMN")
} else {
genenames[n, "N.METH"] <- c("VMN")

}


}


}
}
}#endfor

genenames[,"N.METH"] <- as.factor(genenames[,"N.METH"])



## for tumors
for(n in 1:dim(genenames)[1]){
if(is.na(genenames[n, "T10"] < 0.25 && genenames[n, "T50"] < 0.25 && genenames[n, "T90"] < 0.25)){
	genenames[n, "T.METH"] <- c("UNK")
} else{

	if(genenames[n, "T90"] < 0.25){

	genenames[n, "T.METH"] <- c("CUT")


} else {
	if(genenames[n, "T10"] > 0.75){

	genenames[n, "T.METH"] <- c("CMT")

} else {

if(genenames[n, "T10"] > 0.25 && genenames[n, "T90"] < 0.75){
	
	genenames[n, "T.METH"] <- c("IMT")
} else {
genenames[n, "T.METH"] <- c("VMT")

}


}


}
}
}#endfor

genenames[,"T.METH"] <- as.factor(genenames[,"T.METH"])



table(genenames[,"N.METH"],genenames[,"T.METH"], genenames[,"correlation_call"])
table(genenames[,"T.METH"],genenames[,"correlation_call"])
table(genenames[,"N.METH"],genenames[,"correlation_call"])
 table(genenames[,"N.METH"],genenames[,"T.METH"])





### 450
## Probe classification by gene level

## first, categorize each gene and probe to three classes
## SNC (strong negative correlation), WNC (weak negative correlation), NNC (no negative correlation).

load(file="CpGProbes_collapsed_to_one_gene_450.rd")
genenames.450$correlation_call <- c("UNK")
for(i in 1:dim(genenames.450)[1]){
	if(genenames.450[i,"corvalue_round"] < c(-0.4)){

	genenames.450[i,"correlation_call"]<-c("SNC")
} else {
if(genenames.450[i,"corvalue_round"] >= c(-0.4) && genenames.450[i,"corvalue_round"] <= c(-0.2) ){
	
	genenames.450[i,"correlation_call"]<-c("WNC")
} else {
genenames.450[i,"correlation_call"]<-c("NNC")
}


}
}

genenames.450[,"correlation_call"] <- as.factor(genenames.450[,"correlation_call"])

dimnames(genenames.450)[[1]] <- genenames.450[,"probeID"]
### categorize each gene to Tumor Probe to 10th, 50th, and 90th percentile and the same for the normals
## tumors
load(file="dataFreeze_12_15_2011_dnamethylation_geneexpression_450K_a.rda")

dimnames(y.drs.450)[[1]] <- y.drs.450[,2]
dat.t.450 <- (y.drs.450[as.character(genenames.450[,"probeID"]),])
tum.per.450 <- apply(dat.t.450[,3:76], 1, quantile, c(.1,.5,.9), na.rm=T)
tum.per.450 <- t(tum.per.450)
dimnames(tum.per.450)[[2]] <- c("T10","T50", "T90")
## normals
load(file="md.450k.n_NOfieldEffect.rda") ## received from Swampna
load(file="KIRC_Normals_for_Houtan.rda") ## file received from Hui on Jan 11, 2012
load(file="LUSC_normals_450.rda") ## file downloaded on Jan 11, 2012
dimnames(brca.450.s)[[1]] <- as.character(brca.450.s[,1])
dimnames(lusc.450)[[1]] <- as.character(lusc.450[,1])
brca.normal <- (brca.450.s[,5:94])
kirc.normal <- (temp.N)
lusc.normal <- (lusc.450[,2:43])
brca <- (brca.normal[as.character(genenames.450[,"probeID"]),sample(1:(dim(brca.normal)[2]), 24, replace=F)])
kirc <- (kirc.normal[as.character(genenames.450[,"probeID"]),sample(1:(dim(kirc.normal)[2]), 24, replace=F)])
lusc <- (lusc.normal[as.character(genenames.450[,"probeID"]),sample(1:(dim(lusc.normal)[2]), 24, replace=F)])

dat.n.450 <- cbind(brca, kirc, lusc)
nor.per.450 <- apply(dat.n.450, 1, quantile, c(.1,.5,.9), na.rm=T)
nor.per.450 <- t(nor.per.450)
dimnames(nor.per.450)[[2]] <- c("N10","N50", "N90")


genenames.450 <- cbind(genenames.450, tum.per.450, nor.per.450)


### now, make conditions in the following manner

genenames.450$N.METH <- c("UNK")
genenames.450$T.METH <- c("UNK")
## for normals:


for(n in 1:dim(genenames.450)[1]){
if(is.na(genenames.450[n, "N10"] < 0.25 && genenames.450[n, "N50"] < 0.25 && genenames.450[n, "N90"] < 0.25)){
	genenames.450[n, "N.METH"] <- c("UNK")
} else{

	if(genenames.450[n, "N90"] < 0.25){

	genenames.450[n, "N.METH"] <- c("CUN")


} else {
	if(genenames.450[n, "N10"] > 0.75){

	genenames.450[n, "N.METH"] <- c("CMN")

} else {

if(genenames.450[n, "N10"] > 0.25 && genenames.450[n, "N90"] < 0.75){
	
	genenames.450[n, "N.METH"] <- c("IMN")
} else {
genenames.450[n, "N.METH"] <- c("VMN")

}


}


}
}
}#endfor

genenames.450[,"N.METH"] <- as.factor(genenames.450[,"N.METH"])



## for tumors
for(n in 1:dim(genenames.450)[1]){
if(is.na(genenames.450[n, "T10"] < 0.25 && genenames.450[n, "T50"] < 0.25 && genenames.450[n, "T90"] < 0.25)){
	genenames.450[n, "T.METH"] <- c("UNK")
} else{

	if(genenames.450[n, "T90"] < 0.25){

	genenames.450[n, "T.METH"] <- c("CUT")


} else {
	if(genenames.450[n, "T10"] > 0.75){

	genenames.450[n, "T.METH"] <- c("CMT")

} else {

if(genenames.450[n, "T10"] > 0.25 && genenames.450[n, "T90"] < 0.75){
	
	genenames.450[n, "T.METH"] <- c("IMT")
} else {
genenames.450[n, "T.METH"] <- c("VMT")

}


}


}
}
}#endfor

genenames.450[,"T.METH"] <- as.factor(genenames.450[,"T.METH"])



table(genenames.450[,"N.METH"],genenames.450[,"T.METH"], genenames.450[,"correlation_call"])
table(genenames.450[,"T.METH"],genenames.450[,"correlation_call"])
table(genenames.450[,"N.METH"],genenames.450[,"correlation_call"])
 table(genenames.450[,"N.METH"],genenames.450[,"T.METH"])







#### plots to make a threshold cutoffs
dimnames(y.drs.450)[[1]] <- y.drs.450[,2]
dat.t.450 <- (y.drs.450[as.character(genenames.450[,"probeID"]),])

dimnames(y.drs)[[1]] <- y.drs[,2]
dat.t <- (y.drs[as.character(genenames[,"probeID"]),])


t.27 <- dat.t[,2:281]
t.450 <- dat.t.450[,2:76]
#save(t.450,t.27, genenames,genenames.450, file="27k_450k_complete_data.rda")
load(file="27k_450k_complete_data.rda")

t.27 <- cbind(t.27,genenames[,"T.METH"])
t.450 <- cbind(t.450,genenames.450[,"T.METH"])



names(t.27)[281]<-"T.METH"
names(t.450)[76] <- "T.METH"
t.27<-t.27[,-1]
t.450<-t.450[,-1]

library(ggplot2)
t.27.m <- melt(t.27, id.vars = c("T.METH"))

t.450.m <- melt(t.450, id.vars = c("T.METH"))
t.27.m$type <- c("27K")
t.450.m$type <- c("450K")
dat <- rbind(t.27.m,t.450.m)


m <- ggplot(dat, aes(x=value, colour=T.METH)) + facet_grid(T.METH ~ type)
 x11();m + geom_density()


## Author: Houtan Noushmehr, PhD
## houtana@gmail.com
## 310.570.2362
## date: Jan 12, 2012



