## continued from part 1
## Author: Houtan Noushmehr, PhD
## houtana@gmail.com
## 310.570.2362
## date: Jan 12, 2012
## make patient centric calls for each cgProbe that matches the available gene expression by platform (27k, 450k)
## only need rules and dnamethylation data in probe x sample format, takes less than 1 min to run >350 samples, with over 200,000 probes.

## to run, just source this file. source(file="make_gene_calls_part2.R")


### load the following to start

load(file="27k_450k_complete_data.rda")
rules <- read.delim(file="Rules_Patient_Centric_Calls.txt",sep="\t",header=T)
rules[,"Label"] <- paste(rules[,"Correlation"], rules[,"Normal.Meth"], rules[,"Tumor.Meth"], sep=".")
### end load common data

#### 27 k analysis  #####
genenames$Label <- paste(genenames[,"correlation_call"], genenames[,"N.METH"], genenames[,"T.METH"], sep=".")
## because both rows are identical this can happen:: identical(dimnames(t.27)[[1]],dimnames(genenames)[[1]])
t.27m <- cbind(t.27,genenames[,c("genenames","Label")])
rules.all.27 <- merge(genenames[,c("probeID","Label")], rules[,c(1,5:10)],by.x="Label", by.y="Label",all.x=T)
write.table(rules.all.27, file="rules_27.txt",sep="\t",quote=F,row.names=F)
### in the terminal, ran this: sed s/"NA"/"UNK_in_Normal"/g rules_27.txt > rules_27_a.txt
rules.all.27 <- read.delim(file="rules_27_a.txt",sep="\t",header=T)
dimnames(rules.all.27)[[1]] <- rules.all.27[,"probeID"]
rules.all.27 <- rules.all.27[,-2]
rules.all.27 <- rules.all.27[dimnames(genenames)[[1]],] ## put rules.all.27 in the same order as genenames
patient.call.27k <- t.27m[,2:280]
patient.score.27k <- t.27m[,2:280]

tmp.lt.25.27 <- patient.call.27k < 0.25
tmp.gt.75.27 <- patient.call.27k > 0.75
tmp.ge.25.27 <- patient.call.27k >= 0.25 
tmp.le.75.27 <- patient.call.27k <= 0.75
tmp.27 <- tmp.ge.25.27 + tmp.le.75.27
tmp.ge.25.le.75.27 <- tmp.27 == 2

for(samp in 1:dim(patient.call.27k)[2]){

## convert NAs to FALSE
    tmp.lt.25.27[is.na(tmp.lt.25.27[,samp]),samp] <- FALSE
        tmp.gt.75.27[is.na(tmp.gt.75.27[,samp]),samp] <- FALSE
        tmp.ge.25.le.75.27[is.na(tmp.ge.25.le.75.27[,samp]),samp] <- FALSE


        patient.call.27k[tmp.lt.25.27[,samp],samp] <- as.character(rules.all.27[tmp.lt.25.27[,samp],"Call.TM.LT.0.25"])
        patient.call.27k[tmp.gt.75.27[,samp],samp] <- as.character(rules.all.27[tmp.gt.75.27[,samp],"Call.TM.GT.0.75"])
        patient.call.27k[tmp.ge.25.le.75.27[,samp],samp] <- as.character(rules.all.27[tmp.ge.25.le.75.27[,samp],"Call.0.25.GE.TM.LE.0.75"])

        patient.score.27k[tmp.lt.25.27[,samp],samp] <- as.character(rules.all.27[tmp.lt.25.27[,samp],"Confidence.Score.TM.LT.0.25"])
        patient.score.27k[tmp.gt.75.27[,samp],samp] <- as.character(rules.all.27[tmp.gt.75.27[,samp],"Confidence.Score.TM.GT.0.75"])
        patient.score.27k[tmp.ge.25.le.75.27[,samp],samp] <- as.character(rules.all.27[tmp.ge.25.le.75.27[,samp],"Confidence.Score.0.25.GE.TM.LE.0.75"])

## make each column a factor for easy summary and plot
        patient.call.27k[,samp] <- as.factor(patient.call.27k[,samp])
        patient.score.27k[,samp] <- as.factor(patient.score.27k[,samp])
        cat(date(), " :: Finished Sample : ",samp,"out of ", dim(patient.call.27k)[2],"\n")

}


### 450K analysis ####
genenames.450$Label <- paste(genenames.450[,"correlation_call"], genenames.450[,"N.METH"], genenames.450[,"T.METH"], sep=".")
## because both rows are identical this can happen:: identical(dimnames(t.450)[[1]],dimnames(genenames.450)[[1]])
t.450m <- cbind(t.450,genenames.450[,c("genenames.450","Label")])
rules.all.450 <- merge(genenames.450[,c("probeID","Label")], rules[,c(1,5:10)],by.x="Label", by.y="Label",all.x=T)
write.table(rules.all.450, file="rules_450.txt",sep="\t",quote=F,row.names=F)
### in the terminal, ran this: sed s/"NA"/"UNK_in_Normal"/g rules_450.txt > rules_450_a.txt
rules.all.450 <- read.delim(file="rules_450_a.txt",sep="\t",header=T)
dimnames(rules.all.450)[[1]] <- rules.all.450[,"probeID"]
rules.all.450 <- rules.all.450[,-2]
rules.all.450 <- rules.all.450[dimnames(genenames.450)[[1]],] ## put rules.all.450 in the same order as genenames
patient.call.450k <- t.450m[,2:75]
patient.score.450k <- t.450m[,2:75]

tmp.lt.25.450 <- patient.call.450k < 0.25
tmp.gt.75.450 <- patient.call.450k > 0.75
tmp.ge.25.450 <- patient.call.450k >= 0.25 
tmp.le.75.450 <- patient.call.450k <= 0.75
tmp.450 <- tmp.ge.25.450 + tmp.le.75.450
tmp.ge.25.le.75.450 <- tmp.450 == 2

for(samp in 1:dim(patient.call.450k)[2]){

## convert NAs to FALSE
    tmp.lt.25.450[is.na(tmp.lt.25.450[,samp]),samp] <- FALSE
        tmp.gt.75.450[is.na(tmp.gt.75.450[,samp]),samp] <- FALSE
        tmp.ge.25.le.75.450[is.na(tmp.ge.25.le.75.450[,samp]),samp] <- FALSE


        patient.call.450k[tmp.lt.25.450[,samp],samp] <- as.character(rules.all.450[tmp.lt.25.450[,samp],"Call.TM.LT.0.25"])
        patient.call.450k[tmp.gt.75.450[,samp],samp] <- as.character(rules.all.450[tmp.gt.75.450[,samp],"Call.TM.GT.0.75"])
        patient.call.450k[tmp.ge.25.le.75.450[,samp],samp] <- as.character(rules.all.450[tmp.ge.25.le.75.450[,samp],"Call.0.25.GE.TM.LE.0.75"])

        patient.score.450k[tmp.lt.25.450[,samp],samp] <- as.character(rules.all.450[tmp.lt.25.450[,samp],"Confidence.Score.TM.LT.0.25"])
        patient.score.450k[tmp.gt.75.450[,samp],samp] <- as.character(rules.all.450[tmp.gt.75.450[,samp],"Confidence.Score.TM.GT.0.75"])
        patient.score.450k[tmp.ge.25.le.75.450[,samp],samp] <- as.character(rules.all.450[tmp.ge.25.le.75.450[,samp],"Confidence.Score.0.25.GE.TM.LE.0.75"])

## make each column a factor for easy summary and plot
        patient.call.450k[,samp] <- as.factor(patient.call.450k[,samp])
        patient.score.450k[,samp] <- as.factor(patient.score.450k[,samp])
        cat(date(), " :: Finished Sample : ",samp,"out of ", dim(patient.call.450k)[2],"\n")
}

save(patient.call.27k, patient.score.27k, genenames, patient.call.450k, patient.score.450k, genenames.450, file = "Patient_Centric_Calls_Scores.rda")
## Author: Houtan Noushmehr, PhD
## houtana@gmail.com
## 310.570.2362
## date: Jan 12, 2012
