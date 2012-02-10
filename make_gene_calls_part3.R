## continued from Part 2
## Author: Houtan Noushmehr, PhD
## houtana@gmail.com
## 310.570.2362
## date: Jan 12, 2012
## need to summarize and plot data.

### load the following to start

load(file = "Patient_Centric_Calls_Scores.rda")


identical(dimnames(patient.call.450k)[[1]],dimnames(genenames.450)[[1]])
identical(dimnames(patient.call.27k)[[1]],dimnames(genenames)[[1]])

patient.call.450k <- cbind(patient.call.450k, genenames.450[,c("probeID","genenames.450")]); patient.call.450k$type <- c("call")
patient.score.450k <- cbind(patient.score.450k, genenames.450[,c("probeID","genenames.450")]); patient.score.450k$type <- c("score")
patient.call.27k <- cbind(patient.call.27k, genenames[,c("probeID","genenames")]); patient.call.27k$type <- c("call")
patient.score.27k <- cbind(patient.score.27k, genenames[,c("probeID","genenames")]); patient.score.27k$type <- c("score")

library(ggplot2)
dat.27 <- rbind(patient.call.27k, patient.score.27k); dat.27$platform <- c("27k")
dat.27s <- dat.27[,c(1:279,282:283)]
dat.27m <- melt(dat.27s, id.vars =c("type", "platform"))

dat.450 <- rbind(patient.call.450k, patient.score.450k); dat.450$platform <- c("450k")
dat.450s <- dat.450[,c(1:74,77:78)]
dat.450m <- melt(dat.450s, id.vars =c("type", "platform"))


dat <- rbind(dat.27m, dat.450m)
dat[,1] <- as.factor(dat[,1])
dat[,2] <- as.factor(dat[,2])

save(dat, dat.27, dat.27s, dat.27m, dat.450, dat.450s, dat.450m, file="data_contains_summary_data.rda")

###
###Summary begins
###

library(ggplot2)
load(file="data_contains_summary_data.rda")
## merge 27k, 450k and submit to AWG

dimnames(dat.450)[[2]] <- paste(dimnames(dat.450)[[2]], ".450K", sep="")
dimnames(dat.27)[[2]] <- paste(dimnames(dat.27)[[2]], ".27K", sep="")
dimnames(dat.450)[[2]][76] <- "genenames.450K"
c27 <- subset(dat.27, type.27K=="call")
c450 <- subset(dat.450, type.450K=="call")
c27.450 <- merge(c27, c450, by.x = "genenames.27K", by.y = "genenames.450K", all = T)
c27.450 <- c27.450[,c(2:280,284:357,1,281,358,282:283,359:360)]
c27.450[,354] <- as.factor(c27.450[,354])
c27.450[,355] <- as.factor(c27.450[,355])
c27.450[,356] <- as.factor(c27.450[,356])
c27.450[,357] <- as.factor(c27.450[,357])
c27.450[,358] <- as.factor(c27.450[,358])
c27.450[,359] <- as.factor(c27.450[,359])
c27.450[,360] <- as.factor(c27.450[,360])

s27 <- subset(dat.27, type.27K=="score")
s450 <- subset(dat.450, type.450K=="score")
s27.450 <- merge(s27, s450, by.x = "genenames.27K", by.y = "genenames.450K", all = T)
s27.450 <- s27.450[,c(2:280,284:357,1,281,358,282:283,359:360)]
s27.450[,354] <- as.factor(s27.450[,354])
s27.450[,355] <- as.factor(s27.450[,355])
s27.450[,356] <- as.factor(s27.450[,356])
s27.450[,357] <- as.factor(s27.450[,357])
s27.450[,358] <- as.factor(s27.450[,358])
s27.450[,359] <- as.factor(s27.450[,359])
s27.450[,360] <- as.factor(s27.450[,360])


tcga.gbm.dnameth.scores <- s27.450
tcga.gbm.dnameth.calls <- c27.450

save(tcga.gbm.dnameth.calls, tcga.gbm.dnameth.scores, file="TCGA_GBM_DNAMETHYLATION_CALLS_SCORES_20120112_Noushmehr.rda") ## added to TCGA GBM AWG Dropbox Folder '~/Dropbox/TCGA_GBM_MS_Writing_Group/Working_Group_SampleManifest/'
#####
#####
#####
### summarize and identify genes
###

library(ggplot2)
load(file="data_contains_summary_data.rda")
genedata.27 <- dat.27[,c("probeID","genenames","type","platform")]; genedata.27c <- subset(genedata.27, type=="call"); genedata.27s <- subset(genedata.27, type=="score")
genedata.450 <- dat.450[,c("probeID","genenames.450","type","platform")]; genedata.450c <- subset(genedata.450, type=="call"); genedata.450s <- subset(genedata.450, type=="score")

tmp.27.c <- subset(dat.27, type=="call")
tmp.27.s <- subset(dat.27, type=="score")
tmp.450.c <- subset(dat.450, type=="call")
tmp.450.s <- subset(dat.450, type=="score")

ut.call.450 <- tmp.450.c == "UT"; sumoUT.450 <- rowSums(ut.call.450, na.rm=T);genedata.450c$UT <- sumoUT.450
ml.call.450 <- tmp.450.c == "ML"; sumoML.450 <- rowSums(ml.call.450, na.rm=T);genedata.450c$ML <- sumoML.450
mg.call.450 <- tmp.450.c == "MG"; sumoMG.450 <- rowSums(mg.call.450, na.rm=T);genedata.450c$MG <- sumoMG.450
es.call.450 <- tmp.450.c == "ES"; sumoES.450 <- rowSums(es.call.450, na.rm=T);genedata.450c$ES <- sumoES.450
mt.call.450 <- tmp.450.c == "MT"; sumoMT.450 <- rowSums(mt.call.450, na.rm=T);genedata.450c$MT <- sumoMT.450
uc.call.450 <- tmp.450.c == "UC"; sumoUC.450 <- rowSums(uc.call.450, na.rm=T);genedata.450c$UC <- sumoUC.450

ut.call.27  <- tmp.27.c == "UT"; sumoUT.27 <- rowSums(ut.call.27, na.rm=T);genedata.27c$UT <- sumoUT.27
ml.call.27  <- tmp.27.c == "ML"; sumoML.27 <- rowSums(ml.call.27, na.rm=T);genedata.27c$ML <- sumoML.27
mg.call.27  <- tmp.27.c == "MG"; sumoMG.27 <- rowSums(mg.call.27, na.rm=T);genedata.27c$MG <- sumoMG.27
es.call.27  <- tmp.27.c == "ES"; sumoES.27 <- rowSums(es.call.27, na.rm=T);genedata.27c$ES <- sumoES.27
mt.call.27  <- tmp.27.c == "MT"; sumoMT.27 <- rowSums(mt.call.27, na.rm=T);genedata.27c$MT <- sumoMT.27
uc.call.27  <- tmp.27.c == "UC"; sumoUC.27 <- rowSums(uc.call.27, na.rm=T);genedata.27c$UC <- sumoUC.27

ut.score.450 <- tmp.450.s == "0"; sumo0.450 <- rowSums(ut.score.450, na.rm=T);genedata.450s$s0 <- sumo0.450
ml.score.450 <- tmp.450.s == "1"; sumo1.450 <- rowSums(ml.score.450, na.rm=T);genedata.450s$s1 <- sumo1.450
mg.score.450 <- tmp.450.s == "2"; sumo2.450 <- rowSums(mg.score.450, na.rm=T);genedata.450s$s2 <- sumo2.450
es.score.450 <- tmp.450.s == "3"; sumo3.450 <- rowSums(es.score.450, na.rm=T);genedata.450s$s3 <- sumo3.450
mt.score.450 <- tmp.450.s == "4"; sumo4.450 <- rowSums(mt.score.450, na.rm=T);genedata.450s$s4 <- sumo4.450

ut.score.27  <- tmp.27.s == "0"; sumo0.27 <- rowSums(ut.score.27, na.rm=T);genedata.27s$s0 <- sumo0.27
ml.score.27  <- tmp.27.s == "1"; sumo1.27 <- rowSums(ml.score.27, na.rm=T);genedata.27s$s1 <- sumo1.27
mg.score.27  <- tmp.27.s == "2"; sumo2.27 <- rowSums(mg.score.27, na.rm=T);genedata.27s$s2 <- sumo2.27
es.score.27  <- tmp.27.s == "3"; sumo3.27 <- rowSums(es.score.27, na.rm=T);genedata.27s$s3 <- sumo3.27
mt.score.27  <- tmp.27.s == "4"; sumo4.27 <- rowSums(mt.score.27, na.rm=T);genedata.27s$s4 <- sumo4.27

genedata.27 <- cbind(genedata.27c, genedata.27s[,5:9])
genedata.450 <- cbind(genedata.450c, genedata.450s[,5:9])
names(genedata.450) <- names(genedata.27)
genedata <- rbind(genedata.27, genedata.450); genedata[,"type"] <- as.factor(genedata[,"type"]); genedata[,"platform"] <- as.factor(genedata[,"platform"]); genedata <- genedata[,-3]
genedata.m <- melt(genedata, measure.vars = c("UT", "ML", "MG", "ES", "MT", "UC", "s0", "s1", "s2", "s3", "s4"))



tmp.27 <- dat.27 == "ES"
sumofES.27 <- rowSums(tmp.27, na.rm=T)
dat.27$SumES <- sumofES.27

tmp.450 <- dat.450 == "ES"
sumofES.450 <- rowSums(tmp.450, na.rm=T)
dat.450$SumES <- sumofES.450

length(subset(dat.450, SumES==range(sumofES.450)[2])[,"genenames.450"])
length(subset(dat.27, SumES==range(sumofES.27)[2])[,"genenames"])

length(intersect((subset(dat.450, SumES>c(.5*dim(dat.450)[2]))[,"genenames.450"]), (subset(dat.27, SumES>c(.5*dim(dat.27)[2]))[,"genenames"])))

dat.call <- subset(dat, type=="call"); dimnames(dat.call)[[2]] <- paste(dimnames(dat.call)[[2]], ".call", sep="")
dat.score <- subset(dat, type=="score"); dimnames(dat.score)[[2]] <- paste(dimnames(dat.score)[[2]], ".score", sep="")
dat.t <- cbind(dat.call, dat.score)



library(ggplot2)
theme_white <- function() {

 theme_update (
 plot.background = theme_blank(),
 panel.background=theme_rect(colour="black", size=1),
 axis.text.x= theme_text(colour="black",vjust= 1, size=12),
 axis.text.y= theme_text(colour="black",hjust=1, size=12),
 axis.title.x =theme_text(colour="black",face="bold", size=12),
 axis.title.y =theme_text(colour="black",face="bold", angle = 90, size=12)
 )
}
theme_white()

###  log10 scale
m <- ggplot(dat, aes(x=value, fill=type))
m + geom_bar() + 
	opts(axis.text.x = theme_text(angle=90), title = "Calls & Score\n 27K: 279 samp X 9453 (x2) =  5,274,774\n 450K: 74 samp X 11040 (x2) = 1,633,920") + 
	#scale_y_continuous("Raw Counts", formatter = "comma") + 
	scale_y_log10("Raw Counts - log10") + 
	scale_x_discrete("Call or Score", limits=c("UT","ML","MG","ES","MT","UC","0","1","2","3","4")) +
	facet_grid(type ~ platform)
#ggsave(file="summary_plot_1.pdf") ## raw counts with comma separator
ggsave(file="summary_plot_2.pdf") ## raw counts on a log10 scale


## raw scale with comma separator
m <- ggplot(dat, aes(x=value, fill=type))
m + geom_bar() + 
	opts(axis.text.x = theme_text(angle=90), title = "Calls & Score\n 27K: 279 samp X 9453 (x2) =  5,274,774\n 450K: 74 samp X 11040 (x2) = 1,633,920") + 
	scale_y_continuous("Raw Counts", formatter = "comma") + 
	#scale_y_log10("Raw Counts - log10") + 
	scale_x_discrete("Call or Score", limits=c("UT","ML","MG","ES","MT","UC","0","1","2","3","4")) +
	facet_grid(type ~ platform)
ggsave(file="summary_plot_1.pdf") ## raw counts with comma separator
#ggsave(file="summary_plot_2.pdf") ## raw counts on a log10 scale


 dat.t0 <- subset(dat.t, value.score=="0"); dat.t0$value.score <- c("Se.0")
 dat.t1 <- subset(dat.t, value.score=="1"); dat.t1$value.score <- c("Sd.1")
 dat.t2 <- subset(dat.t, value.score=="2"); dat.t2$value.score <- c("Sc.2")
 dat.t3 <- subset(dat.t, value.score=="3"); dat.t3$value.score <- c("Sb.3")
 dat.t4 <- subset(dat.t, value.score=="4"); dat.t4$value.score <- c("Sa.4")
dat.t.new <- rbind(dat.t0,dat.t1,dat.t2,dat.t3,dat.t4)

## stacked bar chart # raw count
m <- ggplot(dat.t.new, aes(x=value.call, fill=value.score))
m + geom_bar() + 
	opts(axis.text.x = theme_text(angle=90), title = "Calls & Score\n 27K: 279 samp X 9453 (x2) =  5,274,774\n 450K: 74 samp X 11040 (x2) = 1,633,920") + 
	scale_y_continuous("Raw Counts", formatter = "comma") + 
	#scale_y_log10("Raw Counts - log10") + 
	scale_x_discrete("Call", limits=c("UT","ML","MG","ES","MT","UC")) +
	scale_fill_manual(values=c("Sa.4" = "Black", "Sb.3" = "gray25","Sc.2" = "#808080", "Sd.1" = "gray75","Se.0" = "#d3d3d3"), breaks = c("Se.0", "Sd.1", "Sc.2", "Sb.3", "Sa.4"), labels = c("0", "1", "2", "3", "4")) + 
	facet_grid(~ platform.call)
ggsave(file="summary_plot_3.pdf") ## raw counts with comma separator
#ggsave(file="summary_plot_4.pdf") ## raw counts on a log10 scale



## stacked bar chart # raw count
m <- ggplot(dat.t.new, aes(x=value.call, fill=value.score))
m + geom_bar(position="fill") + 
	opts(axis.text.x = theme_text(angle=90), title = "Calls & Score\n 27K: 279 samp X 9453 (x2) =  5,274,774\n 450K: 74 samp X 11040 (x2) = 1,633,920") + 
	scale_y_continuous("Percent Total") + 
	#scale_y_log10("Raw Counts - log10") + 
	scale_x_discrete("Call", limits=c("UT","ML","MG","ES","MT","UC")) +
	scale_fill_manual(values=c("Sa.4" = "Black", "Sb.3" = "gray25","Sc.2" = "#808080", "Sd.1" = "gray75","Se.0" = "#d3d3d3"), breaks = c("Se.0", "Sd.1", "Sc.2", "Sb.3", "Sa.4"), labels = c("0", "1", "2", "3", "4")) + 
	facet_grid(~ platform.call)
#ggsave(file="summary_plot_3.pdf") ## raw counts with comma separator
ggsave(file="summary_plot_4.pdf") ## percent total




## stacked bar chart # raw count
m <- ggplot(genedata.m, aes(x=factor(variable),value))
m + geom_boxplot() + 
	opts(axis.text.x = theme_text(angle=90), title = "Calls & Score\n 27K: 279 samp X 9453 (x2) =  5,274,774\n 450K: 74 samp X 11040 (x2) = 1,633,920") + 
	scale_y_continuous("Raw Counts", formatter = "comma") + 
	#scale_y_log10("Raw Counts - log10") + 
	scale_x_discrete("Call", limits=c("UT","ML","MG","ES","MT","UC")) +
	scale_fill_manual(values=c("4" = "Black", "3" = "gray25","2" = "#808080", "1" = "gray75","0" = "#d3d3d3"), breaks = c("0", "1", "2", "3", "4")) + 
	facet_grid(~ platform.call)
ggsave(file="summary_plot_3.pdf") ## raw counts with comma separator
#ggsave(file="summary_plot_4.pdf") ## raw counts on a log10 scale

## Author: Houtan Noushmehr, PhD
## houtana@gmail.com
## 310.570.2362
## date: Jan 12, 2012



