##Manual checks of correlations; messy script
#refer to more systematic analysis in TAD_BW_analysis.R and distinctions.R

library(tidyverse)
library(reshape2)
library(GenomicRanges)
library(ggsignif)
library(ggpubr)
library(ggplot2)
library(plotly)
library(data.table)

loopdepthCat <- function(x){
  if(x<2){
    return("topLoop")
  }
  if(x<10){
    return("middleLoop")
  }else{
    return("nestedLoop")
  }}

activity <- function(x){
  if(x==1){
    return("active")
  }else{
    return("non-active")
  }}

ctcfCat <- function(x){
  #hist(FullStatObj$vp.Hap1_CTCF_chip,breaks=seq(0,1e5,500), xlim=c(0,2e4))
  if(is.na(x)){
    return(NA)
  }
  if(x<11.695){
    return("lowCTCF")
  }
  if(x< 28.905){
    return("mediumCTCF")
    }
  else{
    return("highCTCF")
  }}

df2GR <- function(statObj4){
  
  df <- data.frame(chr=statObj4$chr, vp=statObj4$vp_X1, anch = statObj4$maxV4CscorePos)
  df[df$vp > df$anch,] <- df[df$vp > df$anch, c(1,3,2)]
  
  ranges <-  GRanges (seqnames=df$chr, IRanges(start = df$vp, end = df$anch) )
  
  return(ranges)
}

setwd('/delaat/group/iwan/peakHiC/scripts/')
FullStatObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/FullStatObj.rds") #statObj/avg_chip/FullStatObj.rds

head(FullStatObj)

#alternative ChIP
# newChips <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/NewChip.rds")
# FullStatObj[, names(newChips[,48:55])] <- newChips[,48:55]

#####add TPM#####
TPM <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/Expression_levels/normHTseq/TPMcounts.rds")
TPM <- cbind(rownames(TPM), data.frame(TPM, row.names=NULL))
TPM[,'WaplKO_1.14'] <- NULL

names(TPM) <- c('ensembl', "vp.TPM.Hap1", "vp.TPM.WaplKO","vp.TPM.SCC4KO","vp.TPM.DKO")
FullStatObj <- left_join(FullStatObj, TPM, by=c("vp.Ens"='ensembl'))

names(TPM) <- c('ensembl', "anch.TPM.Hap1", "anch.TPM.WaplKO","anch.TPM.SCC4KO","anch.TPM.DKO")
FullStatObj <- left_join(FullStatObj, TPM, by=c("anch.Ens"='ensembl'))

FullStatObj$Hap1_TPM <- as.numeric(dplyr::coalesce( FullStatObj$vp.TPM.Hap1,  FullStatObj$anch.TPM.Hap1))
FullStatObj$WaplKO_TPM <- as.numeric(dplyr::coalesce( FullStatObj$vp.TPM.WaplKO,  FullStatObj$anch.TPM.WaplKO))
FullStatObj$SCC4KO_TPM <- as.numeric(dplyr::coalesce( FullStatObj$vp.TPM.SCC4KO,  FullStatObj$anch.TPM.SCC4KO))
FullStatObj$DKO_TPM <- as.numeric(dplyr::coalesce( FullStatObj$vp.TPM.DKO,  FullStatObj$anch.TPM.DKO))
#saveRDS(FullStatObj, "/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/FullStatObj_TPM.rds")

#####filter out enhancers hubs / prom hubs
counts <- table(FullStatObj$id)
counts <- table(FullStatObj$anchor.vpID)
names <- names(counts[counts>7])
manyLoops <- FullStatObj[FullStatObj$id %in% names,]

barplot(table(counts), xlab = 'n.o. contacts', ylab = 'Frequency')
sum(counts[counts>1])


#########loopdistances per condition boxplot#################
conds = c("DKO", "Hap1", "WaplKO_3.3", "SCC4KO", "WaplKO_1.14") 
measurevars <- colnames(FullStatObj)[ (colnames(FullStatObj) %in% paste0("lfc",conds)) ]

melted.stat <- melt(FullStatObj, id.vars = "loopDistCat", measure.vars = measurevars )

boxplot_distcat <- ggplot(data = melted.stat, aes(x=loopDistCat, y=value, fill=loopDistCat)) +
  geom_boxplot() +
  facet_grid(~variable) +
  theme_pubr(border = T)+
  fill_palette(palette = "Reds", direction=-1)+
  coord_cartesian(ylim=c(-2,2))+
  title('Looplength versus Read Coverage')+
  xlab('Looplength Category')+
  ylab('log2(WAPL/WT) Read Coverage')+
  theme(axis.text.x = element_text(angle=30,size = 16, vjust=0.8), axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20), legend.position = "none",
        strip.text.x = element_text(size = 12))
#geom_signif(aes(group=variable), comparisons = list(c("short", "long")), test='t.test') 
boxplot_distcat
#ggsave(filename="/delaat/group/iwan/Hap1_peakHiC/plots/boxplot_distcat.pdf", boxplot_distcat)

dev.off()


#####Take only TSS genes (expr differences comparisons)####
#statObj.tss2 <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/statObj_GeneExpr.rds")
#statObj.tss <- FullStatObj[FullStatObj$anch.TSS | FullStatObj$vp.TSS >= 1,]
statObj.tss <- FullStatObj[ !is.na(FullStatObj$vp.Hap1.RNAcount)  | !is.na(FullStatObj$anch.Hap1.RNAcount), ]
statObj.tss[,grepl('DE',names(statObj.tss))] <- apply(statObj.tss[,grepl('DE',names(statObj.tss))], 2, as.numeric)

######categories from continuous CTCF chip versus DE#######
statObj.tss$Wapl_CTCF.cat <- factor( sapply( apply(statObj.tss[,c("vp.Wapl3_3_CTCF_chip", "anch.Wapl3_3_CTCF_chip")], 1, max), ctcfCat), levels=c("lowCTCF", "mediumCTCF", "highCTCF") )
statObj.tss$Hap_CTCF.cat <- factor( sapply( apply(statObj.tss[,c("vp.Hap1_CTCF_chip","anch.Hap1_CTCF_chip" )], 1, max), ctcfCat), levels=c("lowCTCF", "mediumCTCF", "highCTCF") )
#statObj.tss$Wapl_CTCF.cat <- factor( sapply( apply(statObj.tss[,63:64], 1, median), ctcfCat), levels=c("lowCTCF", "mediumCTCF", "highCTCF") )

###combine anch and vp DE
statObj.tss$WaplKO_DE <- as.numeric(dplyr::coalesce( statObj.tss$vp.DE.WaplKO_vs_Hap1,  statObj.tss$anch.DE.WaplKO_vs_Hap1))
statObj.tss$DKO_DE <- as.numeric(dplyr::coalesce( statObj.tss$vp.DE.DKO_vs_Hap1,  statObj.tss$anch.DE.DKO_vs_Hap1))
statObj.tss$SCC4KO_DE <- as.numeric(dplyr::coalesce( statObj.tss$vp.DE.SCC4_vs_Hap1,  statObj.tss$anch.DE.SCC4_vs_Hap1))



########## H3K27ac at promotor versus TSS expression levels ###### 
#sort of both but much stronger with narrowpeak count; fix chip
statObj.tss$enhanced.promotor <- 'Not enhanced'
#statObj.tss$enhanced.promotor[ (statObj.tss$vp.Hap1_H3k27ac_chip >=7 & statObj.tss$vp.TSS>=1 ) | 
#                                (statObj.tss$anch.Hap1_H3k27ac_chip >=7 & statObj.tss$anch.TSS>=1 ) ] <- 'enhanced'
statObj.tss[ (statObj.tss$vp.H3K27ac >=1 & statObj.tss$vp.TSS>=1 ) | 
               (statObj.tss$anch.H3K27ac >=1 & statObj.tss$anch.TSS>=1 ), 'enhanced.promotor'] <- 'enhanced'

t.test(statObj.tss[statObj.tss$enhanced.promotor == 'enhanced', 'avg.Hap1.TPM'], statObj.tss[statObj.tss$enhanced.promotor == 'Not enhanced', 'avg.Hap1.TPM'])
#median(statObj.tss[statObj.tss$enhanced.promotor == 'enhanced', 'Hap1_TPM'])

conds = c("DKO", "Hap1", "WaplKO_3.3", "SCC4KO", "WaplKO_1.14") 
measurevars <- colnames(statObj.tss)[ (colnames(statObj.tss) %in% paste0(conds,"_TPM")) ] #.RNAcount
melted.stat <- melt(statObj.tss, id.vars = "enhanced.promotor", measure.vars = measurevars )

boxplot_enhancedGene <- ggplot(data = statObj.tss, aes(x=enhanced.promotor, y=avg.Hap1.TPM)) +
  geom_boxplot() +
  theme_pubr(border = T) +
  coord_cartesian(ylim=c(0,100)) +
  xlab('Promoter Status based on pooled counts') +
  ylab("TPM")+
  theme(axis.text.x = element_text(size = 18), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))
  
boxplot_enhancedGene


###############active.genes and gene expression correlation##########
#artifact that active.gene wasnt nicely defined
statObj.tss$active.gene <- 0
statObj.tss[statObj.tss$anch.TSS >= 1 & statObj.tss$anch.RNAPol >=1, ]$active.gene <- 1
statObj.tss[statObj.tss$vp.TSS >= 1 & statObj.tss$vp.RNAPol >=1, ]$active.gene <- 1

statObj.tss$active.gene <- factor( sapply(statObj.tss$active.gene, activity), levels=c("active", "non-active") )

conds = c("DKO", "Hap1", "WaplKO_3.3", "SCC4KO", "WaplKO_1.14") 
for(cond in conds){
  statObj.tss[,paste0(cond,".RNAcount")] <- dplyr::coalesce(statObj.tss[,paste0('vp.',cond,".RNAcount")], statObj.tss[,paste0('anch.',cond,".RNAcount")])
}

#saveRDS(statObj.tss,'/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/statObjTSS.rds')
statObj.tss <- readRDS('/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/statObjTSS.rds')
measurevars <- colnames(statObj.tss)[ (colnames(statObj.tss) %in% paste0(conds,"_TPM")) ] #.RNAcount

melted.stat <- melt(statObj.tss, id.vars = "active.gene", measure.vars = measurevars )

boxplot_activeGene <- ggplot(data = melted.stat, aes(x=active.gene, y=value)) +
  geom_boxplot() +
  facet_grid(~variable) +
  theme_pubr(border = T) +
  coord_cartesian(ylim=c(0,100)) +#2e3
  xlab('RNApolII colocalization with promoter')+
  ylab('TPM')+
  theme(axis.text.x = element_text(size = 18), axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))
#geom_signif(aes(group=variable), comparisons = list(c("short", "long")), test='t.test') 
boxplot_activeGene
#ggsave(filename="/delaat/group/iwan/Hap1_peakHiC/plots/boxplot_activeGene.pdf", boxplot_activeGene)

dev.off()


#####new correlations?#####
#####CTCF as boundary to limit transcription#### 
#see outliers
plot(statObj.tss$vp.TPM.Hap1, statObj.tss$vp.CTCF_WT, type='h') #vp.Hap1_CTCF_chip


###CTCF vs DE TSS genes; no difference####
low <- statObj.tss$WaplKO_DE[statObj.tss$Wapl_CTCF.cat == "lowCTCF"]
medium <- statObj.tss$WaplKO_DE[statObj.tss$Wapl_CTCF.cat ==  "mediumCTCF"]
high <-statObj.tss$WaplKO_DE[statObj.tss$Wapl_CTCF.cat == "highCTCF"]
wilcox.test(low,medium)
wilcox.test(medium,high)
wilcox.test(low,high)
boxplot(Hap1CTCF~WaplKO_DE, data=statObj.tss, main="DE vs CTCF", xlab="CTCF category", ylab="Differential Expression lfc") 
#saveRDS(statObj.tss, "/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/statObjTSS2.rds")

#check for CTCF dependence in DE
statObj.tss$Hap1CTCF <- apply(statObj.tss[,c("vp.Hap1_CTCF_chip", "anch.Hap1_CTCF_chip")], 1, max)
statObj.tss$is.DE.wapl <- 0; statObj.tss$is.DE.wapl[!is.na(statObj.tss$WaplKO_DE)] <- 1
statObj.tss$is.DE.scc4 <- 0; statObj.tss$is.DE.scc4[!is.na(statObj.tss$SCC4KO_DE)] <- 1
statObj.tss$is.DE.dko <- 0; statObj.tss$is.DE.dko[!is.na(statObj.tss$DKO_DE)] <- 1

#wilcox.test(statObj.tss$Hap1CTCF[is.na(statObj.tss$WaplKO_DE)], statObj.tss$Hap1CTCF[!is.na(statObj.tss$WaplKO_DE)])
boxplot(Hap1CTCF~is.DE.wapl, data=statObj.tss, main="CTCF vs DE", xlab="DE no (0) or yes (1)", ylab="CTCF")


##########check for H3K27ac dependence#####
statObj.tss$avg.HAPctcf <- apply(statObj.tss[, c("vp.Hap1_CTCF_chip","anch.Hap1_CTCF_chip" )], 1, max)
statObj.tss$avg.WAPLctcf <- apply(statObj.tss[, c("vp.Wapl3_3_CTCF_chip" ,"anch.Wapl3_3_CTCF_chip" )], 1, max)

plot(statObj.tss$vp.Hap1_CTCF_chip, statObj.tss$vp.TPM.Hap1, type='h', main="Expression vs CTCF", xlab="Hap1 CTCF")
#plot(statObj.tss$avg.WAPLctcf, statObj.tss$WaplKO_TPM, type='h', main="Expression vs CTCF", xlab="Hap1 CTCF")
plot(log(statObj.tss$vp.Wapl3_3_CTCF_chip/statObj.tss$vp.Hap1_CTCF_chip), statObj.tss$WaplKO_DE, type='h', main="Expression vs CTCF",  xlab="lfc CTCF", ylab='WaplKO DE')
plot(log(statObj.tss$vp.Wapl3_3_SMC1_chip/statObj.tss$vp.Hap1_SMC1_chip), statObj.tss$WaplKO_DE, type='h', main="Expression vs CTCF",  xlab="lfc CTCF", ylab='WaplKO DE')
#ctcf chip is less and expression is more (left upper corner)
#for the above x=(statObj.tss$avg.HAPctcf-statObj.tss$avg.WAPLctcf)   xlab='delta CTCF'

#search for the three groups and their differences
#first in the DE>0 and lfcCTCT>-0.15
upDE_downCTCF <- statObj.tss[statObj.tss$WaplKO_DE>0 & log(statObj.tss$vp.Wapl3_3_CTCF_chip/statObj.tss$vp.Hap1_CTCF_chip)<.15, ]
upDE_downCTCF <-upDE_downCTCF[!is.na(upDE_downCTCF$loopID),]

upDE_upCTCF <- statObj.tss[ statObj.tss$WaplKO_DE>0 & log(statObj.tss$vp.Wapl3_3_CTCF_chip/statObj.tss$vp.Hap1_CTCF_chip)>.15 , ]
upDE_upCTCF <-upDE_upCTCF[!is.na(upDE_upCTCF$loopID),]

downDE_normalCTCF <- statObj.tss[ statObj.tss$WaplKO_DE &  statObj.tss$WaplKO_DE<0, ]
downDE_normalCTCF <-downDE_normalCTCF[!is.na(downDE_normalCTCF$loopID),]

reference <- statObj.tss[!is.na( statObj.tss$WaplKO_DE), ]



avgCTCFinGene <- function(df){
  df$CTCFinGene <- 0
  df[ df$vp.CTCFinGene==1 | df$anch.CTCFinGene==1 , ]$CTCFinGene <- 1
  return(df)}

reference1 <- setdiff(reference, upDE_downCTCF)
upDE_downCTCF <- avgCTCFinGene(upDE_downCTCF)
reference1 <- avgCTCFinGene(reference1)
wilcox.test(upDE_downCTCF$CTCFinGene, reference1$CTCFinGene)
t.test(upDE_downCTCF$vp.H3K56ac, reference1$vp.H3K4me3)
t.test(upDE_downCTCF$loopDist, reference1$loopDist)
t.test(is.finite(upDE_downCTCF$lfcWaplKO_3.3), is.finite(reference1$lfcWaplKO_3.3))


reference2 <- setdiff(reference, upDE_upCTCF)
t.test(upDE_upCTCF$prom.enh, reference2$prom.enh) #a tiny bit but barely
t.test(upDE_upCTCF$vp.H3K27ac, reference2$vp.H3K27ac)
t.test(is.finite(upDE_downCTCF$lfcWaplKO_3.3), is.finite(reference2$lfcWaplKO_3.3))
t.test(upDE_upCTCF$loopDist, reference2$loopDist)

upDE_upCTCF[order(-upDE_upCTCF$loopDist),]

#ENSG00000127152 
statObj.tss[which(statObj.tss$anch.Ens == 'ENSG00000127152'),]


reference3 <- setdiff(reference, downDE_normalCTCF) 
t.test(reference3$prom.enh, downDE_normalCTCF$prom.enh) #tiny bit; doesnt seem to work this analysis
t.test(upDE_downCTCF$vp.H3K27ac, reference1$vp.H3K27ac)
t.test(is.finite(upDE_downCTCF$lfcWaplKO_3.3), is.finite(reference$lfcWaplKO_3.3))


statObj.tss$avg.H3K27ac <- apply(statObj.tss[, c('vp.Hap1_H3k27ac_chip', 'anch.Hap1_H3k27ac_chip')], 1,max)
plot(statObj.tss$avg.HAPctcf, statObj.tss$avg.H3K27ac, type='h', main="H3K27ac vs CTCF", xlab="Hap1 CTCF")

########compare DE between prom.enh loops and non-prom.enh loops############# ; failw with continuous
hist(statObj.tss$avg.H3K27ac, breaks=150, xlab = 'average H3K27ac max score')
hist(statObj.tss$vp.Hap1_H3k27ac_chip, breaks=150, xlab = 'vp H3K27ac max score')
hist(statObj.tss$vp.H3K27ac, xlab = 'vp H3K27ac count')

#For max chip
enh.statObj <- statObj.tss[(statObj.tss$vp.TSS >= 1 & statObj.tss$anch.Hap1_H3k27ac_chip >= 7)|
                             (statObj.tss$anch.TSS >= 1 & statObj.tss$vp.Hap1_H3k27ac_chip >= 7), ]

#For summed chip
enh.statObj <-  statObj.tss[(statObj.tss$vp.TSS >= 1 & statObj.tss$anch.Hap1_H3k27ac_chip >= 1.5e4) |
                                (statObj.tss$anch.TSS >= 1 & statObj.tss$vp.Hap1_H3k27ac_chip >= 1.5e4), ]

#For narrowpeak count
enh.statObj <- statObj.tss[statObj.tss$vp.TSS >= 1 & statObj.tss$anch.H3K27ac >= 1 | 
                             statObj.tss$anch.TSS >= 1 & statObj.tss$vp.H3K27ac >= 1, ]

non.enh.statObj <- setdiff(statObj.tss, enh.statObj)

########compare DE between prom.enh loops and non-prom.enh loops#############  #less good at discerning H3K27 levels; great at discerning DE
###test if prom.enh looped genes do really have more H3K27ac
Hap1.H3K27Expr.enh <- apply( data.frame(enh.statObj$vp.Hap1_H3k27ac_chip, enh.statObj$anch.Hap1_H3k27ac_chip), 1, max) 
Hap1.H3K27Expr.non.enh <- apply( data.frame(non.enh.statObj$vp.Hap1_H3k27ac_chip, non.enh.statObj$anch.Hap1_H3k27ac_chip), 1, max) 
t.test(Hap1.H3K27Expr.enh, Hap1.H3K27Expr.non.enh)

#test expression difference
t.test(enh.statObj$Hap1_TPM, non.enh.statObj$Hap1_TPM)
t.test(enh.statObj$WaplKO_DE[enh.statObj$WaplKO_DE<0], non.enh.statObj$WaplKO_DE[non.enh.statObj$WaplKO_DE<0])
median(enh.statObj$WaplKO_DE[enh.statObj$WaplKO_DE<0], na.rm=T); median(non.enh.statObj$WaplKO_DE[non.enh.statObj$WaplKO_DE<0], na.rm=T)
t.test(enh.statObj[enh.statObj$WaplKO_DE>0,]$loopdepth, non.enh.statObj[non.enh.statObj$WaplKO_DE>0,]$loopdepth)


#enh object positive DE have less loops changed than non.enh obj positive DE, and, less CTCF (from narrowpeak)
t.test(is.finite(enh.statObj[enh.statObj$WaplKO_DE>0,]$lfcWaplKO_3.3), is.finite(non.enh.statObj[non.enh.statObj$WaplKO_DE>0,]$lfcWaplKO_3.3)) #avg.HAPctcf
#less ctcf not the case for all general loops (not DE distinguished) loops (from narrowpeak); (yes from bw)
t.test(enh.statObj$avg.HAPctcf, non.enh.statObj$avg.HAPctcf)


#make NA zero; 
for( name in grep('DE', names(statObj.tss), value=T) ){
  enh.statObj[,name][is.na(enh.statObj[,name])] <- 0
  non.enh.statObj[,name][is.na(non.enh.statObj[,name])] <- 0
}

#take absolute values
enh.statObj.abs <- enh.statObj
non.enh.statObj.abs <- non.enh.statObj
for( name in grep('DE', names(statObj.tss), value=T) ){
  enh.statObj.abs[,name] <- abs(as.numeric(enh.statObj[,name]))
  non.enh.statObj.abs[,name] <- abs(as.numeric(non.enh.statObj[,name]))
}

#Not if we use statObj.tss$vp.H3K27ac >= 2 then other way around or max chip; yes if we use sum
t.test(enh.statObj.abs$vp.DE.WaplKO_vs_Hap1, non.enh.statObj.abs$vp.DE.WaplKO_vs_Hap1)
t.test(enh.statObj.abs$vp.DE.SCC4_vs_Hap1, non.enh.statObj.abs$vp.DE.SCC4_vs_Hap1)
t.test(enh.statObj.abs$vp.DE.DKO_vs_Hap1, non.enh.statObj.abs$vp.DE.DKO_vs_Hap1)

#stronger; 
t.test(enh.statObj.abs$WaplKO_DE, non.enh.statObj.abs$WaplKO_DE)


#WAPLKO############
#split in up and down regulated: of the DE loops; 
enh.up.wapl <- enh.statObj[enh.statObj$WaplKO_DE>0,]$WaplKO_DE; non.enh.up.wapl <- non.enh.statObj[non.enh.statObj$WaplKO_DE>0,]$WaplKO_DE
#non prom-enh loops are more upregulated than prom-enh loops
t.test(enh.up.wapl,  non.enh.up.wapl)
hist(enh.up.wapl, breaks = 30) #both roughly 400 loops; some more extreme (>4) values in non.enh.up
hist(non.enh.up.wapl, breaks=30)

WaplGenes <- non.enh.statObj[non.enh.statObj$WaplKO_DE>4,] #all prom.ctcf; all non-prom.polII. Genes: PDE3A; ANO3; NCAM2; GDNF; TLL1

enh.down.wapl <- enh.statObj[enh.statObj$WaplKO_DE<0,]$WaplKO_DE; non.enh.down.wapl <- non.enh.statObj[non.enh.statObj$WaplKO_DE<0,]$WaplKO_DE
t.test( enh.down.wapl, non.enh.down.wapl )

#SCC4###################
#split in up and down regulated: of the DE loops; 
enh.up.SCC4 <- enh.statObj[enh.statObj$SCC4KO_DE>0,]$SCC4KO_DE; non.enh.up.SCC4 <- non.enh.statObj[non.enh.statObj$SCC4KO_DE>0,]$SCC4KO_DE
#non prom-enh loops are more upregulated than prom-enh loops
t.test(enh.up.SCC4,  non.enh.up.SCC4)
hist(enh.up.SCC4, breaks = 30) #both roughly 400 loops; some more extreme (>4) values in non.enh.up
hist(non.enh.up.SCC4, breaks=30)

SCC4Genes <- non.enh.statObj[non.enh.statObj$SCC4KO_DE>4,] #all prom.ctcf; all non-prom.polII. Genes: "WNT10A" "FEV" "UBE2QL1" "VGLL2" "KCNQ5" "GPR85" "COL5A1" "TMEM132E" 

enh.down.SCC4 <- enh.statObj[enh.statObj$SCC4KO_DE<0,]$SCC4KO_DE; non.enh.down.SCC4 <- non.enh.statObj[non.enh.statObj$SCC4KO_DE<0,]$SCC4KO_DE
t.test( enh.down.SCC4, non.enh.down.SCC4 )

t.test( abs(enh.statObj$lfcWaplKO_3.3[is.finite(enh.statObj$lfcWaplKO_3.3)]), abs(non.enh.statObj$lfcWaplKO_3.3[is.finite(non.enh.statObj$lfcWaplKO_3.3)]))

table(sign( enh.statObj$lfcWaplKO_3.3[is.finite(enh.statObj$lfcWaplKO_3.3)] ))
table(sign( non.enh.statObj$lfcWaplKO_3.3[is.finite(non.enh.statObj$lfcWaplKO_3.3)]) )


#Haarhuis et al find that for WAPLKO DE genes are more different in their lfc HiC reads
statObj.tss2 <- statObj.tss[is.finite(statObj.tss$lfcWaplKO_3.3),]
DE.statObj <- statObj.tss2[!is.na(statObj.tss2$WaplKO_DE),]
non.DE.statObj <-  statObj.tss2[is.na(statObj.tss2$WaplKO_DE),]

#other way around with lfc vs 
t.test(abs(DE.statObj$lfcWaplKO_3.3), abs(non.DE.statObj$lfcWaplKO_3.3))
t.test(DE.statObj$lfcWaplKO_3.3[DE.statObj$WaplKO_DE>0], non.DE.statObj$lfcWaplKO_3.3[DE.statObj$WaplKO_DE>0])
t.test(DE.statObj$lfcWaplKO_3.3[DE.statObj$WaplKO_DE<0], non.DE.statObj$lfcWaplKO_3.3[DE.statObj$WaplKO_DE<0]) #not




##some inbetween stats; ctcf loops are longer than non ctcf loops
noctcf <- FullStatObj[FullStatObj$ctcf.ctfc==0,]$loopDist
ctcf <- FullStatObj[FullStatObj$ctcf.ctfc==1,]$loopDist
t.test(noctcf,ctcf)

t.test(FullStatObj[FullStatObj$loopDist<3.35e5,]$vp.Hap1_CTCF_chip, FullStatObj[FullStatObj$loopDist>=3.35e5,]$vp.Hap1_CTCF_chip )
t.test(FullStatObj[FullStatObj$loopDist<3.35e5,]$vp.CTCF_WT, FullStatObj[FullStatObj$loopDist>=3.35e5,]$vp.CTCF_WT )
plot(FullStatObj$loopDist, FullStatObj$vp.CTCF_WT, type='h') #xlim=c(2e5,4e5)
plot(FullStatObj$loopDist, FullStatObj$vp.Hap1_CTCF_chip, type='h')


#ctcf loops are less nested
noctcf <- FullStatObj[FullStatObj$ctcf.ctfc==0,]$loopdepth
ctcf <- FullStatObj[FullStatObj$ctcf.ctfc==1,]$loopdepth
t.test(noctcf,ctcf)


#################TAD Analysis#########################
#source('TAD_analysis.R')
#load in relevant objects
TADcovList <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/TADs/TADcov_norm2.rds")
TADlocs <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/TADs/TADs.rds")
statObj.tss <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/statObjTSS.rds")

#alternative TAD definition based on loopdepth cutoff
# TADfromloops <- statObj.tss[statObj.tss$loopdepth<=1, ]
# TADlocs2 <- df2GR(TADfromloops)
# TADlocs2$TADid <- paste0('TAD.', 1:length(start(TADlocs)))
# 
# #find overlaps of original TADloc and this one and see which dont coincide 
# begin <- findOverlaps(TADlocs, TADlocs2, type = 'start', maxgap = 1e4)
# end <- findOverlaps(TADlocs, TADlocs2, type = 'end', maxgap = 1e4)
# sum(begin%in%end) #only 10 overlapping....


#turn into GRanges and find which loops are in which TAD
loopsforTAD <- df2GR(statObj.tss)
whichTAD <- findOverlaps(TADlocs, resize(loopsforTAD, width=1, fix = 'center'))

statObj.TAD <- statObj.tss[whichTAD@to,]
statObj.TAD$TADid <- whichTAD@from

#add DE for the loops in the TADs
TAD_loops_DE <- right_join(TADcovList[[1]], statObj.TAD[,c(ncol(statObj.TAD),grep('_DE',names(statObj.TAD)))])
TAD_loops_DE$TADid <- as.factor(TAD_loops_DE$TADid)
#sum DE per TAD so we get a DE score per TAD
TAD_DE <- setDT(TAD_loops_DE[,c(1,7:9)])[, lapply(.SD, max, na.rm=T), by = TADid] #sum #max

#for dealing with -Inf when using max
TAD_DE <- as.data.frame(TAD_DE)
TAD_DE[TAD_DE==-Inf]<-0

#Link TAD coverage object with TAD DE object
TADcovList[[1]] <- right_join(TADcovList[[1]], TAD_DE)

#add DE class like in Haarhuis
DE_factor <- function(x){
  if(x==0){return('non_DE')}
  if(x>=0){return('up')}
  if(x<=0){return('down')}
}
TADcovList[[1]]$WaplKO_DEfactor <- relevel(as.factor( sapply(TADcovList[[1]]$WaplKO_DE, DE_factor) ), 'non_DE')
TADcovList[[1]]$DKO_DEfactor <- relevel(as.factor( sapply(TADcovList[[1]]$DKO_DE, DE_factor) ), 'non_DE' )
TADcovList[[1]]$SCC4KO_DEfactor <- relevel(as.factor( sapply(TADcovList[[1]]$SCC4KO_DE, DE_factor) ), 'non_DE' )

#join same table also for n+1 and n+2
TADcovList[[2]] <-  right_join(TADcovList[[2]], TADcovList[[1]][,-(2:6)])
TADcovList[[3]] <-  right_join(TADcovList[[3]], TADcovList[[1]][,-(2:6)])

#add factor for neighbour level
for(i in 1:3){
  TADcovList[[i]]$neighbour <- as.factor(i-1)
  i=i+1
}
TADcov <- rbindlist(TADcovList)

#boxplot for all neighbours the expression category vs wapl/hap coverage of TADs
boxplot_TADexpr <- ggplot(data = TADcov, aes(x=WaplKO_DEfactor, y=log(WaplKO_3.3/Hap1))) + #SCC4KO
  geom_boxplot() +
  facet_grid(~neighbour) +
  theme_pubr(border = T) +
  coord_cartesian(ylim=c(-2,2)) 
boxplot_TADexpr 

non_DE <- TADcovList[[1]][ TADcovList[[1]]$WaplKO_DEfactor == 'non_DE', ][,c(3,4)]
up <- TADcovList[[1]][ TADcovList[[1]]$WaplKO_DEfactor == 'up', ][,c(3,4)]

LFC_nonDE <- log(non_DE$WaplKO_3.3/non_DE$Hap1)
LFC_up <-  log(up$WaplKO_3.3/up$Hap1) 
t.test(LFC_nonDE[is.finite(LFC_nonDE)], LFC_up[is.finite( LFC_up ) ] )


#order on expression change of TADs to compare loops within
TADcov <- TADcovList[[1]]
hist(TADcov$WaplKO_DE, xlim = c(0,5), breaks=30); summary(TADcov$WaplKO_DE) #> or < 2

big_down <- TADcov[TADcov$WaplKO_DE<=-1,]$TADid
big_up <- TADcov[TADcov$WaplKO_DE>=1,]$TADid #2 for sum and 1 for max
big_cov_up <- TADcov[ log(TADcov$WaplKO_3.3/TADcov$Hap1) >=1,]$TADid #summary(log(TADcov$WaplKO_3.3/TADcov$Hap1))
big_cov_down <- TADcov[ log(TADcov$WaplKO_3.3/TADcov$Hap1) <=-0.7,]$TADid
big_cov_up <- big_cov_up[!is.na(big_cov_up)];  big_cov_down <- big_cov_down[!is.na(big_cov_down)]

#big_down and big_up seem unuseable as group_1 because near-zero big-number effects
group1 <- group1[is.finite(group1$WaplKO_DE),]
non_group1 <- non_group1[is.finite(non_group1$WaplKO_DE),]
TPMtester <- data.frame(expression = c(group1$Hap1.RNAcount, non_group1$Hap1.RNAcount), 
                        group = c(rep('group 1', length(group1$Hap1_TPM)), rep('group 2', length(non_group1$Hap1_TPM)) ))
boxplot(expression~group, TPMtester, ylim=c(0,500))

##CHANGES
group1 <- statObj.TAD[statObj.TAD$TADid%in%big_cov_down,]
non_group1 <- statObj.TAD[!statObj.TAD$TADid%in%big_cov_down,]
  
t.test(group1$Hap1_TPM, non_group1$Hap1_TPM) #higher TPM in group1; much lower TPM in down; no diff in TPM in high lfc coverage TADs; 
t.test(group1$Hap1.RNAcount, non_group1$Hap1.RNAcount)
t.test(group1$WaplKO_3.3.RNAcount, non_group1$WaplKO_3.3.RNAcount)

t.test(group1$lfcWaplKO_3.3[is.finite(group1$lfcWaplKO_3.3)], non_group1$lfcWaplKO_3.3[is.finite(non_group1$lfcWaplKO_3.3)]) #difference (though tiny), in high lfc coverage TADs; not in high lfc TADs
t.test( apply(group1[, grep('H3K27ac', names(group1))],1, max ),apply(non_group1[, grep('H3K27ac', names(non_group1))],1, max ) ) #not different in H3K27ac sum up; 
#slightly with max up and max down, espc with stricter cutoff; bit lower in low lfc coverage TADs

t.test( apply(group1[, grep('3_CTCF', names(group1))],1, max ),apply(non_group1[, grep('3_CTCF', names(non_group1))],1, max ) ) #not different in WAPL CTCF or HAP1 1_CTCF
t.test( apply(group1[, grep('[.]CTCF', names(group1))],1, max ),apply(non_group1[, grep('[.]CTCF', names(non_group1))],1, max ) ) #higher CTCF from narrowpeaks in up with max
#lower ctcf from narrowpeaks in down with max strict cutoff; not with lenient cutoff; bit higher in low lfc coverage TADs... contradiction....




########
####look to couple TAD bridging to expression change/readcov/hist marks#####
#turn into GRanges and find which loops are in which TAD
loopsforTAD; TADlocs
TADlocs$TADid <- 1:length(TADlocs$TADid)

TADs.vp <- findOverlaps(resize(loopsforTAD, width=1e4, fix = 'start'), TADlocs, type='any', select = 'first')
TADs.anch <- findOverlaps(resize(loopsforTAD, width=1e4, fix = 'end'), TADlocs, type='any', select='first')

TADstretch <- TADs.anch-TADs.vp
#statObj.tss <- cbind(statObj.tss, TADstretch)
statObj.TAD2 <- cbind(statObj.tss, TADstretch)
statObj.TAD2 <- statObj.TAD2[!is.na(statObj.TAD2$TADstretch),]
barplot(table(statObj.TAD2$TADstretch))

boxplot(lfcWaplKO_3.3~TADstretch, data=statObj.TAD2, xlab='n.o. TADs spanned', ylab='Read Coverage LFC WaplKO')
boxplot(lfcSCC4KO~TADstretch, data=statObj.TAD2, xlab='n.o. TADs spanned', ylab='Read Coverage LFC SCC4KO')

boxplot(WaplKO_DE~TADstretch, data=statObj.TAD2, xlab='n.o. TADs spanned', ylab='DE in WaplKO')
boxplot(vp.Hap1_H3k27ac_chip~TADstretch, data=statObj.TAD2, xlab='n.o. TADs spanned', ylab='Enhancer Max value') #avg.H3K27ac

boxplot(vp.H3K27ac~TADstretch, data=statObj.TAD2, xlab='n.o. TADs spanned', ylab='Enhancer count value') #avg.H3K27ac
t.test(statObj.TAD2$vp.H3K27ac[statObj.TAD2$TADstretch==0], statObj.TAD2$vp.H3K27ac[statObj.TAD2$TADstretch>0])
t.test(statObj.TAD2$vp.H3K4me3[statObj.TAD2$TADstretch==0], statObj.TAD2$vp.H3K4me3[statObj.TAD2$TADstretch>0])
t.test(statObj.TAD2$vp.H3K56ac[statObj.TAD2$TADstretch==0], statObj.TAD2$vp.H3K56ac[statObj.TAD2$TADstretch>0])

t.test(statObj.TAD2$avg.H3K27ac[statObj.TAD2$TADstretch==0], statObj.TAD2$avg.H3K27ac[statObj.TAD2$TADstretch>0])

boxplot(`ctcf.ctcf`~TADstretch, data=statObj.TAD2, xlab='n.o. TADs spanned', ylab='TPM in WT')
t.test(statObj.TAD2$vp.CTCF_WT[statObj.TAD2$TADstretch==0], statObj.TAD2$vp.CTCF_WT[statObj.TAD2$TADstretch>0])

statObj.TAD2$lfc.CTCF <- log(statObj.TAD2$vp.Wapl3_3_CTCF_chip/statObj.TAD2$vp.Hap1_CTCF_chip)
boxplot(statObj.TAD2$lfc.CTCF[is.finite( statObj.TAD2$lfc.CTCF )] ~ TADstretch[is.finite( statObj.TAD2$lfc.CTCF )], data=statObj.TAD2, xlab='n.o. TADs spanned', ylab='Enhancer Max value') #avg.H3K27ac

t.test(statObj.TAD2$loopDist[statObj.TAD2$TADstretch==0], statObj.TAD2$loopDist[statObj.TAD2$TADstretch>0])


#############Gene expression looped vs non-looped###############
#find gene expression for nonlooped genes
 reads <- readRDS( "/delaat/group/iwan/Hap1_peakHiC/objects/Expression_levels/normHTseq/normHTseq.rds")
 reads <- tibble::rownames_to_column(as.data.frame(reads))
FullStatObj$ensembl <- dplyr::coalesce(FullStatObj$vp.Ens, FullStatObj$anch.Ens)
noloop_genename <- setdiff(TPM$ensembl, FullStatObj$ensembl)
nonlooped_genes <- reads[TPM$ensembl %in% noloop_genename,]

length(nonlooped_genes$rowname)

melted.stat <- melt(nonlooped_genes[,-1])
melted.stat$looping <- rep("non-looping", length(melted.stat$variable))

#get gene expression for looped genes
looped_genes1 <- statObj.tss[,61:65]
looped_genes2 <- statObj.tss[,66:70]
names(looped_genes1) <- c("Hap1", "WaplKO_3.3", "SCC4KO", "DKO", "WaplKO_1.14")
names(looped_genes2) <- c("Hap1", "WaplKO_3.3", "SCC4KO", "DKO", "WaplKO_1.14")

looped_genes <- rbind( looped_genes1, looped_genes2 )

length(unique(looped_genes$Hap1))

melted.stat2 <- melt(looped_genes, na.rm=T)
melted.stat2$looping <- rep("looping", length(melted.stat2$variable))

#combine gene expression data for both looped and non-looped genes
gene_expr <- rbind(melted.stat,melted.stat2)

boxplot_geneExpr <- ggplot(data = gene_expr, aes(x=looping, y=value)) +
  geom_boxplot() +
  facet_grid(~variable) +
  theme_pubr(border = T) +
  coord_cartesian(ylim=c(0,2e3)) 
#geom_signif(aes(group=variable), comparisons = list(c("short", "long")), test='t.test') 
boxplot_geneExpr
#ggsave(filename="/delaat/group/iwan/Hap1_peakHiC/plots/boxplot_geneExpr.pdf", boxplot_geneExpr)

dev.off()


###########analysis on newer Object#################
FullstatObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/FullStatObj3.rds") #statObj/avg_chip/FullStatObj.rds

#compare looplength over different looptypes
FullstatObj$looptype <- 'other'
FullstatObj[FullstatObj['ctcf.ctcf']==1,]$looptype <- "ctcf:ctcf"
FullstatObj[FullstatObj['prom.enh'] ==1,]$looptype <- "prom:enh"
FullstatObj[FullstatObj['prom.ctcf']==1,]$looptype <- "prom:ctcf"
FullstatObj[FullstatObj['prom.ctcf']==1 & FullstatObj['ctcf.ctcf']==1,]$looptype <- "prom.ctcf:ctcf"
FullstatObj[FullstatObj['prom.ctcf']==1 & FullstatObj['prom.enh'] ==1,]$looptype <- "prom:enh.ctcf"
FullstatObj[FullstatObj['ctcf.ctcf']==1  & FullstatObj['prom.enh'] ==1,]$looptype <- "prom.ctcf:enh.ctcf"

FullstatObj$looptype <- factor(FullstatObj$looptype, levels = c("ctcf:ctcf","prom.ctcf:ctcf","prom.ctcf:enh.ctcf","prom:enh.ctcf","prom:enh","prom:ctcf","other"))
theme_set(theme_bw())
ggplot(data = FullstatObj, aes(x=looptype)) +
  geom_bar(aes(fill='tomato3'))+
  xlab('Loop Type')+
  ylab('Frequency')+
  #theme_pubr()+
  theme(axis.text.x = element_text(angle=30,size = 18, vjust=0.5), axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20), legend.position = "none")


boxplot_loopdist <- ggplot(data = FullstatObj, aes(x=looptype, y=loopDist/1e3, fill=looptype)) +
  geom_boxplot() +
  fill_palette(palette = "GnBu")+
  xlab('Looptype (grouping A)')+
  ylab('looplength (kb)')+
  coord_cartesian(ylim = c(0,1e3))+
  theme_pubr(border = T) +
  scale_y_continuous()+
  theme(axis.text.x = element_text(angle=30,size = 16, vjust=0.6), axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20), legend.position = "none",
        strip.text.x = element_text(size = 12))
#geom_signif(aes(group=variable), comparisons = list(c("short", "long")), test='t.test') 
boxplot_loopdist


#####find expression correlations######
statObj.tss <- FullstatObj[FullstatObj$vp.TSS==1 | FullstatObj$anch.TSS==1,]

statObj.tss$DEcat <- 'no change'
statObj.tss[!is.na(statObj.tss$avg.WaplKO_DE) & statObj.tss$avg.WaplKO_DE>0.5,]$DEcat <- 'up in WaplKO' #0
statObj.tss[!is.na(statObj.tss$avg.WaplKO_DE) & statObj.tss$avg.WaplKO_DE<(-0.5),]$DEcat <- 'down in WaplKO' #0
statObj.tss[!is.na(statObj.tss$avg.WaplKO_DE) & statObj.tss$avg.WaplKO_DE>1,]$DEcat <- 'high up in WaplKO'
statObj.tss[!is.na(statObj.tss$avg.WaplKO_DE) & statObj.tss$avg.WaplKO_DE<(-1),]$DEcat <- 'high down in WaplKO'

statObj.tss$DEcat <- factor(statObj.tss$DEcat, levels = c("no change","up in WaplKO","down in WaplKO","high up in WaplKO","high down in WaplKO"))


boxplot_DE_looptyes <- ggplot(data = statObj.tss, aes(DEcat, fill=looptype)) +
  geom_bar(position = "fill") +
  theme_pubr(border = T)  +
  theme(axis.text.y=element_blank(), axis.text.x = element_text(angle = 20, hjust = 1, size = 8))
boxplot_DE_looptyes

#looplength for each expression category
boxplot_DElength <- ggplot(data=statObj.tss, aes(x=DEcat, y=loopDist)) +
  geom_boxplot() +
  theme_pubr(border = T) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))
boxplot_DElength

#within looptypes the looplength
boxplot3 <- ggplot(data=statObj.tss, aes(x=looptype, y=loopDist)) +
  geom_boxplot() +
  facet_grid(.~DEcat)+
  theme_pubr(border = T) +
  coord_cartesian(ylim=c(0,1e6)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))
boxplot3



#compare looplength over different looptypes
statObj.tss$looptype2 <- 'other'

statObj.tss[( (statObj.tss$vp.TSS>=1) & (statObj.tss$anch.CTCF_WT>=1) ) |
              ((statObj.tss$anch.TSS>=1) & (statObj.tss$vp.CTCF_WT>=1)) ,]$looptype2 <- "prom:ctcf"

statObj.tss[ (statObj.tss$vp.TSS>=1 & statObj.tss$anch.H3K4me3>=1) |
              (statObj.tss$anch.TSS>=1 & statObj.tss$vp.H3K4me3 >=1) ,]$looptype2 <- "prom:enh2"

statObj.tss[( (statObj.tss$vp.TSS>=1) & (statObj.tss$anch.CTCF_WT>=1 & statObj.tss$anch.H3K4me3>=1) ) |
          ((statObj.tss$anch.TSS>=1) & (statObj.tss$vp.CTCF_WT>=1 & statObj.tss$vp.H3K4me3 >=1 )) ,]$looptype2 <- "prom:ctcf.enh2"

statObj.tss[( (statObj.tss$vp.TSS>=1) & (statObj.tss$anch.CTCF_WT>=1 & statObj.tss$vp.H3K4me3>=1) ) |
              ((statObj.tss$anch.TSS>=1) & (statObj.tss$vp.CTCF_WT>=1 & statObj.tss$anch.H3K4me3 >=1 )) ,]$looptype2 <- "prom.enh2:ctcf"

statObj.tss[ (statObj.tss$vp.TSS>=1 & statObj.tss$vp.H3K4me3>=1 & statObj.tss$anch.CTCF_WT>=1 & statObj.tss$anch.H3K4me3>=1) |
               (statObj.tss$anch.TSS>=1 & statObj.tss$anch.H3K4me3>=1 & statObj.tss$vp.CTCF_WT>=1 & statObj.tss$vp.H3K4me3>=1) ,]$looptype2 <- "prom.enh2:ctcf.enh2"

statObj.tss$looptype2 <- factor(statObj.tss$looptype2, levels = c("prom:enh2","prom:ctcf","prom:ctcf.enh2","prom.enh2:ctcf","prom.enh2:ctcf.enh2","other"))

statObj.tss2 <- statObj.tss[statObj.tss$DEcat=='high up in WaplKO',]

ggplot(data = statObj.tss2, aes(x=looptype2)) +
  geom_bar(aes(fill='tomato3'))+
  xlab('Loop Type')+
  ylab('Frequency')+
  #theme_pubr()+
  theme(axis.text.x = element_text(angle=30,size = 18, vjust=0.6), axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20), legend.position = "none")


#within looptypes the looplength
boxplot4 <- ggplot(data=statObj.tss, aes(x=looptype2, y=loopDist)) +
  geom_boxplot() +
  facet_grid(.~DEcat)+
  theme_pubr(border = T) +
  coord_cartesian(ylim=c(0,1e6)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))
boxplot4

barplot(table(statObj.tss$looptype2))
barplot(table(statObj.tss[statObj.tss$DEcat=='high up in WaplKO',]$looptype2))







##poised genes for plotting
bivalent_genes_up_in_WaplKO <- arrange( statObj.tss[statObj.tss$DEcat=='high up in WaplKO',], desc(window.Hap1_H3K27me3_chip))[c('geneIDvp','geneIDanch')]
bivalent_genes_up_in_WaplKO <- coalesce(bivalent_genes_up_in_WaplKO[,1], bivalent_genes_up_in_WaplKO[,2])[1:20]

bivalent_genes_down_in_WaplKO <- arrange( statObj.tss[statObj.tss$DEcat=='high down in WaplKO',], desc(window.Hap1_H3K27me3_chip))[c('geneIDvp','geneIDanch')]
bivalent_genes_down_in_WaplKO <- coalesce(bivalent_genes_down_in_WaplKO[,1], bivalent_genes_down_in_WaplKO[,2])[1:20]

bivalent_genes_sligthly_down_in_WaplKO <- arrange( statObj.tss[statObj.tss$DEcat=='down in WaplKO',], desc(window.Hap1_H3K27me3_chip))[c('vp.Ens','anch.Ens')]
bivalent_genes_sligthly_down_in_WaplKO <- coalesce(bivalent_genes_sligthly_down_in_WaplKO[,1], bivalent_genes_sligthly_down_in_WaplKO[,2])[1:20]

bivalent_genes_sligthly_up_in_WaplKO <- arrange( statObj.tss[statObj.tss$DEcat=='up in WaplKO',], desc(window.Hap1_H3K27me3_chip))[c('vp.Ens','anch.Ens')]
bivalent_genes_sligthly_up_in_WaplKO <- coalesce(bivalent_genes_sligthly_up_in_WaplKO[,1], bivalent_genes_sligthly_up_in_WaplKO[,2])[1:20]



########general stats plots######
statObj$CTCF <- 'non-CTCF loop'
statObj$CTCF[statObj$ctcf.ctcf==1] <- 'CTCF loop'
statObj$CTCF <- as.factor(statObj$CTCF)

ggplot(data = melted.stat, aes(x=CTCF, y=value)) +
  geom_boxplot() +
  facet_grid(~variable) +
  theme_pubr(border = T)+
  coord_cartesian(ylim=c(-2,2))





##################older############################################
#################Loops of only TSS in VP and/or anchor##############################
statObjDE <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/DE_big_statObj.rds")
statObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/big_statObj_noDoubles.rds")
statObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/FullStatObj3.rds")
#############functional annotation##################
#david; GSA; GO; functional genome annotation
# check for looping structure differences
statObjDE2 <- statObjDE[statObjDE$loopdepth>30,]
statObj2 <- statObj[statObj$loopdepth<=1,]

DE.david <- list( c( statObjDE2$vp.Ens, statObjDE2$anch.Ens[statObjDE2$anch.Ens!="NA"] ))
write.table(DE.david,"/delaat/group/iwan/Hap1_peakHiC/statistics/DAVID/DEdavid_30deep.txt",sep="\n", row.names = F, col.names = F, quote = FALSE)

stat.david <- list( c( statObj2$vp.Ens[statObj2$vp.Ens!="NA"], statObj2$anch.Ens[statObj2$anch.Ens!="NA"] ))
write.table(DE.david,"/delaat/group/iwan/Hap1_peakHiC/statistics/DAVID/statDavid_1deep.txt",sep="\n", row.names = F, col.names = F, quote = FALSE)


statObjDE3 <- statObjDE[ order(as.numeric(statObjDE$DE.WaplKO_vs_Hap1), decreasing=T), ][1:500,]
dif.david <- list( c( statObjDE3$vp.Ens[as.character(statObjDE3$vp.Ens)!="NA"], statObjDE3$anch.Ens[as.character(statObjDE3$anch.Ens)!="NA"] ))
write.table(dif.david,"/delaat/group/iwan/Hap1_peakHiC/statistics/DAVID/DEdavid_500highest_WAPL.txt",sep="\n", row.names = F, col.names = F, quote = FALSE)

##############loopdepth correlations##############
 


statObj$loopdepthCat <- factor( sapply(statObj$loopdepth, loopdepthCat), levels=c("nestedLoop", "middleLoop", "topLoop") )

conds = c("DKO", "Hap1", "WaplKO_3.3", "SCC4KO", "WaplKO_1.14") 
measurevars <- colnames(statObj)[ (colnames(statObj) %in% paste0("lfc",conds)) ]

melted.stat <- melt(statObj, id.vars = "loopdepthCat", measure.vars = measurevars )
  
boxplot_depthcat <- ggplot(data = melted.stat, aes(x=loopdepthCat, y=value)) +
  geom_boxplot() +
  facet_grid(~variable) +
  theme_pubr(border = T)+
  coord_cartesian(ylim=c(-2,2)) 
#geom_signif(aes(group=variable), comparisons = list(c("short", "long")), test='t.test') 
boxplot_depthcat
#ggsave(filename="/delaat/group/iwan/Hap1_peakHiC/plots/boxplot_depthCat.pdf", boxplot_depthcat)

dev.off()


#########loopdepth and expression correlate################
#sum of ctcf scores vs loopdepthCat
statObj$loopdepthCat <- factor( sapply(statObj$loopdepth, loopdepthCat), levels=c("nestedLoop", "middleLoop", "topLoop") )

exprObj <- data.frame( loopdepthCat=statObj$loopdepthCat, Hap1_CTCF=rowMeans(statObj[,59:60]), Wapl3_3_CTCF=rowMeans(statObj[,63:64]), 
                          Hap1_SMC1=rowMeans(statObj[,61:62]), Wapl3_3_SMC1=rowMeans(statObj[,65:66]) )

#exprObj <- data.frame( loopdepthCat=statObj$loopdepthCat, Hap1_CTCF=statObj[,60], Wapl3_3_CTCF=statObj[,64], 
#                      Hap1_SMC1=statObj[,62], Wapl3_3_SMC1=statObj[,66] )

conds = c("Hap1", "Wapl3_3") 
measurevars <- colnames(exprObj)[ (colnames(exprObj) %in% paste0(conds,"_CTCF")) ] #SMC1

melted.stat <- melt(exprObj, id.vars = "loopdepthCat", measure.vars = measurevars )

boxplot_depth_expr <- ggplot(data = melted.stat, aes(x=loopdepthCat, y=value)) +
  geom_boxplot() +
  facet_grid(~variable) +
  theme_pubr(border = T)+
  coord_cartesian(ylim=c(0,4e3)) 
#geom_signif(aes(group=variable), comparisons = list(c("short", "long")), test='t.test') 
boxplot_depth_expr
#ggsave(filename="/delaat/group/iwan/Hap1_peakHiC/plots/boxplot_depth_expr.pdf", boxplot_depth_expr)

dev.off()


#########loopdepth and differential expression correlation################
statObjDE <- statObj.tss
statObjDE$ctcf.ctcf <- sapply(statObjDE$ctcf.ctfc, function(x){if(x==1){return("CTCFloop")}else{return("nonCTCFloop")}})
statObjDE$ctcf.ctcf <- factor(statObjDE$ctcf.ctcf, levels=c("CTCFloop", "nonCTCFloop"))

conds = c("DKO", "WaplKO", "SCC4KO") 
measurevars <- colnames(statObjDE)[ (colnames(statObjDE) %in% paste0(conds,"_DE")) ]
for(var in measurevars){ statObjDE[,var] <- abs(as.numeric(statObjDE[,var])) }

melted.stat <- melt(statObjDE, id.vars = "ctcf.ctcf", measure.vars = measurevars )

boxplot_ctcfExpr <- ggplot(data = statObj.tss, aes(x=ctcf, y=avg.Hap1.TPM)) +
  geom_boxplot() +
  facet_grid(~variable) +
  theme_pubr(border = T)+
  coord_cartesian(ylim=c(0,2)) +
  labs(y="Diff expression")
#  geom_signif(aes(group=variable), comparisons = list(c("CTCFloop", "nonCTCFloop")), test='t.test') 
boxplot_ctcfExpr

#ggsave(filename="/delaat/group/iwan/Hap1_peakHiC/plots/boxplot_activeGene.pdf", boxplot_activeGene)
dev.off()

####################differential expression vs CTCF level##################
#Hap1_CTCF=rowMeans(statObjDE[,59:60]), Wapl3_3_CTCF=rowMeans(statObjDE[,63:64]), Hap1_SMC1=rowMeans(statObjDE[,61:62]), Wapl3_3_SMC1=rowMeans(statObjDE[,65:66]) )

statObjDE$Wapl_CTCF.cat <- factor( sapply(rowMeans(statObjDE[,63:64], na.rm=T), ctcfCat), levels=c("lowCTCF", "mediumCTCF", "highCTCF") )
statObjDE$Hap_CTCF.cat <- factor( sapply(rowMeans(statObjDE[,59:60], na.rm=T), ctcfCat), levels=c("lowCTCF", "mediumCTCF", "highCTCF") )

#statObjDE$Wapl_CTCF.cat <- factor( sapply( apply(statObjDE[,63:64], 1, median), ctcfCat), levels=c("lowCTCF", "mediumCTCF", "highCTCF") )

cond <- "WaplKO"
measurevars <- colnames(statObjDE)[ (colnames(statObjDE) %in% paste0("DE.",cond,"_vs_Hap1")) ]
for(var in measurevars){ statObjDE[,var] <- abs(as.numeric(statObjDE[,var])) }

boxplot(DE.WaplKO_vs_Hap1~Wapl_CTCF.cat, data=statObjDE, main="DE vs CTCF",
        xlab="CTCF category", ylab="Differential Expression lfc", ylim=c(0,0.1)) 

##########continuous version of the above#########
exprObjDE <- data.frame( diff.expr=abs(as.numeric(statObjDE$DE.WaplKO_vs_Hap1)), WaplKO_CTCF= rowMeans(statObjDE[,62:63], na.rm=T) )
plot(exprObjDE$WaplKO_CTCF, exprObjDE$diff.expr, type='h', xlim=c(0,5000))

#now including also non-DE genes
statObj6 <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/big_statObj_noDoubles.rds")
statObj7 <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/DE_big_statObj.rds")
statObj6$geneIDvp <- unlist(statObj6$geneIDvp)
statObj7$geneIDvp <- unlist(statObj7$geneIDvp)
statObj6$geneIDanch <- unlist(statObj6$geneIDanch)
statObj7$geneIDanch <- unlist(statObj7$geneIDanch)
ExprStatObj <- left_join(statObj6, statObj7)
ExprStatObj$DE.WaplKO_vs_Hap1[is.na(ExprStatObj$DE.WaplK