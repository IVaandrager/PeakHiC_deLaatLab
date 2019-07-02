#Possible splits in the data to screen for any type of correlation

#general dataframe
statObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/FullStatObj3.rds")
statObj$loopDist <- statObj$loopDist/1e3
statObj[,grep('DE',names(statObj),value = T)] <- 
  apply(statObj[,grep('DE',names(statObj),value = T)], 2, as.numeric)



########The divides###########
#BASED ON BIGWIGS
#Divide in groups #different each time in this case it is for prom.enh loop based on max chip
plot.title <- "H3K27ac_chip_distinction"
distinction.text <- c("Non-promoter enhancer loop (0)", "Promoter enhancer loop (1)")
statObj$distinction <- 0
statObj[(statObj$vp.TSS >= 1 & statObj$anch.Hap1_H3k27ac_chip >= 7)|
          (statObj$anch.TSS >= 1 & statObj$vp.Hap1_H3k27ac_chip >= 7), ]$distinction <- 1
statObj$distinction <- as.factor(statObj$distinction)

plot.title <- "ATAC_chip_distinction"
distinction.text <- c("Promoter in closed region (0)", "Promoter in open region (1)")
statObj$distinction <- 0
statObj[(statObj$vp.TSS >= 1 & statObj$vp.Hap1_ATAC_chip >= 8.5)|
          (statObj$anch.TSS >= 1 & statObj$anch.Hap1_ATAC_chip >= 8.5), ]$distinction <- 1
statObj$distinction <- as.factor(statObj$distinction)

plot.title <- "ATAC_chip_afar_distinction"
distinction.text <- c("Promoter in closed region (0)", "Promoter in open region (1)")
statObj$distinction <- 0
statObj[(statObj$vp.TSS >= 1 & statObj$anch.Hap1_ATAC_chip >= 8.5)|
          (statObj$anch.TSS >= 1 & statObj$vp.Hap1_ATAC_chip >= 8.5), ]$distinction <- 1
statObj$distinction <- as.factor(statObj$distinction)

plot.title <- "CTCF_chip_afar_distinction"
distinction.text <- c("Promoter looped to CTCF (0)", "Promoter looped to CTCF (1)")
statObj$distinction <- 0
statObj[(statObj$vp.TSS >= 1 & statObj$anch.Hap1_CTCF_chip >= 8)|
          (statObj$anch.TSS >= 1 & statObj$vp.Hap1_CTCF_chip >= 8), ]$distinction <- 1
statObj$distinction <- as.factor(statObj$distinction)


plot.title <- "repressive_chip_distinction"
distinction.text <- c("non-repressed loop (0)", "repressed loop (1)")
statObj <- statObj[!is.na(statObj$window.Hap1_H3K27me3_chip), ]
statObj$distinction <- 0
statObj[ statObj$window.Hap1_H3K27me3_chip > 1 ,]$distinction <- 1
statObj$distinction <- as.factor(statObj$distinction)


hist(statObj$vp.Hap1_CTCF_chip, 200) 
abline(v=8)



##BASED ON COUNTS
#Divide in groups #different each time in this case it is for prom.enh loop based on enhancer counts
plot.title <- "CTCF-promoter_distinction"
distinction.text <- c("Promoter-nonCTCF loop (0)", "Promoter-CTCF loop (1)")
statObj = statObj[statObj$vp.TSS==1 | statObj$anch.TSS==1,]
statObj$distinction <- 0
statObj[(statObj$vp.TSS >= 1 & statObj$anch.CTCF_WT >= 1)|
          (statObj$anch.TSS >= 1 & statObj$anch.CTCF_WT >= 1), ]$distinction <- 1
statObj$distinction <- as.factor(statObj$distinction)

noctcf <- statObj[statObj$distinction==0,]
ctcf <- statObj[statObj$distinction==1,]       
count(noctcf[is.finite(noctcf$avg.WaplKO_DE),])/length(noctcf$loopID)
count(ctcf[is.finite(ctcf$avg.WaplKO_DE),])/length(noctcf$loopID)

#relabelled
plot.title <- "Multi_looptype_count_distinction"
distinction.text <- c("prom.ctcf.ctcf (0),  prom.enh.ctcf.ctcf(1), prom.enh.ctcf (2), prom.enh (3), prom.ctcf (4)","other (5)")
statObj = statObj[statObj$vp.TSS==1 | statObj$anch.TSS==1,]
statObj$distinction <- 5
statObj[statObj$prom.enh==1,]$distinction <- 3
statObj[statObj$prom.ctcf==1,]$distinction <- 4
statObj[statObj$prom.ctcf==1 & statObj$ctcf.ctcf==1,]$distinction <- 0
statObj[statObj$prom.enh==1 & statObj$prom.ctcf==1,]$distinction <- 2
statObj[statObj$prom.enh==1 & statObj$ctcf.ctcf==1,]$distinction <- 1
statObj$distinction <- as.factor(statObj$distinction)
table(statObj$distinction)
statObj2 <- statObj[is.finite(statObj$avg.WaplKO_DE),]
#statObj2 <- statObj[is.finite(statObj$avg.SCC4KO_DE),]

# 0    1    2    3    4    5 
# 7235 2369 4909 4380 5073 3657
#No strong indication for enh.hub vs enh indication
#t.test(statObj[statObj$distinction==2,]$avg.Hap1.TPM, statObj[statObj$distinction==3,]$avg.Hap1.TPM)
# group_by(statObj,distinction) %>%
#    summarise(sum((avg.WaplKO_DE<0), na.rm = T))
# group_by(statObj,distinction) %>%
#    summarise(sum((avg.WaplKO_DE>0), na.rm = T))

#use the above; now further classsify
plot.title <- "Looptype_and_DE_distinction"
distinction.text <- c("prom.enh.ctcf.up.DE(0), prom.enh.ctcf.down.DE (1), prom.enh.up.DE (2), prom.enh.down.DE (3)","other(4)")
statObj2$distinction2 <- 4  
statObj2[statObj2$prom.enh==1 & (statObj2$avg.WaplKO_DE<0),]$distinction2 <- 3
statObj2[statObj2$prom.enh==1 & (statObj2$avg.WaplKO_DE>0),]$distinction2 <- 2
statObj2[statObj2$prom.enh==1 & statObj2$max.CTCF_WT>0 & (statObj2$avg.WaplKO_DE<0),]$distinction2 <- 1
statObj2[statObj2$prom.enh==1 & statObj2$max.CTCF_WT>0  & (statObj2$avg.WaplKO_DE>0),]$distinction2 <- 0
# statObj2$distinction <- as.factor(statObj2$distinction2)
# statObj2$distinction2 <- NULL
# statObj <- statObj2
table(statObj2$distinction2)

#use the above; now further classsify
plot.title <- "Looptype_and_DE_distinction_SCC4KO"
distinction.text <- c("SCC4:prom.enh.ctcf.up.DE(0), prom.enh.ctcf.down.DE (1), prom.enh.up.DE (2), prom.enh.down.DE (3)","other(4)")
statObj2$distinction2 <- 4  
statObj2[statObj2$prom.enh==1 & (statObj2$avg.SCC4KO_DE<0),]$distinction2 <- 3
statObj2[statObj2$prom.enh==1 & (statObj2$avg.SCC4KO_DE>0),]$distinction2 <- 2
statObj2[statObj2$prom.enh==1 & statObj2$max.CTCF_WT>0 & (statObj2$avg.SCC4KO_DE<0),]$distinction2 <- 1
statObj2[statObj2$prom.enh==1 & statObj2$max.CTCF_WT>0  & (statObj2$avg.SCC4KO_DE>0),]$distinction2 <- 0
statObj2$distinction <- as.factor(statObj2$distinction2)
statObj2$distinction2 <- NULL
statObj <- statObj2
table(statObj$distinction)


plot.title <- "prom.enh_length_and_upDE_distinction3"
distinction.text <- c("long prom.enh nonDE(0) long prom.enh upDE(1) short prom.enh nonDE(2)", "short prom.enh upDE(3), other(4)")
statObj = statObj[statObj$vp.TSS==1 | statObj$anch.TSS==1,]
statObj = statObj[statObj$prom.enh==1,]
statObj$distinction <- 4
statObj[statObj$loopDist>=40 & !is.finite(statObj$avg.WaplKO_DE),]$distinction <- 0
statObj[statObj$loopDist>=40 & is.finite(statObj$avg.WaplKO_DE) & statObj$avg.WaplKO_DE>0,]$distinction <- 1
statObj[statObj$loopDist<=40 & !is.finite(statObj$avg.WaplKO_DE),]$distinction <- 2
statObj[statObj$loopDist<=40 & is.finite(statObj$avg.WaplKO_DE) & statObj$avg.WaplKO_DE>0,]$distinction <- 3
statObj$distinction <- as.factor(statObj$distinction)

#sum(statObj$distinction==0)/sum(statObj$distinction==1); sum(statObj$distinction==2)/sum(statObj$distinction==3)
wilcox.test(statObj[statObj$distinction==0 & statObj$avg.Hap1.TPM>1,]$avg.Hap1.TPM, statObj[statObj$distinction==1 & statObj$avg.Hap1.TPM>1,]$avg.Hap1.TPM)
#summary(statObj[statObj$distinction==1,]$max.H3K27ac)
#summary(statObj[statObj$distinction==1,]$avg.Hap1.TPM)

geneset <- statObj[statObj$prom.enh>=1 & is.finite(statObj$avg.WaplKO_DE) & statObj$avg.WaplKO_DE<0,]
geneset <- geneset[order(geneset$avg.WaplKO_DE),]
genes <- geneset[1:20,c('geneIDvp', 'geneIDanch')]
genes <- coalesce(genes[[1]],genes[[2]])


plot.title <- "H3K27ac_count_distinction"
distinction.text <- c("Promoter-non.enhancer loop (0)", "Promoter-enhancer loop (1)")
statObj = statObj[statObj$vp.TSS==1 | statObj$anch.TSS==1,]
statObj$distinction <- 0
statObj[(statObj$vp.TSS >= 1 & statObj$anch.CTCF_WT >= 1)|
          (statObj$anch.TSS >= 1 & statObj$vp.H3K27ac >= 1), ]$distinction <- 1
statObj$distinction <- as.factor(statObj$distinction)



boxplot_loopdist <- ggplot(data = statObj, aes(x=distinction, y=avg.Hap1.TPM, fill=distinction)) +
  geom_boxplot() +
  fill_palette(palette = "GnBu")+
  xlab('Looptype')+
  ylab('TPM')+
  coord_cartesian(ylim = c(0,100))+
  theme_pubr(border = T) +
  scale_x_discrete(labels=c("0" = "Promoter:non.enhancer","1" = "Promoter:Enhancer"))+
  theme(axis.text.x = element_text(,size = 16), axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20), legend.position = "none",
        strip.text.x = element_text(size = 12))
#geom_signif(aes(group=variable), comparisons = list(c("short", "long")), test='t.test') 
boxplot_loopdist

t.test(statObj[statObj$distinction==0,]$avg.Hap1.TPM,statObj[statObj$distinction==1,]$avg.Hap1.TPM)


plot.title <- "CTCF abundance distinction"
distinction.text <- c("no CTCF (0),  1 CTCF:TSS to 1 CTCF (1), 1 CTCF:TSS to many CTCF (2), many CTCF:TSS to 1 CTCF (3), many CTCF:TSS to many CTCF (4)")
statObj = statObj[statObj$vp.TSS==1 | statObj$anch.TSS==1,]
statObj$distinction <- 0
statObj[( (statObj$vp.CTCF_WT==1 & statObj$vp.TSS>=1) & (statObj$anch.CTCF_WT == 1 & statObj$anch.TSS ==0 ) ) |
          (statObj$vp.CTCF_WT==1 & statObj$vp.TSS==0) & (statObj$anch.CTCF_WT == 1 & statObj$anch.TSS >=1 ) ,]$distinction <- 1

statObj[( (statObj$vp.CTCF_WT==1 & statObj$vp.TSS>=1) & (statObj$anch.CTCF_WT > 1 & statObj$anch.TSS ==0 ) ) |
          (statObj$vp.CTCF_WT > 1 & statObj$vp.TSS==0) & (statObj$anch.CTCF_WT == 1 & statObj$anch.TSS >=1 ) ,]$distinction <- 2

statObj[( (statObj$vp.CTCF_WT>1 & statObj$vp.TSS>=1) & (statObj$anch.CTCF_WT == 1 & statObj$anch.TSS ==0 ) ) |
          (statObj$vp.CTCF_WT == 1 & statObj$vp.TSS==0) & (statObj$anch.CTCF_WT>1 & statObj$anch.TSS >=1 ) ,]$distinction <- 3

statObj[( (statObj$vp.CTCF_WT>1 & statObj$vp.TSS>=1) & (statObj$anch.CTCF_WT>1) ) |
          (statObj$vp.CTCF_WT>1) & (statObj$anch.CTCF_WT>1 & statObj$anch.TSS >=1 ) ,]$distinction <- 4

statObj$distinction <- as.factor(statObj$distinction)

noctcf <- statObj[statObj$distinction==0,]
ctcf <- statObj[statObj$distinction==2,]       
count(noctcf[is.finite(noctcf$avg.WaplKO_DE),])/length(noctcf$loopID)
count(ctcf[is.finite(ctcf$avg.WaplKO_DE),])/length(ctcf$loopID)

table(statObj$distinction)
# 0     1     2     3     4 
# 22570  3380   702   669   302 


#BASED ON READS, DISTANCE and (diff) EXPRESSION
plot.title <- "LFC HiC reads distinction"
distinction.text <- c("no change in HiC reads (0), decreased LFC HiC reads in WaplKO (0)", "increased LFC HiC reads in WaplKO (1)")
statObj <- statObj[!is.na(statObj$lfcWaplKO_3.3),]
statObj$distinction <- 0
statObj[statObj$lfcWaplKO_3.3<(-0.5), ]$distinction <- 1
statObj[statObj$lfcWaplKO_3.3>0.5, ]$distinction <- 2
statObj$distinction <- as.factor(statObj$distinction)


plot.title <- "Distance distinction"
distinction.text <- c("Short loops (<50kb) (0)", "Long loops (>50kb) (1)")
statObj$distinction <- 0
statObj[statObj$loopDist>50, ]$distinction <- 1
statObj$distinction <- as.factor(statObj$distinction)


plot.title <- "Differential Expression distinction2_from0"
distinction.text <- c("no change in expression (0), increased expression in WaplKO (1), decreased expression in WaplKO (2), much increased (3)",
                      "much decreased (4)")
statObj <- statObj[statObj$vp.TSS==1 | statObj$anch.TSS==1,]
#statObj <- statObj[!is.na(statObj$avg.WaplKO_DE),]
statObj$distinction <- 0
statObj[!is.na(statObj$avg.WaplKO_DE) & statObj$avg.WaplKO_DE>0,]$distinction <- 1
statObj[!is.na(statObj$avg.WaplKO_DE) & statObj$avg.WaplKO_DE<0,]$distinction <- 2
statObj[!is.na(statObj$avg.WaplKO_DE) & statObj$avg.WaplKO_DE>1,]$distinction <- 3
statObj[!is.na(statObj$avg.WaplKO_DE) & statObj$avg.WaplKO_DE<(-1),]$distinction <- 4
statObj$distinction <- as.factor(statObj$distinction)

table(statObj$distinction)
# 0     1     2     3     4 
# 25582   259  1351   219   212 

plot.title <- "Differential Expression distinction DeltaSCC4"
distinction.text <- c("no change in expression (0), increased expression in WaplKO (1), decreased expression in WaplKO (2), much increased (3)",
                      "much decreased (4)")
statObj <- statObj[statObj$vp.TSS==1 | statObj$anch.TSS==1,]
#statObj <- statObj[!is.na(statObj$avg.WaplKO_DE),]
statObj$distinction <- 0
statObj[!is.na(statObj$avg.SCC4KO_DE) & statObj$avg.SCC4KO_DE>0,]$distinction <- 1
statObj[!is.na(statObj$avg.SCC4KO_DE) & statObj$avg.SCC4KO_DE<0,]$distinction <- 2
statObj[!is.na(statObj$avg.SCC4KO_DE) & statObj$avg.SCC4KO_DE>1,]$distinction <- 3
statObj[!is.na(statObj$avg.SCC4KO_DE) & statObj$avg.SCC4KO_DE<(-1),]$distinction <- 4
statObj$distinction <- as.factor(statObj$distinction)
table(statObj$distinction)


plot.title <- "Abs Expression distinction"
distinction.text <- c("Low expression loops (TPM<40) (0)", "High expression loops (TPM>40) (1)")
statObj <- statObj[statObj$vp.TSS==1 | statObj$anch.TSS==1,]
statObj <- statObj[!is.na(statObj$avg.Hap1.TPM),]
statObj$distinction <- 0
statObj[statObj$avg.Hap1.TPM>40, ]$distinction <- 1
statObj$distinction <- as.factor(statObj$distinction)






