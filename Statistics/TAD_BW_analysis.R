library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(plotly)

significance.cat <- function(significance){
  if(is.na(significance) | significance==0){return('not sig.')} else
  if(significance>0.05){return('not sig.')} else
    if(significance>0.005){return('p<0.05')} else
      if(significance>0.0005){return('p<0.005')} else{
        return('p<0.0005')}
}

plot_violins <- function(statObj, outFolder, plot.title, distinction.text){
  #general column selection based on perspecitve
  #basic columns in every comparison
  basic.cols <- c(1,13,20,22:25,80,ncol(statObj))
  
  #cols specific for vp, anchor, avg, or windown
  vp.cols <- c(basic.cols, grep('^vp.', names(statObj))[-c(1,2,21)])
  anch.cols <-  c(basic.cols, grep('^anch[.]', names(statObj))[-c(19)])
  avg.cols <- c(basic.cols, sort( c(grep('^avg[.]', names(statObj)), grep('^max[.]', names(statObj))) ) )
  avg.cols <- avg.cols[c(1:9,28:35,10:27)] #to keep ordering as above
  cols.list <- list(vp.cols,anch.cols,avg.cols)
  
  for(perspective in cols.list){
    statObj.cat <- statObj[,perspective]
    
    ##add significance for colouring##
    statObj1 <- statObj.cat[statObj.cat$distinction==1,]
    statObj2 <- statObj.cat[statObj.cat$distinction==0,]
    sig <- c()
    
    for(column in names(statObj1)[-c(1,9)]){
    if(sum(!is.na(statObj1[[column]]))>0 & sum(!is.na(statObj2[[column]]))>0){
        p <- wilcox.test(as.numeric(statObj1[[column]]), as.numeric(statObj2[[column]]))
        sig <- c(sig,p$p.value)} else{ sig <- c(sig, NA)}
    }
    
    significance <- p.adjust(sig, method = "hochberg")
    significance <- sapply(as.vector(significance),significance.cat)
    ##     ##
    
    #to remove outliers for nicer plotting
    statObj.cat[,c(10:17)][ statObj.cat[,c(10:17)] >6 ] <- NA
    statObj.cat[,c(24,25)][ statObj.cat[,c(24,25)] >15 ] <- NA
    statObj.cat[,c(28:32)][ statObj.cat[,c(28:32)] >500 ] <- NA
    
    #melt into data frame that has values for all loops with extra column indicating
    #which mark and one indicating which group (0,1, etc)
    melted1 <- melt(statObj.cat[,c(1,2,4:17)], id.vars = c('loopID','distinction'))
    names(melted1) <- c('loop','type','mark','difference')
    melted1$significance <- factor( rep(significance[c(1,3:15)], each= table(melted1$mark)[[1]]),
                                    levels = c("not sig.","p<0.05","p<0.005","p<0.0005"))
    
    melted2 <- melt(statObj.cat[,c(1,9,18:27)], id.vars = c('loopID','distinction'))
    names(melted2) <- c('loop','type','mark','difference')
    melted2$significance <- factor( rep(significance[c(16:25)], each= table(melted2$mark)[[1]]),
                                    levels = c("not sig.","p<0.05","p<0.005","p<0.0005"))
    
    melted3 <- melt(statObj.cat[,c(1,3,9,28:35)], id.vars = c('loopID','distinction'))
    names(melted3) <- c('loop','type','mark','difference')
    melted3$significance <- factor( rep(significance[c(2,26:33)], each= table(melted3$mark)[[1]]),
                                    levels = c("not sig.","p<0.05","p<0.005","p<0.0005"))
    
    
    plot1 <- ggplot(melted1, aes(x=type, y=difference)) +
      geom_violin( aes(fill=significance) ) + #scale = "area"
      facet_wrap(.~mark, scales = 'free_y', nrow=1, ncol=length(unique(melted1$mark))) +
      theme_pubr(border = F)+
      scale_fill_brewer(palette = "OrRd",drop=FALSE) +
      #scale_y_continuous(expand = c(0,0)) +
      #coord_cartesian(ylim=c(-2,5)) +
      theme_pubr() +
      theme( strip.text.x = element_text(size = 6.5))+
      labs(
        title=sprintf('group %s vs group %s',distinction.text[1],distinction.text[2]),
        y= "Values") 
    
    plot2 <- ggplot(melted2, aes(x=type, y=difference)) +
      geom_violin(aes(fill=significance)) + #scale = "area"
      facet_wrap(.~mark, scales = 'free_y', nrow=1, ncol=length(unique(melted2$mark))) +
      theme_pubr(border = F)+
      scale_fill_brewer(palette = "OrRd", drop=FALSE) +
      #scale_y_continuous(expand = c(0,0)) +
      #coord_cartesian(ylim=c(0,15)) +
      theme_pubr() +
      theme(legend.position = "none",strip.text.x = element_text(size = 6.5))+
      labs(
        y= "Values") 
    
    plot3 <- ggplot(melted3, aes(x=type, y=difference)) +
      geom_violin(aes(fill=significance)) + #scale = "area"
      facet_wrap(.~mark, scales = 'free_y', nrow=1, ncol=length(unique(melted3$mark))) +
      theme_pubr(border = F)+
      scale_fill_brewer(palette = "OrRd",drop=FALSE) +
      #scale_y_continuous(expand = c(0,0)) +
      #coord_cartesian(ylim=c(-1e3,1e3)) +
      theme_pubr() +
      theme(legend.position = "none",strip.text.x = element_text(size = 6.5))+
      labs(
        x= 'Differentiable factors',
        y= "Values") 
    
    triple_plot <- ggarrange(plot1,plot2,plot3, nrow=3, align = "v", heights = c(3,2,2))
    
    print(triple_plot)
  }
  
  
  
  
  #plot 4 voor de window view; is alleen de chip waardes
  window.cols <- c(1, ncol(statObj), grep('^window.', names(statObj)) )
  statObj.window <- statObj[,window.cols]
  
  ##add significance for colouring##
  statObj1 <- statObj.window[statObj.window$distinction==1,]
  statObj2 <- statObj.window[statObj.window$distinction==0,]
  sig <- c()
  for(column in names(statObj1)[-c(1,2)]){
    if(length(statObj1[[column]])>0 & length(statObj2[[column]])>0){
    p <- wilcox.test(as.numeric(statObj1[[column]]), as.numeric(statObj2[[column]]))
    sig <- c(sig,p$p.value)} else{ sig <- c(sig, NA)}
  }
  significance <- p.adjust(sig, method = "hochberg")
  significance <- sapply(as.vector(significance),significance.cat)
  ##     ##
  
  #to remove outliers for nicer plotting
  #statObj.cat[,c(10:17)][ statObj.cat[,c(10:17)] >6 ] <- NA
  
  #melt into data frame that has values for all loops with extra column indicating
  #which mark and one indicating which group (0,1, etc)
  melted.window <- melt(statObj.window, id.vars = c('loopID','distinction'))
  names(melted.window) <- c('loop','type','mark','difference')
  melted.window$significance <- factor( rep(significance, each= table(melted.window$mark)[[1]]),
                                        levels = c("not sig.","p<0.05","p<0.005","p<0.0005"))
  
  plot.window <- ggplot(melted.window, aes(x=type, y=difference)) +
    geom_violin(aes(fill=significance)) + #scale = "area"
    facet_wrap(.~mark, scales = 'free_y', nrow=1, ncol=length(unique(melted2$mark))) +
    theme_pubr(border = F)+
    scale_fill_brewer(palette = "OrRd", drop=FALSE) +
    #scale_y_continuous(expand = c(0,0)) +
    #coord_cartesian(ylim=c(0,15)) +
    theme_pubr() +
    theme(legend.position = "none",strip.text.x = element_text(size = 6.5))+
    labs(
      title=sprintf('Window comparison of group %s vs group %s',distinction.text[1],distinction.text[2]),
      y= "Chip value") 
  
  print(plot.window)
}

#Stat framework:
#FullStatObj or statObj.tss and split into two categories based on some criteria:
#Dctcf Dsmc1, Dh3k27ac, Dreads, distance, Dexpr, TAD-based, etc. (any 2 categories)
#then we will make violin plots taken from bigwig for either VP, Anch, Window, (Avg)
#for these (3 bubble plots) we will compare marks LFC based on category division
#bubble height will be degree of change

setwd("/delaat/group/iwan/peakHiC/scripts")
outFolder <- '/delaat/group/iwan/Hap1_peakHiC/plots/violins/'

#take different statObjs from distinctions.R script

pdf(paste0(outFolder,plot.title,".pdf"), width = 20, height = 7)
plot_violins(statObj, outFolder, plot.title, distinction.text)
dev.off()






#Geert's analyse; TAD verschil in overeenkomst met intra-TAD verschillen
genes.info <- readRDS("hg38_genes_with_TADid.rds")
tads <- readRDS("Wapl_KO_vs_WT_sequential_TAD_contacts_tagMatrix.rds")

tads_up <- tads[order(tads$KOvsWT,decreasing=T),] #tads up in the KO versus WT
tads_up_M0 <- tads_up[tads_up$M==0,][1:10,]
tads_up_M1 <- tads_up[tads_up$M==1,][1:10,]
tads_up_M2 <- tads_up[tads_up$M==2,][1:10,]
tadz_up <- c(unique(tads_up_M0$TADid), unique(tads_up_M1$TADid))

tads_down <- tads[order(tads$KOvsWT),]
tads_down_M0 <- tads_down[tads_down$M==0,][1:10,]
tads_down_M1 <- tads_down[tads_down$M==1,][1:10,]
tads_down_M2 <- tads_down[tads_down$M==2,][1:10,]
tadz_down <- c(unique(tads_down_M0$TADid), unique(tads_down_M1$TADid))

TADgenes <- genes.info[genes.info$TAD%in%tads_down_M0,] #c(tadz_down,tadz_up) left_join(TADgenes, tads[,9:11], by=c("TAD"="TADid"))
Genes_down_M0 <- FullstatObj[FullstatObj$ensembl%in%TADgenes$ENSEMBL,]

#Then on all variants do stat.


TADgenePlot <- genes.info[genes.info$TAD%in%c(tadz_down,tadz_up),] 
gene_vps <- left_join(TADgenes, tads[,9:11], by=c("TAD"="TADid"))$gene_name
#plot and feed into DAVID





#export.bed(FullstatObj, './loops.bed')

# computeMatrix reference-point -S $Pol_II -R $loops.bed -a 5000 -b 5000 \
# -out heatmapdata -bs 250 -p 6 --missingDataAsZero
# 
# plotHeatmap -m heatmapdata --colorList white black --heatmapHeight 25 \
# --heatmapWidth 3 --samplesLabel "Pol II" -out heatmap_polII.png \
# --whatToShow 'heatmap and colorbar' --sortUsing max






 