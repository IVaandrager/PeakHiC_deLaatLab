library(peakC)
library("ggplot2")
library("reshape2")
library(rtracklayer)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(plotly)

#paste genes annotation under plot
addGenes <- function(genesGR,plotGR) {
  
  genesGR <- genesGR[unique(findOverlaps(genesGR,plotGR)@from)]
  chr <- as.vector(seqnames(plotGR))[1]
  
  if( !is.null( genesGR ) ) {
    
    minPos <- min( ranges( genesGR )@start )
    maxPos <- max( ranges( genesGR )@start + ranges( genesGR )@width - 1 )
    genesStart <- genesGR[ which( ranges( genesGR )@start >= minPos & ranges(genesGR)@start <= maxPos ), ]$geneName
    genesEnd <- genesGR [ which( ( ranges(genesGR)@start + ranges(genesGR)@width - 1 ) >= minPos & ( ranges(genesGR)@start + ranges(genesGR)@width - 1 ) <= maxPos ), ]$geneName
    genesWhole <- genesGR[ which( ranges( genesGR )@start < minPos & ( ranges( genesGR )@start + ranges( genesGR )@width - 1 )  > maxPos ), ]$geneName
    genesInRange <- unique( c( genesStart, genesEnd, genesWhole ) ) 
    genesGR <- genesGR[ which( genesGR$geneName %in% genesInRange ), ]
    genesGR <- genesGR[ which( genesGR@seqnames == chr ), ]
    
    posX1genes <- ranges( genesGR )@start
    posX2genes <- ranges( genesGR )@start + ranges( genesGR )@width - 1
    posX1genes[ which( posX1genes < minPos ) ] <- minPos
    posX2genes[ which( posX2genes > maxPos ) ] <- maxPos
    genesCenter <- ( posX1genes + posX2genes ) / 2
    
    genesGR.df <- data.frame(start=start(genesGR),end=end(genesGR),gene=genesGR$geneName, middle = genesCenter)
    return(genesGR.df)
  }}   
addHist <- function(genesGR,plotGR, mutant) {
  
  genesGR <- genesGR[unique(findOverlaps(genesGR,plotGR)@from)]
  chr <- as.vector(seqnames(plotGR))[1]
  
  if( !is.null( genesGR ) ) {
    
    minPos <- min( ranges( genesGR )@start )
    maxPos <- max( ranges( genesGR )@start + ranges( genesGR )@width - 1 )
    genesStart <- genesGR[ which( ranges( genesGR )@start >= minPos & ranges(genesGR)@start <= maxPos ), ]$vpID
    genesEnd <- genesGR [ which( ( ranges(genesGR)@start + ranges(genesGR)@width - 1 ) >= minPos & ( ranges(genesGR)@start + ranges(genesGR)@width - 1 ) <= maxPos ), ]$vpID
    genesWhole <- genesGR[ which( ranges( genesGR )@start < minPos & ( ranges( genesGR )@start + ranges( genesGR )@width - 1 )  > maxPos ), ]$vpID
    genesInRange <- unique( c( genesStart, genesEnd, genesWhole ) ) 
    genesGR <- genesGR[ which( genesGR$vpID %in% genesInRange ), ]
    genesGR <- genesGR[ which( genesGR@seqnames == chr ), ]
    
    if(length(grep(genesGR$type, pattern = "CTCF|SMC1"))>0){
      mutant.idx <- grep(genesGR$type, pattern = mutant)
      genesGR <- genesGR[mutant.idx]
    }
    
    posX1genes <- ranges( genesGR )@start
    posX2genes <- ranges( genesGR )@start + ranges( genesGR )@width - 1
    posX1genes[ which( posX1genes < minPos ) ] <- minPos
    posX2genes[ which( posX2genes > maxPos ) ] <- maxPos
    genesCenter <- ( posX1genes + posX2genes ) / 2
    
    genesGR.df <- data.frame(start=start(genesGR),end=end(genesGR), middle = genesCenter, type = genesGR$type)
    return(genesGR.df)
  }}    

#Call peakC on V4C for all conditions
max_in_GRange <- function(i, peakrange, resPeakC){
  peaks_in_peak <- queryHits( findOverlaps( IRanges( start=resPeakC$dbR[,1],end=resPeakC$dbR[,1]) , peakrange[i]))
  max_peak_rows <- resPeakC$dbR[peaks_in_peak,]
  max_peak_frag <-  max_peak_rows[,1][rowSums(max_peak_rows[,-1]) == max(rowSums(max_peak_rows[,-1]))][[1]]
}
reduce_peaks <- function(resPeakC, peakrange = 1000, gapwidth=250){
  peaks_r  <-IRanges(resPeakC$peak-peakrange, resPeakC$peak + peakrange)
  peaks_melted <- IRanges::reduce(peaks_r, min.gapwidth=gapwidth) 
  peaks_reduced <- lapply(FUN = max_in_GRange, seq_along(peaks_melted),peakrange=peaks_melted, resPeakC = resPeakC)
  ###############probleem voor lapply dat peaks_melted geen lijst is
  return(unlist(peaks_reduced))
}

#call the peaks
peakC_call <- function(V4C.list, slices, cond, vpS, wsize, qW, peakrange, gapwidth){
  resPeakC.list <- list()
  
  for(condition in cond){
    slice <- slices[condition]
    num.exp=length(slice[[1]])
    resPeakC <- combined.analysis(V4C.list[slice[[1]]], vp.pos=vpS, wSize=wsize, qW= qW, num.exp=num.exp)
    
    if(! length(resPeakC$peak)==0 ){
      resPeakC$peaks_reduced.x <- reduce_peaks(resPeakC, peakrange, gapwidth)
      resPeakC$peaks_reduced.y <- rowSums(resPeakC$dbR[,-1][match(resPeakC$peaks_reduced.x, resPeakC$dbR[,1]),])/(length(resPeakC$dbR[,-1]))
    }
    
    resPeakC.list[[names(slice)]] <- resPeakC
  }
  
  return(resPeakC.list)
}

#plot function
plot_V4C <- function(vpID,mutant,inFolder,loops, genomeObj){
  #first calculations to make V4C
  GR <- genomeObj$vpsGR[genomeObj$vpsGR$vpID==vpID]
  
  partID <- tryCatch(GR$partID[[1]],
                      error=function(cond) { stop(paste0("vpID does not seem to exist: ", vpID), call.=F) } )

  fragID <- GR$fragID[[1]]
  vpS <- start(GR)[[1]]
  chr <- as.character(seqnames(GR))[[1]]
  
  if(length(loops[loops$id==vpID | loops$anchor.vpID==vpID, 'loopID'])>=1){
  loopID<- loops[loops$id==vpID | loops$anchor.vpID==vpID, 'loopID'][[1]]
  anchor <- if(length(loops[loops$id==vpID,]$id)>=1){data.frame( pos=as.integer(loops[loops$id==vpID,'maxV4CscorePos'][[1]]), 
                                                             score= as.numeric(loops[loops$id==vpID,'maxV4Cscore'][[1]]) ) } else if(
                length(loops[loops$anchor.vpID==vpID,]$id)>=1)  {data.frame( pos=as.integer(loops[loops$anchor.vpID==vpID,'vp_X1'][[1]]), 
                                                                          score= as.numeric(loops[loops$anchor.vpID==vpID,'Hap1'][[1]]) )} }else{
      loopID <- 'noLoop';anchor <- data.frame(pos=0,score=0)}
  
  V4C.folder <- paste0(inFolder, "Hap1_V4C/")
  V4C.folder <- "/delaat/group/iwan/Hap1_peakHiC/analysis/V4Cs/"
  
  V4C.list <-  tryCatch(readRDS(paste0(V4C.folder,"vpReads_",partID,".rds"))[vpID],
                        error= function(cond) { stop(paste0("File does not seem to exist: ", paste0(V4C.folder,"vpReads_",partID,".rds"), "
                                                            Add the V4Cs partition file for ", partID), call.=F) })
  
  #V4C.list <- readRDS('/delaat/group/iwan/Hap1_peakHiC/analysis/Longer_V4C_examples/BCL11B_V4C')
  V4C.list <- lapply(V4C.list[[1]][[1]], function(GR) {data.frame(GR$GR$pos,GR$GR$reads)})
  
  slices <- list(DKO=1:2,WT=3:4,WaplKO=c(5:6,9),SCC4KO=7:8,All=1:9) 
  cond = c("WT", mutant)
  
  resPeakC.list <- peakC_call(V4C.list, slices=slices, cond=cond, vpS=vpS, wsize=21, qW=1, peakrange=1e3, gapwidth=250)
  
  
  ##########from here on for the plot
  y.ave <- list()
  y.ave[["pos"]] <- resPeakC.list[["WT"]]$dbR[,1]
  for(i in cond){
    y.ave[[i]] <- apply(resPeakC.list[[i]]$dbR[,2:(resPeakC.list[[i]]$num.exp+1)], 1, median)
  }
  
  df <- data.frame(pos =y.ave[["pos"]], WT= y.ave[["WT"]], mutant = y.ave[[mutant]])
  names(df)[[3]] <- mutant
  
  ymax <- max(df[-1])/2.5
  
  #plotting the overshoots in different colors
  winner <- apply(df[c("WT",mutant)],1,which.max)
  max.type <- ifelse(winner==1,"WT",mutant)
  max.values <- apply(df[c("WT",mutant)],1,max)
  
  maxipad <- data.frame(pos =y.ave[["pos"]], value = max.values, type = max.type)
  minipad <- data.frame(pos =y.ave[["pos"]], value =apply(df[c("WT",mutant)],1,min))
  
  overlay.df <- df[c('pos', 'WT', mutant)]
  overlay.melt <- melt(overlay.df, id.vars = 'pos')
  names(overlay.melt) <- c('pos', 'type', 'reads')
  
  ######annotation########
  #gene annotation
  genesGR <- readRDS(paste0(inFolder, "hg19_genesGR.rds"))
  plotGR <- GRanges(as.character(chr),IRanges(min(maxipad$pos),max(maxipad$pos)))
  genelist <- addGenes(genesGR,plotGR)
  
  #Differential expression annotation;  DKO_ SCC4_  Wapl2_  WaplKO_ WT2_   vs_Hap1
  diff.exp.files <- grep("WaplKO|SCC4|DKO", list.files(paste0(inFolder, "Diff_expression"), full.names = TRUE), value=T)
  diff.exp.list <-  lapply( diff.exp.files , readRDS)
  names(diff.exp.list) <- c("DKO","SCC4KO","WaplKO")
  
  diff.exp <- diff.exp.list[[mutant]][seqnames(diff.exp.list[[mutant]])==chr]
  expr.df <- data.frame(start=start(diff.exp), end=(end(diff.exp)), lfc=diff.exp$log2FoldChange)
  expr.df <- expr.df[expr.df$start>=maxipad$pos[1] & expr.df$end<=maxipad$pos[nrow(maxipad)],]
  if(length(expr.df[,1])==0){expr.df[1, ] <- c(0,0,0)}
  
  #Chip annotation CTCF, SMC, Histonemarks etc
  type.annot <- lapply(as.list(list.files(paste0(inFolder, "type_annotation"), full.names = TRUE)), readRDS)
  mutantCHIP <- grep(strsplit(mutant,split="")[[1]][1], c("WAPL","SCC4","DKO"), value=T)
  annotlist <- lapply(type.annot, addHist, plotGR=plotGR, mutant=mutantCHIP)
  annot.df <- data.frame()
  for( j in 1:2){partdf <- as.data.frame(annotlist[[j]])
  annot.df <- rbind(annot.df,partdf)}
  
  chip.file <- list.files(inFolder, pattern = '*[.]bw$', full.names = T)
  if(length(chip.file>0)){
    Ranges <- GRanges(seqnames=chr, IRanges(start=maxipad$pos[1], end=maxipad$pos[nrow(maxipad)]))
    chip.ctcf <- import.bw(chip.file, selection=Ranges)
    chip.ctcf <- as.data.frame(chip.ctcf)
  }else{chip.ctcf <- data.frame(start=NULL, end=NULL, score=NULL)}
  
  ######the actual ggplot######
  maxipad$type <- relevel(maxipad$type, "WT")
  
  readplot <- ggplot(maxipad, aes(x=pos, y=value, color=type, xend=pos, yend=0)) +
    geom_segment() +
    geom_segment(data=minipad, aes(color=NULL), color = 'grey80') +
    #geom_segment(data = called_peaks, aes(x=pos, y=10, xend=pos, yend=-10, color=type), inherit.aes = F, alpha = 0.3) +
    #geom_text(data = called_peaks, inherit.aes = F, aes(label = sprintf('%.0f kb', pos/1000), x= pos, y=2), angle = 90) +
    scale_color_brewer(palette = "Dark2") +
    geom_segment(data=anchor, inherit.aes=F, aes(x=pos-100, xend=pos+100, y=0, yend=score+0.2*score), color='red') +
    scale_x_continuous(labels = function(x){sprintf("%.1f Mb", x/1e6)}, limits = c(min(maxipad$pos), max(maxipad$pos))) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(ylim=c(0,ymax)) +
    theme_pubr() +
    labs(
      title = sprintf("Loop %s: V4C for viewpoint %s - %s", loopID, vpID, paste(unique(maxipad$type), collapse =" vs ")),
      x= sprintf("chromosome position %s", chr),
      y= "reads") +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      legend.title = element_blank(),
      legend.position = c(0.8,0.8),
      legend.direction = "vertical")
  
  annotplot <- ggplot(data = annot.df, aes(x=start, y= as.numeric(type)*-0.05-0.45, shape=type, color=type)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_rect(data = chip.ctcf, inherit.aes=F, aes(xmin = start, xmax = end,
                                                   ymin = -1.68, ymax = score*0.01-1.68)) +
    geom_rect(data = expr.df, inherit.aes = F, aes(xmin = start, xmax = end,
                                                   ymin = -1.15, ymax = lfc*0.05-1.15), fill=ifelse(expr.df$lfc>0,"darkgreen","darkred")) + #upregulation of mutant wrt WT is green
    geom_hline(yintercept = -1.15, size = 0.2, color = 'grey50') +
    geom_segment(data = genelist, inherit.aes = F, 
                 aes(x=start, xend=end, y= -0.4, yend =-0.4), arrow=arrow(length = unit(1,"mm"), ends="both", type="closed")) + 
    geom_text_repel(data = genelist,  inherit.aes = F, aes(label = gene, x= middle, y=-0.1), 
                    direction = "y", ylim = c(-0.45, 0), size=2, box.padding=0.1, segment.size=0.3) +
    scale_shape_manual(values = c(15:18,8:(length(levels(annot.df$type))-4))) +
    #geom_segment(data = called_peaks, aes(x=pos, y=10, xend=pos, yend=-10), inherit.aes = F, color = 'red', alpha = 0.3) +
    #scale_color_manual(values= c('#1b9e77','#d95f02')) + c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00')
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(labels = function(x){sprintf("%.1f Mb", x/1e6)}, limits = c(min(maxipad$pos), max(maxipad$pos))) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(ylim=c(-1.7,0)) +
    theme_pubr() +
    labs(
      title = "",
      x= sprintf("chromosome position %s", chr),
      y= "reads") +
    #labs(shape="marks", colour="celltype") +
    # guides(color = guide_legend(order=1),
    #        shape = guide_legend(order=2)) +
    theme(
      plot.title = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.spacing.y = unit(0.01,"pt"),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin=unit(c(-0.2,1,1,1), "cm"))
  
  double_plot <- ggarrange(readplot, annotplot, nrow = 2, align = "v", heights=c(1.2,1))
  double_plot
  return(double_plot)
}
