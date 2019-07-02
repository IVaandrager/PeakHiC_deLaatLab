library(dplyr)
library(reshape2)
library(GenomicRanges)
library(rtracklayer)

#functions
normaliseCov <- function(loopCov){
  median <- median(colSums(loopCov[-1]))
  normfact <- colSums(loopCov[-1])/median
  
  loopCov[-1] <- loopCov[-1] / normfact[col(loopCov[-1])]
  
  return(loopCov)
}

addCov <- function(loopCov, conds){
  sumloopCov=data.frame(temp=rep(0,length(loopCov[,1])))
  
  for(cond in conds){
    cond.ID <-  grep(names(loopCov), pattern=cond)
    sumloopCov$temp <- rowSums(loopCov[cond.ID])
    names(sumloopCov)[names(sumloopCov) == "temp"] <- cond
  }
  
  sumloopCov <- cbind(loopID = loopCov$loopID, sumloopCov)
  
  return(sumloopCov)
}

categoriseDist <- function(dist){
  if(dist<3.5e4){
    return("tiny")
  }
  if(dist<3.35e5){
    return("short")
  }
  if(dist<6.35e5){
    return("medium")
  }
  if(dist<9.35e5){
    return("long")
  }  else{
    return("huge")
  }
}

LFC <- function(statObj){
  statObj$lfcDKO <- log2(statObj$DKO/statObj$Hap1) ########if mut  more expr than WT then we get positive lfc
  statObj$lfcSCC4KO <- log2(statObj$SCC4KO/statObj$Hap1)
  statObj$lfcWaplKO_3.3 <- log2(statObj$WaplKO_3.3/statObj$Hap1)
  if(!is.null(statObj$WaplKO_1.14)){statObj$lfcWaplKO_1.14 <- log2(statObj$WaplKO_1.14/statObj$Hap1)}
     return(statObj)
}

specify_loops <- function(statObj3){
  statObj3$vp.CTCFinGene <- 0
  statObj3$anch.CTCFinGene <- rep(0, length(statObj3$loopID))
  statObj3[ grepl('CTCFinGene', statObj3$id), 'vp.CTCFinGene'] <- 1
  statObj3[ grepl('CTCFinGene', statObj3$anchor.vpID), 'anch.CTCFinGene' ] <- 1
  
  statObj3$ctcf.ctcf <- 0
  statObj3[statObj3$vp.CTCF_WT >= 1 & statObj3$anch.CTCF_WT >=1, ]$ctcf.ctcf <- 1 
  
  statObj3$prom.enh <- 0
  statObj3[statObj3$vp.TSS >= 1 & statObj3$anch.H3K27ac >=1, ]$prom.enh <- 1
  statObj3[statObj3$anch.TSS >= 1 & statObj3$vp.H3K27ac >=1, ]$prom.enh <- 1
  
  statObj3$prom.ctcf <- 0
  statObj3[statObj3$vp.CTCF_WT >= 1 & statObj3$anch.TSS >=1, ]$prom.ctcf <- 1
  statObj3[statObj3$anch.CTCF_WT >= 1 & statObj3$vp.TSS >=1, ]$prom.ctcf <- 1
  
  statObj3$active.gene <- rep(0, length(statObj$loopID))
  statObj3[statObj3$vp.TSS >= 1 & statObj3$vp.RNAPol >=1, ]$active.gene <- 1
  statObj3[statObj3$anch.TSS >= 1 & statObj3$anch.RNAPol >=1, ]$active.gene <- 1
  
  return(statObj3)
}

simplify_marker <- function(statObj){
  if(length(grep(statObj$type, pattern = "CTCF|SMC1"))>0){
  
  simple.idx <- grep(statObj$type, pattern = "CTCF")
  statObj$type[simple.idx] <-  "CTCF"
  
  simple.idx <- grep(statObj$type, pattern = "SMC1")
  statObj$type[simple.idx] <-  "SMC1"
  
  return(statObj)
  
} else{print("No CTCF or SMC1 marks found")}
}

looppartners <- function(statObj){
  anchor_as_VP <- findOverlaps(IRanges(start = statObj$maxV4CscorePos-5e3, statObj$maxV4CscorePos+5e3),
                               IRanges(start = statObj$vp_X1, end = statObj$vp_X2))
  
  return(recips)
}

protein_occurence <- function(vpID, type){ 
  occurence <- grepl(pattern = type, vpID[[2]])
  sum(occurence)
}

add_prot_count_cols <- function(prot.check, vpID_prot.list, statObj){
  
  for(prot in prot.check){
    prot_count <- lapply(vpID_prot.list, protein_occurence, type=prot)
    
    prot_count_df <- data.frame(id = names(prot_count), count = unlist(prot_count, use.names = F) )
    
    statObj <- left_join(statObj, prot_count_df, by=c("id" = "id"))
    
    colnames( statObj )[ which(colnames(statObj) == 'count' ) ]  <- paste0("vp.", prot) 
    
    statObj <- left_join(statObj, prot_count_df, by=c("anchor.vpID" = "id"))
    colnames(statObj)[which(colnames(statObj) == 'count' ) ] <- paste0("anch.", prot)
  }
  
  return(statObj)
  
}

getEnsembl <- function(ID, genemap){
  # print(x$geneIDvp)
  id <- which(as.character(genemap$gene_name) == as.character(ID))
  
  if(length(id)>0){
    ensembl <- as.character(genemap$ensembl[id][[1]]) } 
  else{
    ensembl <- NA
  }
  
  return(ensembl)
}

bwscores <- function(statObj, bw.files, takemax=T){
  #per chip file
  for(file in bw.files){
    name <- grep(strsplit(file,split='/')[[1]], pattern="*.bw", value=T)
    name <- sub("fc.bw","", name)
    print(name)
    
    #for VP and anch
    for(id in c("id", 'anchor.vpID')){
      #lookup once only for unique values of IDs
      uniqueObj <- statObj[!duplicated(statObj[id]),]
      RangesObj <- if(id=='id'){GRanges(seqnames = uniqueObj$chr, IRanges(start = uniqueObj$vp_X1-5e3, end = uniqueObj$vp_X1+5e3))} else{ 
              GRanges(seqnames = uniqueObj$chr, IRanges(start = uniqueObj$maxV4CscorePos-5e3, end = uniqueObj$maxV4CscorePos+5e3)) }
      
      bwcount <- import.bw(file, selection=RangesObj, as='NumericList')
      
      if(takemax){
        scores <- lapply(bwcount, max)}else{
        scores <- lapply(bwcount, function(x){ sum( x[c(TRUE, !diff(x) == 0)] ) } ) }
      
      bwcount_df <- data.frame(id=unique(statObj[id]), count=as.numeric(unlist(scores)) )
      
      statObj <- left_join(statObj, bwcount_df)
      names(statObj)[length(statObj)] <- ifelse(id=='id', paste0('vp.', name, "chip"), paste0('anch.',name, "chip"))
    }
  }
  return(statObj)
}

df2GR <- function(statObj4){
  
  df <- data.frame(chr=statObj4$chr, vp=statObj4$vp_X1, anch = statObj4$maxV4CscorePos)
  df[df$vp > df$anch,] <- df[df$vp > df$anch, c(1,3,2)]
  
  ranges <-  GRanges (seqnames=df$chr, IRanges(start = df$vp, end = df$anch) )
  
  return(ranges)
}

add_sig_DE <- function(row, cond.DE, location){
  #print(as.character(row['loopID']))
  cond <- names(cond.DE)
  cond.DE <- cond.DE[[1]]
  
  genenames <- row[location][[1]]
  
  DE.ids <- vector()
  enriched.gene <- vector()
  
  #if theres multiple genes, check for all if they occur in DE file
  splitnames <-  as.list(unlist(strsplit(genenames, split = ';')))
  for(genename in splitnames){
    DE.idx <- which(genename == cond.DE$transcript)
    
    if(length(DE.idx)>0){
      enriched.gene <- c(enriched.gene,genename)
      DE.ids <- c(DE.ids, DE.idx)}
  }
  
  #return all DE values
  if(length(DE.ids)>0){
    lfcs <- paste(cond.DE[DE.idx]$log2FoldChange, collapse = ';')
  } else{ lfcs <- NA }
  
  return(lfcs)
  
}

averageBW <- function(statObj4){
  
  vpcols <- grep('^vp..*chip$',names(statObj4), value = T)
 
  for(col in vpcols){
    new.colname <- sub("vp.","avg.", col)
    col.idx <- which(names(statObj4)==col)
    statObj4[,new.colname] <- apply(statObj4[, c(col.idx, col.idx+1) ], 1, max)
  }
  
  return(statObj4)
}

windowBW <- function(statObj4){} #still need to do

BW.LFC <- function(statObj4){
  #give LFC in chipmark between WaplKO and WT, for CTCF and SMC1, for vp, anch, max and window
  
  for(protein in c("CTCF","SMC1")){
    prot.cols <- grep(paste0(protein,'_chip$'),names(statObj4),value=T)
    
    for(loc in c('vp.','anch.','max.')){ #'window.'
      loc.col <- grep(loc,prot.cols, value=T)
      new.colname <- paste0(loc,"LFC_",protein,"_chip")

      statObj4[,new.colname] <- log2(statObj4[,loc.col[2]] / statObj4[,loc.col[1]])
    }
  }
  return(statObj4)
}

bw.window <- function(bw.file, bw.col.names, FullStatObj){
  #reads a binned bigwig file and extract the binvalues for ranges between loops, divide by n.o. bins
  bw_df <- read.delim(bw.file,header=T)
  
  names(bw_df) <- bw.col.names
  bw_GR <- makeGRangesFromDataFrame(bw_df,keep.extra.columns=T)
  
  window_ranges_GR <- df2GR(FullStatObj)
  
  ovlp <- findOverlaps(window_ranges_GR, bw_GR)
  norm.coef <- countOverlaps(window_ranges_GR,bw_GR)
  
  ovlp.df <- as.data.frame(ovlp) %>%
    group_by(queryHits) %>% 
    summarise(start = min(subjectHits), end = max(subjectHits))
  
  ovlp.ranges <- mapply(":", ovlp.df$start, ovlp.df$end)
  
  sums.list <-  lapply(ovlp.ranges, function(x,bw_df){ colSums(bw_df[x,-c(1:3)]) }, bw_df) #takes long
  bw.sums <- rbind_list(sums.list)
  
  return(bw.sums/norm.coef)
}



infile_cov <- "/delaat/group/iwan/Hap1_peakHiC/analysis/loop_coverage/loopCoverage.rds"
infile_loops <- "/delaat/group/iwan/Hap1_peakHiC/analysis/peaks/GW_nReps_9_peakHiC_wSize_41_qWr_1_alphaFDR_0.1_loops_reduced.rds"
outFolder <- '/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/'
genomeObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/Hap1_CTCF_TSS_enhancers_reduced2.rds")

#bigwig and vp prots/marks counts, TPM, and TAD filenames are in the code below


getStatObj <- function(infile_cov,infile_loops, genomeObj, geneCount){
  #######Initial loop object with loopcoverage##########
  conds = c("DKO", "Hap1", "WaplKO_3.3", "SCC4KO", "WaplKO_1.14")  #conds = c("DKO", "Hap1", "WaplKO_3.3", "SCC4KO")
  
  #get the rds object with all called loops
  loops <- readRDS(infile_loops)
  loops$id <- as.character(loops$id)
  loops$anchor.vpID <- as.character(loops$anchor.vpID)
  
  #get rds object with coverage for all called loops
  loopCov <- as.data.frame(readRDS(infile_cov))
  
  #add count per condition and normalise loopCoverage 
  sumloopCov <- addCov(loopCov, conds)
  normloopCov <- normaliseCov(sumloopCov)
  
  #create statistics object
  statObj <- right_join(normloopCov, loops, by='loopID')
  
  statObj$loopDist <- abs(statObj$maxV4CscorePos - statObj$vp_X1)
  #min(statObj$loopDist)
  #hist(statObj$loopDist,breaks=500)
  
  statObj$loopDistCat <- factor(sapply(statObj$loopDist, categoriseDist), levels = c("tiny", "short", "medium", "long", "huge"))
  
  #add HiC reads per loop viewpoint/anchor log fold change wrt wild type
  statObj <- LFC(statObj)
  lfcs <- c('lfcDKO', 'lfcSCC4KO', 'lfcWaplKO_3.3', 'lfcWaplKO_1.14')
  
  head(statObj)
  #saveRDS(statObj,  paste0(outFolder,"small_statObj.rds"))
  #readRDS(paste0(outFolder,"small_statObj.rds"))
  
  ######further extension of StatObj; left join with genomeObj info########
  genomeObj$vpsGR$vpID <- as.character(genomeObj$vpsGR$vpID)
  statObj$id <- as.character(statObj$id)
  genomeObj$vpsGR <- genomeObj$vpsGR[!duplicated(genomeObj$vpsGR$vpID)]
  
  statObj2 <- left_join(statObj, as.data.frame(elementMetadata(genomeObj$vpsGR)), by = c("id" = "vpID") )
  
  
  ######add protein marks at start and end of loop, only discrete counts here#########
  vpID_prot.list <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/vpID_to_marks.rds")
  
  prot.check <- list('TSS', 'CTCF_WT', 'SMC1_WT', 'H3K27ac','H3K4me3', 'H3K56ac','RNAPol')
  
  #add chip info, narrowpeak counts
  statObj3 <- add_prot_count_cols(prot.check, vpID_prot.list, statObj2)
  
  #mark some looptypes eg prom.enh or ctfc.ctfc based on chip counts
  statObj3 <- specify_loops(statObj3)
  
  #saveRDS(statObj3,  paste0(outFolder,"medium_statObj.rds"))
  #statObj3 <-readRDS(paste0(outFolder,"medium_statObj.rds"))
  
  #######bw CTCF and SMC1 and H3K27ac###########
  bw.files1 <- list.files('/delaat/group/geert/WAPL_HAP1/chipseq/MACS2', pattern = '*fc[.]bw$', full.names = T)
  bw.files1 <- bw.files1[-c(3,5)]
  
  bw.files2 <- list.files('/delaat/group/iwan/peakHiC/chip_for_vp/histone_marks2', pattern = '[.]bw$', full.names = T)
  bw.files <- c(bw.files1,bw.files2)

  statObj4 <- bwscores(statObj3, bw.files, takemax=T)
  statObj4 <- averageBW(statObj4)
  #statObj4 <- windowBW(statObj4) still need to do
  statObj4 <- BW.LFC(statObj4)
  
  #saveRDS(statObj4,paste0(outFolder,"statObj_allbw_noWindow.rds"))
  #statObj4 <- readRDS(paste0(outFolder,"statObj_allbw_noWindow.rds"))
  
  #####Add loop level########
  rangesObj <- df2GR(statObj4)
  loopdepth <- countOverlaps(rangesObj, type = "within")
  statObj4$loopdepth <- loopdepth
  #statObj4 <- readRDS(paste0(outFolder,"NewChip2.rds"))
  
  
  #####add TAD info#######
  TADcovList <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/TADs/TADcov_norm2.rds")
  TADlocs <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/TADs/TADs.rds")
  
  loopsforTAD <- df2GR(statObj4)
  whichTAD <- findOverlaps(TADlocs, resize(loopsforTAD, width=1, fix = 'center'))
  
  statObj.TAD <- statObj4[whichTAD@to,]
  statObj.TAD$TADid <- whichTAD@from
  statObj4 <- left_join(statObj4,statObj.TAD)
  
  TADs.vp <- findOverlaps(resize(loopsforTAD, width=1e4, fix = 'start'), TADlocs, type='any', select = 'first')
  TADs.anch <- findOverlaps(resize(loopsforTAD, width=1e4, fix = 'end'), TADlocs, type='any', select='first')
  
  #add the number of TADs a loop crosses
  statObj4$TADstretch <- TADs.anch-TADs.vp
  
  #saveRDS(statObj4,paste0(outFolder,"statObj_allbw_noWindow_TADs.rds"))
  #statObj4 <- readRDS(paste0(outFolder,"statObj_allbw_noWindow_TADs.rds"))
  #statObj4 <- statObj4[,-c(64:66)]
  
  #######from here on work with subset of statObj which contains TSS; later to be joined to original table#########
  #keep only entries with TSS in anchor or VP (relevant for expression)  
  statObj.tss <- statObj4[statObj4$vp.TSS | statObj4$anch.TSS >= 1,]
  
  #add which gene
  statObj.tss$geneIDvp <- sapply( statObj.tss$id, function(x){if(!grepl('CTCFinGene|vpID',x)){x}else{NA}} )
  statObj.tss$geneIDanch <- sapply( statObj.tss$anchor.vpID, function(x){if(!grepl('CTCFinGene|vpID',x)){x}else{NA}} )
  
  
  #load normalised reads
  #geneCount <- "/delaat/group/iwan/Hap1_peakHiC/objects/Expression_levels/normHTseq/normHTseq.rds"
  TPM <- "/delaat/group/iwan/Hap1_peakHiC/objects/Expression_levels/normHTseq/TPMcounts.rds"
  reads <- readRDS(TPM) #readRDS(geneCount)
  reads <- tibble::rownames_to_column(as.data.frame(reads))
  names(reads)[3] <- "WaplKO"
    
  #add Ensembl ID
  genemap <- readRDS("/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_ENSEMBL_Genes.rds")
  statObj.tss$vp.Ens <- unname(sapply(statObj.tss$geneIDvp, getEnsembl, genemap=genemap))
  statObj.tss$anch.Ens <- unname(sapply(statObj.tss$geneIDanch, getEnsembl, genemap=genemap))
  
  statObj.tss <- left_join(statObj.tss, reads, by = c("vp.Ens" = "rowname"))
  statObj.tss <- left_join(statObj.tss, reads, by = c("anch.Ens" = "rowname"))
  
  names(statObj.tss)[(length(statObj.tss)-9):length(statObj.tss)] <- paste0(rep(c('vp.','anch.'), each=5), names(reads[-1]), ".TPM")
  
  statObj.tss <- statObj.tss[ rowSums(!is.na(statObj.tss[,c(length(statObj.tss)-9):length(statObj.tss)]))>0, ]
 
  statObj.tss$avg.Hap1.TPM <- as.numeric(dplyr::coalesce( statObj.tss$vp.Hap1.TPM,  statObj.tss$anch.Hap1.TPM))
  statObj.tss$avg.WaplKO.TPM <- as.numeric(dplyr::coalesce( statObj.tss$vp.WaplKO.TPM,  statObj.tss$anch.WaplKO.TPM))
  statObj.tss$avg.SCC4KO.TPM <- as.numeric(dplyr::coalesce( statObj.tss$vp.SCC4KO.TPM,  statObj.tss$anch.SCC4KO.TPM))
  statObj.tss$avg.DKO.TPM <- as.numeric(dplyr::coalesce( statObj.tss$vp.DKO.TPM,  statObj.tss$anch.DKO.TPM))
  statObj.tss$avg.WaplKO_1.14.TPM <- as.numeric(dplyr::coalesce( statObj.tss$vp.WaplKO_1.14.TPM,  statObj.tss$anch.WaplKO_1.14.TPM))
  
  head(statObj.tss)
  
  
  ############Add Differential xpression Data################
  diff_exp <-  lapply(list.files("/delaat/group/iwan/Hap1_peakHiC/objects/Expression_levels/diff_expression/", full.names = T)[c(1,2,4)],
                      FUN = readRDS)
  
  diff_comp <- c('DKO_vs_Hap1', 'SCC4_vs_Hap1', 'WaplKO_vs_Hap1')
  names(diff_exp) <- diff_comp
  
  #add DE
  for(location in c('geneIDvp','geneIDanch')){
    for(cond in diff_comp){
      cond.DE <- diff_exp[cond]
      lfcs <- apply(statObj.tss, 1, add_sig_DE, cond.DE=cond.DE, location=location)
      statObj.tss[, paste0( ifelse(grepl('vp',location),'vp','anch'), '.DE.', cond ) ] <- lfcs
    }
  }
  
  ###combine anch and vp DE
  statObj.tss$avg.WaplKO_DE <- as.numeric(dplyr::coalesce( statObj.tss$vp.DE.WaplKO_vs_Hap1,  statObj.tss$anch.DE.WaplKO_vs_Hap1))
  statObj.tss$avg.DKO_DE <- as.numeric(dplyr::coalesce( statObj.tss$vp.DE.DKO_vs_Hap1,  statObj.tss$anch.DE.DKO_vs_Hap1))
  statObj.tss$avg.SCC4KO_DE <- as.numeric(dplyr::coalesce( statObj.tss$vp.DE.SCC4_vs_Hap1,  statObj.tss$anch.DE.SCC4_vs_Hap1))

  saveRDS(statObj.tss, paste0(outFolder,"statObj_TSS2.rds"))
  
  
  #######create full statObj##########
  FullStatObj <- left_join(statObj4, statObj.tss[,-(2:6)])
  #FullStatObj <- FullStatObj[,-c(71:75)]
  
  #Forgot to average counts as well
  for(i in seq(28,43,2)){
    FullStatObj[, sub('vp.', 'max.', names(FullStatObj)[i] ) ] <- apply(FullStatObj[, c(i,i+1)], 1, max)
  }
  
  
  #Finally, obtaining all the window-view values of the chip
  binned_bw_folder <- '/delaat/group/iwan/peakHiC/chip_for_vp/histone_marks2/'
  bw.files <- paste0(binned_bw_folder, c("scores_per_bin_ctcf.tab", "scores_per_bin.tab"))
  bw.col.name.list <- list(c("chr",'start','end','Hap1_CTCF',"Hap1_H3k27ac","Hap1_SMC1", "Wapl3_3_CTCF","Wapl3_3_SMC1"),
                           c("chr",'start','end','Hap1_ATAC',"Hap1_H2AK119ub1","Hap1_H3K27me3"))
  
  #get the binned reads from the bigwig of the interloop region. Take average value. 
  FullStatObj[, paste0("window.",bw.col.name.list[[1]][-c(1:3)],"_chip") ] <- 
    bw.window(bw.file=bw.files[[1]], bw.col.names=bw.col.name.list[[1]], FullStatObj) 
  
  FullStatObj[, paste0("window.",bw.col.name.list[[2]][-c(1:3)], "_chip") ] <- 
    bw.window(bw.file=bw.files[[2]], bw.col.names=bw.col.name.list[[2]], FullStatObj)
  
  FullStatObj <- BW.LFC(FullStatObj) #function interior adapted to window
  
  #saveRDS(FullStatObj, paste0(outFolder,"FullStatObj3.rds"))
}


##combinatorial prom enh ctcf categories
add_categories <- function(statObj){
  for(pers in c('vp','anch')){
    statObj[,paste0(pers,".prom")] <- 0
    statObj[,paste0(pers,".prom")][statObj[,paste0(pers,".TSS")]>0] <- 1
    
    statObj[,paste0(pers,".ctcf")] <- 0
    statObj[,paste0(pers,".ctcf")][statObj[,paste0(pers,".CTCF_WT")]>0] <- 1
    
    statObj[,paste0(pers,".enh")] <- 0
    statObj[,paste0(pers,".enh")][statObj[,paste0(pers,".H3K27ac")]>0] <- 1
    
    statObj[,paste0(pers,".prom.enh")] <- 0
    statObj[,paste0(pers,".prom.enh")][statObj[,paste0(pers,".TSS")]>0 & statObj[,paste0(pers,".H3K27ac")]>0] <- 1
    
    statObj[,paste0(pers,".prom.ctcf")] <- 0
    statObj[,paste0(pers,".prom.ctcf")][statObj[,paste0(pers,".TSS")]>0 & statObj[,paste0(pers,".CTCF_WT")]>0] <- 1
    
    statObj[,paste0(pers,".ctcf.enh")] <- 0
    statObj[,paste0(pers,".ctcf.enh")][statObj[,paste0(pers,".H3K27ac")]>0 & statObj[,paste0(pers,".CTCF_WT")]>0] <- 1
    
    statObj[,paste0(pers,".prom.ctcf.enh")] <- 0
    statObj[,paste0(pers,".prom.ctcf.enh")][statObj[,paste0(pers,".ctcf.enh")]>0 & statObj[,paste0(pers,".TSS")]>0] <- 1
  }
  
  vp.loopstr <-  statObj[,c((ncol(statObj)-13):(ncol(statObj)-7))]
  vp.loopstr <- apply(vp.loopstr,1,paste0,collapse="")
  vp.loopstr <- strtoi(vp.loopstr, base=2)
  
  anch.loopstr <- statObj[,c((ncol(statObj)-6):ncol(statObj))]
  anch.loopstr <- apply(anch.loopstr,1,paste,collapse="")
  anch.loopstr <- strtoi(anch.loopstr, base=2)
  
  strtoi("0000001", base = 2)
  c("0","P","C","E","P.C","P.E","C.E","P.C.E")
  
  return(statObj)
}





###For Ruiqi####
FullStatObj <- readRDS(paste0(outFolder,"FullStatObj.rds"))
statObj.ruiqi <- FullStatObj[ FullStatObj$vp.CTCFinGene | FullStatObj$anch.CTCFinGene >= 1,]

#add info if RNApol within 1kb of GeneinCTCF viewpoint region
enhancer.dir <- "/delaat/group/iwan/peakHiC/chip_for_vp/histone_marks/"
marks <- as.list(rep(c('RNAPolII'), each = 2))

enhfiles <- list.files(path = enhancer.dir, pattern = "POL2.*[.]cod$",full.names = TRUE)
hg19_files <- lapply(enhfiles, read.delim) 
hg19_GR <- lapply(hg19_files, GRanges)
hg19_GR_ext <- mapply(GR = hg19_GR, mark = marks, function(GR,mark){ GR$type <- mark; GR }  )
hg19_enh <- do.call("c", hg19_GR_ext)

#make hg19_enh genomeObj 
enh.genomeObj <- GRanges(seqnames = seqnames(hg19_enh), ranges = resize(ranges(hg19_enh), width=1, fix="center"),
                         strand = strand(hg19_enh), type = hg19_enh$type)

vpOverlap <- countOverlaps(GRanges(seqnames=statObj.ruiqi$chr, IRanges(start=statObj.ruiqi$vp_X1-5e2,end=statObj.ruiqi$vp_X1+5e2)),
                           enh.genomeObj)

anchOverlap <- countOverlaps(GRanges(seqnames=statObj.ruiqi$chr, 
                                    IRanges(start=statObj.ruiqi$maxV4CscorePos-5e2,end=statObj.ruiqi$maxV4CscorePos+5e2)), enh.genomeObj)

statObj.ruiqi$vp.PolII.1kb <- vpOverlap
statObj.ruiqi$anch.PolII.1kb <- anchOverlap

#saveRDS(statObj.ruiqi, "/delaat/group/iwan/Hap1_peakHiC/ruiqi/newer/manyLoops.rds")

statObj.ruiqi2 <- statObj.ruiqi[ statObj.ruiqi$vp.PolII.1kb > 0 & grepl("CTCFinGene", statObj.ruiqi$id) | 
                                  statObj.ruiqi$anch.PolII.1kb  > 0 & grepl("CTCFinGene", statObj.ruiqi$anchor.vpID), ]


#saveRDS(statObj.ruiqi2, "/delaat/group/iwan/Hap1_peakHiC/ruiqi/newer/CTCFinGene_PolII1kb_loops_withCHIP.rds")
#write.csv( statObj.ruiqi2, "/delaat/group/iwan/Hap1_peakHiC/ruiqi/newer/CTCFinGene_PolII1kb_loops_withCHIP.txt")






################junk##############
# onces <- count2[[1]][c(TRUE, !diff(count[[1]]) == 0)] 
# 
# count1=import.bw(file, selection=GRanges(seqnames='chr1',IRanges(start=990074, end= 1000074)))
# count2=import.bw(file, selection=GRanges(seqnames='chr1',IRanges(start=71330674, end=71340674)), as='NumericList')
# sum ( count1$score )

# bwscores <- function(row, file){
#   
#   bigwig <- import.bw(file, selection = row)
#   bwscore <- as.integer(sum(bigwig$score))
#   
#   return(bwscore)
# }
# #file <- '/delaat/group/geert/WAPL_HAP1/chipseq/MACS2/Hap1_CTCF_fc.bw'
# for(file in bw.files){
#   name <- grep(strsplit(file,split='/')[[1]], pattern="*.bw", value=T)
#   name <- sub("fc.bw","", name)
#   print(name)
#   uniqueObj <- statObj6[unique(statObj6$id),]
#   RangesObj <- GRanges(seqnames = uniqueObj$chr, IRanges(start = uniqueObj$vp_X1-5e3, end = uniqueObj$vp_X1+5e3))
#   bwcount <- import.bw(file, selection=RangesObj, as='NumericList')
#   statObj6 <- left_join(statObj6, bwcount_df)
#   names(statObj6)[length(statObj6)] <- paste0('vp.', name, "chip")
#   
#   
#   import.bw(file, selection=RangesObj, as='NumericList')
#   bwcount <- lapply(as.list(as(statObject, "GRangesList")), bwscores, file=file)
#   bwcount_df <- data.frame(id=as.character(unique(statObj6$id)), count=unlist(bwcount))
#   
#   
#   
#   bigwig <- import.bw(file, selection=as.list(as(statObject, "GRangesList"))[[1]] )
#   bwscore <- as.integer(sum(bigwig$score))
#   
#   
#   uniqueObj <- statObj6[!duplicated(statObj6$anchor.vpID),]
#   statObject <- GRanges(seqnames = uniqueObj$chr, IRanges(start = uniqueObj$maxV4CsorePos-5e3, end = uniqueObj$maxV4CsorePos+5e3))
#   bwcount <- lapply(as.list(as(statObject, "GRangesList")), bwscores, file=file)
#   bwcount_df <- data.frame(anchor.vpID=as.character(unique(statObj6$anchor.vpID)), count=unlist(bwcount))
#   
#   statObj6 <- left_join(statObj6, bwcount_df)
#   names(statObj6)[length(statObj6)] <- paste0('anch.',name, "chip")
# }




# bwscores <- function(row, file, anchor=F){
#   if(anchor){
#     anch.avg <- (row$anchor_X1+row$anchor_X2)/2
#     range <- GRanges(seqnames = row$chr,ranges = IRanges(anch.avg-5e3, anch.avg+5e3))
#   }else{
#     range <- GRanges(seqnames = row$chr,ranges = IRanges(row$vp_X1-5e3, row$vp_X2+5e3))}
#   
#   selection <- BigWigSelection(range)
#   bigwig <- import.bw(file, selection = selection)
#   
#   bwscore <- sum(bigwig$score)
#   
#   return(bwscore)
# }

# tshirt <- c("M","S","S","M","XL")
# tshirt_factor <- factor(tshirt, ordered = T, levels = c("S","M","XL"))
# #ordered = T for calling things as factor[1]<factor[2]
# 
# ############ Junk
# type.annot <- lapply(as.list(list.files(paste0(annot.dir,"type_annotation"), full.names = TRUE)), readRDS)
# genesGR <- readRDS("/delaat/group/iwan/peakHiC/rds/genomeObj/hg19_genesGR.rds")
# 
# 
###find over all bigwig regions (for genes) the nearest one (possbibly within X bp window)
# for(chrom in paste0('chr', c(1:22,"X"))){
#   
#   statObj_chr <- statObj6[statObj6$chr == chrom,]
#   max(statObj6[statObj6$chr=='chr1',]$vp_X1)
#   
#   rl <- RangesList()
#   rl[[chrom]] <- IRanges(1, max(statObj_chr$vp_X2))
#   
#   selection <- BigWigSelection(rl)
#   bigwig_chr <- import.bw(genecounts[[3]], selection=selection)
#   
#   #findnearestneighbour
# }
# 
# 
# for(chrom in paste0('chr', c(1:22,"X"))){
#   
#   statObj_chr <- statObj6[statObj6$chr == chrom,]
#   genemap_chr <- genemap[seqnames(genemap) == chrom]
#   
#   chr_list <-  lapply(unlist(statObj_chr$geneIDvp), get_gene_ranges, genemap=genemap_chr)
#   chr_vector <- do.call(c, unlist(chr_list))
#   
#   rl[[chrom]] <- chr_vector
# }
# 
# selection <- BigWigSelection(rl)
# 
# 
# 
# get_gene_ranges <- function(ID, genemap){
#   # print(x$geneIDvp)
#   id <- which(as.character(genemap$gene_name) == as.character(ID))
#   
#   if(length(id)>0){
#     ranges <- ranges(genemap[id,0]) } 
#   else{
#     return()
#   }
#   
#   return(ranges)
# }
# 
# 
# 
# 
# library(broom) 
# dist.type <- group_by(statObj3, ctcf.ctfc) %>%
#   summary(meanLoopDist = mean(loopDist))
# 
# #summarise(avg= t(ctcf.ctfc) )
# 
# ggplot(data = melted.stat, aes(x=loopDistCat, y=value)) +
#   geom_boxplot() +
#   facet_grid(~variable) +
#   theme_pubr(border = T)+
#   coord_cartesian(ylim=c(-2,2)) 
#
#old functions

#onegene <- function(ID){
# if(grepl(";", ID)){
#   less <- strsplit(ID, split = ';')[[1]][[1]]
# }else{
#   less <- ID
# }
# return(less)
# }

## statObj.tss$geneIDvp = statObj.tss$id #for all id that grep(non ctcfingene and vpID_XXX)
#Is it necessary or can we just use $id directly for function below?
# statObj.tss <-addGeneIDs(statObj.tss, vpID_prot.list)
# 
# addGeneIDs <- function(statObj.tss, vpID_prot.list){
#   vpID_prot.list.tss <- vpID_prot.list[! grepl( "vpID*|CTCFinGene*", names(vpID_prot.list)) ]
#   
#   for(vpID in vpID_prot.list.tss){
#     genename <-as.character(vpID_prot.list[[vpID]]$marks)[1]
#     
#     statObj.tss$geneIDvp[statObj.tss$id==vpID] <- genename
#     statObj.tss$geneIDanch[statObj.tss$anchor.vpID==vpID] <- genename
#   }
# }
######
# looptypes <- function(statObj){
#   
#   typeObj <- data.frame()
#   
#   for(chr in paste0("chr", c(1:22, 'X', "Y"))){
#     statObj_chr <- statObj[statObj$chr==chr,]
#     partners <- looppartners(statObj_chr)
#     
#     statObj_chr$looptype <- paste(unique(statObj_chr[partners@from,])$type, unique(statObj_chr[partners@to,])$type, sep = '-')
#                                   
#     typeObj <- rbind(typeObj, unique(statObj_chr))
#   }
#   
#   return(typeObj)
# }

# 
# for(vpID in vpID_prot.list){
#   vpID_prot.list.tss <- vpID_prot.list[! grepl( "vpID*|CTCFinGene*", names(vpID_prot.list)) ]
#   gene.idx <- which(as.character(vpID_prot.list[[vpID]]$type) =="TSS")
#   
#   marks <-as.character(vpID_prot.list[[vpID]]$marks)[gene.idx]
#   marks <- paste(marks, collapse=";")
#   statObj.tss$geneIDvp[statObj.tss$id==vpID] <- marks
#   statObj.tss$geneIDanch[statObj.tss$anchor.vpID==vpID] <- marks
# 
# }



