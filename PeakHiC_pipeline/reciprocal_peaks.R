library(dplyr)
library(data.table)
library(GenomicRanges)

#filter peaks to retain only the reciprocal ones
check_reciprocal <- function(inFolder){
  
  infile <- list.files(inFolder, full.names = T, include.dirs = F)
  in.idx <- grep(infile, pattern = "*loops.txt")
  loops_orig <- read.table(infile[in.idx], sep = '\t')
  
  names(loops_orig) <- c('chr','vp_X1', 'vp_X2', 'anchor_X1', 'anchor_X2', 'anchorSize', 'maxV4Cscore','maxV4CscorePos', 'delta', 'ratio', 'minPval','id')
  
  #create reciprocal loops object
  reciprocals <- setNames(data.table(matrix(ncol = 13, nrow = 0)), 
                          c('chr','vp_X1', 'vp_X2', 'anchor_X1', 'anchor_X2', 'anchorSize', 'maxV4Cscore','maxV4CscorePos', 'delta', 'ratio', 'minPval','id', 'anchor.vpID'))
  
  
  #old neccesity where overlapping viewpoints werent deleted and TSS or CTCFinGene were given preference
  redundant.vpIDs <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/redudant_vpIDs")
  loops_orig <- loops_orig[ !(loops_orig$id %in% redundant.vpIDs), ]
  
  #loop per chromosome
  chrom <- paste0("chr", c(1:22, 'X'))
  for(chr in chrom){
    
    loops_chr <- filter(loops_orig, chr==!!chr) #bangbang
    
    #find if the anchor also occurs as viewpoint(s), and what these viewpoints are
    #check anchor1 exists as viewpoint2
    anchor_also_VP <- findOverlaps(IRanges(start = loops_chr$maxV4CscorePos-5e3,loops_chr$maxV4CscorePos+5e3),
                                 IRanges(start = loops_chr$vp_X1, end = loops_chr$vp_X2))
    
    #check if the anchor, beloning to the viewpoint (viewpoint which originally was the anchor), falls within the original viewpoint
    #check if viewpoint1 is within anchor2
    VP_also_anchor <- findOverlaps( IRanges(start = loops_chr$vp_X1, end = loops_chr$vp_X2),
                                     IRanges(start = loops_chr$maxV4CscorePos-5e3,loops_chr$maxV4CscorePos+5e3))
    
    #take the peaks that call each other
    recips <- anchor_also_VP[which(anchor_also_VP %in% VP_also_anchor)]
   
    #some loops are called twice (eg chr X index 42 anchor exists as 138 and 139, and both loops can be validated.) However, upon selection by
    #recips@from this means a duplicate of loop 42, which we dont want cause it is still a single loop. The double will come back in that it
    #does allow both 138 and 139 to be valid loops, but will not allow loop 42 to double exist. For the mapping from anchor to vp for protein
    #retrieval in the statistics part we also only want to keep the mappings to the unique indexes.
    for( i in 1:length(recips)){
      recips <- recips[!recips@from==recips[i]@to]
      
      if(i>=length(recips)){
        break }
    }
    
    recips <- recips[ !duplicated(recips@from) ]
    #if loop 133 is reciprocal to 142; then all 142 are deleted; However, if 142 loops 643 and 563 then those recips (643:142 and 563:142) arent deleted (chr5)
    recips <- recips[ !recips@from>recips@to]

    recipLoops <- loops_chr[recips@from,]
    #get the vpID that corresponds to the anchor; for prot retrieval later (5kb region so fine for counts)
    recipLoops$anchor.vpID <- loops_chr[recips@to,'id']
    
    reciprocals <-  rbind(reciprocals, recipLoops)
  
  }
  
  reciprocals$loopID = 1:length(reciprocals$id)
  
  #sum(lengths(looppartners[-1])) == length(reciprocals$loopID)
  saveRDS(reciprocals, file=gsub("loops[.]txt","loops_reduced.rds", infile[in.idx]) )
  #saveRDS(looppartners, file = gsub("loops[.]txt","loop_partners.rds", infile[in.idx]) )
}






# #remove loops that are called twice, due to reciprocal requirement, as not all were. Some are literal duplicates due to VPs, statObj[c(45780,45781,45829,45830),]
# #Some were called double and others werent due to unique() step in reciprocal over grouped anchor; see reciprocal_peaks for thorough explanation
# IDobj <- FullStatObj[,18:19]
# names(IDobj) <- NULL
# rownames(IDobj) <- NULL
# IDobj.list <- split(IDobj, seq(nrow(IDobj)))
# 
# doubled <- lapply(IDobj.list, FUN=function(x){which(FullStatObj[,19]==x[[1]] & FullStatObj[,18]==x[[2]])})
# doubled <- doubled[lapply(doubled,length)>0] 
# 
# doubled.df <- data.frame(row = as.integer( rep(names(doubled), sapply(doubled, length))),
#                          recip = unlist(doubled))
# 
# sorted.doubled <- t(apply(doubled.df, 1, sort))
# duplicates <- which(duplicated(sorted.doubled))
# 
# tripled <- unlist( doubled[lapply(doubled,length)>1] )
# duplicates.idx <- c(unlist(tripled),doubled.df[duplicates,]$recip)
# 
# #filter away duplicate loops for which id == vpID and vpID == id
# FullStatObj_nodouble <- FullStatObj[-duplicates.idx,]
#
#
# ######alternative syntax leading to same result#######
# reciprocals2 <- setNames(data.table(matrix(ncol = 12, nrow = 0)),
#                          c('chr','vp_X1', 'vp_X2', 'anchor_X1', 'anchor_X2', 'anchorSize', 'maxV4Cscore','maxV4CscorePos', 'delta', 'ratio', 'minPval','id'))
# 
# #take original anchor location and the anchor as a viewpoint
# #take anchor1 and viewpoint2
# orig_anchor <- loops_chr[anchor_also_VP@from,]
# vpanchor_as_vp <- loops_chr[anchor_also_VP@to,]
# 
# overlaps <- orig_anchor$vp_X1 >= vpanchor_as_vp$maxV4CscorePos-5e3 & orig_anchor$vp_X1 <= vpanchor_as_vp$maxV4CscorePos+5e3
# 
# 
# looppartners[[chr]] <- anchor_also_VP[overlaps]
# 
# 
# reciprocals2 <- rbind(reciprocals2, orig_anchor[overlaps,] )
# 
# ############################################################












