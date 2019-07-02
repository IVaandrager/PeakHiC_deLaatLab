library(data.table)
library(GenomicRanges)
library(parallel)

#### get PE valid Hi-C pairs (READS) from a specific Hi-C "track" (cond/rep)
getPEReads <- function(partID,trackID,genomeObj, readsFldr){
  
  designMat <- genomeObj$hic$design
  partGR <- genomeObj$partition$partGR

  IDX.READS <- match(trackID,designMat$trackID)
  IDX.PART <- match(partID,partGR$partID)

  if(is.na(IDX.READS)|is.na(IDX.PART)) {
    
    out <- NULL
    
  } else {
    
    gR <- partGR[IDX.PART]
    chr <- as.vector(seqnames(gR))
    fRDS <- paste0(readsFldr,designMat$condID[IDX.READS],"/",trackID,"/",chr,".rds")

    if(file.exists(fRDS)){
      
      out <- readRDS(fRDS)
      PE1.ovl <- findOverlaps(out$PE1,gR)
      PE2.ovl <- findOverlaps(out$PE2,gR)
      IDX.OVL <- intersect(PE1.ovl@from,PE2.ovl@from)
      out <- list(PE1=out$PE1[IDX.OVL],PE2=out$PE2[IDX.OVL])
      
    } else{
      
      out <- NULL
      
    }
  }
  
  return(out)
  
}

#  [1]     chr1 [122500001, 127500000]      * |     part.50 gR


### Merge Hi-C pairs (READS) from multiple reps / conditions
reads_per_tracks <- function(partID, trackIDs, genomeObj, readsFldr){
  tracks <- unique(genomeObj$hic$design$trackID)[trackIDs]
  trackReads <- list()

  for(i in 1:length(tracks)){
    trackReads[[i]] <- getPEReads(partID=partID,trackID=tracks[i],genomeObj=genomeObj,readsFldr= readsFldr)
  }
  names(trackReads) <- tracks

  reads <- trackReads[[1]] #firs rep/cond for init
  if(length(trackReads)>1){
    for(i in 2:length(tracks)) { #combine all reps/conds
      reads$PE1 <- c(reads$PE1,trackReads[[i]]$PE1)
      reads$PE2 <- c(reads$PE2,trackReads[[i]]$PE2)
    }
  }
  
  return(reads)
}

readsloop <- function(partID, trackIDs, genomeObj, part.dirs, readsFldr){
  part.dir <- paste0(part.dirs,partID)
  tracknomer <- paste0("/track_",paste(lapply(trackIDs, FUN = paste, collapse = ""), collapse = ":"),".rds")
  
  if(!file.exists(paste0(part.dir,tracknomer))){
    dir.create(part.dir)
    reads <- lapply(trackIDs, reads_per_tracks, partID=partID, genomeObj=genomeObj, readsFldr= readsFldr)
    saveRDS(reads, file = paste0(part.dir, tracknomer))
    print(partID)
}}

##########function call###############
partition_reads <- function(inFolder, outFolder, genomeObj, tracks){
  
  vpsGR <- genomeObj$vpsGR
  partIDs <- unique(vpsGR$partID)
  
  readsFldr <- inFolder       #"/delaat/group/iwan/peakHiC/READS/"
  part.dirs <- outFolder      #"/delaat/group/iwan/peakHiC/rds/partition_reads/test/"
  trackIDs <- tracks       #as.list(1:9)  #list(c(1,2),c(3,4),c(5,6),c(7,8),9) 
  # [1] "DKO_3.3-A"     "DKO_3.3-B"     "Hap1-A"        "Hap1-B"       
  # [5] "WaplKO_3.3-A"  "WaplKO_3.3-B"  "SCC4KO-A"      "SCC4KO-B"     
  # [9] "WaplKO_1.14-A"
  
  
  cl = makeCluster(6, outfile="")
  clusterEvalQ(cl, c(library("GenomicRanges"), library("data.table")))
  clusterExport(cl, c('getPEReads', 'reads_per_tracks'))
  
  parLapply(cl, partIDs, fun = readsloop, trackIDs = trackIDs, genomeObj = genomeObj, part.dirs=part.dirs, readsFldr=readsFldr)
  #not every partition has a viewpoint
  
  on.exit(stopCluster(cl))

}


