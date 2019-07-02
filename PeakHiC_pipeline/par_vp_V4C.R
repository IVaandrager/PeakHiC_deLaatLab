#########functions#################
#extract reads per partition per vp per rep
extractReads <- function( fragID, vpS, vpZoom=c(1e6,1e6), reads, frags, k=31, vpSize=31 ){
  
  fragsChr <- frags
  # first we search for the fragment containing the viewpoint
  vpIdx <- match(fragID,fragsChr$fragID)
  
  # then we select all the fragments around the viewpoint in the range provided
  vpF <- fragsChr[ queryHits( findOverlaps( ranges( fragsChr ), IRanges( start=vpS-vpZoom[[1]], end=vpS+vpZoom[[2]] ) ) ) ]
  
  # fragments pooled as viewpoint
  metaFrags <- (vpIdx-floor(vpSize/2)):(vpIdx+floor(vpSize/2))
  metaFrags <- metaFrags[ which( metaFrags > 0 ) ]
  
  if ( length( metaFrags ) != vpSize ){
    # we check if we are at the start or the end of a chromosome
    if ( max(metaFrags) < vpSize ) {
      metaFrags <- seq( 1, vpSize, 1 )
    }
  }
  if ( max(metaFrags) > fragsChr$fragID[length(fragsChr)] ){
    metaFrags <- seq( fragsChr$fragID[length(fragsChr)]-vpSize+1, fragsChr$fragID[length(fragsChr)], 1 )
  } 
  
  # take the viewpoint extremes for computing vp coverage
  vpGR <- reduce(fragsChr[ metaFrags ])
  minVPpos <- min( start(vpGR) ) - 100000
  maxVPpos <- max( end( vpGR )) + 100000
  
  if ( minVPpos < 1 ) {
    minVPpos <- 1
  }
  
  if ( maxVPpos > max( end( fragsChr ) ) ){
    maxVPpos <- max( end( fragsChr ) )
  }
  
  # selection of the reads that share the viewpoint fragments
  idxPE1.c1 <- queryHits( findOverlaps(reads$PE1,fragsChr[metaFrags]) )
  idxPE2.c1 <- queryHits( findOverlaps(reads$PE2,fragsChr[metaFrags]) )
  
  # the position of the fragment is set to be in the center.
  vpF$pos <- start( resize( vpF, width=1, fix="center") )
  
  # now we sum up all the reads found in each position in range that share one of the ends with the viewpoint fragment
  if ( ( length( vpF ) > 0 ) & ( length( idxPE1.c1 ) > 0 ) & ( length( idxPE2.c1 ) > 0 ) ) {
    vpF$reads <- countOverlaps( vpF,reads$PE2[ idxPE1.c1 ] ) + countOverlaps( vpF,reads$PE1[ idxPE2.c1 ] )
    if ( length( which( vpF$reads > 0 ) ) > 2 ) {
      vpF <- normV4C( vpF )
      vpF$normV4C <- rollmean( x=vpF$normReads, k=k, fill="extend" )
    } else {
      vpF$normReads <- 0
      vpF$normV4C <- 0
    }
  } else {
    vpF$reads <- 0
    vpF$normReads <- 0
    vpF$normV4C <- 0
  }
  
  vpCovGR <- vpF[ which( start( vpF ) >= minVPpos & end( vpF ) <= maxVPpos ) ]
  vp_cov <- 100*( length( which( vpCovGR$reads > 0 ) ) / length( vpCovGR ) )
  
  nReads <- length( which( vpF$reads > 0 ) )
  
  if ( nReads != 0 ){
    plotCov <- 100*( nReads / length( vpF ) )
  } else {
    plotCov <- 0
  }
  
  return( list( GR=vpF, vpGR=vpGR, vpCov=vp_cov, plotCov=plotCov ) )
}

#get vps from genomeObj per parition
getVPs <- function(partID,genomeObj) {
  
  if(is.null(genomeObj[["vpsGR"]])) {
    
    out <- NULL
    
  } else {
    
    vpsGR <- genomeObj[["vpsGR"]]
    out <- vpsGR[vpsGR$partID==partID]
    
  }
  
  return(out)
}

#READ rds FILE PARTITION READS
getPeakHiCData <- function(partID,frags,genomeObj,hicCond, partReadFldr, tracks, wSize=21,vpSize=31,viewSize=1e6) {

  vps <- getVPs(partID,genomeObj)
  #vps <- 'BCL11B'
  vpReads <- list()

  part.reads <- list.files(paste0(partReadFldr,partID), full.names = T)
  
  
  if(file.exists(part.reads)){
    reads <- readRDS(part.reads)
    #reads <- readRDS("/delaat/group/iwan/Hap1_peakHiC/analysis1/partition_reads/part.896/track_1:2:3:4:5:6:7:8:9.rds")
    
    if(length(reads)>0){
      for(vpID in vps$vpID) {
          vpReads[[vpID]] <- list()
          vpReads[[vpID]][[hicCond]] <- list()
          fragID <- vps$fragID[match(vpID,vps$vpID)]
          vpPos <- start(vps[match(vpID,vps$vpID)])
          
          #extract the actual reads for the V4C
          for(trackID in tracks) {
            vpReads[[vpID]][[hicCond]][[trackID]] <- extractReads(fragID=fragID,vpS=vpPos,vpZoom=c(viewSize,viewSize),reads=reads[[trackID]],frags=frags,k=wSize,vpSize=vpSize)
            }
      }
    }else{vpReads <- list()} } else{vpReads <- list() } #if the reads file doesnt exist or is empty; return an empty list
  
  return(vpReads)
}

normV4C <- function( readsGR, nReads=10e3, nTop=1 ) {
  readsGR$normReads <- 0
  sumTop <- sum( -sort( -readsGR$reads )[ 1:nTop] )
  wNorm <- nReads/( sum( readsGR$reads )-sumTop )
  readsGR$normReads <- wNorm*readsGR$reads
  return(readsGR)
}

vp_reads <- function(partID, genomeObj, frags, hicCond, V4CFldr, partReadFldr, tracks){
  message(partID)
  
  if(!file.exists(paste0(V4CFldr,"vpReads_",partID,".rds"))){
    if(length(genomeObj$vpsGR[genomeObj$vpsGR$partID==partID])>0){
      tryCatch({
        chr <- as.character(unique(seqnames(genomeObj$vpsGR[genomeObj$vpsGR$partID==partID])))
        frags <- frags[[chr]]
        vpReads <- getPeakHiCData(partID,frags,genomeObj,hicCond, partReadFldr, tracks)
      
      if(length(vpReads>0)){
        fRDS <- paste0(V4CFldr,"vpReads_",partID,".rds")
        saveRDS(vpReads,file=fRDS)}
      
      }, error=function(e) { message(paste0(partID," failed")) })
    }
  }
}



####function calls####
c(library("parallel"), library("GenomicRanges"), library("data.table"), library("zoo"), library('peakC'))

VP_V4Cs <- function(partReadsFldr, outFolder, genomeObj, hicCond, tracks, frags){

  V4CFldr <- outFolder
  ids <- as.list(genomeObj$partition$partGR$partID)

  cl = makeCluster(6, outfile="")
  clusterEvalQ(cl, c(library("GenomicRanges"), library("data.table"), library("zoo"), library('peakC')))
  clusterExport(cl, c('extractReads', 'getVPs', 'getPeakHiCData', 'normV4C'))

  parLapply(cl, ids, fun = vp_reads, genomeObj=genomeObj, frags=frags, hicCond=hicCond, V4CFldr=V4CFldr, partReadFldr=partReadsFldr, tracks=tracks)

  stopCluster(cl)

  #non-parallel test
  #lapply(ids, FUN = vp_reads, genomeObj=genomeObj, frags=frags, hicCond=hicCond, V4CFldr=V4CFldr, partReadFldr=partReadsFldr, tracks=tracks)

}


