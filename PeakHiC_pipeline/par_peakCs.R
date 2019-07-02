library("GenomicRanges")
library("peakC")
library("zoo")

righttailgamma = function(r,k,n) 1 - pgamma(-log(r/(n+1)^k),k,scale=1)
getThreshold <- function(resids,qW=5) {
  q75 <- quantile(resids,probs=0.75) #75% quantile of the residuals
  qd50 <- diff(quantile(resids,probs=c(0.25,0.75))) #the range between the 25% and 75% quantiles
  threshold <- q75 + qW*qd50
  return(threshold)
}
multi.seq <- function( start, end ){
  x <- rep(start, end-start+1)->x
  df <- diff(x)
  df <- df + 1
  low <- which(df > 1)
  df[low] <- -diff(c(0,low))+1
  add <- c(0,cumsum(df))
  x + add
}

non.zero.quantile <- function( x, probs ){
  quantile(x[x > 0], probs)
}
rank.product.p <- function( data, num.exp,method="diff"){
  if(method=="diff") {
    stats <- data[,2:(num.exp+1)]-data[,(2:(num.exp+1))+num.exp]
  } else {
    stats <- data[,2:(num.exp+1)]/data[,(2:(num.exp+1))+num.exp]
  }
  rp <- nrow(data)-apply(stats,2,rank)+1
  rp <- apply(rp,1,prod)
  p <- righttailgamma(rp,num.exp,length(rp))
}
significant.fragments <- function( p.value, pos, window = 21, FDR = 0.01 ){
  #correct the nominal p-value for multiple hypothesis testing
  p.combined <- p.adjust(p.value, method="fdr")
  #determine the significant windows and select the fragments therein
  sig.i <- which(p.combined < FDR)
  if(length(sig.i)>0) {
    sig.i.start <- sig.i-floor(window/2); sig.i.end <- sig.i+floor(window/2)
    sig.i <- unique(multi.seq(sig.i.start,sig.i.end))
    sig.i <- sig.i[sig.i >= 1 & sig.i <= length(pos)]
    sigFrags <- pos[sig.i]
    return(sigFrags)
  } else {
    return(NULL)
  }
}
thresholdFrags <- function(resids,frags,wSize=21,qW=5) {
  
  qMax <- getThreshold(resids=resids,qW=qW)
  
  sel.i <- which(resids > qMax)
  if(length(sel.i)>0) {
    sel.i.start <- sel.i-floor(wSize/2); sel.i.end <- sel.i+floor(wSize/2)
    sel.i <- unique(multi.seq(sel.i.start,sel.i.end))
    sel.i <- sel.i[sel.i >= 1 & sel.i <= length(frags)]
    selFrags <- frags[sel.i]
    return(selFrags)
  }else{
    return(NULL)
  }
  
  
}
collapseFrags <- function(frags, peakFrags) {
  
  N <- length(frags)
  
  if(length(peakFrags)>0){
    
    fragRanges <- IRanges(frags,frags)
    end(fragRanges)[1:(N-1)] <- start(fragRanges)[2:N]-1
    
    binaryFrags <- ifelse(frags%in%peakFrags,1,0)
    
    return(reduce(split(fragRanges,binaryFrags)[["1"]],min.gapwidth=1))
    
  } else {
    
    return(NULL)
  }
  
}

peakAnalysis <- function (data, vp.pos, wSize, alphaFDR, qWr, minDist) {
  
  num.exp <- data$num.exp

  if (length(vp.pos) == 1) {
    vp.pos <- c(vp.pos, vp.pos)
  }
  
  vp.pos <- sort(vp.pos)
  
  #combine conditions into one dataframe, for most reliable peak calling; use this as loopset under which to count coverage later
  db <- tryCatch({ combine.experiments(data, num.exp, vp.pos) }, error=function(e) { return (NULL) } )
  
  #combined.analysis()
  if ( !is.null( db ) ){
    dbR <- db
    dbR[, 2:(num.exp + 1)] <- apply(db[, 2:(num.exp + 1)], 2, 
                                    caTools::runmean, k = wSize, endrule = "mean")
    dbR[, 2:(num.exp + 1) + num.exp] <- apply(db[, 2:(num.exp + 
                                                        1) + num.exp], 2, caTools::runmean, k = 5, endrule = "mean")
    pseudoCount <- apply(db[, 2:(num.exp + 1)], 2, non.zero.quantile, 
                         probs = 0.05)
    if(sum(is.na(pseudoCount))==num.exp){message(sprintf("too few reads for vp.pos %.0f", vp.pos[[1]] )); return() } else{
    pseudoCount <- sum(pseudoCount, na.rm = T)/sum(!is.na(pseudoCount)) }
    ratio <- cbind(db[, 1], (dbR[, 2:(num.exp + 1)] + pseudoCount)/(dbR[, 
                                                                        (2:(num.exp + 1)) + num.exp] + pseudoCount))
    delta <- cbind(db[, 1], dbR[, 2:(num.exp + 1)] - dbR[, (2:(num.exp + 
                                                                 1)) + num.exp])
    p.val <- rank.product.p(data = dbR, num.exp = num.exp, method = "diff")
    sfr <- significant.fragments(p.value = p.val, pos = db[, 
                                                           1], window = wSize, FDR = alphaFDR)
    sel.frag <- db[which((db[, 1] < vp.pos[1] & vp.pos[1] - db[, 
                                                               1] > minDist) | (db[, 1] > vp.pos[2] & db[, 1] - vp.pos[2] > 
                                                                                  minDist)), 1]
    idx <- delta[, 1] %in% sel.frag
    tfr <- thresholdFrags(resids = apply(ratio[idx, 2:(num.exp + 
                                                         1)], 1, mean, na.rm = T), frags = ratio[idx, 1], wSize = wSize, ####na.rm = T
                          qW = qWr)
    sfr <- intersect(sfr, tfr)
    list(dbR = dbR, peak = sfr, num.exp = num.exp, p.value = p.val, 
         ratio = apply(ratio[, 2:(num.exp + 1)], 1, mean), delta = apply(delta[, 2:(num.exp + 1)], 1, mean), sel = sel.frag)
  }
}

getPeakCPeaksWithReps <- function( vpID, vpReads, vpInfo, wSize=41, qWr=1, alphaFDR=0.1,  minDist=30e3, nReps) {
  
  vpPos <- start(vpInfo[match(vpID,vpInfo$vpID)])
  vpChr <- as.character(seqnames(vpInfo)[1])
  
  num.exp <- nReps
  
  peakCReads <- list()
  
  #turn into data frame
  for(i in 1:num.exp) {
    dat <- data.frame(pos=vpReads[[i]]$GR$pos,reads=vpReads[[i]]$GR$reads)
    colnames(dat) <- c("pos",paste0("reads.R",i))
    peakCReads[[i]] <- dat
  }
  
  peakCReads$num.exp <- num.exp
  
  #call peaks
  peakCRes <- peakAnalysis( data=peakCReads, vp.pos=vpPos, wSize=wSize, minDist=minDist, qWr=qWr, alphaFDR=alphaFDR )
  
  #give specific form, necessary for loopCov script later
  if( length( peakCRes$peak ) > 0 ) {
    
    peakGR <- sort(collapseFrags(frags=peakCRes$dbR[,1],peakFrags=peakCRes$peak))
    y.ave <- apply(peakCRes$dbR[, 2:(nReps + 1)], 1, median)
    gR <- IRanges(peakCRes$dbR$pos,peakCRes$dbR$pos)
    finalGR <- GRanges( unique( vpChr ), gR, normV4C=y.ave )
    ovl <- findOverlaps(gR,peakGR)
    maxV4C <- as.vector(tapply(y.ave[ovl@from],ovl@to,max))
    minPval <- as.vector(tapply(peakCRes$p.value[ovl@from],ovl@to,min))
    maxRatio <- as.vector(tapply(peakCRes$ratio[ovl@from],ovl@to,max))
    maxDelta <- as.vector(tapply(peakCRes$delta[ovl@from],ovl@to,max))
    ttt <- data.frame( gR, y.ave )[queryHits(ovl),]
    ttt$loop <- subjectHits( ovl )
    ttt <- split( ttt, ttt$loop )
    maxV4CscorePos <- sapply( 1:length( ttt ), function(x) ttt[[x]]$start[ which.max( ttt[[x]]$y.ave ) ] )
    finalDF <- data.frame( 
      chr=vpChr
      , vp_X1=vpPos
      , vp_X2=vpPos
      , anchor_X1=start( peakGR )
      , anchor_X2=end( peakGR )
      , anchorSize=width( peakGR )
      , maxV4Cscore=maxV4C
      , maxV4CscorePos=maxV4CscorePos
      , delta=maxDelta
      , ratio=maxRatio
      , minPval=minPval
      , id=vpID
      , stringsAsFactors=FALSE
    )
  } else {
    
    peakGR <- NULL
    maxV4C <- NULL
    minPval <- NULL
    maxRatio <- NULL
    maxDelta <- NULL
    finalGR <- NULL
    maxV4CscorePos <- NULL
    finalDF <- NULL
    
  }  
  
  return( list( gR=finalGR, df=finalDF ) )
  
}

getPartitionPeaks <- function(partID, genomeObj, wSize=41, qWr=1, alphaFDR=0.1,  minDist=30e3 , nReps=NULL, folders) {
  
  message(partID)
  
  rdsFldr <- folders$rdsFldr
  loopsFldr <- folders$loopsFldr
  loopFile <- paste0(loopsFldr,"GW_nReps_",nReps,"_peakHiC_wSize_",wSize,"_qWr_",qWr,"_alphaFDR_",alphaFDR,"_loops.txt")	
    
  fRDS <- paste0(rdsFldr,"vpReads_",partID,".rds")
  vpReads <- readRDS(fRDS) 
  
  vpInfo <- genomeObj$vpsGR[genomeObj$vpsGR$partID==partID]
  
  ids <- names(vpReads)
  
  V4Cs <- GRangesList()
  
  #call per vpID
  for(id in ids) {
    id_vpReads <- vpReads[[id]][[1]]
    peakRes <- getPeakCPeaksWithReps(vpID=id, vpReads=id_vpReads, vpInfo=vpInfo, wSize=wSize, qWr=qWr, alphaFDR=alphaFDR,
                                   minDist=minDist, nReps=nReps)
    
    #save all the loops in a loops file
    if(length(peakRes$df)>0) {
          
      write.table(peakRes$df,file=loopFile,append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
      
      V4Cs[[id]] <- peakRes$gR
    }
    
  }
  
  saveRDS(V4Cs, paste0(loopsFldr,"normV4Cs/",partID,"_normV4Cs.rds" ))
}


############variables############
###this script takes all the virtual V4C and combines the data from all replicate conditions to most reliably, with most cover
###calls peaks. These peak locations are then later, per condition, checked for coverage
call_peaks <- function(inFolder, outFolder, genomeObj){
  
  folders <- list(rdsFldr = inFolder, loopsFldr = outFolder)
  
  wSize=41
  qWr=1
  alphaFDR=0.1
  minDist=30e3
  nReps <- length(readRDS(list.files(folders$rdsFldr, full.names = T)[1])[[1]][[1]])
  
  #list partitions on the basis of V4C files, saved per partition
  parts <- lapply(list.files(folders$rdsFldr), function(x){sub("^.*_", "", gsub("[.]rds","",x))})
  
  #call per partition
  #############call functions###########
  library(parallel)
  nThreads=6
  cl = makeCluster(nThreads, outfile="")
  clusterEvalQ(cl, c(library("GenomicRanges"), library("data.table"), library("zoo"), library('peakC')))
  clusterExport(cl, c('collapseFrags','getPartitionPeaks','getPeakCPeaksWithReps','getThreshold','multi.seq','non.zero.quantile','peakAnalysis',
                      'rank.product.p','righttailgamma','significant.fragments','thresholdFrags'))


  parLapply(cl, parts, fun = getPartitionPeaks, genomeObj=genomeObj, wSize=wSize, qWr=qWr,
            alphaFDR=alphaFDR, folders=folders, minDist=minDist, nReps=nReps)

  on.exit(stopCluster(cl))
  
  
  ##### non-parallel test#########
#   for(i in 1:length(parts)){
#     print(parts[i])
#     getPartitionPeaks(partID=parts[i], genomeObj=genomeObj, wSize=wSize, qWr=qWr, alphaFDR=alphaFDR,  minDist=minDist, nReps=nReps, folders=folders)
#   }


  }


# call_peaks(inFolder, outFolder, genomeObj)
# inFolder <-  "/delaat/group/iwan/peakHiC/rds/V4C_VPs_Hap1/test/"
# outFolder <- "/delaat/group/iwan/peakHiC/rds/V4C_peaks/test/"
# call_peaks(inFolder,outFolder)


