subChr <- function(gR,chr) {

idx <- as.vector(seqnames(gR))==chr
return(gR[idx])

}

extractReads <- function( fragID, vpS, vpZoom=c(1e6,1e6), reads, frags, k=31, vpSize=31 ){
	
	fragsChr <- frags
	# first we search for the fragment containing the viewpoint
	vpIdx <- match(fragID,fragsChr$fragID)

	# then we select all the fragments around the viewopint in the range provided
	vpF <- fragsChr[ queryHits( findOverlaps( ranges( fragsChr ), IRanges( start=vpS-vpZoom[1], end=vpS+vpZoom[2] ) ) ) ]

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

normV4C <- function( readsGR, nReads=10e3, nTop=1 ) {
	readsGR$normReads <- 0
	sumTop <- sum( -sort( -readsGR$reads )[ 1:nTop] )
	wNorm <- nReads/( sum( readsGR$reads )-sumTop )
	readsGR$normReads <- wNorm*readsGR$reads
	return(readsGR)
}

doPeakCCall <- function(vpID,genomeObj,wSize=21,alphaFDR=0.1,qWr=2.0,minDist=20e3) {

rdsFldr <- "/delaat/group/geert/peakHiC/ANALYSES/hg19_WAPL_peakHiC/rds/"

hicCond <- "WAPL_HiC"

vpsGR <- genomeObj[["vpsGR"]]
partID <- vpsGR$partID[match(vpID,vpsGR$vpID)]
vpPos <- start(vpsGR[match(vpID,vpsGR$vpID)])

fRDS <- paste0(rdsFldr,"vpReads_",partID,".rds")

if(file.exists(fRDS)) {

	vpReads <- readRDS(fRDS)
	nextVPReads <- vpReads[[vpID]][[hicCond]]

	num.exp <- length(nextVPReads)

	db <- list()
	db$num.exp <- num.exp

	for(i in 1:num.exp) {
		
		db[[i]] <- data.frame(pos=nextVPReads[[i]]$GR$pos,reads=nextVPReads[[i]]$GR$reads)
	
	}

	resPeakC <- combined.analysis(data=db,num.exp=3,vp.pos=vpPos,wSize=wSize,alphaFDR=alphaFDR,qWr=qWr,minDist=minDist)
	return(resPeakC)

} else {

}

}


combine.reps <- function (data) {

	num.exp <- length(data)

    data.m <- data[[1]]
    
    for (i in 2:num.exp) {
        data.m <- merge(data.m, data[[i]], by = 1)
    }
    
    sumReads <- apply(data.m[,2:(num.exp+1)],1,sum)
    out <- data.frame(pos=data.m[,1],reads=sumReads)

    return(out)
}

peakAnalysis <- function (data, num.exp = 3, vp.pos, wSize = 21, alphaFDR = 0.1,qWr = 1, minDist = 25e3) {
    
    if (num.exp == 0) {
        num.exp = data$num.exp
    }
    if (length(vp.pos) == 1) {
        vp.pos <- c(vp.pos, vp.pos)
    }
    vp.pos <- sort(vp.pos)
    db <- tryCatch({ combine.experiments(data, num.exp, vp.pos) }, error=function(e) { return (NULL) } )
    if ( !is.null( db ) ){
      dbR <- db
      dbR[, 2:(num.exp + 1)] <- apply(db[, 2:(num.exp + 1)], 2, 
          caTools::runmean, k = wSize, endrule = "mean")
      dbR[, 2:(num.exp + 1) + num.exp] <- apply(db[, 2:(num.exp + 
          1) + num.exp], 2, caTools::runmean, k = 5, endrule = "mean")
      pseudoCount <- apply(db[, 2:(num.exp + 1)], 2, non.zero.quantile, 
          probs = 0.05)
      # print(pseudoCount)
      pseudoCount <- sum(pseudoCount)/num.exp
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
          1)], 1, mean), frags = ratio[idx, 1], wSize = wSize, 
          qW = qWr)
      sfr <- intersect(sfr, tfr)
      list(dbR = dbR, peak = sfr, num.exp = num.exp, p.value = p.val, 
          ratio = apply(ratio[, 2:(num.exp + 1)], 1, mean), delta = apply(delta[, 2:(num.exp + 1)], 1, mean), sel = sel.frag)
    }
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

righttailgamma = function(r,k,n) 1 - pgamma(-log(r/(n+1)^k),k,scale=1)

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

thresholdFrags <- function(resids,frags,wSize=21,qW=5) {

	qMax <- getThreshold(resids=resids,qW=qW)

	if(is.na(qMax)) {

		return(NULL)

	} 

	else{
		
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
}

non.zero.quantile <- function( x, probs ){
  quantile(x[x > 0], probs)
}

rem <- function(a, n ){
  half.window <- floor(n/2)
  head(tail(a, -half.window),-half.window)
}

#quick way of generating a vector with the required indexes
multi.seq <- function( start, end ){
  x <- rep(start, end-start+1)->x
  df <- diff(x)
  df <- df + 1
  low <- which(df > 1)
  df[low] <- -diff(c(0,low))+1
  add <- c(0,cumsum(df))
  x + add
}

getThreshold <- function(resids,qW=5) {

	q75 <- quantile(na.omit(resids),probs=0.75) #75% quantile of the residuals
	qd50 <- diff(quantile(na.omit(resids),probs=c(0.25,0.75))) #the range between the 25% and 75% quantiles
	threshold <- q75 + qW*qd50
	return(threshold)
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

getPartitionPeaks <- function(partID, genomeObj, hicCond="WAPL_HiC", wSize=41, qWr=1, alphaFDR=0.1,  minDist=30e3 , nReps=NULL) {

	if(is.null(nReps)) {

		nReps <- sum(as.vector(genomeObj[["hic"]][["design"]][["condID"]])==hicCond)

	}

	rdsFldr <- paste0("/delaat/group/geert/peakHiC/ANALYSES/mESC_OCT4_GFP/mm10_lncRNA_peakHiC/rds/")
	profilesFldr <- paste0("/delaat/group/geert/peakHiC/ANALYSES/mESC_OCT4_GFP/mm10_lncRNA_peakHiC/rds/profiles/")
	loopsFldr <- paste0("/delaat/group/geert/peakHiC/ANALYSES/mESC_OCT4_GFP/mm10_lncRNA_peakHiC/loops/")

	vpsGR <- genomeObj[["vpsGR"]]
	partVPs <- vpsGR[vpsGR$partID==partID]

	loopFile <- paste0(loopsFldr,"GW_nReps_",nReps,"_peakHiC_wSize_",wSize,"_qWr_",qWr,"_alphaFDR_",alphaFDR,"_loops.txt")	
	profilesFile <- paste0(profilesFldr,partID,".rds")

	if(file.exists(profilesFile)) {

		V4Cs <- readRDS(profilesFile)
	
	} 
	
	else {
	
			V4Cs <- list()
	
	}
	
	V4Cs[[hicCond]] <- list()


	if(length(partVPs)>0) {
		
		fRDS <- paste0(rdsFldr,"vpReads_",partID,".rds")

		if(file.exists(fRDS)) {
			
			ids <- partVPs$vpID
			vpReads <- readRDS(fRDS)

			for(id in ids) {
				
				peakRes <- getPeakCPeaksWithReps(vpID=id,vpReads=vpReads, partVPs=partVPs, hicCond=hicCond, wSize=wSize, qWr=qWr, alphaFDR=alphaFDR,  minDist=minDist)
			
				if(length(peakRes$df)>0) {

					write.table(peakRes$df,file=loopFile,append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
				}

				V4Cs[[hicCond]][[id]] <- peakRes$gR

			}

			saveRDS(V4Cs,file=profilesFile)

		}
	}

}

getChrPeaks <- function(chr, genomeObj, hicCond="WAPL_HiC", wSize=41, qWr=1, alphaFDR=0.1,  minDist=30e3 , nReps=NULL, nThreads=8) {
	
	suppressPackageStartupMessages(require(doParallel))
	registerDoParallel(cores=nThreads)

	vpsGR <- genomeObj[["vpsGR"]]
	ids <- unique(subChr(vpsGR,chr)$partID)

	foreach(i = 1:length(ids)) %dopar% {
	
		getPartitionPeaks(partID=ids[i], genomeObj=genomeObj, hicCond=hicCond, wSize=wSize, qWr=qWr, alphaFDR=alphaFDR,  minDist=minDist, nReps=nReps)
	
	}

}

getPeakCPeaksWithReps <- function( vpID, vpReads, partVPs, hicCond="WAPL_HiC", wSize=41, qWr=1, alphaFDR=0.1,  minDist=30e3 ) {

	vpPos <- start(partVPs[match(vpID,partVPs$vpID)])
	vpChr <- as.vector(seqnames(partVPs[match(vpID,partVPs$vpID)]))
	
	nextVPReads <- vpReads[[vpID]][[1]]

	num.exp <- length(nextVPReads)

  	peakCReads <- list()
		
	for(i in 1:num.exp) {
		
		dat <- data.frame(pos=nextVPReads[[i]]$GR$pos,reads=nextVPReads[[i]]$GR$reads)
		colnames(dat) <- c("pos",paste0("reads.R",i))
		peakCReads[[i]] <- dat
	}
	
	peakCReads$num.exp <- num.exp
	nReps <- num.exp

	peakCRes <- peakAnalysis( data=peakCReads, num.exp=nReps, vp.pos=vpPos, wSize=wSize, minDist=minDist, qWr=qWr, alphaFDR=alphaFDR )

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

doPeakCPlot <- function(vpID,genomeObj,wSize=21,alphaFDR=0.1,qWr=2.0,minDist=20e3) {

rdsFldr <- "/delaat/group/geert/peakHiC/ANALYSES/hg19_WAPL_peakHiC/rds/"

hicCond <- "WAPL_HiC"

vpsGR <- genomeObj[["vpsGR"]]
partID <- vpsGR$partID[match(vpID,vpsGR$vpID)]
vpPos <- start(vpsGR[match(vpID,vpsGR$vpID)])

fRDS <- paste0(rdsFldr,"vpReads_",partID,".rds")

if(file.exists(fRDS)) {

	vpReads <- readRDS(fRDS)
	nextVPReads <- vpReads[[vpID]][[hicCond]]

	num.exp <- length(nextVPReads)

	db <- list()
	db$num.exp <- num.exp

	for(i in 1:num.exp) {
		
		db[[i]] <- data.frame(pos=nextVPReads[[i]]$GR$pos,reads=nextVPReads[[i]]$GR$reads)
	
	}

	resPeakC <- combined.analysis(data=db,num.exp=3,vp.pos=vpPos,wSize=wSize,alphaFDR=alphaFDR,qWr=qWr,minDist=minDist)

} else {

}

}


getPeakHiCData <- function(partID,frags,genomeObj,hicCond="mESC_OCT4_GFP",wSize=21,vpSize=31,viewSize=1e6) {
	
	designMat <- genomeObj[["hic"]][["design"]]
	hicTracksByCondition <- split(as.vector(designMat$trackID),designMat$condID)
	
	vps <- getVPs(partID,genomeObj)
	vpReads <- list()

	if(length(vps$vpID)>0) {
		
		reads <- list()
		tracks <- hicTracksByCondition[[hicCond]]

		for(trackID in tracks) {

			reads[[trackID]] <- getPartitionReads(partID=partID,trackID=trackID,genomeObj=genomeObj)

		}
		
		for(vpID in vps$vpID) {
			
			vpReads[[vpID]] <- list()
			vpReads[[vpID]][[hicCond]] <- list()
			fragID <- vps$fragID[match(vpID,vps$vpID)]
			vpPos <- start(vps[match(vpID,vps$vpID)])

			for(trackID in tracks) {

				vpReads[[vpID]][[hicCond]][[trackID]] <- extractReads(fragID=fragID,vpS=vpPos,vpZoom=c(viewSize,viewSize),reads=reads[[trackID]],frags=frags,k=wSize,vpSize=vpSize)

			}
		}
	}

	return(vpReads)

}

getPartitionReads <- function(partID,trackID,genomeObj,baseFldr=NULL) {

require(data.table)
require(GenomicRanges)

baseFldr <- "/delaat/group/geert/peakHiC/READS/mESC_OCT4_GFP/"
pairixBinary <- "/home/geert/localdev/prog/pairix/bin/pairix"

out <- NULL

partIdx <- match(partID,genomeObj$partition$partGR$partID)

if(!is.na(partIdx)){

	qPartition <- genomeObj$partition$partGR[partIdx]
	chr <- as.vector(seqnames(qPartition[1]))
	qRegion <- paste0("\'",chr,":",start(qPartition),"-",end(qPartition),"\'")

	tmpFldr <- tempdir()
	tmpFile <- paste0(tmpFldr,"/pairix_query.txt")
	pairixFile <- paste0(baseFldr,trackID,"/",chr,"_pairix.gz")
	
	if (file.exists(pairixFile)&file.exists(pairixBinary)) {

		cmd <- paste0(pairixBinary," ",pairixFile," ",qRegion," > ",tmpFile)
		system(cmd)
		dat <- fread(file=tmpFile,header=FALSE,sep="\t",stringsAsFactors=FALSE)
		PE1 <- GRanges(seqnames=dat$V2,IRanges(dat$V3,dat$V3))
		PE2 <- GRanges(seqnames=dat$V4,IRanges(dat$V5,dat$V5))
		out <- list(PE1=PE1,PE2=PE2)

	}

}

return(out)

}


getPEReads <- function(partID,trackID,genomeObj){

	if(is.null(genomeObj[["partition"]][["partGIntervals2D"]][[partID]])) {

		out <- NULL

	} else {

		vpInt2d <- genomeObj[["partition"]][["partGIntervals2D"]][[partID]]
		
		points <- gextract(trackID, vpInt2d )
		PE1 <- GRanges(points$chrom1,IRanges(points$start1,points$start1))
		PE2 <- GRanges(points$chrom2,IRanges(points$start2,points$start2))
		out <- list(PE1=PE1,PE2=PE2)
	}

	return(out)

}

getVPs <- function(partID,genomeObj) {

	if(is.null(genomeObj[["vpsGR"]])) {

		out <- NULL

	} else {

		vpsGR <- genomeObj[["vpsGR"]]
		out <- vpsGR[vpsGR$partID==partID]

	}

	return(out)
}

plotChromHMM <- function(gR,plotRanges,yrange=c(0.85,0.9),xdiv=1e6,...) {

ovlIdx <- findOverlaps(gR,plotRanges)@from

if(length(ovlIdx)) {

feat <- as.data.frame(gR[ovlIdx])

rect(feat[["start"]]/xdiv,yrange[1],feat[["end"]]/xdiv,yrange[2],border=NA,col=feat$color,...)


}

}
