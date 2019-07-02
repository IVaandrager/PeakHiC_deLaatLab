#functions 

cpet <- function(i,Sq1,Sq2){return(length(intersect(Sq1[[i]],Sq2[[i]])))}


getVPs <- function(partID,genomeObj) {
  
  if(is.null(genomeObj[["vpsGR"]])) {
    
    out <- NULL
    
  } else {
    
    vpsGR <- genomeObj[["vpsGR"]]
    out <- vpsGR[vpsGR$partID==partID]
    
  }
  
  return(out)
}


tagRegions <- function(PE1,PE2,qR1,qR2) {

olapPE1q1 <- findOverlaps(qR1,PE1)
olapPE2q2 <- findOverlaps(qR2,PE2)

olapPE2q1 <- findOverlaps(qR1,PE2)
olapPE1q2 <- findOverlaps(qR2,PE1)

Sq1 <- split(olapPE1q1@to,factor(olapPE1q1@from,1:length(qR1)))
Sq2 <- split(olapPE2q2@to,factor(olapPE2q2@from,1:length(qR2)))

tags <- sapply(1:length(qR1),cpet,Sq1=Sq1,Sq2=Sq2)

Sq1 <- split(olapPE1q2@to,factor(olapPE1q2@from,1:length(qR1)))
Sq2 <- split(olapPE2q1@to,factor(olapPE2q1@from,1:length(qR2)))

tags <- tags+sapply(1:length(qR1),cpet,Sq1=Sq1,Sq2=Sq2)

return(tags)

}

getLoopCovbyPartition <- function(partID,loops,genomeObj,anchorSize=10e3, tracks, partReadsFldr){
 
   message( paste0( '> analyzing partition ..' , partID) )
  
	if(is.null(loops$loopID)|is.null(loops$id)){stop(paste0("loop data.frame must have loopID (identifier of loop) and id (identifier of viewpoint) columns"))}

	vps <- getVPs(partID,genomeObj)$vpID
	partLoops <- loops[loops$id%in%vps,]
	designMat <- genomeObj$hic$design[unlist(tracks),]
	
	if(nrow(partLoops)>0) {
		
		lVP <- floor((partLoops$vp_X1+partLoops$vp_X2)/2)
		lAnchor <- partLoops$maxV4CscorePos
		xPos <- pmin(lVP,lAnchor)
		yPos <- pmax(lVP,lAnchor)
		lChr <- partLoops$chr
		lx <- resize(GRanges(lChr,IRanges(xPos,xPos)),width=anchorSize,fix="center")
		ly <- resize(GRanges(lChr,IRanges(yPos,yPos)),width=anchorSize,fix="center")

		out <- as.data.frame(matrix(0,nrow=length(lx),ncol=(length(tracks)+1)))
		colnames(out) <- c("loopID", designMat$trackID) #sapply(FUN = tail, strsplit(tracks, split = '/'), n=1)
		out$loopID <- partLoops$loopID
    
		PRreads.all <- readRDS(list.files(paste0(partReadsFldr,partID,"/"), full.names = T))
		
		for(trackID in tracks) {
			PRreads <- PRreads.all[[trackID]]
			loopTags <- tagRegions(PE1=PRreads$PE1,PE2=PRreads$PE2,qR1=lx,qR2=ly)
			out[[trackID+1]] <- loopTags

		}
	}

	return(out)

}



