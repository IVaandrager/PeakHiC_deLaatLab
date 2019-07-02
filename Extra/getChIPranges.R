library(rtracklayer)

getChIPAtGRanges <- function(bwFile,gR,lab="ChIP.tags",q=NULL) {
  
  ids <- gR$ID
  
  if(is.null(ids)) {
    
    ids <- paste0("region_",1:length(gR))
    
  }
  
  gR$score <- 0
  
  bwSel <- BigWigSelection(ranges=gR,colnames="score")
  
  dat <- import(con=bwFile,selection=bwSel)
  olap <- findOverlaps(gR,dat)
  
  if(is.null(q)){
    
    meanScores <- tapply(dat$score[olap@to],olap@from,mean)
    
  } else {
    
    meanScores <- tapply(dat$score[olap@to],olap@from,quantile,probs=q)
    
  }
  gR$score[as.numeric(names(meanScores))] <- meanScores
  
  return(gR)
  
}