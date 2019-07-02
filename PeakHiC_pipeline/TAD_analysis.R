source('peakHiC_functions_getloopCov.R')
tracks <- as.list(1:9)
partReadsFldr <- "/delaat/group/iwan/Hap1_peakHiC/analysis1/partition_reads/"

TADs <- read.table("/delaat/group/iwan/Hap1_peakHiC/objects/TADs/H1-ESC_Dixon2015-TADs.txt", sep=',', header = T)
TADs <- GRanges(TADs)
TADs$X <- NULL
TADs$TADid <- paste0('TAD.',1:length(start(TADs)))
saveRDS(TADs, "/delaat/group/iwan/Hap1_peakHiC/objects/TADs/TADs.rds")

genomeObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/Hap1_CTCF_TSS_enhancers_reduced2.rds")
part.index <- resize(genomeObj$partition$partGR, width=1, fix = "center")
TADs$partID<- paste0("part.",nearest(TADs, part.index))
TADs$vp1 <- start(TADs)

#n+0
TADs$vp2 <- end(TADs)

#n+1
TADs1 <- TADs
TADs1.vp2 <- data.frame(vp=TADs$vp2)
TADs1.vp2 <- transform(TADs1.vp2, vp = c(vp[-1], NA))
TADs1$vp2 <- as.numeric(unlist(TADs1.vp2))
TADs1 <- TADs1[-length(TADs1$TADid)]

#n+2
TADs2 <- TADs
TADs2.vp2 <- data.frame(vp=TADs$vp2)
TADs2.vp2 <- transform(TADs2.vp2, vp = c(vp[-c(1,2)], c(NA,NA)))
TADs2$vp2 <- as.numeric(unlist(TADs2.vp2))
TADs2 <- TADs2[-((length(TADs2$TADid)-1):length(TADs2$TADid))]

getTADCovbyPartition <- function(partID,TADs,genomeObj, tracks, partReadsFldr){
  partLoops=TADs[TADs$partID==partID]

  if( length(partLoops) >0) {
    xPos <- partLoops$vp1
    yPos <- partLoops$vp2
    lChr <- seqnames(partLoops)
    lx <- resize(GRanges(lChr,IRanges(xPos,xPos)),width=10e3,fix="center")
    ly <- resize(GRanges(lChr,IRanges(yPos,yPos)),width=10e3,fix="center")
    
    out <- as.data.frame(matrix(0,nrow=length(lx),ncol=(length(tracks)+1)))
    names(out)[2:10] <- c("DKO_3.3-A","DKO_3.3-B","Hap1-A","Hap1-B","WaplKO_3.3-A","WaplKO_3.3-B", "SCC4KO-A", "SCC4KO-B","WaplKO_1.14-A")
    out$TADid <- partLoops$TADid
    
    PRreads.all <- readRDS(list.files(paste0(partReadsFldr,partID,"/"), full.names = T))
    
    for(trackID in tracks) {
      PRreads <- PRreads.all[[trackID]]
      loopTags <- tagRegions(PE1=PRreads$PE1,PE2=PRreads$PE2,qR1=lx,qR2=ly)
      out[[trackID+1]] <- loopTags
    }
    return(out)
  }
}

partIDs <- unique(TADs$partID)

TADcov <- lapply(as.list(partIDs), getTADCovbyPartition,TADs=TADs[-1],genomeObj=genomeObj, tracks=tracks, partReadsFldr=partReadsFldr)
TADcov<-rbindlist(TADcov)

TADcov1 <- lapply(as.list(partIDs), getTADCovbyPartition,TADs=TADs1[-1],genomeObj=genomeObj, tracks=tracks, partReadsFldr=partReadsFldr)
TADcov1<-rbindlist(TADcov1)

TADcov2 <- lapply(as.list(partIDs), getTADCovbyPartition,TADs=TADs2[-1],genomeObj=genomeObj, tracks=tracks, partReadsFldr=partReadsFldr)
TADcov2<-rbindlist(TADcov2)

TADcovList <- list("n+0"=TADcov, "n+1"=TADcov1, 'n+2'=TADcov2)

saveRDS(TADcovList, "/delaat/group/iwan/Hap1_peakHiC/objects/TADs/TADcov2.rds")
#TADcovList <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/TADs/TADcov2.rds")

#normalise and add counts
normaliseCov <- function(loopCov){
  median <- median(colSums(loopCov[-1]))
  normfact <- colSums(loopCov[-1])/median
  
  loopCov[-1] <- loopCov[-1] / normfact[col(loopCov[-1])]
  
  return(loopCov)
}

addCov <- function(loopCov, conds){
  loopCov <- as.data.frame(loopCov)
  sumloopCov=data.frame(temp=rep(0,length(loopCov[,1])))
  
  for(cond in conds){
    cond.ID <-  grep(names(loopCov), pattern=cond, value = T)
    
    if(length(names(loopCov[,cond.ID]))>0){sumloopCov$temp <- rowSums(loopCov[,cond.ID])} else{
      sumloopCov$temp <- loopCov[,cond.ID]}
    
    names(sumloopCov)[names(sumloopCov) == "temp"] <- cond
  }
  
  sumloopCov <- cbind(TADid = loopCov$TADid, sumloopCov)
  
  return(sumloopCov)
}

conds = c("DKO", "Hap1", "WaplKO_3.3", "SCC4KO", "WaplKO_1.14") 


TADcov.sum <- lapply(TADcovList, addCov, conds=conds)
TADcov.norm <- lapply(TADcov.sum, normaliseCov)
  
saveRDS(TADcov.norm, "/delaat/group/iwan/Hap1_peakHiC/objects/TADs/TADcov_norm2.rds")



#####boxplot for coverage dependency on neighbour########
#we dont find relation, maybe many outerloops exceeding 1MB limit
#something wrong with normalization?

TADcovList <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/TADs/TADcov_norm2.rds")
#add factor for neighbour level
for(i in 1:3){
  TADcovList[[i]]$neighbour <- as.factor(i-1)
  i=i+1
}
TADcov <- rbindlist(TADcovList)
TADcov$logWaplKO_WT <- log2((TADcov$WaplKO_3.3+0.001)/(TADcov$Hap1+0.001))

#boxplot for all neighbours the expression category vs wapl/hap coverage of TADs
boxplot_TADexpr <- ggplot(data = TADcov, aes(x=neighbour, y=logWaplKO_WT)) + #SCC4KO
  geom_boxplot() +
  theme_pubr(border = T) +
  coord_cartesian(ylim=c(-2,2)) 
boxplot_TADexpr 

  
  
