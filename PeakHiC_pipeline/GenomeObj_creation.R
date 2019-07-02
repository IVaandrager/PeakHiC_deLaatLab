library("rtracklayer")
library("dplyr")
library("GenomicRanges")

genomeObj<- readRDS("/delaat/group/iwan/peakHiC/rds/genomeObj/hg19_WAPL_HiC_CTCFs.rds")
part.index <- resize(genomeObj$partition$partGR, width=1, fix = "center")

chrFrags <- readRDS("/delaat/group/iwan/peakHiC/rds/hg19_endonuclease_frags/hg19_MboI_frags_byChr.rds")
genemap <- readRDS("/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_ENSEMBL_Genes.rds")



###########creation of CTCF and SMC object
infile.dir <- "/delaat/group/iwan/peakHiC/chip_for_vp/CTCF_marks/"
rdsfiles <- list.files(path = infile.dir, pattern = "[.]bed.gz$",full.names = TRUE)
types <- list("CTCF_WT","CTCF_WAPL","CTCF_SCC4","CTCF_DKO","SMC1_WT","SMC1_WAPL","SMC1_SCC4","SMC1_DKO")
ctcf.files <- lapply(rdsfiles,import)
ctcf_withtype <- mapply(GR = ctcf.files, mark = types, function(GR,mark){ GR$type <- mark; GR }  )
ctcf.object <- do.call("c",ctcf_withtype)
colnames(mcols(ctcf.object))[1] <- "vpID"
ctcf.object<- ctcf.object[!(seqnames(ctcf.object)=="chrM")]
ctcf.object$vpID <- lapply(ctcf.object$vpID,
                           function(str){return(paste0(unlist(strsplit(str, split="_"))[c(3,4,11)], collapse="_"))})
#saveRDS(ctcf.object,"/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_CTCFs.rds")
ranges(ctcf.object) <- resize(ranges(ctcf.object), width=1, fix = 'center')
ctcf.object$partID <- paste0("part.",nearest(ctcf.object, part.index))



########add enhancers from file
enhancer.dir <- "/delaat/group/iwan/peakHiC/chip_for_vp/histone_marks/"
marks <- as.list(rep(c('H3K27ac', 'H3K4me3', 'H3K56ac', 'RNAPolII'), each = 2))

enhfiles <- list.files(path = enhancer.dir, pattern = "[.]cod$",full.names = TRUE)
hg19_files <- lapply(enhfiles, read.delim) 
hg19_GR <- lapply(hg19_files, GRanges)

#check ovelap between replicates; Later addition (after analysis1)
for(i in 1:length(unique(marks))){
  hg19_GR_ovl.idx <- unique( findOverlaps( ranges(hg19_GR[[i]]), ranges(hg19_GR[[i+1]]) )@from )
  hg19_GR[[i]] <- hg19_GR[[i]][hg19_GR_ovl.idx]
  i=i+2
}
hg19_GR <- hg19_GR[c(1,3,5,7)]

#add mark and concatenate
hg19_GR_ext <- mapply(GR = hg19_GR, mark = marks, function(GR,mark){ GR$type <- mark; GR }  )
hg19_enh <- do.call("c", hg19_GR_ext)

hg19_enh <- hg19_enh[-which(seqnames(hg19_enh)=='chrY'|seqnames(hg19_enh)=='chrM')] 
enh.vpID <- paste0(hg19_enh$type,1:length(hg19_enh$type))

#make hg19_enh genomeObj 
enh.genomeObj <- GRanges(seqnames = seqnames(hg19_enh), ranges = ranges(hg19_enh), strand = strand(hg19_enh),
                         vpID= enh.vpID, score= NA, type = hg19_enh$type)
#saveRDS(enh.genomeObj2,"/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_enhancers2.rds")

ranges(enh.genomeObj) <- resize(ranges(hg19_enh), width=1, fix = 'center')
enh.genomeObj$partID <- paste0("part.",nearest(hg19_enh, part.index))




#### reduce peaks and add fragIDs
all_vps <- c(ctcf.object,enh.genomeObj)
#saveRDS(all_vps, "/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/Hap1_CTCF_enhancers.rds")

#initialise objects
genomeObj.ranges <- GRanges()
site_tracker <- list(help = "This object contains a mapping from all original sites to the reduced VPs used to call chips. The hitsObj maps from origObj to final genomeObj per chr. In this way one can retrieve all sites associated with a VP, to characterise loops later on.")

#reduce and add fragID
for( chr in paste0("chr",c(1:22,"X"))){ #######chrY missing from chrFrags !!########

  vps_per_chr <- all_vps[seqnames(all_vps)==chr]

  #combine nearby viewpoints
  ranged_GO  <-IRanges(start(vps_per_chr)-2500, end(vps_per_chr) + 2500)
  ranged_GO2 <- IRanges::reduce(ranged_GO)

  hitsObj <- findOverlaps(ranges(vps_per_chr), ranges(ranged_GO2))
  site_tracker[chr] <- list(list(origSites = vps_per_chr, hitsObj=hitsObj))

  reduced_vps <- vps_per_chr[distanceToNearest(resize(ranged_GO2,width=1,fix="center"),ranges(vps_per_chr))@to]
  
  reduced_vps$type <- NULL
  reduced_vps$vpID <- NULL
  reduced_vps$score <- NULL
  
  #add fragIDs
  reduced_vps$fragID <-
    chrFrags[[chr]][ findOverlaps(ranges( chrFrags[[chr]] ), ranges(reduced_vps), type = "any" )@from ]$fragID
  
  genomeObj.ranges <- c(genomeObj.ranges,reduced_vps)
  }

#add vpID
genomeObj.ranges$vpID <- paste0("vpID_",1:length(genomeObj.ranges$fragID))

#attach to old genomeObj
genomeObj_backbone <- readRDS("/delaat/group/iwan/peakHiC/rds/genomeObj/old_TSS/hg19_WAPL_HiC_CTCF_TSS_enhancer_reduced.rds")
genomeObj_backbone$vpsGR <- genomeObj.ranges

genomeObj <- genomeObj_backbone



#######make vpID_to_prot Object for previous vps
#map VPs to proteins
vpID_to_prots <- function(vpID,newObj,hits,origObj){
  newidx <- which(newObj$vpID==vpID)
  oldidx <- which(hits@to == newidx)
  prots <- data.frame(marks=unlist(origObj[oldidx]$vpID), type=origObj[oldidx]$type)
  return(prots)
}

#map VPs to proteins
#site_tracker <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/genomeObj_sitetracker.rds")

vpID_to_prot <- list()
for(chr in paste0("chr",c(1:22,"X"))){
  newObj <- genomeObj$vpsGR[seqnames(genomeObj$vpsGR)==chr]
  origObj <- site_tracker[[chr]][['origSites']]
  hits <- site_tracker[[chr]][['hitsObj']]
  
  vpIDs <- as.list(newObj$vpID)
  
  protein.list <- lapply(vpIDs, vpID_to_prots, newObj, hits, origObj)
  
  names(protein.list) <- vpIDs
  
  vpID_to_prot <- c(vpID_to_prot, protein.list)
}


##################For TSS and geneCTCF different (no) reduce ################
###################add to genomeObj#########################
#FUNCTIONS
#add fragID
add_fragID <- function(TSS.genomeObj, chrFrags){
  TSS.whole <- GRanges()
  for( chr in names(chrFrags) ){
    Frags.chr <- chrFrags[[chr]]
    TSS.chr <- TSS.genomeObj[seqnames(TSS.genomeObj)==chr]
    
    idx <- findOverlaps(ranges(TSS.chr), ranges(Frags.chr))@to
    
    TSS.chr$fragID <- Frags.chr[idx]$fragID
    
    TSS.whole <- c(TSS.whole, TSS.chr)
  }
  return(TSS.whole)
}

#add all proteins in 5kb region of TSS vpIDs
get_prots <- function(row, all_vps){
  
  all_vps.chr <- all_vps[ seqnames(all_vps) == as.character(row[[1]]) ] 
  idx <- findOverlaps(IRanges(as.integer(row[[2]])-2500, as.integer(row[[2]])+2500), ranges(all_vps.chr) )
  
  marks2 <-  unlist( all_vps.chr[idx@to]$vpID )
  types <- unlist( all_vps.chr[idx@to]$type )
  
  df <- data.frame(marks=c(as.character(row[[6]]),marks2), type=c("TSS",types))
  return(df)}

#add all ctfc sites under, and the relevant gene to, the viewpoints of the type CTCF in Gene in the vpID_to_prot Object
geneCTCFs <- function(vpID, reduced_vps2, hitsObj2, ctcfsites, reduced_genes){
  old.index <- reduced_vps2[reduced_vps2$vpID==vpID]$index
  vpOrigins <-  which(hitsObj2@to==old.index)
  prots <- data.frame(marks=c(reduced_genes[old.index]$gene_name, unlist(ctcfsites[vpOrigins]$vpID) ), 
                      type= c("TSS", ctcfsites[vpOrigins]$type) )
  return(prots)
}



###########creation of proper TSS object
TSS <- readRDS("/delaat/group/iwan/peakHiC/rds/genomeObj/hg19_TSS_GM12878_single.rds")
TSS <- TSS[-which(seqnames(TSS)=='chrY'|seqnames(TSS)=='chrM')] 

TSS.genomeObj <- GRanges(seqnames = seqnames(TSS), ranges = ranges(TSS), strand = strand(TSS),
                         vpID= as.character(TSS$gene_symbol), type = "TSS")
TSS.genomeObj$partID <- paste0("part.",nearest(TSS, part.index))

TSS.genomeObj <- add_fragID(TSS.genomeObj, chrFrags)

#update vpID_to_prot list 
vpID_to_prot[TSS.genomeObj$vpID] <- apply(as.data.frame(TSS.genomeObj), 1, all_vps=all_vps, get_prots)
TSS.genomeObj$type <- NULL
TSS.genomeObj <- TSS.genomeObj[,c(2,3,1)]

#add TSS to genomeObj
#saveRDS(TSS.genomeObj,"/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/hg19_TSS.rds")
genomeObj$vpsGR <- c(genomeObj$vpsGR, TSS.genomeObj)



#ctcf.obj.chr[ctcf.obj.chr$vpID == 'CTCF_WT_30890']

###########creation of proper geneCTFC object
#do it per chr
#maybe use TSS object as it has part and fragID already
geneCTCF.obj <- GRanges()
partID.counter <- 1
for(chr in names(chrFrags)){
  
  genemap.chr <- genemap[seqnames(genemap)==chr]
  ctcf.obj.chr <- ctcf.object[seqnames(ctcf.object)==chr]
  
  #find ctcf sites within genes
  overlap <- findOverlaps(ranges(genemap.chr), ranges(ctcf.obj.chr))
  
  #the genes containing ctcf
  ctcfgenes <- genemap.chr[overlap@from]
  
  #the ctcf sites
  ctcfsites <- ctcf.obj.chr[overlap@to]
  
  wide  <-IRanges(start(ctcfsites)-2500, end(ctcfsites) + 2500)
  reduced <- IRanges::reduce(wide)
  
  hitsObj2 <- findOverlaps(ranges(ctcfsites), ranges(reduced))
  
  reduced_vps <- ctcfsites[distanceToNearest(resize(reduced,width=1,fix="center"),ranges(ctcfsites))@to]
  reduced_genes <- ctcfgenes[distanceToNearest(resize(reduced,width=1,fix="center"),ranges(ctcfgenes))@to]
  
  reduced_vps$index <- 1:length(reduced_vps)
  reduced_genes$index <- 1:length(reduced_genes)
  
  #take only the ctcfs that are 5000bp from TSS or end of gene
  inside.idx <- which(start(reduced_vps) > (start(reduced_genes)+5000) & end(reduced_vps) < (end(reduced_genes)-5000))
  reduced_vps2 <- reduced_vps[inside.idx]
  reduced_genes2 <- reduced_genes[inside.idx]
  
  #add fragIDs
  reduced_vps2$fragID <-
    chrFrags[[chr]][ findOverlaps(ranges( chrFrags[[chr]] ), ranges(reduced_vps2), type = "any" )@from ]$fragID
  
  #add vpID
  reduced_vps2$score <- NULL
  reduced_vps2$vpID <- paste0('CTCFinGene', partID.counter:(partID.counter+length(reduced_vps2$index)-1))
  partID.counter <- partID.counter+length(reduced_vps2$index)
  
  
  #add all proteins and relevant gene under reduced viewpoint (whats under the umbrella)
  vpIDs <- reduced_vps2$vpID
  protein.list <- lapply(vpIDs, geneCTCFs, reduced_vps2, hitsObj2, ctcfsites, reduced_genes)
  vpID_to_prot[vpIDs] <- protein.list
  
  reduced_vps2$index <- NULL
  reduced_vps2$type <- NULL
  reduced_vps2 <- reduced_vps2[,c(2,3,1)]
  
  geneCTCF.obj <-  c(geneCTCF.obj,reduced_vps2)
}

#saveRDS(geneCTCF.obj, "/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/Hap1_CTCFinGenes.rds")

genomeObj$vpsGR <- c(genomeObj$vpsGR, geneCTCF.obj)
# saveRDS(vpID_to_prot, "/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/vpID_to_marks.rds")


#Delete overlapping regions; give preference to TSS and CTCFinGene as they contain as much info as underlying regions
#genomeObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/Hap1_CTCF_TSS_enhancers_reduced.rds")
ObjCTCFs <- genomeObj$vpsGR[ grepl(pattern="vpID", genomeObj$vpsGR$vpID) ]
RangesCTCFs <- GRanges(seqnames=seqnames(ObjCTCFs), IRanges(start=start(ObjCTCFs)-5e3, end=end(ObjCTCFs)+5e3))

ObjCTCFinGene <- genomeObj$vpsGR[ grepl(pattern="CTCFinGene", genomeObj$vpsGR$vpID) ]
RangesCTCFinGene <- GRanges(seqnames=seqnames(ObjCTCFinGene), IRanges(start=start(ObjCTCFinGene)-5e3, end=end(ObjCTCFinGene)+5e3))

ObjGene <- genomeObj$vpsGR[ !grepl(pattern="CTCFinGene|vpID", genomeObj$vpsGR$vpID) ]
RangesGenes <- GRanges(seqnames=seqnames(ObjGene), IRanges(start=start(ObjGene)-5e3, end=end(ObjGene)+5e3))


CTCFalsoTSS <- unique( findOverlaps(RangesCTCFs, RangesGenes, minoverlap=2.5e3)@from )
CTCFalsoCTCFinGene <- unique( findOverlaps(RangesCTCFs, RangesCTCFinGene, minoverlap=2.5e3)@from )
redundant.idx <- union(CTCFalsoCTCFinGene,CTCFalsoTSS)

#save redundant vpIDs to adept old genomeObj to modernity
#saveRDS(ObjCTCFs[redundant.idx]$vpID, "/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/redudant_vpIDs")
ObjCTCFs <- ObjCTCFs[-redundant.idx]

genomeObj$vpsGR <- c(ObjCTCFs, ObjCTCFinGene, ObjGene)

#saveRDS(genomeObj, file = "/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/Hap1_CTCF_TSS_enhancers_reduced2.rds")






#function to change designMat
change.obj.design <- function(genomeObj){
  names(genomeObj$hic) <- "reps_as_reps"
  onerep.idx <- grep(genomeObj$hic$reps_as_reps$trackID, pattern = 'A$')
  genomeObj$hic$design <- genomeObj$hic$reps_as_reps[onerep.idx,]
  genomeObj$hic$design$trackID <-substr(genomeObj$hic$conds_as_reps$trackID,1,nchar(genomeObj$hic$conds_as_reps$trackID)-2)
  genomeObj$hic$design$sampleID <- paste0("sample.", 1:length(onerep.idx))
  return(genomeObj)
}










####creation of TSS object OLDER VERSIONS

#sort to remove non-unique fragIDs; then sort on chr position
# sorted_vps <- reduced_vps[order(reduced_vps$score, decreasing = TRUE),]
# unique_reduced_vps <- sorted_vps[!duplicated(sorted_vps$fragID)]
# unique_reduced_vps <-unique_reduced_vps[order(start(unique_reduced_vps)),]


#  hg19_ensembl <- readRDS("/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_ENSEMBL_Genes.rds")
#  TSS <- hg19_ensembl[levels(hg19_ensembl$type)=="start_codon"]
#  TSS <- TSS[-which(seqnames(TSS)=='chrY'|seqnames(TSS)=='chrM')] 

#  #TSS2 <- GRanges(inner_join( as.data.frame(TSS),
#                          data.frame(gene_symbol=genemap$gene_name,ensembl=genemap$ensembl, gene_type=genemap$gene_type)))  
#here we lose 1500 genes
# 
# saveRDS(TSS2,"/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_TSS_ensembl.rds")

#  #make TSS genomeObj 
#  TSS.genomeObj <- GRanges(seqnames = seqnames(TSS), ranges = ranges(TSS), strand = strand(TSS),
#                           vpID= as.character(TSS$ensembl), score= NA, type = "TSS")
# saveRDS(TSS.genomeObj,"/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_TSS_old.rds")
# 
# ranges(TSS.genomeObj) <- resize(ranges(TSS), width=1, fix = 'center')
# TSS.genomeObj$partID <- paste0("part.",nearest(TSS, part.index))


# newObj=genomeObj$vpsGR[seqnames(genomeObj$vpsGR)=='chr4']
# origObj=site_tracker$chr4$origSites
# hits=site_tracker$chr4$hitsObj
# newidx=which(newObj$vpID=='vpID_15333')
# oldidx=which(hits@to== newidx)
# prots=origObj[oldidx]



#genomeObj reduction rationale########
#For calling all peaks not so much a problem (30kb non-call range)
# sum(width(ranged_GO2)>=3e4) #chr='chr2'
# [1] 6
#For calling all marks in the area, more of a problem (300/5000=6%)
# > sum(width(ranged_GO2)>=1e4)
# [1] 291
#All reduced regions longer than 10kb, split into two. This has the upside that we later retrieve all proteins (in 10kb) of viewpoint
#so will give a more whole looping profile. The downside is that vp called in two seperate (minimally 5kb seperated) VPs, will give
# garantueed double loops (10kb res) to other sites. The question is whether this is biologically realistic (indeed 2 loops), or also
#how otherwise to choose the valid one and select relevant loop proteins/genes. I find it most likely these loops do occur, and I would
#want to retreive all binding partners to get valid loop characterisation. Solution: keep track which proteins went where, to also not
#overrepresent loops, but do correctly characterise them. (the splitting of regions larger than 20kb already is negligible (1%))


####add enhancers from second file; JUST USE BOTH REPLICATES, NO OVERLAP CHECK
# enhancer.dir <- "/delaat/group/iwan/peakHiC/chip_for_vp/histone_marks2/"
# marks2 <- as.list(rep(c('H3K27me3', 'H2AK119ub1'), each = 2))
# enhfiles2 <- list.files(path = enhancer.dir, pattern = "[.]cod$",full.names = TRUE)
# hg19_GR2 <- lapply(enhfiles2, import.bw)
# hg19_GR2 <- mapply(GR = hg19_GR2, mark = marks2, function(GR,mark){ GR$type <- mark; GR }  )
# hg19_enh2 <- do.call("c", hg19_GR2)
# 
# hg19_enh2 <- hg19_enh2[-which(seqnames(hg19_enh2)=='chrY'|seqnames(hg19_enh2)=='chrM')] 
# enh.vpID2 <- paste0(hg19_enh2$type,1:length(hg19_enh2$type))
# hg19_enh2$vpID <- enh.vpID2
# 
# #make hg19_enh2 genomeObj 
# saveRDS(hg19_enh2,"/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_enhancers2.rds")



#####combine Enhancers (Enh.genomeObj), TSS (TSS.genomeObj) and CTCF (ctcf.object) objects
# all_vps0 <- c(ctcf.object,TSS.genomeObj,enh.genomeObj)
# 
# all_vps <- GRangesList()
# for( chr in paste0("chr",c(1:22,"X"))){
#   vps_chr <- all_vps0[seqnames(all_vps0)==chr]
#   vps_chr$fragID <-
#   chrFrags[[chr]][ findOverlaps(ranges( chrFrags[[chr]] ), ranges(vps_chr), type = "any" )@from ]$fragID
#   all_vps[[chr]] <-  vps_chr
# }
# all_vps <- unlist(all_vps)
#
#saveRDS(all_vps, "/delaat/group/iwan/peakHiC/rds/genomeObj/hg19_WAPL_HiC_CTCF_TSS_enhancer.rds")

# reduce peaks and add fragIDs
# genomeObj.ranges <- GRanges()
# all_vps <- c(ctcf.object,TSS.genomeObj,enh.genomeObj)
# saveRDS(all_vps, "/delaat/group/iwan/peakHiC/rds/genomeObj/new_all/hg19_WAPL_HiC_CTCF_TSS_enhancer.rds")



######################junk########################
#genomeObj$partition$partGR <- c(genomeObj$partition$partGR,GRanges(seqnames="chrY",ranges=IRanges(seq(from = 0, to = 55e6, by = 2.5e6),seq(from = 5e6, to = 60e6, by = 2.5e6)),
#                                                                 strand="*",partID=paste0("part.",1174:1196)))
# RUN ONCE
# for (file in peak.files){
#   print(file)
#   system(paste0("zcat ",file, "| cut -f1-5 | gzip > ", outfile.dir,gsub("[.]gz",".bed.gz",basename(file)))) 
#   importfile <-import(paste0(outfile.dir, gsub("[.]gz",".bed.gz",basename(file)))) 
#   GR.peakfile <-resize(importfile, width = 1, fix = "center")
#   colnames(mcols(GR.peakfile))[1]<-"vpID"
#   GR.peakfile <- GR.peakfile[!(seqnames(GR.peakfile)=="chrM")]
#   GR.peakfile$type <- paste(strsplit(basename(file),"_")[[1]][4],strsplit(basename(file),"_")[[1]][5],sep="_")
#   GR.peakfile$partID <- paste0("part.",nearest(GR.peakfile, part.index))
#   saveRDS(GR.peakfile, file = paste0(outfile.dir,gsub("[.]gz",".rds",basename(file))))
# }

######combining seperate chip peak files into a single file for the genomeObj##########
# rdsfiles <- list.files(path = outfile.dir, pattern = "[.]rds$",full.names = TRUE)
# readrds <- lapply(rdsfiles,readRDS)
#CTCF.genomeObj <- do.call("c", readrds)


#Genes <- readRDS('/delaat/group/iwan/peakHiC/rds/genomeObj/hg19_genesGR.rds')

# all_types <- function(count, i, j, genomeObj){
#   genomeObj$vpsGR$type[j] <- paste(unique(genomeObj$vpsGR$type[seq(i,i+count-1)]), collapse=".")
#   i=i+count
#   j=j+1
#   return(genomeObj)
# }
# genomeObj <- lapply(as.list(count_repeats), all_types, i=1, j=1, genomeObj=big.genomeObj)







