#need to have:
#infolder/...
#loopsObj "loops.rds"
#MacsObj "Hap1_CTCF.bw"
#hg19 genes "hg19_genesGR.rds"
#Diff_expression "/Diff_expression/*WaplKO|SCC4|DKO*"
#type_annotation "/type_annotation/*"
#vpReads partID obj paste0("V4Cs/",vpReads_",partID,".rds")
#
#find corresponding partID for vpID to find which file to put
#loops <- readRDS('/delaat/group/iwan/Hap1_peakHiC/ruiqi/newer/CTCFinGene_PolII1kb_loops_withCHIP.rds')
#partID <- loops[loops$id==vpID,'partID'][[1]]
#loops <- readRDS('/delaat/group/iwan/Hap1_peakHiC/ruiqi/newer/manyLoops.rds')


#input
source("/delaat/group/iwan/peakHiC/scripts/plot_V4C_vpID_functions.R")
inFolder <- "/delaat/group/iwan/Hap1_peakHiC/ruiqi/Objects/"
outFolder <- "/delaat/group/iwan/Hap1_peakHiC/plots/"
loops <-  readRDS('/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/FullStatObj.rds')
genomeObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/Hap1_CTCF_TSS_enhancers_reduced2.rds")
#all loops; not filtered on RNApolII at TSS
#loops2 <-  readRDS('/delaat/group/iwan/Hap1_peakHiC/ruiqi/newer/manyLoops.rds')
#loops2[ loops2$id == vpID | loops2$anchor.vpID == vpID, ]

#####Single input#####
vpID <- 'MEGF11'#'BCL11B' #'FOXO1' #'vpID_43172' #'BCL11B' #loopsWaplKO[2,'id'] #'FAM69B' #'COL5A1' 'SOX6'  vpID_16983 'ANO3''EBF1'
mutant <- 'WaplKO'     #"DKO", WaplKO", "SCC4KO", "All"

pdf( paste0(outFolder, mutant, "_", vpID, ".pdf"), width = 14, height = 7 )
plot_V4C(vpID,mutant,inFolder,loops,genomeObj)
dev.off()


#####multiple input#####
vpIDs <- c( 'DNAH14', 'SIPA1', 'GJB2', 'CNIH2', 'LEFTY1')  #loopsWaplKO[1:10,'id'] # c('SOX6', 'CTCFinGene27')  
vpIDs <- genes; #bivalent_genes_down_in_WaplKO
mutant <- 'WaplKO'     #"DKO", WaplKO", "SCC4KO", "All"

pdf(paste0(outFolder,"downDE_prom_enh_Loops.pdf"))  #'name.pdf' paste0(outFolder, mutant, "_", paste(vpIDs, collapse = "_"))
for(vpID in vpIDs){ print( plot_V4C(vpID,mutant,inFolder,loops, genomeObj) ) }
dev.off()


### some analysis ###
######grand analysis######
#combine DE of vp and anch
loops2$DKO_DE <- as.numeric(dplyr::coalesce( loops2$vp.DE.DKO_vs_Hap1,  loops2$anch.DE.DKO_vs_Hap1))
loops2$WaplKO_DE <- as.numeric(dplyr::coalesce( loops2$vp.DE.WaplKO_vs_Hap1,  loops2$anch.DE.WaplKO_vs_Hap1))
loops2$SCC4KO_DE <- as.numeric(dplyr::coalesce( loops2$vp.DE.SCC4_vs_Hap1,  loops2$anch.DE.SCC4_vs_Hap1))
names(loops2)[[4]] <- "WaplKO"

## see if there is differential expression in either VP or anchor. If gene expression is higher wrt wildtype if then the peak strength is lower and vice versa
#H3K56ac=T; if there is no H3K56ac indicative of cancer
#Dist=T; if not loop tiny so that we know its a proper loop
#ordering='ctcf' then if ordered by CTCF strength; ordering='expression' most change in gene expression
grand_filter <- function(loops, mutant, H3K56ac=F, Dist=F, ordering='ctcf'){
  index <- grep(paste0(mutant,"_DE"), names(loops))
  filtered <- loops[!is.na(loops[,index]),]
  filtered <- filtered[ (( (filtered[,mutant] > (filtered[,'Hap1'] +5)) & (filtered[,index]  < 0) ) | 
                           ( (filtered[,'Hap1'] > (filtered[,mutant]+5)) & (filtered[,index]  > 0) )), ]
  
  if(H3K56ac){filtered=filtered[filtered$vp.H3K56ac==0,]} #to avoid cancer effects
  if(Dist){filtered=filtered[filtered$loopDistCat!='tiny',]} #to ensure looping effect
  
  if(ordering=='ctcf'){ #strongest ctcf 
  filtered$Hap_CTCF.avg<-rowMeans(filtered[,c("vp.Hap1_CTCF_chip","anch.Hap1_CTCF_chip" )], na.rm=T)
  filtered=filtered[order(filtered$Hap_CTCF.avg,decreasing=T),] } 
  
  if(ordering=='expression'){ filtered=loops_Wapl[order(abs(filtered[,index]),decreasing=T),] }  #most changes expression
  
  return(filtered)
}

ordering <- 'ctcf' #'expression'
loopsWaplKO <- grand_filter(loops=loops2, mutant='WaplKO', H3K56ac=T, Dist=T, ordering=ordering)
loopsDKO <- grand_filter(loops=loops2, mutant='DKO', H3K56ac=T, Dist=T, ordering=ordering)
loopsSCC4KO <- grand_filter(loops=loops2, mutant='SCC4KO', H3K56ac=T, Dist=T, ordering=ordering)


###possible interesting sites
#loop 7066;  CTCFinGene5509.   Gene downregulation upon formation of loop in WaplKO - ctcfinGene loop lowers gene expression; downside, its far away
#loop 9291; CTCFinGene7116. Same as above
#loop 16224; CTCFinGene12063; Not the found anchor, possibly in the middle of gene ADK where we see upregulation upon removal of loop in WaplKO;
#downside unkown if recip loop, !!loopID 22748!! looper to anch = CTCFinGene16354 to vp=CTCFinGene16351. Interesting: gene bending on itself
#loopID 22748; see above
#loopID 31040; downregulated upon increased looping in WaplKO; two times a strong CTCF
#loopID 7662; clean loop; upregulated upon more looping.. wrong way around (several of this 'reversed' type)
#loopID 31856 22427  6640; nice and clean CTCFinGene to CTCFinGene loop

#check further in objects
#find info of a gene
genesGR <- readRDS(paste0(inFolder, "hg19_genesGR.rds"))

#find CTCFs within a gene and if loops exist at these sites
ADKgene <- GRanges(seqnames='chr10', IRanges(75910943, 76469061))
idx <- findOverlaps(GRanges(seqnames='chr10', IRanges(76.2e6, 76.25e6)), GRanges(seqnames='chr10', IRanges(loops2$vp_X1, loops2$vp_X1)))
loops2[idx@to,]



