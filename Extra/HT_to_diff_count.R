library("rtracklayer")
library("GenomicRanges")
library("dplyr")  
library("reshape2")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")

#systemctl status rstudio-server.service 

infile.dir <- "/delaat/group/iwan/peakHiC/HTseq_for_diff_expr"
outfile.dir <- "/delaat/group/iwan/peakHiC/HTseq_for_diff_expr/rds/"

#transfrom txt file to rds file
HT.files <- list.files(path = infile.dir, full.names=TRUE, pattern = "GSM")
for (file in HT.files){
  HT.peakfile <- read.table(gzfile(file))
  colnames(HT.peakfile)<-c("ensID","counts")
  saveRDS(HT.peakfile, file = paste0(outfile.dir,gsub("[.]txt[.]gz",".rds",basename(file))))
}

#read data
rdsfiles <- list.files(path = outfile.dir, pattern = "[.]rds$",full.names = TRUE)
readrds <- lapply(rdsfiles,readRDS)
big_HT <- data.frame(readrds)
all_HT <- big_HT[,grepl("(counts)",colnames(big_HT))]
colnames(all_HT) <- paste(rep(c("Hap1","WaplKO","SCC4","DKO","WT2","Wapl2"), each=3), 1:3, sep = "_")
rownames(all_HT) <- big_HT[,1]

#add coldata
coldata <- data.frame(condition = rep(c("Hap1","WaplKO","SCC4","DKO","WT2","Wapl2"), each=3))
rownames(coldata) <- colnames(all_HT)

#########
all_HT2 <- data.frame(Hap1=rowMeans(all_HT[,c(1:3, 13:14)]), WaplKO_3.3=rowMeans(all_HT[,4:6]), SCC4KO=rowMeans(all_HT[,7:9]),
                      DKO= rowMeans(all_HT[,10:12]), WaplKO_1.14=rowMeans(all_HT[,15:17]))

# #filter out low reads
keep <- rowSums(all_HT2) >= 10
all_HT2 <- all_HT2[keep,]

cds <- estimateSizeFactorsForMatrix(all_HT2)
normalizedCounts <- t(t(all_HT2)/cds)

saveRDS(normalizedCounts, "/delaat/group/iwan/Hap1_peakHiC/objects/Expression_levels/normHTseq/normHTseq.rds")

#########
#create TPM ]:  RPKM = ( C/(10^-6 * colsums) )/L ; TPM = (C/L)/(colsum*1e-6)
genemap <- readRDS("/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_ENSEMBL_Genes.rds")
#transcript ID vs gene ID
genemap <- left_join( data.frame(ensembl=rownames(all_HT2)) , as.data.frame(genemap))

GeneLength <- data.frame(ensembl=genemap$ensembl, length=(genemap$width*1e-3))
ScalingFactor <- colSums(all_HT2/GeneLength$length, na.rm = T)*1e-6

TPM <- t(t(( all_HT2/GeneLength$length )) / as.numeric(ScalingFactor))

saveRDS(TPM, "/delaat/group/iwan/Hap1_peakHiC/objects/Expression_levels/normHTseq/TPMcounts.rds")
#statObj.TSS2 contains the new TPM scores



#create DESeq data object
dds <- DESeqDataSetFromMatrix(countData=all_HT,
                       colData = coldata,
                       design = ~ condition
)

#filter out low reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#indicate control group
dds$condition <- relevel(dds$condition, ref = "Hap1")

#run DESeq
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)

resnames <- resultsNames(dds)
############per condition comparison visualisation below#####################
#resnames[[i]]
# [1] "Intercept"                "condition_DKO_vs_Hap1"   
# [3] "condition_SCC4_vs_Hap1"   "condition_Wapl2_vs_Hap1" 
# [5] "condition_WaplKO_vs_Hap1" "condition_WT2_vs_Hap1" 


#perform shrinkage on LogFoldChange values for visualisation  [[5]] for waplKO_vs_Hap
resLFC <- lfcShrink(dds, coef=resnames[[5]], type="apeglm") 
summary(resLFC)
sum(resLFC$padj < 0.05, na.rm=TRUE)
resLFC.ordered <- resLFC[order(resLFC$pvalue),]

#select rows with an adjusted p-value of 0.05 or lower
sig.resLFC <- resLFC[resLFC$padj< 0.05 & !is.na(resLFC$padj),]


plotMA(resLFC, ylim=c(-2,2))
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]


#independent hypothesis weighting
#install.packages('IHW') 
resIHW <- results(dds, coef=resnames[[5]], filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)

#heatmap for quality assessment
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[4:43]
df <- as.data.frame(colData(dds)[,"condition"])
rownames(df) <- colnames(dds); colnames(df) <- "condition"
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


#clustering of high p-value genes for quality assessment
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#weird to see that not neccesarily cluster on DKO+WT and WaplKO and SCC4, but rather SCC4+someWT and the rest
#then the rest subdivides into WT+WaplKO_onetype and DKO+Wapl_othertype

#PCA
pdf(file=paste0("/delaat/group/iwan/peakHiC/plots/PCA/RNA_PCA.pdf"))
plotPCA(vsd, intgroup="condition")
dev.off()

########individual gene all conditions visualisation######
#see per gene how it relates between different conds; plot the one with lowest p-value between WaplKO and Hap1
plotCounts(dds, gene=which.min(resLFC$padj), intgroup="condition")



#####################Save rds files of conditions, used apeglm shrinkage and alpha=0.05 for shrinkage#############
# dds <- readRDS("/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_DESeqDataSet_Hap1.rds")
# dds <- DESeq(dds)

resnames <- resultsNames(dds)
for(i in c(2,3,5)){
  #resnames[[i]]
  # [1] "Intercept"                "condition_DKO_vs_Hap1"
  # [3] "condition_SCC4_vs_Hap1"   "condition_Wapl2_vs_Hap1"
  # [5] "condition_WaplKO_vs_Hap1" "condition_WT2_vs_Hap1"
  library(apeglm)
  resLFC <- lfcShrink(dds, coef=resnames[[i]], type="apeglm")
  
  #select rows with an adjusted p-value of 0.05 or lower
  sig.resLFC <- resLFC[resLFC$padj< 0.05 & !is.na(resLFC$padj),]
  
  genemap <- readRDS("/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_ENSEMBL_Genes.rds")
  
  mapped.row <- lapply(FUN = function(gene){return(genemap[as.integer(which(genemap$ensembl == gene)),])},
                       rownames(sig.resLFC))
  ens.GR <- do.call("c",mapped.row)
  
  diff.exp <- GRanges(seqnames = seqnames(ens.GR), ranges = ranges(ens.GR), strand = strand(ens.GR), type = ens.GR$type,
                      transcript = ens.GR$transcript_name)
  mcols(diff.exp) <- c(mcols(diff.exp), sig.resLFC)
  
  saveRDS(diff.exp, paste0("/delaat/group/iwan/peakHiC/rds/plot_annotation/Diff_expression/",resnames[[i]],".rds"))
}

