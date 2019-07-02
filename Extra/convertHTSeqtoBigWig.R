library(BSgenome.Hsapiens.UCSC.hg19)
# genome <- BSgenome.Hsapiens.UCSC.hg19
# bins <- tileGenome(seqlengths=seqlengths(genome),tilewidth=2e3)
# bins <- unlist(bins)

# genes <- readRDS(file="/delaat/group/iwan/peakHiC/rds/plot_annotation/Diff_expression/condition_WaplKO_vs_Hap1.rds")
# bwBins <- bins[ovl@to]

reads <- readRDS("/delaat/group/iwan/peakHiC/HTseq_for_diff_expr/rds/GSM2493895_3844_10_DKO_1_TAGCTTA_S28_L003_topHat2-standed.HTSeq.rds")
reads2 <- readRDS("/delaat/group/iwan/peakHiC/HTseq_for_diff_expr/rds/GSM2493896_3844_11_DKO_2_GGCTACA_S29_L003_topHat2-standed.HTSeq.rds")
reads3 <- readRDS("/delaat/group/iwan/peakHiC/HTseq_for_diff_expr/rds/GSM2493897_3844_12_DKO_3_CTTGTAA_S30_L003_topHat2-standed.HTSeq.rds")
reads$counts <- (reads$counts+reads2$counts+reads3$counts)/3
  
genesGR <- readRDS("/delaat/group/iwan/peakHiC/rds/plot_annotation/hg19_ENSEMBL_Genes.rds")
genesGR <- genesGR[!is.na(match(genesGR$ensembl,reads$ensID))]

bwGenes <- GRanges(seqnames=seqnames(genesGR),ranges=ranges(genesGR),seqinfo=seqinfo(genome))
bwGenes$score <- reads$counts[match(genesGR$ensembl,reads$ensID)]

bwGenes <- unique(bwGenes)

#then reduce object because for a bigwig it cannot be overlapping
bwOut <- GenomicRanges::reduce(bwGenes)
ovl <- findOverlaps(bwOut,bwGenes)
bwOut$score <- tapply(bwGenes$score[ovl@to],ovl@from,sum)

export.bw(object=bwOut,con="/delaat/group/iwan/peakHiC/HTseq_for_diff_expr/Ruiqi/Hap1_DKO_GeneCount.bw")


