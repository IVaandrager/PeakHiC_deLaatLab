#PeakHiC pipeline

setwd("/delaat/group/iwan/peakHiC/scripts")

#genomeObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/Hap1_CTCF_TSS_enhancers_reduced.rds")
genomeObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/genomeObj/Hap1_CTCF_TSS_enhancers_reduced2.rds")
#made in GenomeObj_creation.R
#change $hic$design if different grouping of conditions; function in GenomeObj_creation, change.obj.design

frags <- readRDS("/delaat/group/iwan/Hap1_peakHiC/objects/hg19_endonuclease_frags/hg19_MboI_frags_byChr.rds")
tracks <- as.list(1:9)
hicCond <- "Hap1_9reps"

###############create partition reads################
# source("par_partition_reads.R")
# 
# inFolder <- "/delaat/group/iwan/peakHiC/READS/"
# outFolder <- "/delaat/group/iwan/Hap1_peakHiC/analysis/partition_reads/"
# 
# partition_reads(inFolder, outFolder, genomeObj, tracks)
# 
# message("partition reads created")

###############create V4Cs for all VPs################
source("par_vp_V4C.R")

partReadsFldr <- "/delaat/group/iwan/Hap1_peakHiC/analysis1/partition_reads/"
outFolder <- "/delaat/group/iwan/Hap1_peakHiC/analysis/V4Cs/"

VP_V4Cs(partReadsFldr, outFolder, genomeObj, hicCond, tracks, frags)

message("V4Cs created")


################call peaks on combined replicates##########
source("par_peakCs.R")

inFolder <- outFolder
outFolder <- "/delaat/group/iwan/Hap1_peakHiC/analysis/peaks/"

call_peaks(inFolder, outFolder, genomeObj)

message("peaks called")

############verify called peaks by checking reciprocality###########
source("reciprocal_peaks.R")

inFolder <- outFolder

check_reciprocal(inFolder)

message("peaks validated")

##########retrieve loop coverage per condition, per peaks##########
source("peakHiC_getLoopCov.R")

outFolder <- "/delaat/group/iwan/Hap1_peakHiC/analysis/loop_coverage/"
  
loopCov(inFolder, outFolder, genomeObj, tracks, partReadsFldr)

message("loop coverage calculated")


##########create object to perform statistics on#########
#needs to be run from Drakoon because folder of geert not properly mounted
source("statObj_creation.R")

infile_cov <- "/delaat/group/iwan/Hap1_peakHiC/analysis/loop_coverage/loopCoverage.rds"
infile_loops <- "/delaat/group/iwan/Hap1_peakHiC/analysis/peaks/GW_nReps_9_peakHiC_wSize_41_qWr_1_alphaFDR_0.1_loops_reduced.rds"
outFolder <- "/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/"
geneCount <- "/delaat/group/iwan/Hap1_peakHiC/objects/Expression_levels/normHTseq/normHTseq.rds"
#DE.expr, vpID_to_prot, genemap for EnsemblID and chip_files are hardcoded

getStatObj(infile_cov,infile_loops,outFolder,genomeObj,geneCount)

message("statObj created")


#Go to loop_statistics.R


