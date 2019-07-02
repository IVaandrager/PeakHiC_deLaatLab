library("rtracklayer")
library("dplyr")
library("GenomicRanges")


hist.dir <- "/delaat/group/iwan/peakHiC/chip_for_vp/histone_marks2/"
histfiles <- list.files(path = enhancer.dir, pattern = "[.]bed$",full.names = TRUE)

marks <- as.list(rep(c('H3K27me3', 'H2AK119ub1'), each = 2))

