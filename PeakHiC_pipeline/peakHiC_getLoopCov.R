library("GenomicRanges")
library( "parallel")
library( "data.table")
source("peakHiC_functions_getloopCov.R")

loopCov <- function(inFolder, outFolder, genomeObj, tracks, partReadsFldr){

  #reduced 
  inFiles <- list.files(inFolder, full.names = T, include.dirs = F)
  in.idx <- grep(inFiles, pattern = "*reduced.rds")
  loops <- readRDS(inFiles[in.idx])
  
  partIDs <- unique(genomeObj$vpsGR$partID[match(unique(loops$id),genomeObj$vpsGR$vpID)])
  
  #run parallel
  nThreads <- 6
  cl <- makeCluster(nThreads, outFile="")
  clusterEvalQ(cl, c(library("GenomicRanges"), library("data.table"),source("peakHiC_functions_getloopCov.R")))

  loopTags <- parLapply(cl, as.list(partIDs), getLoopCovbyPartition,loops=loops,genomeObj=genomeObj,
                        anchorSize=10e3, tracks=tracks, partReadsFldr=partReadsFldr)

  loopTags <-rbindlist(loopTags)

  saveRDS(loopTags,file=paste0(outFolder,"loopCoverage.rds"))

  on.exit(stopCluster(cl))
}



#non-parallel test
# testerlist <- list()
# for(i in 1:length(partIDs)){
#     	tryCatch({
#         print(i)
#     		nextTags <- getLoopCovbyPartition(partID=partIDs[i],loops=loops,genomeObj=genomeObj,anchorSize=10e3, tracks=tracks, partReadsFldr=partReadsFldr)
#     		testerlist[[i]] <- nextTags
# 
#     	}, error=function(e) { message(paste0(partIDs[i]," failed")) })
# }
# 
# wat <- rbindlist(testerlist)
# saveRDS(wat,file=paste0(outFolder,"loopCoverage.rds"))


# VERSION <- '1.00'
# 
# get_script_path <- function(path=NULL) {
#   if( is.null( path ) ){
#     cmdArgs = commandArgs(trailingOnly = FALSE)
#     needle = "--file="
#     match = grep(needle, cmdArgs)
#     if (length(match) > 0) {
#         # Rscript
#         return(normalizePath(sub(needle, "", cmdArgs[match])))
#     } else {
#         ls_vars = ls(sys.frames()[[1]])
#         if ("fileName" %in% ls_vars) {
#             # Source'd via RStudio
#             return(normalizePath(sys.frames()[[1]]$fileName))
#         } else {
#             # Source'd via R console
#             return(normalizePath(sys.frames()[[1]]$ofile))
#         }
#     }
#   } else {
#     return(path)
#   }
# }
# 
# #################################################################################################################
# ### PARSING THE INPUT ###########################################################################################
# #################################################################################################################
# if( !suppressMessages(require( "argparse", character.only = TRUE ) ) ) stop( "Package not found: argparse" )
# 
# parser <- ArgumentParser()
# 
# parser$add_argument('-genomeObjFile', help='path to peakHi-C genomeObj file', metavar='/path/to/vp_file', required=TRUE )
# parser$add_argument('-loopsFile', help='path to peakHi-C loopFile (.rds)', metavar='/path/to/loop_file', required=TRUE )
# parser$add_argument('-outFolder', help='name of the output file [default %(default)s]', metavar='/path/to/output_folder/',required=TRUE  )
# parser$add_argument('-cores', type='integer', help='number of cores for parallelization [default %(default)s]', default=1, metavar='XX' )
# 
# argsL <- parser$parse_args()

