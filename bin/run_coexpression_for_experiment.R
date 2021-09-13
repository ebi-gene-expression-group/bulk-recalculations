#!/usr/bin/env Rscript
# Create a gene coexpression matrix for a baseline Expression Atlas experiment
# Input: undecorated tsv file, first column gene identifier, subsequent columns

# Load BiocParallel package
suppressMessages(library(BiocParallel))
suppressMessages(library(clusterSeq))

print( paste0('node: ', Sys.info()[["nodename"]] )  )
options(ports=sample(1024:49151,1)  )

# Get the commandline arguments
args <- commandArgs(TRUE)
if (length(args) == 3) {
  # Input file
  expressionsFile <- args[1]

  # Destination (if not finishing with .gz, clusterSeq will append .gz to the
  # end)
  outputPath <- args[2]
  wdir <- as.character(args[3])
  source(paste0(wdir,'/bin/kcluster_parallel.R'))
} else {
  # Print a usage message and exit.
  stop("\nUsage:
        \n\trun_coexpression_for_experiment.R
        <expressions file: gene identifier with quartiles> <output path .gz> <workflow.basedir>\n")
}

# read file
exp <- read.delim(expressionsFile)

# Make sure there are at least three columns (gene ID, expression columns) --
# doesn't make sense to try coexpression with only one or two columns of
# expression data.
if (ncol(exp) < 3) {
  warning("Fewer than three columns in total. Cannot do coexpression on only one data column.")
  # Quit without exit code.
  q(save = "no")
}

exp[, 2:ncol(exp)] <- sapply(exp[2:ncol(exp)], function(x) sub("(^[^,]+[,][^,]+[,])([^,]+)(,.+$)",
  "\\2", x))  #get the middle value for each gene/tissue
expL <- sapply(exp[, 2:ncol(exp)], as.numeric)  # make sure the values are numeric
rownames(expL) <- exp[, 1]
expL <- log(expL)  # get the natural logarithm

expL[is.na(expL)] <- 0  # turn any NAs to 0

cD <- expL[rowSums(expL) > ncol(expL), ]  # Filter out non-expressed genes

# run kCluster function to create coexpression matrices, set to use 8 cores.
# It can be changed to use less or more

cores <- 8
max_avail_workers <- multicoreWorkers()

use_cores <- max(1, min(cores, max_avail_workers -2) )

# up to 5 retries in case of error/warning
for (i in 1:5){
  print( paste0( 'kClust will run with ' , use_cores  , ' cores... Starting try number: ', i)  )
  out <- tryCatch(
      {
        message("This is the 'try' part")
        kCluster(cD, ncores= use_cores ,  matrixFile = outputPath)
      },
      error=function(cond) {
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
      },
      warning=function(cond) {
        message("Here's the original warning message:")
        message(cond)
        # Choose a return value in case of warning
        return(NA)
      },
      finally={
        message(paste("Finished try number:", i))
      }
      )

  if ( is.matrix(out) == TRUE  ) {
    print(paste0('kClust finished successfully after ', i, ' iterations.') )
    break
  }
}


