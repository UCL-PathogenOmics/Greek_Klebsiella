#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(BactDating)
  library(seqinr)
  library(lubridate)
})
options(stringsAsFactors = FALSE)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript run_bactdate.R <fastapath> <tree_path> <date_path> <out_nwk>")
}

fastapath <- args[1]
tree_path <- args[2]
date_path <- args[3]
out_nwk   <- args[4]

# Read FASTA file used to make tree
fasta <- read.fasta(fastapath)
site_length <- length(fasta[[1]])

# Read Tree 
tre <- read.tree(tree_path)
tre$edge.length <- tre$edge.length * site_length # Rescale branches by sites length

# Read dates file
dates  <- read.table(date_path, sep = "\t",
                     col.names = c("taxon","date"),
                     header = FALSE, stringsAsFactors = FALSE)
dates <- dates[order(match(dates$taxon, tre$tip.label)),]
dates <- decimal_date(as.Date(dates$date,format = "%Y/%m/%d"))

## Run BactDate
treR <- initRoot(tre, dates) # initial tree
bd <- bactdate(treR, dates, # Main function
               nbIts = 10000,     
               useRec = FALSE,
               updateRoot = FALSE,
               showProgress = FALSE)

write.tree(bd$tree, file = out_nwk)
