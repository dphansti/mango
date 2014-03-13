#' aligns reads using bowtie
#'
#' This function aligns reads using bowtie
#'  
#' @param fastq full path to fastq file
#' @param bowtiepath full path to bowtie
#' @param bowtieref full path to bowtie reference
#' @param arguments arguments to give to bowtie
#' @param verbose boolean whether or not to print command
#' @param nlines the number of lines to look at to determine the scoring method
#' @export
#' 
alignBowtie <- function(fastq,bowtiepath,bowtieref,
                        arguments="v 2 -k 1 --sam-nohead --mapq 40 --best --strata -m 1 -",
                        verbose=TRUE)
{
  # determine the illumina score encoding
  illuminascore = findScore(fastq,nlines =10000) 
  
  # establish bowtie settings
  bowtievar = "v 2 -k 1 --sam-nohead --mapq 40 --best --strata -m 1 -"
  
  # form full command
  bowtiecommand = paste (bowtiepath, bowtieref, illuminascore, bowtievar)
  
  # print command if desirec
  if (verbose ==TRUE)
  {
    print (bowtiecommand)
  }
  
  # execute command
  system(bowtiecommand)
}