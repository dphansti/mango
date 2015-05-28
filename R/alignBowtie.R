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
alignBowtie <- function(fastq,output,bowtiepath,bowtieref,
                        shortreads,verbose=TRUE)
{
  
  # choose alignment parameters
  # note- "-m 1" ensures that only uniquely mapped reads are reported.
  bowtievar="-S -v 0 -k 1 --chunkmbs 500--sam-nohead --mapq 40 -m 1"
  if (shortreads == FALSE)
  {
    #bowtievar="-S -v 2 -k 1 --sam-nohead --mapq 40 -m 1 --best"
    bowtievar="-S -n 2 -l 50 -k 1 --chunkmbs 500 --sam-nohead --mapq 40 -m 1 --best"
  }
  
  # determine the illumina score encoding
  illuminascore = findScore(fastq,nlines =10000) 
  
  # form full command
  bowtiecommand = paste (bowtiepath, bowtieref, fastq, output, illuminascore, bowtievar)
  
  # print command if desirec
  if (verbose ==TRUE)
  {
    print (bowtiecommand)
  }
  
  # execute command
  system(bowtiecommand)
}

