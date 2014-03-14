#' determines the illumina scoring system of the fastq
#'
#' This function reads in the first nlines of a fastq and makes a best guess at the illumina scoring system used
#'  
#' @param fastq full path to fastq file
#' @param nlines the number of lines to look at to determine the scoring method
#' @export
#' 
#fastq = "~/Desktop/mango2014test/NH.K562_RAD21_K562_std_2.1_1.head.fastq"
findScore <- function(fastq,nlines =10000) 
{
  # initialize
  scoretype = "unknown"
  
  # open connection
  con  <- file(fastq, open = "r")
  
  # read first lines
  lines <- readLines(con, n = nlines*4, warn = FALSE)
  
  # get every 4th line
  lines = paste(lines[seq(1, length(lines), 4)],collapse="")
  
  # check score
  if (grepl("0", lines) == TRUE)
  {
    # Sanger or Illumina 1.8
    if (grepl("J", lines) == TRUE)
    {
      # Illumina 1.8
      print (paste("score type detected:","Illumina 1.8"))
      scoretype = "--phred33-quals"
    }
    if (grepl("J", lines) == FALSE)
    {
      # Sanger
      scoretype = "--phred33-quals"
      print (paste("score type detected:","Sanger"))
    }
  }
  if (grepl("0", lines) == FALSE)
  {
    #  Solexa, Illumina 1.3, Illumina 1.5 
    if (grepl("=", lines) == TRUE)
    {
      # Solexa
      scoretype = "--solexa-quals"
      print (paste("score type detected:","Solexa"))
    }
    if (grepl("A", lines) == TRUE)
    {
      # Illumina 1.3
      scoretype = "--phred64-quals"
      print (paste("score type detected:","Illumina 1.3"))
    }
    if (grepl("A", lines) == FALSE & grepl("=", lines) == FALSE)
    {
      # Illumina 1.5
      scoretype = "--phred64-quals"
      print (paste("score type detected:","Illumina 1.5"))
    }
  }  
  
  # close connection
  close(con)
  
  # return the score
  return (scoretype)
}