# Define a function that filters files (either bed or bedpe) for overlap with peaks
filterPeakOverlap <- function(chroms,outname,bedtoolspath,peaksfileslop,
                              bedpe=FALSE,peakslop=0,verbose=FALSE)
{
  
  if (bedpe==FALSE)
  {
    #Overlap with peak files
    for (chrom in chroms)
    {
      readsfile    = paste(outname,"." ,chrom, ".bed",sep="")
      overlapfile  = paste(outname,"." ,chrom, ".peakfilt.bed",sep="")
      if (file.exists(readsfile) == TRUE)
      {
        command = paste(bedtoolspath , " intersect -u -a " ,readsfile ,
                        " -b ", peaksfileslop, " > ",overlapfile,sep="")
        if (verbose == TRUE)
        {
          print (command)
        }
        system(command)
      }
    }
  }
  
  if (bedpe==TRUE)
  {
    #Overlap with peak files
    for (chrom in chroms)
    {
      petsfile     = paste(outname,"." ,chrom, ".bedpe",sep="")
      overlapfile  = paste(outname,"." ,chrom, ".peakfilt.bedpe",sep="")
      tmpfile      = paste(outname,"." ,chrom, ".peakfilt.bedpe.tmp",sep="")
      if (file.exists(petsfile) == TRUE)
      {
        command = paste(bedtoolspath , " pairtobed -type both -a " ,petsfile ,
                        " -b ", peaksfileslop, " > ",tmpfile,sep="")
        if (verbose == TRUE)
        {
          print (command)
        }
        system(command)
        
        # bedtools prints out two lines for each one.  We need to trim for uniq ones
        everyotherline(tmpfile,overlapfile)
        
        # remove tmp file
        if (file.exists(tmpfile)) file.remove(tmpfile)
      }
    }
  }
}
