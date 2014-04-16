# Define a function that finds PET / peak overlaps
peakOverlap <- function(outname,datatype,peaksfile,verbose=FALSE)
{
  # (1) Split files by chromosome
  # files to split
  totreadsfile  = paste(outname,".", datatype,".bed",  sep="")
  totpetsfile   = paste(outname,".", datatype,".bedpe",sep="")
  
  # split reads by chromosome
  readschroms = splitBedbyChrom(totreadsfile,paste(outname, ".",datatype,sep="")) 
  petschroms  =      splitBedpe(totpetsfile, paste(outname, ".",datatype,sep=""), printreads=FALSE)
  
  # (2) Overlap with peak files
  for (chrom in petschroms)
  {
    readsfile    = paste(outname,"." ,datatype,"." ,chrom, ".bed",sep="")
    overlapfile  = paste(outname,"." ,datatype,"." ,chrom, ".bedNpeak",sep="")
    if (file.exists(readsfile) == TRUE)
    {
      command = paste(bedtoolspath , " intersect -wo -a " ,readsfile ,
                      " -b ", peaksfile, " > ",overlapfile,sep="")
      if (verbose == TRUE)
      {
        print (command)
      }
      system(command)
    }
  }
  
  # (3) gather information
  
  
  for (chrom in petschroms)
  {
    interactionfile = paste(outname,".", datatype , ".", chrom,".pairs.bedpe",sep="")
    if (file.exists(interactionfile)) file.remove(interactionfile)
    
    overlapfile   = paste(outname,"." ,datatype,"." ,chrom, ".bedNpeak",sep="")
    petpairsfile  = paste(outname,"." ,datatype,"." ,chrom, ".bedpe",sep="")
    
    findPairs(overlapfile,petpairsfile,interactionfile)
  }
  
  # (4) clean up temp files
  
  
}
