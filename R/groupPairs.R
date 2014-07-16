# Define a function that finds PET / peak overlaps
groupPairs <- function(bedpefilesortrmdup,outname,peaksfile, bedtoolspath,
                     verbose=FALSE)
{

  # split reads by chromosome
  petschroms  =      splitBedpe(bedpefilesortrmdup, outname, printreads=TRUE)
  print(petschroms)

  petschroms = unique(petschroms)

  # (2) Overlap with peak files
  for (chrom in petschroms)
  {
    readsfile    = paste(outname,"." ,chrom, ".bed",sep="")
    overlapfile  = paste(outname,"." ,chrom, ".bedNpeak",sep="")
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

    interactionfile = paste(outname , ".", chrom,".pairs.bedpe",sep="")
    if (file.exists(interactionfile)) file.remove(interactionfile)

    overlapfile   = paste(outname,"." ,chrom, ".bedNpeak",sep="")
    petpairsfile  = paste(outname,"." ,chrom, ".bedpe",sep="")
    peakscount    = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
    findPairs(overlapfile,petpairsfile,interactionfile,peakscount)

  }

  
  # (4) clean up temp files
  for (chrom in  petschroms)
  {
    overlapfile   = paste(outname,"." ,chrom, ".bedNpeak",sep="")
    readsfile     = paste(outname,"." ,chrom, ".bed",sep="")
    petssfile     = paste(outname,"." ,chrom, ".bedpe",sep="")
#     if (file.exists(readsfile)) file.remove(readsfile)
#     if (file.exists(petssfile)) file.remove(petssfile)
#     if (file.exists(overlapfile)) file.remove(overlapfile)
  }
  
  return (petschroms)

}






