# Define a function that finds PET / peak overlaps
groupPairs <- function(bedpefilesortrmdup,outname,peaksfile, bedtoolspath,
                     bedtoolsgenome,extendreads=120,verbose=FALSE)
{

  # split reads by chromosome
  print ("---Splitting PETs by chromosome")
  petschroms  =      splitBedpe(bedpefilesortrmdup, outname, printreads=TRUE)

  petschroms = unique(petschroms)

  # (2) Overlap with peak files
  print ("---Intersecting PETs with peaks")
  for (chrom in petschroms)
  {
    readsfile     = paste(outname,"." ,chrom, ".bed",sep="")
    readexendfile = paste(outname,"." ,chrom, ".extend.bed",sep="")
    overlapfile   = paste(outname,"." ,chrom, ".bedNpeak",sep="")
    if (file.exists(readsfile) == TRUE)
    {
      command = paste(bedtoolspath , " slop -r ",extendreads," -s -i " ,readsfile ,
                      " -g ", bedtoolsgenome, " > ",readexendfile,sep="")
      
      command = paste(bedtoolspath , " intersect -wo -a " ,readexendfile,
                      " -b ", peaksfile, " > ",overlapfile,sep="")
      if (verbose == TRUE)
      {
        print (command)
      }
      system(command)
      
    }
  }
  
  # (3) gather information
  print ("---Building putative interaction set")
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
    readexendfile = paste(outname,"." ,chrom, ".extend.bed",sep="")
    if (file.exists(readexendfile)) file.remove(readexendfile)
#      if (file.exists(readsfile)) file.remove(readsfile)
#      if (file.exists(petssfile)) file.remove(petssfile)
#      if (file.exists(overlapfile)) file.remove(overlapfile)
  }
  
  return (petschroms)

}






