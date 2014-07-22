

# Define a function that sorts bedpe file for duplicate removal
sortbedpe <- function(bedpefile,outname,bedpefilesort)
{
  # (1) split by chromosome
  chroms = splitBedpe(bedpefile, outname, printreads = FALSE , printpets = TRUE, skipstars=FALSE, skipinter=FALSE)
  chroms = unique(chroms)
  
  # (2) sort
  sortedchromfiles = c()
  for (chrom in chroms)
  {
    chromfile       = paste(outname, ".", chrom, ".bedpe",sep="")
    if (file.exists(chromfile))
    {
      sortedchromfile = paste(outname, ".", chrom, ".sort.bedpe",sep="")
      external_sort(chromfile, sortedchromfile)
      sortedchromfiles = c(sortedchromfiles,sortedchromfile)
    }
  }
  
  # (3) join sorted files
  if (file.exists(bedpefilesort)){file.remove(bedpefilesort)}
  joinchromfiles(sortedchromfiles,bedpefilesort)

  # (4) remove temp file
  for (chrom in chroms)
  {
    chromfile       = paste(outname, ".", chrom, ".bedpe",sep="")
    sortedchromfile = paste(outname, ".", chrom, ".sort.bedpe",sep="")
    if (file.exists(chromfile)){file.remove(chromfile)}
    if (file.exists(sortedchromfile)){file.remove(sortedchromfile)}
  }
}
