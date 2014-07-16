
# Define a function that establishes bin borders for use with findInterval later
combineputativepairs <- function (chromosomes,outname)
{
  putpairs = c()
  for (chrom in chromosomes)
  {
    pairsfile = paste(outname,".",chrom,".pairs.bedpe",sep="")
    putpairs = rbind(putpairs,read.table(pairsfile,header=FALSE,sep="\t"))
  }
  # return
  return (putpairs)
}