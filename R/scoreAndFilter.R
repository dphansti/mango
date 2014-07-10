
# Define a function that calculates P-values of interactions
calcP <- function(v)
{
  P = 1 - pbinom(q=v[1]-1,size=v[2],prob=v[3])
  return(P)
}

# Define a function that calculates P for each pair
scoreAndFilter <-function(chromosomes,outname ,mindist,maxdist,averageDepth,spline,N,corMeth="BY",normPmeth )
{
  # keep track of PETs filtered for size
  smallPETs = 0
  mediumPETs = 0
  longPETs = 0
  
  mindist = log10(mindist)
  maxdist = log10(maxdist)
  
  # make bins
  allpairs = c()
  for (chrom in chromosomes)
  { 
    # get total intra IAB
    pairsfile = paste(outname,".",chrom,".pairs.bedpe",sep="")
    pairs = read.table(pairsfile,header=FALSE,sep="\t")
    pairs = pairs[which(as.character(pairs$V8) != as.character(pairs$V9)),]
    pairs$dist = log10(abs( (pairs[,3]+pairs[,2]/2) - (pairs[,6]+pairs[,5]/2) ) )
    pairs = pairs[which(pairs$dist>mindist & pairs$dist<maxdist),]
    
    if (nrow(pairs) == 0)
    {
      next
    }
    
    # calculate p(success)
    if (normPmeth == "sum")
    {
      # calculate D (depth)
      pairs$D = pairs[,10] + pairs[,11]
      pairs$psuccess = predict(spline,pairs$dist)$y * (pairs$D / averageDepth)
    }
    
    # calculate p(success)
    if (normPmeth == "product")
    {
      # calculate D (depth)
      pairs$D = pairs[,10] * pairs[,11]
      pairs$psuccess = predict(spline,pairs$dist)$y * (pairs$D / averageDepth)
    }
    
    # add N
    pairs$N = N
    
    allpairs=rbind(allpairs,pairs)
  }
  
  # do the actual P-value calculations
  allpairs$P = apply(cbind(allpairs$V12,allpairs$N,allpairs$psuccess),1,calcP)    
  allpairs$Q = p.adjust(allpairs$P,method=corMeth,n=N)

  return(allpairs)
}