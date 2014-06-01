
############################# estimateP sub functions #############################

# Define a function that calculates P for each pair
calcIAB <-function(chromosomes,outname,numofbins ,binrange,outliers)
{
  # make bins
  bins = seq(log10(binrange[1]),log10(binrange[2]),length.out=numofbins+1)
  
  # prepare genome wide counters of PETs for P(D|IB) and P(D|!IAB)
  IABS                      = rep(0,numofbins)
  names(IABS)               = (1:numofbins)
  for (chrom in chromosomes)
  {
    print(chrom)

    # get total intra IAB
    pairsfile = paste(outname,".",chrom,".pairs.bedpe",sep="")
    pairs = read.table(pairsfile,header=FALSE,sep="\t")
    pairs = pairs[which(as.character(pairs$V8) != as.character(pairs$V9)),]
    pairs$dist = log10(abs( (pairs[,3]+pairs[,2]/2) - (pairs[,6]+pairs[,5]/2) ) )
    pairs$bin  = findInterval(pairs$dist, bins)
    pairs = pairs[which(pairs$bin>0 & pairs$bin<length(bins)),]
    
    # set outlier interactions to zero
    if (length(which(pairs[,7] %in% outliers)) > 0)
    {
      pairs[which(pairs[,7] %in% outliers),12] = 0
    }
    print(length(pairs[which(as.character(pairs[,7]) %in% as.character(outliers)),7]))
    
    if (nrow(pairs) == 0)
    {
      next
    }
    for (eachbin in  (1:numofbins))
    {
      # sum all IABs
      #IABSbin = length(which(pairs$bin == eachbin))
      IABSbin = sum(pairs[which(pairs$bin == eachbin),12])
      IABS[as.character(eachbin)] = IABS[as.character(eachbin)] + IABSbin
    }
  }
  
  return(IABS)
}


# Define a function that calculates M (# combos) and D (average read depth: CA + CB) from the peak overlap files
calcMandD <- function(chromosomes,outname,numofbins,binrange)
{
  
  # make bins
  bins = seq(log10(binrange[1]),log10(binrange[2]),length.out=numofbins+1)
  length(bins)
  
  # prepare genome wide counters of PETs for P(D|IB) and P(D|!IAB)
  readdepths                = rep(0,numofbins)
  combos                    = rep(0,numofbins)
  distances                 = rep(0,numofbins)
  names(readdepths)         = (1:numofbins)
  names(combos)             = (1:numofbins)
  names(distances)          = (1:numofbins)
  
  for (chrom in chromosomes)
  {
    print (chrom)
    
    peaksfile    = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
    
    # read in data
    if (file.exists(peaksfile)     == FALSE){next}
    if (file.info(peaksfile)$size  == 0    ){next}
    chrpeaks = read.table(peaksfile,header=F,sep="\t")
    names(chrpeaks) = c("chr","start","end","name","score","strand")
    
    # calc center position
    chrpeaks$pos = (chrpeaks$start + chrpeaks$end) /2
    
    for (i in (1:nrow(chrpeaks)))
    {
      #calc dist
      chrpeaks$dist  = log10(abs(chrpeaks$pos - chrpeaks$pos[i]))
      
      #calc depth
      chrpeaks$depth = chrpeaks$score + chrpeaks$score[i]
      
      # assign to bins
      chrpeaks$bin  = findInterval(chrpeaks$dist, bins)
      for (eachbin in  (1:numofbins))
      {
        # add to M
        combosinbin = length(which(chrpeaks$bin == eachbin))
        combos[as.character(eachbin)] = combos[as.character(eachbin)] + combosinbin
        
        # add to D
        readdepthsbin = sum(chrpeaks$depth[(which(chrpeaks$bin == eachbin))])
        readdepths[as.character(eachbin)] = readdepths[as.character(eachbin)] + readdepthsbin
        
        # add to distances
        distsbin = sum(10^chrpeaks$dist[(which(chrpeaks$bin == eachbin))])
        distances[as.character(eachbin)] = distances[as.character(eachbin)] + distsbin
      } 
    }
  }
  
  # convert distances to averages
  distances = log10(distances/combos)
  
  return(cbind(distances,combos,readdepths))
}


############################# estimateP function #############################

# Define a function that estimates p of success
estimateP <- function(chromosomes,outname,numofbins ,binrange,outliers)
{
  
  MDD    = calcMandD(chromosomes,outname,numofbins,binrange)
  IAB    = calcIAB(chromosomes,outname,numofbins,binrange,outliers)
  p_tab = data.frame(cbind(MDD,IAB))
  names(p_tab) = c("dist","M","D","IAB")
  p_tab$AvIAB = (p_tab$IAB/p_tab$M)
  averageDepth = sum(p_tab$D) / sum(p_tab$M)
  N = sum(p_tab$IAB)
  p_tab$p = (p_tab$AvIAB/N)
  spline = smooth.spline(p_tab$dist,p_tab$p,spar=.4)
  
  pEstimates = list( p_tab,spline,averageDepth,N) 
  names(pEstimates) = c("p_table","spline","averageDepth","N")
  
  return(pEstimates)
}
  
  
  


