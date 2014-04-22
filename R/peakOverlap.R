# Define a function that finds PET / peak overlaps
groupPairs <- function(bedpefilesortrmdup,outname,datatype,peaksfile,  min_distance = 4,verbose=FALSE)
{
  # (1) Split files by chromosome
  # files to split
  #totreadsfile  = paste(outname,".", datatype,".bed",  sep="")
  #totpetsfile   = paste(outname,".", datatype,".bedpe",sep="")

  # split reads by chromosome
  #readschroms = splitBedbyChrom(totreadsfile,paste(outname, ".",datatype,sep="")) 
  petschroms  =      splitBedpe(bedpefilesortrmdup, outname, printreads=TRUE)[2]
  
  # (2) Overlap with peak files
  for (chrom in petschroms[[1]])
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
  for (chrom in petschroms[[1]])
  {
    interactionfile = paste(outname , ".", chrom,".pairs.bedpe",sep="")
    if (file.exists(interactionfile)) file.remove(interactionfile)
    
    overlapfile   = paste(outname,"." ,chrom, ".bedNpeak",sep="")
    petpairsfile  = paste(outname,"." ,chrom, ".bedpe",sep="")
    peakscount    = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
    findPairs(overlapfile,petpairsfile,interactionfile,peakscount)
  }
  
  # (4) clean up temp files
  for (chrom in petschroms[[1]])
  {
    overlapfile   = paste(outname,"." ,chrom, ".bedNpeak",sep="")
    readsfile     = paste(outname,"." ,chrom, ".bed",sep="")
    petssfile     = paste(outname,"." ,chrom, ".bedpe",sep="")
    if (file.exists(readsfile)) file.remove(readsfile)
    if (file.exists(petssfile)) file.remove(petssfile)
    if (file.exists(overlapfile)) file.remove(overlapfile)
  }
  
  # (5) score interactions
  allpairs = c()
  
  # make bins
  numofbins = 25
  bins = seq(log10(1000),8,length.out=numofbins)
  finalframe = c()
  
  for (chrom in petschroms[[1]])
  {
    peakscount    = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
    
    # read in data
    chrpeaks = read.table(peakscount,header=F,sep="\t")
    
    names(chrpeaks) = c("chr","start","end","name","score","strand")
    
    # calc center position
    chrpeaks$pos = (chrpeaks$start + chrpeaks$end) /2

    # prepare counters
    tots = rep(0,numofbins+1)
    lnks = rep(0,numofbins+1)
    names(tots) = (0:numofbins)
    names(lnks) = (0:numofbins)
    
    for (i in (1:nrow(chrpeaks)))
    {
      #calc dist
      chrpeaks$dist = log10(abs(chrpeaks$pos - chrpeaks$pos[i]))
      
      # assign to bins
      chrpeaks$bin  = findInterval(chrpeaks$dist, bins)
      
      # sum over bins
      summedvalues =  tapply(chrpeaks$score,chrpeaks$bin,sum)
      tots[names(summedvalues)] = tots[names(summedvalues)] + summedvalues
      
      # now you need add the peaks own value
      numberineachbin = table(chrpeaks$bin) * chrpeaks$score[i]
      tots[names(numberineachbin)] = tots[names(numberineachbin)] + numberineachbin
    }
    
    ######################## gather linkages (and adjust totals) ########################
    
    pairsfile    = paste(outname,"." ,chrom,".pairs.bedpe",sep="")
    
    chrpairs = read.table(pairsfile,header=FALSE,sep="\t")
    
    names(chrpairs) = c("chr1","start1","end1","chr2","start2","end2",
                        "name","p1name","p2name","CA","CB","IAB","dist")
    
    chrpairs$pos1 = (chrpairs$start1 + chrpairs$end1) / 2
    chrpairs$pos2 = (chrpairs$start2 + chrpairs$end2) / 2
    chrpairs$dist = log10(abs(chrpairs$pos1 - chrpairs$pos2))
    
    # assign to bins
    chrpairs$bin  = findInterval(chrpairs$dist, bins)
    
    # keep pairs
    allpairs = rbind(allpairs,chrpairs)
    
    # sum over bins
    summedvalues =  tapply(chrpairs$IAB,chrpairs$bin,sum)
    lnks[names(summedvalues)] = lnks[names(summedvalues)] + summedvalues
    
    # make dataframe
    pairframe = data.frame(chrom=chrom,binnumber=names(lnks),binlabel=c(NA,bins),links=lnks,totals=tots)
    pairframe$probs = pairframe$links / (pairframe$totals - pairframe$links)
    
    finalframe = rbind(finalframe,pairframe)
    
    #plot(pairframe$binlabel, pairframe$probs,type="l",col="blue",main=chrom)
  }
  
  # combine data from different chromosomes
  links  = tapply(finalframe$links,finalframe$binnumber,sum)
  totals = tapply(finalframe$totals,finalframe$binnumber,sum)
  combined = data.frame(binnumber= (0:numofbins),binlabel=c(NA,bins))
  row.names(combined) = combined$binnumber
  combined$links = 0
  combined$totals = 0
  combined[names(links),]$links  = links
  combined[names(totals),]$totals = totals
  combined$probs = combined$links / (combined$totals - combined$links)
  #plot(combined$binlabel, combined$probs,type="l",col="blue",ylim=c(0,0.001))
  
  ################ score interaction #################
  
  # filter for distance
  obsall = allpairs[which(allpairs$dist >= min_distance),]
  
  # look up probability
  obsall$prob = combined$probs[match(obsall$bin, combined$binnumber)]
  
  # calculate pvalue
  obsall$q    = obsall$IAB
  obsall$size = obsall$CA + obsall$CB - obsall$IAB
  
  # Define a function that calculates P
  calcP <- function(v)
  {
    P = 1 - pbinom(q=v[1]-1,size=v[2],prob=v[3])
    return(P)
  }
  
  obsall$P = apply(cbind(obsall$q,obsall$size,obsall$prob),1,calcP)
  obsall$Q = p.adjust(obsall$P,method="bonferroni")
  
  
}






