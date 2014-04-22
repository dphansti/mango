
# Define a function that scores interactions
scorePairs <- function(chromosomes,outname,min_distance = 1000,maxPval=0.01,
                       numofbins = 25,binrange=c(1000,250000000),
                       corrMethod="bonferroni",verbose=FALSE)
{
  logmaxP = -log10(maxPval)
  
  # Define a function that calculates P
  calcP <- function(v)
  {
    P = 1 - pbinom(q=v[1]-1,size=v[2],prob=v[3])
    return(P)
  }
  
  interactionfiletemp  = paste(outname ,".interactions.P.bedpe",sep="")
  interactionfilefinal = paste(outname ,".interactions.bedpe",sep="")
  
  # make bins
  numofbins = numofbins
  bins = seq(log10(binrange[1]),log10(binrange[2]),length.out=numofbins)
  finalframe = c()

  # prepare genome wide counters
  tots_gw = rep(0,numofbins+1)
  lnks_gw = rep(0,numofbins+1)
  names(tots_gw) = (0:numofbins)
  names(lnks_gw) = (0:numofbins)
  
  # update progess
  if (verbose == TRUE )
  {
    print ("Deterimining interaction probabilities")
  }
  
  for (chrom in chromosomes)
  {
    if (verbose == TRUE )
    {
      print (chrom)
    }
    
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
    
    # add to genome wide totals
    tots_gw = tots_gw + tots
    
    ######################## gather linkages (and adjust totals) ########################
    
    pairsfile    = paste(outname,"." ,chrom,".pairs.bedpe",sep="")
    
    chrpairs = read.table(pairsfile,header=FALSE,sep="\t")
    
    names(chrpairs) = c("chr1","start1","end1","chr2","start2","end2",
                        "name","p1name","p2name","CA","CB","IAB","dist")

    # get distances
    chrpairs$pos1 = (chrpairs$start1 + chrpairs$end1) / 2
    chrpairs$pos2 = (chrpairs$start2 + chrpairs$end2) / 2
    chrpairs$dist = log10(abs(chrpairs$pos1 - chrpairs$pos2))
    
    # assign to bins
    chrpairs$bin  = findInterval(chrpairs$dist, bins)
    
    # keep pairs
    #allpairs = rbind(allpairs,chrpairs)
    
    # sum over bins
    summedvalues =  tapply(chrpairs$IAB,chrpairs$bin,sum)
    lnks[names(summedvalues)] = lnks[names(summedvalues)] + summedvalues
    
    # add to genome wide totals
    lnks_gw = lnks_gw + lnks
    
    # make dataframe
    pairframe = data.frame(chrom=chrom,binnumber=names(lnks),binlabel=c(NA,bins),links=lnks,totals=tots)
    pairframe$probs = pairframe$links / (pairframe$totals - pairframe$links)
    
    finalframe = rbind(finalframe,pairframe)
  }
  
  # make a data frame of genome wide probabilities
  combined        = data.frame(binnumber= (0:numofbins),binlabel=c(NA,bins))
  combined$links  = lnks_gw
  combined$totals = tots_gw
  combined$probs  = combined$links / (combined$totals - combined$links)
  
  # keep track of the minimum bin for plotting purposes
  minbin = numofbins
  
  if (verbose == TRUE )
  {
    print ("Calculating P-values")
  }
  
  # Calculate the P-values
  Pvalues = c()
  firstchrom = TRUE;
  for (chrom in chromosomes)
  {
    if (verbose == TRUE )
    {
      print (chrom)
    }
    
    pairsfile    = paste(outname,"." ,chrom,".pairs.bedpe",sep="")
    chrpairs = read.table(pairsfile,header=FALSE,sep="\t")
    names(chrpairs) = c("chr1","start1","end1","chr2","start2","end2",
                        "name","p1name","p2name","CA","CB","IAB","dist")
    
    # get distances
    chrpairs$pos1 = (chrpairs$start1 + chrpairs$end1) / 2
    chrpairs$pos2 = (chrpairs$start2 + chrpairs$end2) / 2
    chrpairs$dist = log10(abs(chrpairs$pos1 - chrpairs$pos2))
    
    # assign to bins
    chrpairs$bin  = findInterval(chrpairs$dist, bins)
    
    # filter for distance
    chrpairs = chrpairs[which(chrpairs$dist >= log10(min_distance)),]
    
    minbin = min(minbin,min(chrpairs$bin))
    
    # look up probability
    chrpairs$prob = combined$probs[match(chrpairs$bin, combined$binnumber)]
    
    # calculate pvalue
    chrpairs$q    = chrpairs$IAB
    chrpairs$size = chrpairs$CA + chrpairs$CB - chrpairs$IAB
    chrpairs$P    = apply(cbind(chrpairs$q,chrpairs$size,chrpairs$prob),1,calcP)
    Pvalues       = c(Pvalues,chrpairs$P)
    chrpairs$P    = -log10(chrpairs$P)
    chrpairs$P[which(is.infinite(chrpairs$P ) == TRUE)] = 16
    
    # print to output
    if (firstchrom == FALSE)
    {
      write.table(chrpairs, file=interactionfiletemp,quote = FALSE, sep = "\t",row.names = FALSE,
                  col.names = FALSE,append=TRUE)
    }
    if (firstchrom == TRUE)
    {
      write.table(chrpairs, file=interactionfiletemp,quote = FALSE, sep = "\t",row.names = FALSE,
                  col.names = TRUE)
      firstchrom = FALSE
    }
  }
  
  # Adjust the P-values
  Q = -log10(p.adjust(Pvalues,method=corrMethod))
  Q[which(is.infinite(Q) == TRUE)] = 16
  
  numsig = length(which(Q > logmaxP))
  
  # print number of signicant interactions
  print (paste("Number of significant interactions (at P <= ",10^-logmaxP,") = ",numsig,sep="") )

  # Add Q-values to file
  AddQvals(interactionfiletemp,interactionfilefinal,Q,logmaxP)
  
  # remove the zero bin
  combined = combined[-1,]
  
  # make plots of probabilities
  pdfofprobs = paste(outname, ".probabilities.pdf",sep="")
  pdf(pdfofprobs)
  par(mgp = c(3,.3, 0),mar=c(5,4,3,4))
  plot(10^combined$binlabel, combined$probs,type="l",col="blue",xlab="",ylab="",log="x",yaxt="n",xaxt="n")
  axis(side=1,las=2,tck=0.015, labels=T)
  axis(side=2,las=2,tck=0.015, labels=T)
  mtext("probability",side=2,line=3,font=2)
  mtext("interaction distance",side=1,line=3,font=2)
  
  par(fig=c(0.4,.9,.4,0.9), new=T,mar=c(2,2,1,1))
  subcombined = combined[minbin:length(combined$binlabel),]
  plot(10^subcombined$binlabel, subcombined$probs,type="l",col="blue",xlab="",ylab="",log="x",yaxt="n",xaxt="n")
  axis(side=1,las=2,tck=0.015, labels=T)
  axis(side=2,las=2,tck=0.015, labels=T)

  mtext("probability",side=2,line=3,font=2)
  mtext("interaction distance",side=1,line=3,font=2)
  
  par(fig=c(0,1,0,1))
  dev.off()
}



# nrow(sig)
# head(sig)
# sig = res[which(res$adjP > 2 & res$q >1),]
# plot(density(sig$dist))
# res =  read.table(interactionfilefinal,header=TRUE,sep="\t")
# res = res[which(res$dist > log10(min_distance)),]
# sig = res[which(res$"X.log10.adjP." > 2),]
# not = res[which(res$"X.log10.adjP." < 2),]
# nrow(sig)
# nrow(not)
# 
# hist(sig$q[which(sig$q <=30)],breaks=(1:30))
# hist(not$q[which(not$q <=30)],breaks=(2:30))
# 
# 
# 
# head(sig)
# names(res)
# 
# res$otherP = -log10(res$P)
# head(res)
# -log10(0.01)
# plot(density(sig$dist))
# 
# nrow(res[which(res$dist > 7 & res$q == 1),])
# head(res[which(res$dist > 7 & res$q == 1),])



#   # combine data from different chromosomes
# numofbins = 50
#   links  = tapply(finalframe$links,finalframe$binnumber,sum)
#   totals = tapply(finalframe$totals,finalframe$binnumber,sum)
#   combinedNEW = data.frame(binnumber= (0:numofbins),binlabel=c(NA,bins))
#   row.names(combinedNEW) = combinedNEW$binnumber
#   combinedNEW$links = 0
#   combinedNEW$totals = 0
#   combinedNEW[names(links),]$links  = links
#   combinedNEW[names(totals),]$totals = totals
#   combinedNEW$probs = combinedNEW$links / (combinedNEW$totals - combinedNEW$links)
#   combined
#   #plot(combined$binlabel, combinedNEW$probs,type="l",col="blue",ylim=c(0,0.001))
#   
#   ################ score interaction #################
# 
#   # filter for distance
#   obsall = allpairs[which(allpairs$dist >= min_distance),]
#   
#   # look up probability
#   obsall$prob = combined$probs[match(obsall$bin, combined$binnumber)]
#   
#   # calculate pvalue
#   obsall$q    = obsall$IAB
#   obsall$size = obsall$CA + obsall$CB - obsall$IAB
#   
#   # Define a function that calculates P
#   calcP <- function(v)
#   {
#     P = 1 - pbinom(q=v[1]-1,size=v[2],prob=v[3])
#     return(P)
#   }
#   
#   obsall$P = apply(cbind(obsall$q,obsall$size,obsall$prob),1,calcP)
#   obsall$Q = p.adjust(obsall$P,method="bonferroni")
#   
#   1.668131e+02
#   (obsall$P * nrow(obsall))[200:205]
#   (obsall$Q)[200:205]
#   
#   # filter and print out significant
#   sig = obsall[which(obsall$Q <= maxPval),]
#   
#   # print to output
#   interactionfile = paste(outname, ".interactions.bedpe",sep="")
#   write.table(sig, file=interactionfile,quote = FALSE, sep = "\t",row.names = FALSE,
#       col.names = TRUE)
# 
# }