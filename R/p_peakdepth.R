# Define a function that models IAB as a function of peak depth
p_peakdepth <- function(chromosomes,outname,pdepthbins = 50 ,distrange=c(5000,10000000),requireddepth=5000,deptheffectplot)
{
  #------ (1) read in peaks to determine bins -------#
  depths = c()
  for (chrom in chromosomes)
  {    
    print (chrom)
    peaksfile    = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
    if (file.exists(peaksfile)     == FALSE){next}
    if (file.info(peaksfile)$size  == 0    ){next}
    chrpeaks = read.table(peaksfile,header=F,sep="\t")
    names(chrpeaks) = c("chr","start","end","name","score","strand")
    
    # calc center position
    chrpeaks$pos = (chrpeaks$start + chrpeaks$end) /2
    curdepths = c()
    for (i in (1:nrow(chrpeaks)))
    {
      #calc dist
      chrpeaks$dist  = log10(abs(chrpeaks$pos - chrpeaks$pos[i]))
      
      #calc depth
      chrpeaks$depth = as.numeric(chrpeaks$score) * as.numeric(chrpeaks$score[i])
      
      # filter for only pairs within the correct distance
      chrpeakssize = chrpeaks[which(chrpeaks$dist < log10(max(distrange[2])) & chrpeaks$dist > log10(max(distrange[1]))),]
      
      if ( length(which(is.finite(chrpeakssize$depth) ==FALSE)) >0)
      {
        print (chrpeakssize[i,])
        print (head(chrpeakssize[which(is.finite(chrpeakssize$depth) ==FALSE),]))
      }
      # add to curdepths
      depths = c(depths,sample(chrpeakssize$depth,size=min(20 , length(chrpeakssize$depth))))
    }
  }
  
  depths = sort(depths)
  a = seq(log2(min(depths)),log2(max(depths)),length.out=pdepthbins+1) 
  bins  = c()
  for (i in 2:(pdepthbins-1))
  {
    bins =c(bins,2^a[i])
  }
  depthmedianestimate = median(depths)

  #------ (2) get the total values per bin -------#

  # prepare genome wide counters of PETs for P(D|IB) and P(D|!IAB)
  readdepths                = rep(0,pdepthbins+1)
  combos                    = rep(0,pdepthbins+1)
  distances                 = rep(0,pdepthbins+1)
  IABS                      = rep(0,pdepthbins+1)
  names(readdepths)         = (0:pdepthbins)
  names(combos)             = (0:pdepthbins)
  names(distances)          = (0:pdepthbins)
  names(IABS)               = (0:pdepthbins)
  
  # combine peaks
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
      chrpeaks$depth = chrpeaks$score * chrpeaks$score[i]

      # filter for only pairs within the correct distance
      chrpeakssize = chrpeaks[which(chrpeaks$dist < log10(max(distrange[2])) & chrpeaks$dist > log10(max(distrange[1]))),]
      
      # assign to bins
      chrpeakssize$bin  = findInterval(chrpeakssize$depth, bins)
      for (eachbin in  (0:pdepthbins))
      {
        # add to M
        combosinbin = length(which(chrpeakssize$bin == eachbin))
        combos[as.character(eachbin)] = combos[as.character(eachbin)] + combosinbin
        
        # add to D
        readdepthsbin = sum(as.numeric(chrpeakssize$depth[(which(chrpeakssize$bin == eachbin))]))
        readdepths[as.character(eachbin)] = readdepths[as.character(eachbin)] + readdepthsbin
        
        # add to distances
        distsbin = sum(as.numeric(10^chrpeakssize$dist[(which(chrpeakssize$bin == eachbin))]))
        distances[as.character(eachbin)] = distances[as.character(eachbin)] + distsbin
      } 
    }
  }
  
  #------ (3) get the total values per bin -------#
  
  allpairs = c()
  for (chrom in chromosomes)
  {    
    print (chrom)
    pairsfile    = paste(outname,"." ,chrom, ".pairs.bedpe",sep="")
    if (file.exists(pairsfile)     == FALSE){next}
    if (file.info(pairsfile)$size  == 0    ){next}
    chrpairs= read.table(pairsfile,header=F,sep="\t")    
    allpairs = rbind(chrpairs)
  }

  # calc depth  
  allpairs$depth = allpairs[,10] * allpairs[,11]
  
  # assign to bins
  allpairs$bin  = findInterval(allpairs$depth, bins)

  for (eachbin in  (0:pdepthbins))
  {
    # add to D
    IABSbin = sum(as.numeric(allpairs[(which(allpairs$bin == eachbin)),12]))
    IABS[as.character(eachbin)] = IABS[as.character(eachbin)] + IABSbin
  } 
  
  combos=combos[0:pdepthbins]
  readdepths=readdepths[0:pdepthbins]
  distances=distances[0:pdepthbins]
  IABS=IABS[0:pdepthbins]
  
  depth_PETs = data.frame(distances=distances,combos=combos,depths=readdepths,IABS=IABS)
  
  
  #--- plot the model ----#
  
  pdf(deptheffectplot)
  par(mfrow=c(2,1))
  #-- plot peak depth model --# 
  x =  (depth_PETs$depths/depth_PETs$combos)[which(depth_PETs$combos>0)]
  y = (depth_PETs$IABS/depth_PETs$combos)[which(depth_PETs$combos>0)]
  plot(x,y,xaxt='n',yaxt='n',ann=F,pch=19,col="grey",log="xy")
  points(x[which(depth_PETs$combos>requireddepth)],y[which(depth_PETs$combos>requireddepth)],pch=19,col="dodgerblue3")
  axis(side=1,las=2,tcl=0.25)
  axis(side=2,las=2,tcl=0.25)
  mtext("depth",side=1,line=3,font=2)
  mtext("PETs per combo",side=2,line=3,font=2)
  spline_depth = smooth.spline(x[which(depth_PETs$combos>requireddepth)],y[which(depth_PETs$combos>requireddepth)],spar=1.7)
  lines(predict(spline_depth,x))
  
  y2 =combos
  barplot(y2,xaxt='n')
  mtext("bin",side=1,line=3,font=2)
  mtext("Possible pairs of peaks",side=2,line=3,font=2)
  axis(side=1,las=2,tcl=-0.25)
  dev.off()
  
  return(list(depth_PETs,spline_depth,depthmedianestimate)  )
}
