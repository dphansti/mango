
# Define a function that determines distance cutoff
calcDistBias <- function(distancefile,distancecutpdf,range,biascut)
{
  # Define a function that calculates self-ligation bias
  selfligationbias <- function(Sobs,Total)
  {
    bias = -(Sobs/0.5 - Total) / Total
    return (bias)
  }
  
  mindist = range[1]
  maxdist = range[2]
  
  # read in data
  dist = read.table(distancefile,header=FALSE,sep="\t")
  names(dist) = c("distance","orientation")
  
  # filter for min distance
  dist = dist[which(dist$distance > mindist & dist$distance < maxdist),]
  
  # make bins
  numofbins = 50
  bins = seq(log10(mindist),log10(maxdist),length.out=numofbins)
  
  # assign to bins
  dist$bin  = findInterval(log10(dist$distance), bins)
  S = hist(dist$bin[which(dist$orientation == "S")],breaks=(0:numofbins),plot=FALSE)
  D = hist(dist$bin[which(dist$orientation == "D")],breaks=(0:numofbins),plot=FALSE)
  
  
  # calculate the biases
  biases = c()
  for (i in 1:numofbins)
  {
    Sobs =   S$counts[i]
    Total =  D$counts[i] + Sobs
    bias = selfligationbias(Sobs,Total)
    biases= c (biases,bias)
  }
  distancecutoff = signif(10^ (head(bins[which(biases<biascut)],n=1)[1] ),digits=5)
  
  #- plot results -#
  
  pdf(distancecutpdf)
  
  # plot S/D ratio
  par(mgp = c(2,.3, 0),mar=c(5,4,3,4))
  plot(x=10^bins ,y=log2(S$counts/D$counts),log="x",type='b',pch=19,
       ylab="",xlab="",col="dodgerblue2",yaxt='n')
  axis(side=2,las=2,tck=0.01,col="dodgerblue2", labels=F)
  at = axTicks(2)
  mtext(side = 2, text = at, at = at, col = "dodgerblue2", line = 0.2,las=2)
  mtext("log2 (same / dif)",side=2,line=2,font=2,col="dodgerblue2")
  mtext("PET distance",side=1,line=2,font=2)
  
  par(new=T)
  
  # plot calculated bias
  plot(x=10^bins ,y=biases,log="x",type='b',pch=19,col="firebrick2",
       ann=FALSE,xaxt='n',yaxt='n')
  axis(side=4,las=2,tck=0.01,ylab='', labels=F)
  at = axTicks(4)
  mtext(side = 4, text = at, at = at, col = "firebrick2", line = 0.2,las=2)
  mtext("% due to bias",side=4,line=2,font=2,col="firebrick2")
  
  abline(v=distancecutoff)  
  abline(h=biascut,col="firebrick2",lty=2) 
  legend("right",inset=0.05,legend=paste("cutoff = ",distancecutoff))
  
  dev.off()
  
  return (distancecutoff)
}