

# Define a function that plots the distribution of PET distances
plotdistancedistribution <- function(bedpefile,pdffile,npets=1000000)
{

  bedpe = read_tsv(bedpefile,col_names = F,n_max=npets)
  distances = (bedpe[,5]-bedpe[,2])[,1]
  log10distances = log10(distances[which(distances>0)])
  
  pdf(pdffile)
  par(mgp=c(3,.3,0))
  plot(density(log10distances),xaxt='n',ann=F,yaxt="n",col="dodgerblue2")
  at = c(2,4,6,8,10,12)
  axis(side=1,at=at,10^at,tcl=.25)
  axis(side=2,tcl=.25,las=2)
  mtext(side=1,line=2,font=2,text="distance (bp)")
  mtext(side=2,line=2,font=2,text="density")
  mtext(side=3,line=.5,font=2,text="Distribution of PET sizes")
  dev.off()
  
}