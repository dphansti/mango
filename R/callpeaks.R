
# Define a function that calls peaks using macs2
callpeaks <- function(macs2path,tagAlignfile,peaksfile,qvalue=0.05,
                      bedtoolspath,bedtoolsgenome,peakslop=0,MACS_shiftsize="NULL",gsize="hs")
{
  # call peaks
  shiftsize = ""
  if (MACS_shiftsize != "NULL")
  {
    shiftsize =  paste("--nomodel --extsize",MACS_shiftsize,sep=" ")
  }
  if (MACS_shiftsize == "NULL")
  {
    shiftsize = ""
  }

  command = paste(macs2path," callpeak -t ",tagAlignfile, shiftsize ," -g ",gsize,  " -f BED -n ",peaksfile," -q ",qvalue,sep=" ")
  print (command)
  system(command)
  
  # now shorten peak names
  peaks = read.table(paste(peaksfile,"_peaks.narrowPeak",sep=""),header=FALSE,sep="\t")
  peaks[,4] = paste("peak_",1:nrow(peaks),sep="")  
  
  write.table(peaks,file=paste(peaksfile,"_peaks.narrowPeak",sep=""),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=FALSE)
  
  # remove unneccesary files 
  listofsuffixes = c("peaks.xls","summits.bed","model.r","","")
  for (suf in listofsuffixes)
  {
    fname = paste(peaksfile,suf,sep="_")
    if (file.exists(fname) == TRUE)
    { 
        file.remove(fname)
    }
  }
}