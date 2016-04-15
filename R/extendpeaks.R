
# Define a function that extends peaks and merges ones that overlap
extendpeaks <- function(peaksfile=peaksfile,peaksfileslop=peaksfileslop,
                        bedtoolspath=bedtoolspath,bedtoolsgenome=bedtoolsgenome,
                        peakslop=peakslop,blacklist,verbose=FALSE)
{
  
  # make peakfile with slop
  command = paste(bedtoolspath, " slop -i ",peaksfile, " -g ",bedtoolsgenome," -b ",peakslop,
                  " | ", bedtoolspath, " merge -c 4 -o collapse > ",peaksfileslop,sep="")
  if (verbose == TRUE){ print (command) }
  system(command)
  
  if (blacklist != "NULL")
  {
    command = paste(bedtoolspath, " intersect -v -a ",peaksfileslop, " -b ",blacklist,
                    " > tmp.bed ; mv tmp.bed ",peaksfileslop,sep="")
    if (verbose == TRUE){ print (command) }
    system(command)
  }
  
  # clean up the peak names
  mergedpeaks = read_tsv(peaksfileslop, col_names = FALSE)
  mergedpeaks$V4 = paste("speak_",1:nrow(mergedpeaks),sep="")
  write.table(mergedpeaks,file=peaksfileslop,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=FALSE)
  
  # count lines in peak files
  numpeaks = nrow(read_tsv(peaksfile, col_names=FALSE))
  nummergepeaks = nrow(mergedpeaks)
  return(c(numpeaks,nummergepeaks)) 
}