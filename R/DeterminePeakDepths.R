
# Define a function that determines the reads depths for the extended peaks
DeterminePeakDepths <-function(bedtools,bedtoolsgenome,extendreads,tagAlignfile,
                               tagAlignfileExt,peaksfileslop,temppeakoverlap)
{
  # overlap the peak and tagAlign files
  cmd1 = paste(bedtools," slop -l 0 -r ", extendreads ," -s -i ",tagAlignfile, " -g ",bedtoolsgenome, " > ", tagAlignfileExt,sep="")
  system(cmd1)
  cmd2 = paste(bedtools, " intersect -wo -a ",tagAlignfileExt," -b " , peaksfileslop, " > " ,temppeakoverlap ,sep = "")
  system(cmd2)
  
  # make a new peak file with depth info
  DeterminePeakDepthsC(temppeakoverlap,peaksfileslopdepth)
}