#' removes duplicates from bedpe files

removeDups <- function(path,expname="mango")
{
  bedperawfiles = list.files(path = path, pattern = "*raw.bedpe")
  for (bedperawfile in bedperawfiles)
  {
    rawbedpe = read.table(paste(path,bedperawfile,sep=""),header=FALSE,sep="\t")
    bedpe = rawbedpe[!duplicated(rawbedpe[,c('V1','V2','V4','V5')]), ]
    if (nrow(bedpe) > 0)
    {
      write.table(bedpe,file= paste(path,expname,".",bedpe[1,1],".rmdup.bedpe",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t",append=TRUE) 
    }
  }  
}
