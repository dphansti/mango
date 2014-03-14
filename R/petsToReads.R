#' splits bedpe into reads

petsToReads <- function(path=path,expname="mango")
{
  bedpefiles = list.files(path = path, pattern = "*rmdup.bedpe")
  for (bedpefile in bedpefiles)
  {
    bedpe = read.table(paste(path,bedpefile,sep=""),header=FALSE,sep="\t")
    bedtmp = rbind ( bedpe[,c(1,2,3,7,9)] ,   setNames( bedpe[,c(4,5,6,8,10)] , names( bedpe[,c(1,2,3,7,9)] ) ) )
    
    # print to chrom reads files
    for (chrom in names(table(bedtmp[,1])))
    {
      outname = paste(path,expname,".",chrom,".reads.bed",sep="")
      chromdata = bedtmp[which(bedtmp[,1] == chrom),]
      if (nrow(chromdata) > 0)
      {
        write.table(chromdata,file=outname,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t") 
      }
    }  
  }
}
