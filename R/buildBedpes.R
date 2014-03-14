#' combines 2 same files into bepde format for every chromosome
#'
#' This function combines 2 same files into bepde format for every chromosome
#'  
#' @param samfilt1 filtered sam file
#' @param samfilt2 filtered sam file
#' @param bedpebase basename for bedpe file
#' @param linesatatime number o flines to read into memory at a time
#' @export

buildBedpes <- function(samfilt1,samfilt2,bedpebase,linesatatime=10000)
{
  unpaired1 = c()
  con1  <- file(samfilt1, open = "r")
  con2  <- file(samfilt2, open = "r")
  i = 0
  while (length(lines1 <- readLines(con1, n = linesatatime, warn = FALSE)) > 0)
  {
    i = i + 1
    print (i)
    lines2 <- readLines(con2, n = linesatatime, warn = FALSE)
    lines1 = unlist(lapply(lines1,strsplit,split="\t"))
    lines1 = matrix(lines1,ncol=7,byrow=TRUE)
    lines2 = unlist(lapply(lines2,strsplit,split="\t"))
    lines2 = matrix(lines2,ncol=7,byrow=TRUE)
    
    # make dataframe
    df1 = data.frame(chr = lines1[,3],start = as.numeric(lines1[,4]), 
                     stop =  as.numeric(lines1[,4]) + unlist(lapply(lines1[,5],nchar)),
                     strand = lines1[,2])
    df2 = data.frame(chr = lines2[,3],start = as.numeric(lines2[,4]), 
                     stop =  as.numeric(lines2[,4]) + unlist(lapply(lines2[,5],nchar)),
                     strand = lines2[,2])
    
    # intialize new PET info
    read1 = data.frame(chr=c(),start=c(),stop=c(),strand=c())
    read2 = data.frame(chr=c(),start=c(),stop=c(),strand=c())
    
    # reorder interchrom
    read1      = rbind(read1,df1[which(as.character(df1$chr) < as.character(df2$chr)),])
    read2      = rbind(read2,df2[which(as.character(df1$chr) < as.character(df2$chr)),])
    read2      = rbind(read2,df1[which(as.character(df1$chr) > as.character(df2$chr)),])
    read1      = rbind(read1,df2[which(as.character(df1$chr) > as.character(df2$chr)),])
    
    # reorder intrachrom
    read1      = rbind(read1,df1[which(as.character(df1$chr) == as.character(df2$chr) &
                                         as.numeric(df1$start) <= as.numeric(df2$start)),])
    read2      = rbind(read2,df2[which(as.character(df1$chr) == as.character(df2$chr) &
                                         as.numeric(df1$start) <= as.numeric(df2$start)),])
    read2      = rbind(read2,df1[which(as.character(df1$chr) == as.character(df2$chr) &
                                         as.numeric(df1$start) > as.numeric(df2$start)),])
    read1      = rbind(read1,df2[which(as.character(df1$chr) == as.character(df2$chr) &
                                         as.numeric(df1$start) > as.numeric(df2$start)),])
    
    names(read1) = paste(names(read1),"1",sep="")
    names(read2) = paste(names(read2),"2",sep="")
    
    bedpetemp = data.frame(cbind(read1[,1:3],read2[,1:3],".",".",read1[,4],read2[,4]))
    names(bedpetemp)[7:10] = c("name","score","strand1","strand2")
    
    for (chrom in names(table(bedpetemp[,1])))
    {
      outname = paste(bedpebase,".",chrom,".raw.bedpe",sep="")
      chromdata = bedpetemp[which(bedpetemp[,1] == chrom),]
      if (nrow(chromdata) > 0)
      {
        write.table(chromdata,file=outname,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t") 
      }
    }    
  }
  close(con1)
  close(con2)
}