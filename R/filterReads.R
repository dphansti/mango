

filterReads <- function(sam1,sam2,linesatatime=10000)
{
  print ("finding matched pairs")
  # go through files and find aligned pairs
  unpaired1 = c()
  con  <- file(sam1, open = "r")
  while (length(lines <- readLines(con, n = linesatatime, warn = FALSE)) > 0)
  {
    lines = unlist(lapply(lines,strsplit,split="\t"))
    lines = matrix(lines,ncol=7,byrow=TRUE)
    unpaired1 = c(unpaired1,lines[which(lines[,7] == 0),1])
  }
  close(con)
  
  unpaired2 = c()
  con  <- file(sam2, open = "r")
  while (length(lines <- readLines(con, n = linesatatime, warn = FALSE)) > 0)
  {
    lines = unlist(lapply(lines,strsplit,split="\t"))
    lines = matrix(lines,ncol=7,byrow=TRUE)
    unpaired2 = c(unpaired2,lines[which(lines[,7] == 0),1])
  }
  close(con)
  
  # remove end of header line
  unpaired1 = matrix(unlist(lapply(unpaired1,strsplit,split="_")),ncol=2,byrow=TRUE)[,1]
  unpaired2 = matrix(unlist(lapply(unpaired2,strsplit,split="_")),ncol=2,byrow=TRUE)[,1]
  
  # get intersection
  paired = intersect(unpaired1, unpaired2)
  
  # go through files again and print out aligned pairs
  print ("printing matched pairs")
  for (filename in c(sam1,sam2))
  {
    con     <- file(filename, open = "r") 
    outname = paste(strsplit(filename,split="\\.")[[1]][1:length(strsplit(filename,split="\\.")[[1]])-1],collapse=".")
    outname = paste(outname,".filt.sam",sep="") 
    conout  <- file(outname, open = "w")
    while (length(lines <- readLines(con, n = linesatatime, warn = FALSE)) > 0)
    {
      lines = unlist(lapply(lines,strsplit,split="\t"))
      lines = matrix(lines,ncol=7,byrow=TRUE)
      lines[,1] =  matrix(unlist(lapply(lines[,1],strsplit,split="_")),ncol=2,byrow=TRUE)[,1]
      lines = lines[which(lines[,1] %in% paired),]
      write.table(lines,conout,quote = FALSE,sep="\t",row.names = FALSE,col.names = FALSE)
    }
    close(con)
    close(conout)
  }
}