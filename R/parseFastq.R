#' trims linkers
#'
#' This function find and trims linkers from fastq files and output 'same' and 'chim' files
#'  
#' @param fastqsin fastq files to parse
#' @param fastqsout full path of output files
#' @param linkers vector of length 2 containing the linkers
#' @param minlength minlength of read after linker trimming. NULL for no min length.
#' @param maxlegnth maxlength of read after linker trimming. NULL for no max length.
#' @param keepempty boolean whether or not reads with no linker should be kept
#' @param linesatatime number of lines to read at a time
#' @export
#' @examples
#' 
#' fastqs = c("/Users/dougphanstiel/Dropbox/data/NH.K562_RAD21_K562_std_2.1_1.fastq.gz","/Users/dougphanstiel/Dropbox/data/NH.K562_RAD21_K562_std_2.1_2.fastq.gz")
#' parseFastq(fastqsin=fastqs,minlegnth = 15,maxlength = 25,keepempty=FALSE,linesatatime = 10000,linkers = c("GTTGGATAAG","GTTGGAATGT"))
#' 
fastqs = c("/Volumes/HD3/projects/newmango/data/NH.K562_RAD21_K562_std_2.1_1.head.fastq",
           "/Volumes/HD3/projects/newmango/data/NH.K562_RAD21_K562_std_2.1_2.head.fastq")

outdir = NULL
outname="mango"
fastqsin = fastqs
fastqsout = NULL
linkers = c("GTTGGATAAG","GTTGGAATGT")
minlength = 15
maxlength = 25
keepempty=FALSE
linesatatime = 1000

parseFastq <- function(fastqsin,outdir=NULL,outname="mango",
                       linkers = c("GTTGGATAAGATATCGC","GTTGGAATGTATATCGC"),
                       minlegnth = 15,maxlength = 25,
                       keepempty=TRUE,linesatatime = 10000)
{  
  
  # open the connection
  if (tail(strsplit(fastqs[1],split="\\.")[[1]], n=1) == "gz")
  {
    con1  <- gzfile(fastqs[1], open = "rb")
  }
  if (tail(strsplit(fastqs[1],split="\\.")[[1]], n=1) != "gz")
  {
    con1  <- file(fastqs[1], open = "r")
  }
  if (tail(strsplit(fastqs[2],split="\\.")[[1]], n=1) == "gz")
  {
    con2  <- gzfile(fastqs[2], open = "rb")
  }
  if (tail(strsplit(fastqs[2],split="\\.")[[1]], n=1) != "gz")
  {
    con2  <- file(fastqs[2], open = "r")
  }
  
  # open output files
  if (is.null(outdir) ==TRUE)
  {
    outdir = dirname(fastqs[1])
  }
  same1 = file(paste(outdir,"/",outname,"_1.same.fastq",sep=""), open = "w")
  same2 = file(paste(outdir,"/",outname,"_2.same.fastq",sep=""), open = "w")
  chim1 = file(paste(outdir,"/",outname,"_1.chim.fastq",sep=""), open = "w")
  chim2 = file(paste(outdir,"/",outname,"_2.chim.fastq",sep=""), open = "w")
  

  print(Sys.time())  
  i = 0
  while (length(lines1 <- readLines(con1, n = linesatatime, warn = FALSE)) > 0) {
    lines2 <- readLines(con2, n = linesatatime, warn = FALSE)
    
    # increment counter and report lines read
    i = i + 1
    linesread = i * linesatatime
    print (linesread)
    print(Sys.time())

    # convert lines to dataframe
    df1 = data.frame(matrix(lines1,ncol=4,byrow = TRUE))
    lines1 = NULL
    df2 = data.frame(matrix(lines2,ncol=4,byrow = TRUE))
    lines2 = NULL
    
    # add names
    names(df1) = c("info","sequence","third","score")
    names(df2) = names(df1)
  
    # trim sequences
    maxreadlength = maxlength + max(c(nchar(linkers[1]),nchar(linkers[2])))
    df1$sequence = unlist(lapply(df1$sequence,substr,start=1,stop=maxreadlength))
    df2$sequence = unlist(lapply(df2$sequence,substr,start=1,stop=maxreadlength))
    df1$score = unlist(lapply(df1$score,substr,start=1,stop=maxreadlength))
    df2$score = unlist(lapply(df2$score,substr,start=1,stop=maxreadlength))

    # find index of linkers
    df1$link1 = unlist(lapply(df1$sequence,regexpr,pattern=linkers[1]))
    df1$link2 = unlist(lapply(df1$sequence,regexpr,pattern=linkers[2]))
    df2$link1 = unlist(lapply(df2$sequence,regexpr,pattern=linkers[1]))
    df2$link2 = unlist(lapply(df2$sequence,regexpr,pattern=linkers[2]))

    # remove pairs with both linkers
    badones = which((df1$link1 > 0 & df1$link2 > 0) | (df2$link1 > 0 & df2$link2 > 0))
    if (length(badones) > 0)
    {
      df1 = df1[-badones,]
      df2 = df2[-badones,]
    }
    
    # set linker
    df1$linker = 0
    df2$linker = 0
    df1$linker[which(df1$link2 > df1$link1)] = 2
    df2$linker[which(df2$link2 > df2$link1)] = 2
    df1$linker[which(df1$link2 < df1$link1)] = 1
    df2$linker[which(df2$link2 < df2$link1)] = 1

    # set length
    df1$length = lapply(as.character(df1$sequence),nchar)
    df2$length = lapply(as.character(df2$sequence),nchar)
    df1$length[which(df1$link2 > df1$link1)] = df1$link2[which(df1$link2 > df1$link1)] - 1
    df1$length[which(df1$link2 < df1$link1)] = df1$link1[which(df1$link2 < df1$link1)] - 1
    df2$length[which(df2$link2 > df2$link1)] = df2$link2[which(df2$link2 > df2$link1)] - 1
    df2$length[which(df2$link2 < df2$link1)] = df2$link1[which(df2$link2 < df2$link1)] - 1
    
    # remove pairs with reasds that are too short or too long
    badones = which(df1$length < minlength | df1$length > maxlength |
                    df2$length < minlength | df2$length > maxlength)
    if (length(badones) > 0)
    {
      df1 = df1[-badones,]
      df2 = df2[-badones,]
    }
    
    # remove pairs with missing linkers (if desired)
    if (keepempty == FALSE)
    {
      badones = which(df1$linker == 0 | df2$linker == 0)
      if (length(badones) > 0)
      {
        df1 = df1[-badones,]
        df2 = df2[-badones,]
      }
    }
    
    # trim remaining reads
    df1$sequence = substr(df1$sequence, 1,df1$length)
    df2$sequence = substr(df2$sequence, 1,df2$length)
    df1$score = substr(df1$score, 1,df1$length)
    df2$score = substr(df2$score, 1,df2$length)
    
    # annotate as same, chom, or ambi
    df1$pair = "unknown"
    df1$pair[which(df1$linker == df2$linker)] = 'same'
    df1$pair[which(df1$linker != df2$linker)] = 'chim'
    
    # print same and chim
    writeLines(text=as.vector(t(as.matrix(df1[which(df1$pair == 'same'),(1:4)]))), con= same1)
    writeLines(text=as.vector(t(as.matrix(df2[which(df1$pair == 'same'),(1:4)]))), con= same2)
    writeLines(text=as.vector(t(as.matrix(df1[which(df1$pair == 'chim'),(1:4)]))), con= chim1)
    writeLines(text=as.vector(t(as.matrix(df2[which(df1$pair == 'chim'),(1:4)]))), con= chim2)
  } 
  
  # close connections
  close(con1)
  close(con2)
  close(same1)
  close(same2)
  close(chim1)
  close(chim2)
}