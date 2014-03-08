
fastqs = c("/Users/dougphanstiel/Desktop/mango2014test/r1.fastq","/Users/dougphanstiel/Desktop/mango2014test/r2.fastq")

fastqs = c("/Users/dougphanstiel/Desktop/mango2014test/NH.K562_RAD21_K562_std_2.1_1.fastq.gz","/Users/dougphanstiel/Desktop/mango2014test/NH.K562_RAD21_K562_std_2.1_2.fastq.gz")
linkers = c("GTTGGATAAGATATCGC","GTTGGAATGTATATCGC")
range=c(15,25),keepempty=FALSE

parseFastq <- function(fastqsin,fastqsout,
                       linkers = c("GTTGGATAAGATATCGC","GTTGGAATGTATATCGC"),
                       lengths=c(15,200),keepempty=TRUE,linesatatime = 10000)
{  
  
  # open the connection
  if (tail(strsplit(fastqs[1],split="\\.")[[1]], n=1) == "gz")
  {
    con1  <- gzfile(fastqs[1], open = "r")
  }
  else
  {
    con1  <- file(fastqs[1], open = "r")
  }
  if (tail(strsplit(fastqs[2],split="\\.")[[1]], n=1) == "gz")
  {
    con2  <- gzfile(fastqs[2], open = "r")
  }
  else
  {
    con2  <- file(fastqs[2], open = "r")
  }

  
  print(Sys.time())  
  i = 0
  while (length(lines1 <- readLines(con1, n = linesatatime, warn = FALSE)) > 0) {
    
    i = i + 1
    lines2 <- readLines(con2, n = linesatatime, warn = FALSE)

    head(lines1)
    # convert lines to dataframe
    df1 = data.frame(matrix(lines1,ncol=4,byrow = TRUE))
    lines1 = NULL
    df2 = data.frame(matrix(lines2,ncol=4,byrow = TRUE))
    lines2 = NULL
    
    names(df1) = c("info","sequence","third","score")
    names(df2) = names(df1)
    
    # find index of linkers
    df1$link1 = unlist(lapply(df1$sequence,regexpr,pattern=linkers[1]))
    df1$link2 = unlist(lapply(df1$sequence,regexpr,pattern=linkers[2]))
    df2$link1 = unlist(lapply(df2$sequence,regexpr,pattern=linkers[1]))
    df2$link2 = unlist(lapply(df2$sequence,regexpr,pattern=linkers[2]))
    
    linesread = i * linesatatime
    if (linesread  %/% 100000 == 0)
    {
      print (linesread)
      print(Sys.time())
    }
    
    # remove pairs with both linkers
    badones = which((df1$link1 > 0 & df1$link2 > 0) |
                               (df2$link1 > 0 & df2$link2 > 0))
    
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
    badones = which(df1$length < lengths[1] | df1$length > lengths[2] |
                    df2$length < lengths[1] | df2$length > lengths[2])
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
    
    # print same and chim
  } 
  
close(con1)
close(con2)
}



########################## testing section ##########################
# inputFile <- "/Users/dougphanstiel/Dropbox/Sushi/datasets/GWAS/ICBP-summary-Nature.csv"
# con  <- file(inputFile, open = "r")
# while (length(lines1 <- readLines(con, n = 1000000, warn = FALSE)) > 0) {
#   m = matrix(lines1,ncol=4)
#   lines1 = NULL
#   print(m[1:2,1:2])
# } 
# close(con)
# 
# 
# m = matrix(c("a","b","c","d"),ncol=2)
# m[1:2,1:2]
