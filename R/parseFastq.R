




parseFastq <- function(fastqsin,fastqsout,
                       linkers = c("GTTGGATAAGATATCGC","GTTGGAATGTATATCGC"),
                       range=c(15,25),keepempty=FALSE)
{  
  

  con1  <- file(fastqs[1], open = "r")
  con2  <- file(fastqs[2], open = "r")
  
  i = 0
  while (length(lines1 <- readLines(con1, n = 1000000, warn = FALSE)) > 0) {
    i = i + 1
    print(i)
    lines2 <- readLines(con2, n = 1000000, warn = FALSE)
    
    # convert lines to dataframe
    df1 = data.frame(matrix(lines1,ncol=4))
    lines1 = NULL
    df2 = data.frame(matrix(lines2,ncol=4))
    lines2 = NULL
    
    names(df1) = c("info","sequence","third","score")
    names(df2) = names(df1)
    
    # find index of linkers
    
    # trim reads
    
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
