




parseFastq <- function(fastqs,linkers,range)
  
  

con1  <- file(fastqs[1], open = "r")
con2  <- file(fastqs[2], open = "r")

i = 0
while (length(lines1 <- readLines(con1, n = 1000000, warn = FALSE)) > 0) {
  i = i + 1
  print(i)
  lines2 <- readLines(con2, n = 1000000, warn = FALSE))

  # convert lines to matrix
  

} 

close(con1)
close(con2)

########################## testing section ##########################

inputFile <- "/Users/dougphanstiel/Dropbox/Sushi/datasets/GWAS/ICBP-summary-Nature.csv"
con1  <- file(inputFile, open = "r")while (length(lines1 <- readLines(con1, n = 1000000, warn = FALSE)) > 0) {
  i = i + 1
  print(i)
  lines2 <- readLines(con2, n = 1000000, warn = FALSE))

# convert lines to matrix




} 
