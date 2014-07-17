
# Define a function that determines paths to requires programs
DefinePaths <- function(progs = c("bedtools","MACS2","bowtie"))
{
  syspath = strsplit( Sys.getenv(c("PATH")),split=":")[[1]]
  progpaths = c()
  for (prog in progs)
  {
    progpathfound = "notfound"
    for (p in syspath)
    {
      path = file.path(p,prog)
      if (file.exists(path) == TRUE)
      {
        progpathfound = path
        break
      }
    }
    progpaths = c(progpaths,progpathfound)    
  }
  return(progpaths)
}


