

# Define a funtion that establishes all parameters
establishParameters <-function(argscmdline)
{
  args = hash()
  
  # determine if there is a argument file
  argfile = NULL
  for (arg in argscmdline)
  {
    argparts = strsplit(arg,split="=")[[1]]
    if (argparts[1] == "argfile")
    {
      argfile = argparts[2] 
    }
  }
  
  # first read in default parameters
  args = readParameters()
  
  # overwrite with supplied file
  if (is.null(argfile) == FALSE)
  {
    args = readParameters(argfile,args)
  }
  
  # overwrite with command line arguments
  for (arg in argscmdline)
  {  
    argparts = strsplit(arg,split="=")[[1]]
    args[[argparts[1]]] = argparts[2]
  }
  
  # overwrite with command line arguments
  args[["outname"]] = ""
  if (args[["outdir"]] == "NULL")
  {
    args[["outname"]] = file.path(getwd(),args[["prefix"]])
  }
  if (args[["outdir"]] != "NULL")
  {
    args[["outname"]] = file.path( args[["outdir"]] ,args[["prefix"]]) 
  }
  
  # correct stages
  stages = strsplit(args[["stages"]],split=":")[[1]]
  if (length(stages) == 2)
  {
    args[["stages"]] = seq(as.numeric(stages[1]),as.numeric(stages[2]))
  }
  if (length(stages) == 1)
  {
    args[["stages"]] = as.numeric(stages[1])
  }
  
  return (args)

}
