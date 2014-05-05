readParameters <-function(file="mango.default.parameters.txt",args=hash(),verbose=FALSE)
{
  lines = readLines(file)
  
  for (line in lines)
  {
    # remove spaces
    line = gsub(pattern=" ",x=line,replace="")
    
    if (line == "")
    {
      next
    }
    else if (strsplit(line,split="")[[1]][1] == "#" )
    {
      next
    }
    else
    {
      lineinfo = strsplit(line,split="#")[[1]][1]
      arginfo  = strsplit(lineinfo,split="=")[[1]]
      args[[arginfo[1]]] = arginfo[2]
      if (verbose == TRUE)
      {
        print(paste(arginfo[1] , "        " , arginfo[2]))
      }
    }
  }
  return (args)
}
