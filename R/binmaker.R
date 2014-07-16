
# Define a function that establishes bin borders for use with findInterval later
binmaker <- function (vectortobin,binmethod="equalocc",numberbins=10)
{
  
  # sort vector
  sortedvec = sort(vectortobin)
  
  if (binmethod == "equalocc")
  {    
    # determine the number of items per bin
    itemsperbin = length(vectortobin)/numberbins
    
    # bin the data
    splitdata = split(sortedvec, ceiling(seq_along(sortedvec)/itemsperbin))
    
    # get max and min of each
    mins  = c()
    maxes = c()
    for (i in (1:length(splitdata)))
    {
      mins  = c(mins,min(splitdata[[i]]))
      maxes = c(maxes,max(splitdata[[i]]))
    }
    
    # set borders
    borders = (mins[-1] + maxes[-length(maxes)]) /2 
    
    # return
    return (borders)
  }
  
  if (binmethod == "log2")
  {
    log10sortedvec = log2(sortedvec)
    borders = seq(min(log10sortedvec),max(log10sortedvec),length.out=numberbins+1)[2:numberbins]
    
    # return
    return (borders)
  }
  
  if (binmethod == "log10")
  {
    log10sortedvec = log10(sortedvec)
    borders = seq(min(log10sortedvec),max(log10sortedvec),length.out=numberbins+1)[2:numberbins]
    
    # return
    return (borders)
  }
  
  if (binmethod == "equalsize")
  {
    borders = seq(min(sortedvec),max(sortedvec),length.out=numberbins+1)[2:numberbins]
    
    # return
    return (borders)
  }
}
