
# Define a function that calculates distance between regions in a bedpe file
calcDistance <- function(df,type="average")
{

  if (type == "outer")
  {
    pos1 = apply(df[,c(2,3,5,6)],1,min)
    pos2 = apply(df[,c(2,3,5,6)],1,max)
    distance = abs(pos2 -pos1)
  }
  if (type == "average")
  {
    pos1 = apply(df[,c(2,3)],1,mean)
    pos2 = apply(df[,c(5,6)],1,mean)
    distance = abs(pos2 -pos1)
  }
  if (type == "inner")
  {
    # find max start
    maxstart = apply(df[,c(2,5)],1,max)
    minstop  = apply(df[,c(3,6)],1,min)
    distance = maxstart - minstop
  }
  return(distance)
}
