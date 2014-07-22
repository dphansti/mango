
# Define a function that calculates depths of peaks in interactions
calcDepths <- function(df,type="product")
{
  if (type == "product")
  {
    depths = as.numeric(df[,1]) * as.numeric(df[,2])
  }
  if (type == "sum")
  {
    depths = as.numeric(df[,1]) + as.numeric(df[,2])
  }
  return(depths)
}
