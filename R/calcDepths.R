
# Define a function that calculates depths of peaks in interactions
calcDepths <- function(df,type="product")
{
  if (type == "product")
  {
    depths = df[,1] * df[,2]
  }
  if (type == "sum")
  {
    depths = df[,1] + df[,2]
  }
  return(depths)
}
