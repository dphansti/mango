# Define a function that calculates P-values of interactions
calcP <- function(v)
{
  P = 1 - pbinom(q=v[1]-1,size=v[2],prob=v[3])
  return(P)
}
