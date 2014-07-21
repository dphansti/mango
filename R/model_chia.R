
# define a function that calculates values per bin
model_chia <- function(x,y=NA,borders)
{
    sumofy  = c()
    meanofx = c()
    sumofx  = c()
    countofx  = c()
    bin = findInterval(x, borders)
    for (b in 0:length(borders))
    {
      binindexes = which(bin == b)
      if (is.na(y) == TRUE)
      {
        sumofy =  c(sumofy,length(binindexes))
      }
      if (is.na(y) == FALSE)
      {
        sumofy  = c(sumofy,  sum( y[binindexes]))
      }
      meanofx = c(meanofx, mean(as.numeric(x[binindexes])))
      sumofx  = c(sumofx, sum(as.numeric(x[binindexes])))
      countofx = c(countofx,length(binindexes))
    }
    pvals = sumofy / sum(sumofy)
    
    return (cbind(meanofx,sumofy,pvals,sumofx,countofx))
}

a = rep(422664,20000)
b = mean(c(a,NA))
b
