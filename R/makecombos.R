
# Define a function that makes all possible combos of peaks
makecombos <- function(chrom,outname,mindist,maxdist)
{

  peaksfile    = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
  chrpeaks = read.table(peaksfile,header=F,sep="\t")
  names(chrpeaks) = c("chr","start","end","name","score","strand")
  
  # sort by distance
  chrpeaks = chrpeaks[order(chrpeaks[,2]),]
  
  npeaks = nrow(chrpeaks)
  n = npeaks^2

  # intialize dataframe  
  start1=numeric(n)
  end1=numeric(n)
  start2=numeric(n)
  end2=numeric(n)
  score1=numeric(n)
  score2=numeric(n)

  curlength = 0
  i = 1
  
  for (i in (1:npeaks))
  {
    #print (i)
    j = i + 1

    temp_score2=chrpeaks[j:npeaks ,5]
    temp_score1=rep(chrpeaks[i,5],length(temp_score2))
    temp_end2=chrpeaks[j:npeaks ,3]
    temp_end1=rep(chrpeaks[i,3],length(temp_end2))
    temp_start2=chrpeaks[j:npeaks ,2]
    temp_start1=rep(chrpeaks[i,2],length(temp_start2))
    
    dist = abs( (temp_start1 + temp_end1 ) / 2 - (temp_start2 + temp_end2 ) / 2  )
    keepers = which(dist < maxdist & dist> mindist)
  
    if( length(keepers) > 0)
    {
      first = curlength + 1
      curlength = curlength + length(keepers)
      start1[first:curlength]=temp_start1[keepers]
      end1[first:curlength]  =temp_end1[keepers]
      start2[first:curlength]=temp_start2[keepers]
      end2[first:curlength]  =temp_end2[keepers]
      score1[first:curlength]=temp_score1[keepers]
      score2[first:curlength]=temp_score2[keepers]
    }
  }
  
  chrpairs = data.frame(chr1=rep(chrom,curlength),start1=start1[1:curlength],end1=end1[1:curlength],
                        chr2=rep(chrom,curlength),start2=start2[1:curlength],end2=end2[1:curlength],
                        score1=score1[1:curlength],score2=score2[1:curlength])
  
  return(chrpairs)
}

