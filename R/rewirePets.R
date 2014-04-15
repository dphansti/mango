
rewire <- function(readsfile,petsfile,obspetsfile,rwrpetsfile,obsreadsfile,rwrreadsfile,reps = 10,counter=0)
{
  bedpe     = read.table(petsfile,header=FALSE,sep="\t")
  # remove interchroms
  bedpe = bedpe[which(as.character(bedpe[,1]) == as.character(bedpe[,4] )),]
  
  # calculate distance
  bedpe$dis = log10(abs(bedpe[,2] - bedpe[,5]))
  
  # get number of PETs
  numPETS = length(bedpe$dis)
  
  # remove distance = 0
  bedpe     = bedpe[which(bedpe$dis > 0),]
  
  # make the distribution
  distribution = logspline(bedpe$dis,lbound=min(bedpe$dis),ubound=max(bedpe$dis))
  
  # randomly sample from the distribution
  randomvals  =  rlogspline(numPETS,distribution)
  
  # get reads for this chrom
  reads     = read.table(readsfile,header=FALSE,sep="\t")
  
  # count reads
  numreads = nrow(reads)
  
  # reorder
  reads = reads[order(reads[,2]),]
  
  # determine the average separation
  avg = (max(reads[,c(2,3)]) - min(reads[,c(2,3)])) / numreads

  # add an index
  reads$pos1 = (1:numreads)

  # generate random pairs (rep # of times)    
  randpair =c()
  i = 1
  for (i in (1:reps))
  {
    # pick a bunch of reads to starts with 
    randomdata = reads[sample((1:numreads),numPETS),] 

    names(randomdata) = c("chr","start","stop","name","score","strand","pos1")
    randomdata$randdist = sample(randomvals)
    
    # add random pos shift
    randomdata$shift = ceiling((10^randomdata$randdist) / avg)
    
    # find pos2
    randomdata$pos2 = randomdata$shift + randomdata$pos1
    randomdata$pos2[which(randomdata$pos2 > numreads)] =  randomdata$pos1[which(randomdata$pos2 > numreads)] - randomdata$shift[which(randomdata$pos2 > numreads)]
    randomdata$pos2[which(randomdata$pos2 <= 0)] = sample((1:numreads),1)
    
    # make the pair
    randpairsmall = cbind(reads[randomdata$pos1,(1:3)], reads[randomdata$pos2,(1:3)],reads[randomdata$pos1,4],".",reads[randomdata$pos1,6],reads[randomdata$pos2,6],randomdata$randdist )
    randpair = rbind(randpair,randpairsmall)
  }

  # correct names
  names(randpair) = c("chr1","start1","stop1","chr2","start2","stop2","name","score","strand1","strand2","targetdist")
  
  # calculate the dist
  randpair$dist = log10(abs(randpair[,2] - randpair[,5]))
  
  # find the difference in distance
  randpair$dif = abs(randpair$dist - randpair$targetdist)
  
  # sort the data
  randpair.sort = randpair[order(randpair$dif,decreasing=FALSE),]
  
  # remove infinites
  randpair.sort = randpair.sort[which(is.finite(randpair.sort$dif) == TRUE),]
  
  # keep first of each randdist
  randpair.sort.uniq = randpair.sort[!duplicated(randpair.sort$targetdist), ]
  chromosome = sub(pattern = "chr",replacement="",x= randpair.sort.uniq[1,1])

  densities = list()
  densities[[1]] = density(bedpe$dis)
  densities[[2]] = density(randpair.sort.uniq$dist)
  densities[[3]] = chromosome
  print(paste("Chromosome", chromosome))

  randpair.sort.uniq[,7] = paste("rwr",( 1:nrow(randpair.sort.uniq) + counter),sep="_")
  bedpe[,7]              = paste("obs",( 1:nrow(bedpe              ) + counter),sep="_")
  
  # write to new bedpe files
  write.table(randpair.sort.uniq[,1:10],file=rwrpetsfile,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
  write.table(bedpe[,1:10],file=obspetsfile,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
  
  # write to new bed files
  write.table(randpair.sort.uniq[,c(1,2,3,7,8,9)],file=rwrreadsfile,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
  write.table(randpair.sort.uniq[,c(4,5,6,7,8,10)],file=rwrreadsfile,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
  
  write.table(bedpe[,c(1,2,3,7,8,9)],file=obsreadsfile,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
  write.table(bedpe[,c(4,5,6,7,8,10)],file=obsreadsfile,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
  
  densities[[4]] = nrow(bedpe) + counter
  
  return (densities)
}
