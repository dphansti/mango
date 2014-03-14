
rewire <- function(path=path,expname="mango",reps = 10)
{
  # make a pdf to print results
  pdf(paste(path,expname,".rewire.pdf",sep=""),height=10.5,width=10.5)
  par(mfrow=c(3,3))
  
  bedpefiles = list.files(path = path, pattern = "*rmdup.bedpe")
  for (bedpefile in bedpefiles)
  {
    bedpe     = read.table(paste(path,bedpefile,sep=""),header=FALSE,sep="\t")
    
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
    readfile = sub("rmdup.bedpe", "reads.bed", bedpefile)
    reads     = read.table(paste(path,readfile,sep=""),header=FALSE,sep="\t")
    
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
      print (i)
      
      # pick a bunch of reads to starts with 
      randomdata = reads[sample((1:numreads),numPETS),] 
      
      names(randomdata) = c("chr","start","stop","name","strand","pos1")
      randomdata$randdist = sample(randomvals)
      
      # add random pos shift
      randomdata$shift = ceiling((10^randomdata$randdist) / avg)
      
      # find pos2
      randomdata$pos2 = randomdata$shift + randomdata$pos1
      randomdata$pos2[which(randomdata$pos2 > numreads)] =  randomdata$pos1[which(randomdata$pos2 > numreads)] - randomdata$shift[which(randomdata$pos2 > numreads)]
      randomdata$pos2[which(randomdata$pos2 <= 0)] = sample((1:numreads),1)

      # make the pair
      randpairsmall = cbind(reads[randomdata$pos1,(1:3)], reads[randomdata$pos2,(1:3)],".",".",reads[randomdata$pos1,5],reads[randomdata$pos2,5],randomdata$randdist )
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
    head(randpair.sort.uniq)
    
    chromosome = sub(pattern = "chr",replacement="",x= randpair.sort.uniq[1,1])

    par(ann=F)
    plot(density(bedpe$dis),col='red',xaxt='n')
    axis(side=1,at=seq(1,8,by=1),labels = 10^seq(1,8,by=1),las=2)
    points(density(randpair.sort.uniq$dist),col='black',type='l')
    mtext("distance",side=1,line=3.5,font=2)
    mtext(paste("Chromosome", chromosome),side=3,line=1,font=1.0,cex=2)
    legend("topright",inset = .05,legend=c("observed","rewired"),col=c("red","black"),pch=15)
  }
  dev.off()
}
