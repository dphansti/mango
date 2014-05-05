# 
# 
# # Variables
# chromosomes=c("chr1",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
#               "chr17", "chr18", "chr19" ,"chr2",  "chr20","chr21", "chr22", "chr3",  
#               "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9",  "chrX") 
# 
# bedpefilesortrmdup = "~/Desktop/mango_alan/NH.K562_RAD21_K562_std_2.1.sort.rmdup.bedpe"
# peaksfileslop      = "~/Desktop/mango_alan/NH.K562_RAD21_K562_std_2.1_peaks.slopPeak"
# outname            = "~/Desktop/mango_alan/NH.K562_RAD21_K562_std_2.1"
# numofbins          = 100
# mindist            = 20000
# maxdist            = 10000000
# binrange           = c(mindist, maxdist)
# 
# 
# 
# 
# 
# MDD    = calcMandD(chromosomes,outname,numofbins,binrange)
# IAB    = calcIAB(chromosomes,outname,numofbins,binrange)
# p_tab = data.frame(cbind(MDD,IAB))
# names(p_tab) = c("dist","M","D","IAB")
# p_tab$AvIAB = (p_tab$IAB/p_tab$M)
# averageDepth = sum(p_tab$D) / sum(p_tab$M)
# N = sum(p_tab$IAB)
# p_tab$p = (p_tab$AvIAB/N)
# spline = smooth.spline(p_tab$dist,p_tab$p,spar=.7)
# 
# # plot the spline
# plot(p_tab$dist,p_tab$p,pch=19)
# lines(spline, col = "blue")
# legend("topright",inset=0.05,legend=c("binned values","spline fit"),pch=c(19,NA),lwd=c(NA,1),col=c("black","blue"))
# 
# # score the interactions
# allpairs = scoreinteractions(chromosomes,outname,mindist=log10(mindist),maxdist=log10(maxdist),
#                              averageDepth=averageDepth,spline=spline,N=p_tab$N[1])
# 
# # Define a function that calculates P for each pair
# scoreinteractions <-function(chromosomes,outname ,mindist,maxdist,averageDepth,spline,N)
# {
#   # make bins
#   allpairs = c()
#   for (chrom in chromosomes)
#   {
#     print (chrom)
#     # get total intra IAB
#     pairsfile = paste(outname,".",chrom,".pairs.bedpe",sep="")
#     pairs = read.table(pairsfile,header=FALSE,sep="\t")
#     pairs = pairs[which(as.character(pairs$V8) != as.character(pairs$V9)),]
#     pairs$dist = log10(abs( (pairs[,3]+pairs[,2]/2) - (pairs[,6]+pairs[,5]/2) ) )
#     pairs = pairs[which(pairs$dist>mindist & pairs$dist<maxdist),]
# 
#     if (nrow(pairs) == 0)
#     {
#       next
#     }
#     
#     # calculate D (depth)
#     pairs$D = pairs[,10] + pairs[,11]
#     
#     # calculate p(success)
#     pairs$psuccess = predict(spline,pairs$dist)$y * (pairs$D / averageDepth)
# 
#     # add N
#     pairs$N = N
# 
#     allpairs=rbind(allpairs,pairs)
#   }
#   
#   print(head(allpairs))
#   
#   # do the actual P-value calculations
#   allpairs$P = apply(cbind(allpairs$V12,allpairs$N,allpairs$psuccess),1,calcP)    
#   allpairs$Q = p.adjust(allpairs$P,method="BY")
#   
#   return(allpairs)
# }
# 
# 
# # Define a function that calculates P for each pair
# calcIAB <-function(chromosomes,outname,numofbins ,binrange)
# {
#   # make bins
#   bins = seq(log10(binrange[1]),log10(binrange[2]),length.out=numofbins+1)
# 
#   # prepare genome wide counters of PETs for P(D|IB) and P(D|!IAB)
#   IABS                      = rep(0,numofbins)
#   names(IABS)               = (1:numofbins)
#   for (chrom in chromosomes)
#   {
#     print(chrom)
#     
#     # get total intra IAB
#     pairsfile = paste(outname,".",chrom,".pairs.bedpe",sep="")
#     pairs = read.table(pairsfile,header=FALSE,sep="\t")
#     pairs = pairs[which(as.character(pairs$V8) != as.character(pairs$V9)),]
#     pairs$dist = log10(abs( (pairs[,3]+pairs[,2]/2) - (pairs[,6]+pairs[,5]/2) ) )
#     pairs$bin  = findInterval(pairs$dist, bins)
#     pairs = pairs[which(pairs$bin>0 & pairs$bin<length(bins)),]
#     
#     if (nrow(pairs) == 0)
#     {
#       next
#     }
#     
#     for (eachbin in  (1:numofbins))
#     {
#       # add to M
#       IABSbin = length(which(pairs$bin == eachbin))
#       IABS[as.character(eachbin)] = IABS[as.character(eachbin)] + IABSbin
#     }
# 
#   }
#   
#   return(IABS)
# }
# 
# # Define a function that calculates P-values of interactions
# calcP <- function(v)
# {
#   P = 1 - pbinom(q=v[1]-1,size=v[2],prob=v[3])
#   return(P)
# }
# 
# # Define a function that calculates M (# combos) and D (average read depth: CA + CB) from the peak overlap files
# calcMandD <- function(chromosomes,outname,numofbins,binrange)
# {
#   
#   # make bins
#   bins = seq(log10(binrange[1]),log10(binrange[2]),length.out=numofbins+1)
#   length(bins)
#   
#   # prepare genome wide counters of PETs for P(D|IB) and P(D|!IAB)
#   readdepths                = rep(0,numofbins)
#   combos                    = rep(0,numofbins)
#   distances                 = rep(0,numofbins)
#   names(readdepths)         = (1:numofbins)
#   names(combos)             = (1:numofbins)
#   names(distances)          = (1:numofbins)
#   
#   for (chrom in chromosomes)
#   {
#     print (chrom)
#     
#     peaksfile    = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
#     
#     # read in data
#     chrpeaks = read.table(peaksfile,header=F,sep="\t")
#     
#     names(chrpeaks) = c("chr","start","end","name","score","strand")
#     
#     # calc center position
#     chrpeaks$pos = (chrpeaks$start + chrpeaks$end) /2
# 
#     for (i in (1:nrow(chrpeaks)))
#     {
#       #calc dist
#       chrpeaks$dist  = log10(abs(chrpeaks$pos - chrpeaks$pos[i]))
#       
#       #calc depth
#       chrpeaks$depth = chrpeaks$score + chrpeaks$score[i]
#       
#       # assign to bins
#       chrpeaks$bin  = findInterval(chrpeaks$dist, bins)
#       for (eachbin in  (1:numofbins))
#       {
#         # add to M
#         combosinbin = length(which(chrpeaks$bin == eachbin))
#         combos[as.character(eachbin)] = combos[as.character(eachbin)] + combosinbin
#         
#         # add to D
#         readdepthsbin = sum(chrpeaks$depth[(which(chrpeaks$bin == eachbin))])
#         readdepths[as.character(eachbin)] = readdepths[as.character(eachbin)] + readdepthsbin
#         
#         # add to distances
#         distsbin = sum(10^chrpeaks$dist[(which(chrpeaks$bin == eachbin))])
#         distances[as.character(eachbin)] = distances[as.character(eachbin)] + distsbin
#       } 
#     }
#   }
#   
#   # convert distances to averages
#   distances = log10(distances/combos)
#   
#   return(cbind(distances,combos,readdepths))
# }
# 
# 
#   
# # Define a function that calculates N (# linking PETs per bin)  from the pair files
# calcN <- function(chromosomes,outname,numofbins,binrange)
# {
#   # make bins
#   bins = seq(log10(binrange[1]),log10(binrange[2]),length.out=numofbins)
#   
#   # prepare genome wide counters of PETs for P(D|IB) and P(D|!IAB)
#   linkingPETs               = rep(0,numofbins+1)
#   names(linkingPETs)         = (0:numofbins)
#   
#   for (chrom in chromosomes)
#   {
#     print (chrom)
#     
#     # get pairs bins
#     pairs = read.table(paste(paste(outname,"." ,chrom, ".pairs.bedpe",sep="")),header=F,sep="\t")
#     pairs$dist = (pairs[,6]+pairs[,5])/2 - (pairs[,3]+pairs[,2])/2
#     pairsdists =   log10(pairs$dist [which(pairs$dist > 0)])
#     
#     # bin distances
#     pairbins =  findInterval(pairsdists, bins)
#     
#     for (eachbin in  (0:numofbins))
#     {
#       linkingPETs[as.character(eachbin)] = 
#         linkingPETs[as.character(eachbin)] + length(which(pairbins == eachbin))
#     }
#   }
#   return (linkingPETs)
# }
#   
#   
# 
#   
# 
# #allpairs =  allpairs[which(allpairs$dist > log10(20000) & allpairs$dist < log10(10000000) ),]
# allpairs$Q = p.adjust(allpairs$P,method="BY") 
# 
# par(mfrow=c(2,3))
# hist(allpairs$P)
# hist(allpairs$Q)
# #sig= allpairs[which(allpairs$Q < 0.01 & allpairs$V12 > 1),]
# sig= allpairs[which(allpairs$Q < 0.01 & allpairs$V12 > 0),]
# nrow(sig)
# 
# plot(density(sig$dist),main="Q < 0.01",xlab="log10(dist)")
# hist(sig$V12[which(sig$V12<30)],breaks=seq(0,30),main="Q < 0.01",xlab="# PETs")
# 
# color.paletteX = colorRampPalette(c("black","deepskyblue"))
# col=color.paletteX(10)
# plot(density(allpairs[which(allpairs$Q < 0.01 & allpairs$V12 ==2),]$dist),col[1],xlim=c(4,7))
# for (i in (3:10))
# {
#   lines(density(allpairs[which(allpairs$Q < 0.01 & allpairs$V12 ==i),]$dist),col=col[i])
# }
# 
# hist(allpairs$V12[which(allpairs$V12<30 )],breaks=seq(0,30),main=10^bins[bin],xlab="# PETs",col="red",ylim=c(0,5000))
# hist(sig$V12[which(sig$V12<30)],breaks=seq(0,30),main="Q < 0.01",xlab="# PETs",add=TRUE,col="skyblue2")
