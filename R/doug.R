# # Define a function that computes p of success
# psuccess_calc <- function(v,full=FALSE)
# {
#   pD_g_IAB = v[1]
#   pD_n_IAB = v[2]
#   PofIAB   = v[3]
#   PnIAB    = v[4]
#   
#   if (full == FALSE){ p_success = pD_g_IAB * PofIAB}
#   if (full == TRUE){p_success = pD_g_IAB * PofIAB / ( pD_g_IAB * PofIAB + pD_n_IAB*PnIAB )}
#   # if (full == TRUE){p_success = pD_g_IAB * PofIAB / ( pD_g_IAB * PofIAB + (1/50)*PnIAB )}
#   # if (full == TRUE){p_success = pD_g_IAB * PofIAB / (1/50)}
#   #print( paste(pD_g_IAB * PofIAB, "    ", pD_n_IAB*PnIAB ))
#   
#   return(p_success)
# }
# 
# 
# 
# 
# # Define a function that calculates P-values of interactions
# calcP <- function(v)
# {
#   P = 1 - pbinom(q=v[1]-1,size=v[2],prob=v[3])
#   return(P)
# }
# 
# 
# 
# 
# 
# # Define a function that calculates first term
# psuccess <- function(chromosomes,outname,bedpefilesortrmdup,distancefile,numofbins,binrange,peaksfileslop)
# {
#   
#   # split reads by chromosome
#   #petschroms  =      splitBedpe(bedpefilesortrmdup, outname, printreads=FALSE)[2]
#   
#   # make bins
#   bins = seq(log10(binrange[1]),log10(binrange[2]),length.out=numofbins)
#   
#   # prepare genome wide counters of PETs for P(D|IB) and P(D|!IAB)
#   obs_binned_pets           = rep(0,numofbins+1)
#   rwr_binned_pets           = rep(0,numofbins+1)
#   names(obs_binned_pets)    = (0:numofbins)
#   names(rwr_binned_pets)    = (0:numofbins)
#   
#   # prepare genome wide counters  peakcombos for P(IAB|D)
#   obs_binned_peak_combos           = rep(0,numofbins+1)
#   names(obs_binned_peak_combos)    = (0:numofbins)
#   peaks = read.table(peaksfileslop,header=FALSE,sep="\t")
#   
#   IABsbinned  = rep(0,numofbins+1)
#   names(IABsbinned)    = (0:numofbins)
# 
#   for (chrom in chromosomes)
#   {
#     
#     print(chrom)
#     
#     #---- Collect info for P(D|IB) and P(D|!IAB) ----#
#     
#     petsfile    = paste(outname,"." ,chrom, ".bedpe",sep="")
#     
#     # read in pets file and make rewired pets
#     obspets = read.table(petsfile,header=F,sep="\t")
#     
#     rwrpets = rbind( obspets[,1:3] , setNames( obspets[,4:6] , names( obspets[,1:3]) ) )[sample(2*nrow(obspets)),]
#     rwrpets = cbind(rwrpets[1:nrow(obspets),], rwrpets[(nrow(obspets)+1):(nrow(obspets)*2),] )
#     
#     # calc distances
#     obspets$dist = (obspets[,6]+obspets[,5])/2 - (obspets[,3]+obspets[,2])/2
#     rwrpets$dist = (rwrpets[,6]+rwrpets[,5])/2 - (rwrpets[,3]+rwrpets[,2])/2
#     
#     # filter distances
#     obslogdists = log10(obspets$dist [which(obspets$dist > 0)])
#     rwrlogdists = log10(rwrpets$dist [which(rwrpets$dist > 0)])
#     
#     # bin distances
#     obsbincounts =  findInterval(obslogdists, bins)
#     rwrbincounts =  findInterval(rwrlogdists, bins)
#     
#     # collect binned pets info
#     for (eachbin in  (0:numofbins))
#     {
#       obs_binned_pets[as.character(eachbin)] = 
#         obs_binned_pets[as.character(eachbin)] + length(which(obsbincounts == eachbin))
#       rwr_binned_pets[as.character(eachbin)] = 
#         rwr_binned_pets[as.character(eachbin)] + length(which(rwrbincounts == eachbin))
#     }
#     
#     #---- Collect P(IAB|D) ----#
#     
#     # read in data
#     chrpeaks = peaks[which(peaks$V1 == chrom),]
# 
#     names(chrpeaks) = c("chr","start","end","score")
#     
#     # calc center position
#     chrpeaks$pos = (chrpeaks$start + chrpeaks$end) /2
#     
#     for (i in (1:nrow(chrpeaks)))
#     {
#       #calc dist
#       chrpeaks$dist = log10(abs(chrpeaks$pos - chrpeaks$pos[i]))
#       
#       # assign to bins
#       chrpeaks$bin  = findInterval(chrpeaks$dist, bins)
#       for (eachbin in  (0:numofbins))
#       {
#         valsfromeachbin = length(which(chrpeaks$bin == eachbin))
#         obs_binned_peak_combos[as.character(eachbin)] = obs_binned_peak_combos[as.character(eachbin)] + valsfromeachbin
#       } 
#     }
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
#       IABsbinned[as.character(eachbin)] = 
#         IABsbinned[as.character(eachbin)] + length(which(pairbins == eachbin))
#     }
#   } 
#   
#   # now put together table of info
#   p_df = data.frame(bin=1:numofbins,distance=bins,
#                     obspets     = obs_binned_pets[2:length(obs_binned_pets)],
#                     rwrpets     = rwr_binned_pets[2:length(rwr_binned_pets)],
#                     pD          = obs_binned_pets[2:length(obs_binned_pets)] / sum(obs_binned_pets[2:length(obs_binned_pets)]),
#                     peakcombos  = obs_binned_peak_combos[2:length(obs_binned_peak_combos)],
#                     IABsbinned  = IABsbinned[2:length(IABsbinned)]   )
# 
#   return(p_df)
# }
# 
# 
# scoreinteractions <-function(chromosomesless,outname,numofbins ,binrange,ptable)
# {
#   
#   # add numbers to all pairs
#   numofbins = numofbins
#   bins = seq(log10(binrange[1]),log10(binrange[2]),length.out=numofbins)
#   allpairs = c()
#   for (chrom in chromosomesless)
#   {
#     print(chrom)
#     # get total intra IAB
#     pairsfile = paste(outname,".",chrom,".pairs.bedpe",sep="")
#     pairs = read.table(pairsfile,header=FALSE,sep="\t")
#     pairs = pairs[which(as.character(pairs$V8) != as.character(pairs$V9)),]
#     pairs$dist = log10(pairs$V13)
#     pairs$bin  = findInterval(pairs$dist, bins)
#     pairs = pairs[which(pairs$bin>0),]
#     
#     if (nrow(pairs) == 0)
#     {
#       next
#     }
#     
#     # and p(success)
#     pairs$psuccess = ptable$psuccess[pairs$bin]
#     allpairs=rbind(allpairs,pairs)
#   }
#   
#   # do the actual P-value calculations
#   allpairs$size = allpairs$V10 + allpairs$V11 - allpairs$V12
#   allpairs$P = apply(cbind(allpairs$V12,allpairs$size,allpairs$psuccess),1,calcP)    
#   allpairs$Q = p.adjust(allpairs$P,method="BH")
#   
#   return(allpairs)
# }
# 
# 
# ########################################################################
# #                                                                      #
# #                               CODE                                   # 
# #                                                                      #
# ########################################################################
# 
# a = read.table(peaksfileslop,header=F,sep="\t")
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
# numofbins          = 50
# binrange           = c(1.0e+03, 2.5e+08)
# 
# # calculate probabilities of success
# ptable = psuccess(chromosomes,outname=outname,
#                   bedpefilesortrmdup = bedpefilesortrmdup,
#                   numofbins =numofbins,
#                   binrange=binrange,
#                   peaksfileslop=peaksfileslop)
# 
# plot(ptable$distance,ptable$psuccess,type='l')
# 
# ptable$psuccess = ptable$pD * ptable$IABsbinned / ptable$peakcombos
# 
# # calculate pvalues for interactions
# allpairs = scoreinteractions(chromosomes,
#                              outname   = outname,
#                              numofbins = numofbins,
#                              binrange= binrange,
#                              ptable)
# 
# allpairs =  allpairs[which(allpairs$dist > log10(20000)),]
# allpairs$Q = p.adjust(allpairs$P,method="BH") 
# 
# 
# pdf("~/Dropbox/mangoRpackage/extra/probabilitymethod2.pdf")
# par(mfrow=c(2,3))
# hist(allpairs$P)
# hist(allpairs$Q)
# #sig= allpairs[which(allpairs$Q < 0.01 & allpairs$V12 > 1),]
# sig= allpairs[which(allpairs$Q < 0.001 & allpairs$V12 > 0),]
# nrow(sig)
# 
# plot(density(sig$dist),main="Q < 0.001",xlab="log10(dist)")
# hist(sig$V12[which(sig$V12<30)],breaks=seq(0,30),main="Q < 0.001",xlab="# PETs")
# 
# col=color.paletteX(10)
# plot(density(allpairs[which(allpairs$Q < 0.01 & allpairs$V12 ==2),]$dist),col[1],xlim=c(4,7))
# for (i in (3:10))
# {
#   lines(density(allpairs[which(allpairs$Q < 0.01 & allpairs$V12 ==i),]$dist),col=col[i])
# }
# 
# dev.off()
# 
# 
# # filter according to distance cutoff
# subpairs = allpairs[which(allpairs$dist > log10(20000)),]
# minval = min(subpairs$bin)
# bins = seq(log10(binrange[1]),log10(binrange[2]),length.out=numofbins)
# bindist  = signif(10^(bins[minval:(length(bins)-1)]),digits=2)
# 
# par(mfrow=c(3,4))
# 
# for (PETcount in (1:12))
# {
#   breakdown = c()
#   for (i in minval:49)
#   {
#     sub = subpairs[which(subpairs$bin == i),]
#     curbreakdown= c(length(which(sub$V12 == PETcount & sub$Q < 0.01)), length(which(sub$V12 == PETcount & sub$Q >= 0.01)))
#     breakdown = rbind(breakdown,curbreakdown)
#   }
#   total = sum(breakdown)
#   breakdown=breakdown/rowSums(breakdown)
#   
#   #breakdown[,1][which(is.na(breakdown[,1])==TRUE)] = 1
#   #breakdown[,2][which(is.na(breakdown[,2])==TRUE)] = 0
#   bp =barplot(t(breakdown),las=2,col=c("deepskyblue2","dodgerblue4"),xaxt='n',main=total)
#   positions = seq(1,length(bp),by=5)
#   axis(1,at=bp[positions],labels=bindist[positions],las=2)
# }
# 
# 
# 
# 
# 
# lengthsig = c()
# lengthall = c()
# ratio = c()
# for (PETcount in (1:30))
# {
#   lengthsig =  c(lengthsig,length(which(allpairs$V12 == PETcount & allpairs$V12 > 1 & allpairs$Q < 0.01)))
#   lengthall =  c(lengthall,length(which(allpairs$V12 == PETcount )))
# }
# ratio = lengthsig / lengthall
# bp =barplot(ratio,las=2,col="deepskyblue2")
# axis(side=1,at=bp,label=(1:30))
# 
# 
# p_bin    = 1 / 26560
# peak1    = 80
# peak2    = 50
# k        = 5
# N        = 130
# N_bin    = 10000
# sumpeaksavgbin = 100
# 
# p_p1p2 = (130 * 100) * p_bin
# 
# sum(p_p1p2) = sum(p_bin)
# 
# 
# ########################################################################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
