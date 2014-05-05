# require('hash')
# 
# ########################################################################
# #                                                                      #
# #                               FUNCTIONS                              # 
# #                                                                      #
# ########################################################################
# 
# 
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
#   # prepare genome wide counters of IABs and peakcombos for P(IAB) and P(!IAB) 
#   totalIAB = 0
#   all_poss_peak_pairs = 0
#   peaks = read.table(peaksfileslop,header=FALSE,sep="\t")
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
#     #---- Collect info for P(IAB) and P(!IAB) ----#
#     # get total pairs of peaks info
#     n = nrow(peaks[which(peaks$V1 == chrom),])
#     combos =  (n * n-1 ) / 2
#     all_poss_peak_pairs = all_poss_peak_pairs + combos
#     
#     # get total intra IAB
#     pairsfile = paste(outname,".",chrom,".pairs.bedpe",sep="")
#     pairs = read.table(pairsfile,header=FALSE,sep="\t")
#     pairs = pairs[which(as.character(pairs$V8) != as.character(pairs$V9)),]
#     totalIAB = totalIAB + sum(pairs$V12)
#   } 
#   
#   # now put together table of info
#   p_df = data.frame(bin=1:numofbins,distance=bins,
#                     obspets  = obs_binned_pets[2:length(obs_binned_pets)],
#                     rwrpets  = rwr_binned_pets[2:length(rwr_binned_pets)],
#                     pD_g_IAB = obs_binned_pets[2:length(obs_binned_pets)] / sum(obs_binned_pets[2:length(obs_binned_pets)]),
#                     pD_n_IAB = rwr_binned_pets[2:length(rwr_binned_pets)] / sum(rwr_binned_pets[2:length(rwr_binned_pets)]),
#                     pIAB     = totalIAB / all_poss_peak_pairs,
#                     pnIAB    = 1 - totalIAB / all_poss_peak_pairs)
#     
#   # calculate the p(success)
#   p_df$psuccess_part = apply(  p_df[,5:8] , 1 , psuccess_calc,full=FALSE)
#   p_df$psuccess_full = apply(  p_df[,5:8] , 1 , psuccess_calc,full=TRUE)
#     
#   return(p_df)
# }
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
#     pairs$psuccess_part = ptable$psuccess_part[pairs$bin]
#     pairs$psuccess_full = ptable$psuccess_full[pairs$bin]
#     allpairs=rbind(allpairs,pairs)
#   }
# 
#   # do the actual P-value calculations
#   allpairs$size = allpairs$V10 + allpairs$V11 - allpairs$V12
#   allpairs$P_part = apply(cbind(allpairs$V12,allpairs$size,allpairs$psuccess_part),1,calcP)
#   allpairs$P_full = apply(cbind(allpairs$V12,allpairs$size,allpairs$psuccess_full),1,calcP)
#   allpairs$Q_part = p.adjust(allpairs$P_part,method="BH")
#   allpairs$Q_full = p.adjust(allpairs$P_full,method="BH")
#   
#   return(allpairs)
# }
# 
# ########################################################################
# #                                                                      #
# #                               CODE                                   # 
# #                                                                      #
# ########################################################################
# 
# 
# # Variables
# chromosomes=c("chr1",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
#                   "chr17", "chr18", "chr19" ,"chr2",  "chr20","chr21", "chr22", "chr3",  
#                   "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9",  "chrX","chrY","chrM") 
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
# # calculate the p(success)
# ptable$psuccess_part = apply(  ptable[,5:8] , 1 , psuccess_calc,full=FALSE)
# ptable$psuccess_full = apply(  ptable[,5:8] , 1 , psuccess_calc,full=TRUE)
# 
# 
# ptable$psuccess_part = ptable$pD_g_IAB * ptable$psuccess_full
# 
# # calculate pvalues for interactions
# allpairs = scoreinteractions(chromosomes,
#                              outname   = outname,
#                              numofbins = numofbins,
#                              binrange= binrange,
#                              ptable)
# 
# 
# 
# #pdf("~/Dropbox/mangoRpackage/extra/probabilitymethod.pdf")
# par(mfrow=c(2,3))
# 
# plot(ptable$distance,ptable$pD_g_IAB,type='b',col="dodgerblue2",pch=19,xlab="log10(dist)",ylab="prob",ylim=c(0,.10))
# points(ptable$distance,ptable$pD_n_IAB,type='b',col="firebrick2",pch=19)
# legend("top",legend=c("p(D|IAB)","p(D|!IAB)"),bty='n',ncol=2,pch=19,col=c("dodgerblue2","firebrick2"))
# 
# plot(ptable$distance,ptable$psuccess_part,pch=19,xlab="log10(dist)",ylab="prob",main="p(D|IAB) * p(IAB)")
# legend("topright",legend=paste("p(IAB) =",signif(ptable$pIAB[1],digits=5)))
# plot(ptable$distance, ptable$psuccess_full, pch=19,xlab="log10(dist)",ylab="prob",main="p(D|IAB) * p(IAB) / p(D)")
# legend("topright",legend=paste("p(IAB) =",signif(ptable$pIAB[1],digits=5)))
# 
# plot(1, type="n", axes=F, xlab="", ylab="")
# hist(allpairs$P_part,main="p(D|IAB) * p(IAB)")
# hist(allpairs$P_full,main="p(D|IAB) * p(IAB) / p(D)")
# #dev.off()
# 
# 
# 
# suball = allpairs[which(allpairs$dist > log10(20000)),]
# suball$Q_part = p.adjust(suball$P_part,method="BH")
# 
# length(which(suball$Q_part < 0.001))
# sig = suball[(which(suball$Q_part < 0.001)),]
# head (sig)
# plot(density(sig$dist))
# hist(sig$V12[which(sig$V12<30)],breaks=seq(0,30))
# 
# 
# 
# 
