# runs mango chia pet analysis pipeline
library("Rcpp")
library("hash")
library("mango")

##################################### initialization #####################################

argsscript = c(
"bowtiepath=/path/to/bowtie",
"bowtieref=/path/to/indexes/hg19",
"bedtoolspath=/path/to/bedtools",
"bedtoolsgenome=/path/to/human.hg19.genome",
"macs2path=/path/to/macs2"
)

print ("Starting mango ChIA PET analysis tool")
Sys.time()
set.seed(1)

##################################### read in arguments #####################################

# read in command line args
argscmdline <- commandArgs(trailingOnly = TRUE)

# gather all argum
args = establishParameters(argscmdline,argsscript)

print ("using the following parameters:")
print (args)
       
resultshash = hash()


##################################### parse fastqs #####################################

if (1 %in% args[["stages"]])
{
  checkRequired(args,c("fastq1","fastq2","outname","minlength","maxlength","keepempty","linkerA","linkerB"))
  
  # gather arguments
  outname         = args[["outname"]]
  fastq1=args[["fastq1"]]
  fastq2=args[["fastq2"]]
  basename = args[["outname"]]
  minlength = as.numeric(args[["minlength"]])
  maxlength = as.numeric(args[["maxlength"]]) 
  keepempty=eval(parse(text=args[["keepempty"]]))
  linker1=args[["linkerA"]]
  linker2=args[["linkerB"]]
  
  print ("finding linkers")
  parsingresults = parseFastq( fastq1=fastq1,
              fastq2=fastq2,
              basename = basename,
              minlength = minlength,
              maxlength = maxlength, 
              keepempty = keepempty,
              linker1=linker1,
              linker2=linker2)
  
  resultshash[["total PETs"]] = sum(parsingresults)
  resultshash[["same PETs"]] = parsingresults[1]
  resultshash[["chimeric PETs"]] = parsingresults[2]
  resultshash[["ambigious PETs"]] = parsingresults[3]
}
  
###################################### align reads #####################################

if (2 %in% args[["stages"]])
{
  checkRequired(args,c("outname","bowtiepath","bowtieref","shortreads"))
  
  # gather arguments
  outname         = args[["outname"]]
  bowtiepath      = args[["bowtiepath"]]
  bowtieref       = args[["bowtieref"]]
  shortreads      = args[["shortreads"]]
  
  print ("aligning reads")
  # filenames
  fastq1 = paste(outname ,"_1.same.fastq",sep="")
  fastq2 = paste(outname ,"_2.same.fastq",sep="")
  sam1   = paste(outname ,"_1.same.sam",sep="")
  sam2   = paste(outname ,"_2.same.sam",sep="")
  
  # align both ends of each PET
  alignBowtie(fastq=fastq1,output=sam1,bowtiepath=bowtiepath,bowtieref=bowtieref,shortreads)
  alignBowtie(fastq=fastq2,output=sam2,bowtiepath=bowtiepath,bowtieref=bowtieref,shortreads)
}

##################################### filter reads #####################################

if (3 %in% args[["stages"]])
{
  checkRequired(args,c("outname"))
  
  # gather arguments
  outname         = args[["outname"]]
  
  # filenames
  bedpefile          = paste(outname ,".bedpe",sep="")
  bedpefilesort      = paste(outname ,".sort.bedpe",sep="")
  bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")
  sam1 = paste(outname ,"_1.same.sam",sep="")
  sam2 = paste(outname ,"_2.same.sam",sep="")
  
  # build bedpe
  print ("building bedpe")
  if (file.exists(bedpefile)){file.remove(bedpefile)}
  buildBedpe(sam1 =sam1, sam2 = sam2, bedpefile = bedpefile);
  
  # sort bedpe
  print ("sorting bedpe")
  if (file.exists(bedpefilesort)){file.remove(bedpefilesort)}
  external_sort(bedpefile, bedpefilesort)
  
  # filter duplicates
  print ("removing PCR duplicates")
  if (file.exists(bedpefilesortrmdup)){file.remove(bedpefilesortrmdup)}
  rmdupresults = removeDupBedpe(bedpefilesort,bedpefilesortrmdup,renamePets=TRUE);
  resultshash[["duplicate PETs"]] = rmdupresults[1]
  resultshash[["nonduplicate PETs"]] = rmdupresults[2]
  resultshash[["interchromosomal PETs"]] = rmdupresults[3]
  resultshash[["intrachromosomal PETs"]] = rmdupresults[4]
}

##################################### call peaks #####################################

if (4 %in% args[["stages"]])
{
  checkRequired(args,c("outname","macs2path","MACS_pvalue",
                       "bedtoolspath","bedtoolsgenome","peakslop","peakinput"))
  
  # gather arguments
  outname         = args[["outname"]]
  macs2path       = args[["macs2path"]]
  MACS_pvalue     = args[["MACS_pvalue"]]
  bedtoolspath    = args[["bedtoolspath"]]
  bedtoolsgenome  = args[["bedtoolsgenome"]]
  peakslop        = args[["peakslop"]]
  peakinput       = args[["peeakinput"]]
    
  # filenames
  bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")
  tagAlignfile       = paste(outname,".tagAlign",sep="")
  peaksfile          = paste(outname,"_peaks.narrowPeak",sep="")
  peaksfileslop      = paste(outname,"_peaks.slopPeak",sep="")
  
  if (is.null(peakinput) == FALSE)
  {
    peaksfile = peakinput
  }
  
  if (is.null(peakinput) == TRUE)
  {
    print ("building tagAlign file")
    # reverse strands for peak calling
    buildTagAlign(bedpefilesortrmdup ,tagAlignfile )
    
    # call peaks 
    print ("calling peaks")
    callpeaks(macs2path=macs2path,tagAlignfile,outname,pvalue=MACS_pvalue,
              bedtoolspath=bedtoolspath,bedtoolsgenome=bedtoolsgenome,
              peakslop=peakslop)
  }
  
  # extend and merge peaks according to peakslop
  print ("extending peaks")
  peakcounts = extendpeaks(peaksfile,peaksfileslop,bedtoolspath=bedtoolspath,
              bedtoolsgenome=bedtoolsgenome,peakslop=peakslop)
  resultshash[["peaks"]] = peakcounts[1]
  resultshash[["mergedpeaks"]] = peakcounts[2]
}

##################################### group pairs #####################################

if (5 %in% args[["stages"]])
{
  checkRequired(args,c("outname","distcutrangemin","biascut","maxPval",
                       "numofbins","corrMethod","bedtoolspath","bedtoolsgenome",
                       "maxinteractingdist","FDR","minPETS","corrMethod",
                       "chrominclude","chromexclude","reportallpairs"))
  
  # gather arguments
  outname = args[["outname"]]
  distcutrangemin = as.numeric(args[["distcutrangemin"]])
  distcutrangemax = as.numeric(args[["distcutrangemax"]])
  bedtoolspath = args[["bedtoolspath"]]
  bedtoolsgenome    = args[["bedtoolsgenome"]]
  biascut = as.numeric(args[["biascut"]])
  maxinteractingdist = as.numeric(args[["maxinteractingdist"]])
  numofbins = as.numeric(args[["numofbins"]])
  FDR = as.numeric(args[["FDR"]])
  minPETS = as.numeric(args[["minPETS"]])
  corrMethod = args[["corrMethod"]]
  chrominclude      = args[["chrominclude"]]
  chromexclude      = args[["chromexclude"]]
  reportallpairs    = args[["reportallpairs"]]
  
  # filenames
  peaksfile          = paste(outname ,"_peaks.narrowPeak",sep="")
  peaksfileslop      = paste(outname ,"_peaks.slopPeak",sep="")
  bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")
  distancefile       = paste(outname ,".distance",sep="")
  distancecutpdf     = paste(outname ,".distance.pdf",sep="")
  pestimatepdf       = paste(outname ,".pEstimate.pdf",sep="")
  pestimate1txt       = paste(outname ,"pestimate1.txt",sep="")
  pestimate2txt       = paste(outname ,"pestimate2.txt",sep="")
  summarypdf         = paste(outname ,".summary.pdf",sep="")
  allpairsfile       = paste(outname ,".interactions.all.bedpe",sep="")
  fdrpairsfile       = paste(outname ,".interactions.fdr.bedpe",sep="")
  
  # build a file of just distances and same / dif
  print ("deterimining self-ligation distance")
  makeDistanceFile(bedpefilesortrmdup,distancefile,
                   distcutrangemin,
                   distcutrangemax)
  
  # calculate bias and cutoff
  distancecutoff = calcDistBias(distancefile,distancecutpdf=distancecutpdf,
                   range=c(distcutrangemin,distcutrangemax),
                   biascut= biascut)
  
  # group PETs into interactions
  print ("grouping PETs into interactions")
  chromosomes = groupPairs(bedpefilesortrmdup=bedpefilesortrmdup,
                           outname=outname,
                           peaksfile=peaksfileslop,
                           bedtoolspath = bedtoolspath,
                           verbose=FALSE)
  
  print(chromosomes)
  # filter out unwanted chromosomes
  originalchroms = chromosomes
  
  # get chromosomes from bedtools
  bedtoolsgenome = read.table(bedtoolsgenome,header=FALSE,sep="\t")
  chromosomes = bedtoolsgenome[,1]
  chromosomes = chromosomes[grep("_",chromosomes,invert=TRUE)]
  if(chrominclude[1] != "NULL")
  {
    chromosomes = unlist(strsplit(chrominclude,split=","))
  }
  
  if (chromexclude[1] !=  "NULL")
  {
    chromosomestpremove = unlist(strsplit(chromexclude,split=","))
    chromosomes = chromosomes[-which(chromosomes %in% chromosomestpremove)] 
  }  

  # estimate probabilities
  print ("estimating p-values")
  pEstimates = estimateP (chromosomes,outname,numofbins=numofbins ,binrange=c(distancecutoff,maxinteractingdist),outliers=NULL)
  
  # score and filter interactions
  print ("scoring interactions")
  allpairs   = scoreAndFilter(chromosomes = chromosomes,outname = outname,
                              mindist = distancecutoff, maxdist = maxinteractingdist,
                              averageDepth = pEstimates$averageDepth,
                              spline = pEstimates$spline,
                              N      = pEstimates$N,corrMethod)

  # filter out outliers
  outliercut = 1 / sum(pEstimates$p_table$M)
  outliers   = allpairs[which( allpairs$Q < outliercut),7]

  # re-estimate probabilities excluding outliers
  print ("estimating p-values (2nd iteration)")
  pEstimates2 = estimateP (chromosomes,outname,numofbins=numofbins ,binrange=c(distancecutoff,maxinteractingdist),outliers=outliers)


  # score and filter interactions with new estimates
  print ("scoring interactions (2nd iteration)")
  allpairs   = scoreAndFilter(chromosomes = chromosomes,outname = outname,
                              mindist = distancecutoff, maxdist = maxinteractingdist,
                              averageDepth = pEstimates2$averageDepth,
                              spline = pEstimates2$spline,
                              N      = pEstimates2$N,corrMethod)
  
  allpairs = cbind(allpairs[,c(1,2,3,4,5,6)],paste("pair_",(1:nrow(allpairs)),sep=""),allpairs[,c(10,11,12,14,18,19)])
  names(allpairs) = c("chrom1","start1","end1","chrom2","start2","end2","name",
                      "peak1","peak2","PETs","log10distance","P","Q")

  # write results to output
  if (reportallpairs == TRUE)
  {
    write.table(x=allpairs,file=allpairsfile,quote = FALSE, sep = "\t",row.names = FALSE)
  }
  sig= allpairs[which(allpairs$Q < FDR & allpairs$PETs >= minPETS),]
  write.table(x=sig,file=fdrpairsfile,quote = FALSE, sep = "\t",row.names = FALSE)
  
  resultshash[["putative interactions"]]    = nrow(allpairs)
  resultshash[["significant interactions"]] = nrow(sig)
  
  #- plot results -#
  
  # plot the spline
  pdf(pestimatepdf)
    plot(pEstimates$p_table$dist,pEstimates$p_table$p,pch=19,xlab="log10(distance)",ylab="p estimate", col = "firebrick2")
    points(pEstimates2$p_table$dist,pEstimates2$p_table$p,pch=19, col = "blue")
    lines(pEstimates$spline, col = "firebrick2")
    lines(pEstimates2$spline, col = "blue")
    legend("topright",inset=0.05,legend=c("bins 1","bins 2","spline fit 1","spline fit 2"),
           pch=c(19,19,NA,NA),lwd=c(NA,NA,1,1),col=c("firebrick2","blue","firebrick2","blue"))
  dev.off()
  
  # write spline to output
  save(pEstimates, file=pestimate1txt)
  save(pEstimates2,file=pestimate2txt)     

  pdf(summarypdf,height=6,width=10)
    par(mfrow=c(1,2))
    #hist(allpairs$P,main="P-value distribution",xlab="P",col="dodgerblue4")
    #hist(allpairs$Q,main="Q-value distribution",xlab="Q",col="dodgerblue4")

    plot(density(sig$log10distance),main="Size of sig pairs",xlab="log10(distance)",lwd=2,col="dodgerblue4",xlim=log10(c(distancecutoff,maxinteractingdist)))
    legend("topright",inset=0.05,legend=c(paste("Q <",FDR),paste("n =",nrow(sig))),pch=NA)
  
    maxy = max(hist(allpairs$PETs[which(allpairs$PETs<30 )],breaks=seq(0,30),plot=F)$counts)
    hist(allpairs$PETs[which(allpairs$PETs<30 )],breaks=seq(0,30),main="PETs in sig pairs",xlab="# PETs",
         ylab=paste("pairs (trunc ",maxy,")",sep=""),
         col="red",ylim=c(0,1.2*length(which(allpairs$PETs==2 )  )))
    hist(sig$PETs[which(sig$PETs<30)],breaks=seq(0,30),main="Q < 0.01",xlab="# PETs",add=TRUE,col="skyblue2")
    legend("topright",inset=0.05,legend=paste("Q <",FDR),pch=NA)
    legend("right",inset=0.05,legend=c("all observed","significant"),pch=15,col=c("red","skyblue2"),bty='n')
  dev.off()

  # clean up extra files
  if (file.exists(distancefile)) file.remove(distancefile)
  for (chrom in originalchroms)
  {
    print (chrom)
    peaksizecount  = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
    pairsbedpe     = paste(outname,"." ,chrom, ".pairs.bedpe",sep="")
    bedpefile      = paste(outname,"." ,chrom, ".bedpe",sep="")
    bedfile        = paste(outname,"." ,chrom, ".bed",sep="")
    print (bedfile)
    if (file.exists(peaksizecount)) file.remove(peaksizecount)
    if (file.exists(pairsbedpe)) file.remove(pairsbedpe)
    if (file.exists(bedpefile)) file.remove(bedpefile)
    if (file.exists(bedfile)) file.remove(bedfile)
  }
}

##################################### print to logfile #####################################

logfile = paste(args[["outname"]],".mango.log",sep="")
write(Sys.time(),file=logfile,append=TRUE)
write("Analyzed by Mango using the following parameters:",file=logfile,append=TRUE)
for (key in keys(args))
{
  write(paste( key, ":",args[[key]]),file=logfile,append=TRUE)
}
write("",file=logfile,append=TRUE)
write("With the following results:",file=logfile,append=TRUE)
print (resultshash)
for (key in keys(resultshash))
{
  write(paste( key, ":",resultshash[[key]]),file=logfile,append=TRUE)
}

print("done")
