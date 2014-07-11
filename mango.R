# runs mango chia pet analysis pipeline
suppressPackageStartupMessages(library("Rcpp"))
suppressPackageStartupMessages(library("hash"))
suppressPackageStartupMessages(library("mango"))
suppressPackageStartupMessages(library("optparse"))

print ("Starting mango ChIA PET analysis tool")
Sys.time()
set.seed(1)

##################################### check for dependencies #####################################

progs = c("bedtools","MACS2","bowtie")
Paths = DefinePaths(progs = progs)
bedtoolspath  = Paths[1]
macs2path     = Paths[2]
bowtiepath    = Paths[3]

i = 0
pathsOK = T
for (p in Paths)
{
  i = i+ 1
  if (p == "notfound")
  {
    pathsOK = F
    print ("! Configuration Error !")
    print (paste("     Path to ",progs[i]," not in PATH",sep=""))
    print ("     Please add to PATH and try again")
    print ("")
  }
}
if (pathsOK == F)
{
  stop ("Exiting Mango.R pipeline.  Check system PATH")
}  

##################################### read commandline paramters #####################################

# Should be able to delete readParameters and establishParameters if this works

# read in parameters
option_list <- list(
  
  #---------- GENERAL PARAMETERS ----------#
  
  make_option(c("--stages"),  default="1:5",help="stages of the pipeline to execute"),
  make_option(c("--prefix"), default="mango",help="prefix for all output files"),
  make_option(c("--outdir"),  default="NULL",help="out put directory"),
  make_option(c("--argsfile"),  default="NULL",help="optional argument file used to pass in command line paramters"),
  make_option(c("--bowtieref"),   default="NULL",help="genome reference file for bowtie"),
  make_option(c("--bedtoolsgenome"),  default="NULL",help="bedtools genome file"),
  make_option(c("--chrominclude"),  default="NULL",help="comma separated list of chromosomes to use (e.g. chr1,chr2,chr3,...).  Only these chromosomes will be processed"),
  make_option(c("--chromexclude"),  default="NULL",help="comma separated list of chromosomes to exclude (e.g. chrM,chrY)"),
  
  #---------- STAGE 1 PARAMETERS ----------#
  
  make_option(c("--fastq1"),  default="NULL",help="fastq read 1 file"),
  make_option(c("--fastq2"),  default="NULL",help="fastq read 2 file"),
  make_option(c("--linkerA"),  default="GTTGGATAAG",help="linker sequence A to look for"),
  make_option(c("--linkerB"),  default="GTTGGAATGT",help="linker sequence B to look for"),
  make_option(c("--minlength"),  default="15",help="min length of reads after linker trimming"),
  make_option(c("--maxlength"),  default="25",help="max length of reads after linker trimming"),
  make_option(c("--keepempty"),  default="FALSE",help="Should reads with no linker be kept"),
  
  #---------- STAGE 2 PARAMETERS ----------#
  
  make_option(c("--shortreads"),  default="TRUE",help="should bowtie alignments be done using paramter for very short reads (~20 bp)"),

  #---------- STAGE 4 PARAMETERS ----------#
  
  make_option(c("--MACS_pvalue"),  default="0.00001",help="MACS values"),
  make_option(c("--MACS_shiftsize"),  default="NULL",help="MACS shiftize.  NULL allows MACS to determine it"),
  make_option(c("--peakslop"),  default="500",help="Number of basespairs to extend peaks on both sides"),
  make_option(c("--peakinput"),  default="NULL",help="user supplied peaks file"),
  
  #---------- STAGE 5 PARAMETERS ----------#
  
  make_option(c("--distcutrangemin"),  default="1000",help="range in which to look for the self-ligation distance"),
  make_option(c("--distcutrangemax"),  default="100000",help="range in which to look for the self-ligation distance"),
  make_option(c("--biascut"),  default="0.05",help="Self ligation bias cutoff"),
  make_option(c("--maxPval"),  default="0.01",help="P-value cutoff"),
  make_option(c("--numofbins"),  default="30",help="number of bins for probability calculations"),
  make_option(c("--corrMethod"),  default="BY",help="multiple hypothesis tersting correction method"),
  make_option(c("--maxinteractingdist"),  default="10000000",help="maximum disance allowed for an interaction"),
  make_option(c("--FDR"),  default="0.01",help="FDR cutoff for interactions"),
  make_option(c("--minPETS"),  default="2",help="minimum number of PETs required for an interaction (applied after FDR filtering)"),
  make_option(c("--reportallpairs"),  default="FALSE",help="Should all pairs be reported or just significant pairs")  
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

# correct stages
if (grepl( ":",opt$stages))
{
  stages = strsplit(opt$stages,split=":")[[1]]
  if (length(stages) == 2)
  {
    opt$stages = seq(as.numeric(stages[1]),as.numeric(stages[2]))
  }
  if (length(stages) == 1)
  {
    opt$stages = as.numeric(stages[1])
  }
}

# set basename for output
opt["outname"] = ""
if (opt["outdir"] == "NULL")
{
  opt["outname"] = file.path(getwd(),opt["prefix"])
}
if (opt["outdir"] != "NULL")
{
  opt["outname"] = file.path( opt["outdir"] ,opt["prefix"]) 
}


##################################### read in arguments #####################################

lines = c()
if (opt$argsfile != "NULL")
{
  if (file.exists(opt$argsfile) == FALSE)
  {
    stop (paste("Exiting Mango.R pipeline.  No file named ",opt$argsfile,sep=""))
  }
  lines = readLines(argsfile)
}

for (line in lines)
{
  # remove spaces
  line = gsub(pattern=" ",x=line,replace="")
  
  if (line == "")
  {
    next
  }
  else if (strsplit(line,split="")[[1]][1] == "#" )
  {
    next
  }
  else
  {
    lineinfo = strsplit(line,split="#")[[1]][1]
    arginfo  = strsplit(lineinfo,split="=")[[1]]
    opt[arginfo[1]] = arginfo[2]
  }
}
       
resultshash = hash()

##################################### parse fastqs #####################################

if (1 %in% opt$stages)
{
  checkRequired(opt,c("fastq1","fastq2","outname","minlength","maxlength","keepempty","linkerA","linkerB"))

  # gather arguments
  outname         = as.character(opt["outname"])
  fastq1=as.character(opt["fastq1"])
  fastq2=as.character(opt["fastq2"])
  basename = as.character(opt["outname"])
  minlength = as.numeric(as.character(opt["minlength"]))
  maxlength = as.numeric(as.character(opt["maxlength"]) )
  keepempty=eval(parse(text=as.character(opt["keepempty"])))
  linker1=as.character(opt["linkerA"])
  linker2=as.character(opt["linkerB"])
  
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

if (2 %in% opt$stages)
{
  checkRequired(opt,c("outname","bowtieref","shortreads"))
  
  # gather arguments
  outname         = as.character(opt["outname"])
  bowtieref       = as.character(opt["bowtieref"])
  shortreads      = as.character(opt["shortreads"])
  
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

if (3 %in% opt$stages)
{
  checkRequired(opt,c("outname"))
  
  # gather arguments
  outname         =  as.character(opt["outname"])
  
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

if (4 %in% opt$stages)
{
  checkRequired(opt,c("outname","MACS_pvalue",
                       "bedtoolsgenome","peakslop","peakinput"))
  
  # gather arguments
  outname         = as.character(opt["outname"])
  MACS_pvalue     = as.character(opt["MACS_pvalue"])
  bedtoolsgenome  = as.character(opt["bedtoolsgenome"])
  peakslop        = as.character(opt["peakslop"])
  peakinput       = as.character(opt["peakinput"])
  MACS_shiftsize  = as.character(opt["MACS_shiftsize"])
  
  # filenames
  bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")
  tagAlignfile       = paste(outname,".tagAlign",sep="")
  peaksfile          = paste(outname,"_peaks.narrowPeak",sep="")
  peaksfileslop      = paste(outname,"_peaks.slopPeak",sep="")
  
  if (peakinput != "NULL")
  {
    peaksfile = peakinput
  }
  
  if (peakinput == "NULL")
  {
    print ("building tagAlign file")
    # reverse strands for peak calling
    if (file.exists(tagAlignfile)){file.remove(tagAlignfile)}
    buildTagAlign(bedpefilesortrmdup ,tagAlignfile )
    
    # call peaks 
    print ("calling peaks")
   callpeaks(macs2path=macs2path,tagAlignfile,outname,pvalue=MACS_pvalue,
             bedtoolspath=bedtoolspath,bedtoolsgenome=bedtoolsgenome,
             peakslop=peakslop,MACS_shiftsize)
  }
  
  # extend and merge peaks according to peakslop
  print ("extending peaks")
  peakcounts = extendpeaks(peaksfile,peaksfileslop,bedtoolspath=bedtoolspath,
             bedtoolsgenome=bedtoolsgenome,peakslop=peakslop)
  resultshash[["peaks"]] = peakcounts[1]
  resultshash[["mergedpeaks"]] = peakcounts[2]
}

##################################### group pairs #####################################

if (5 %in% opt$stages)
{
  checkRequired(opt,c("outname","distcutrangemin","biascut","maxPval",
                       "numofbins","corrMethod","bedtoolsgenome",
                       "maxinteractingdist","FDR","minPETS","corrMethod",
                       "chrominclude","chromexclude","reportallpairs"))
  
  # gather arguments
  outname = as.character(opt["outname"])
  distcutrangemin = as.numeric(as.character(opt["distcutrangemin"]))
  distcutrangemax = as.numeric(as.character(opt["distcutrangemax"]))
  bedtoolsgenome    = as.character(opt["bedtoolsgenome"])
  biascut = as.numeric(as.character(opt["biascut"]))
  maxinteractingdist = as.numeric(as.character(opt["maxinteractingdist"]))
  numofbins = as.numeric(as.character(opt["numofbins"]))
  FDR = as.numeric(as.character(opt["FDR"]))
  minPETS = as.numeric(as.character(opt["minPETS"]))
  corrMethod = as.character(opt["corrMethod"])
  chrominclude      = as.character(opt["chrominclude"])
  chromexclude      = as.character(opt["chromexclude"])
  reportallpairs    = as.character(opt["reportallpairs"])
  normPmeth          = "product"
  
  # filenames
  peaksfile          = paste(outname ,"_peaks.narrowPeak",sep="")
  peaksfileslop      = paste(outname ,"_peaks.slopPeak",sep="")
  bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")
  distancefile       = paste(outname ,".distance",sep="")
  distancecutpdf     = paste(outname ,".distance.pdf",sep="")
  pestimatepdf       = paste(outname ,".pEstimate.pdf",sep="")
  pestimate1txt      = paste(outname ,".pestimate1.txt",sep="")
  pestimate2txt      = paste(outname ,".pestimate2.txt",sep="")
  summarypdf         = paste(outname ,".summary.pdf",sep="")
  allpairsfile       = paste(outname ,".interactions.all.bedpe",sep="")
  fdrpairsfile       = paste(outname ,".interactions.fdr.bedpe",sep="")
  
  
  # build a file of just distances and same / dif
  print ("determining self-ligation distance")
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
  pEstimates = estimateP (chromosomes,outname,numofbins=numofbins ,binrange=c(distancecutoff,maxinteractingdist),outliers=NULL,normPmeth = normPmeth)
  
  # score and filter interactions
  print ("scoring interactions")
  allpairs   = scoreAndFilter(chromosomes = chromosomes,outname = outname,
                              mindist = distancecutoff, maxdist = maxinteractingdist,
                              averageDepth = pEstimates$averageDepth,
                              spline = pEstimates$spline,
                              N      = pEstimates$N,corrMethod,
                              normPmeth = normPmeth)

  # filter out outliers
  outliercut = 1 / sum(pEstimates$p_table$M)
  outliers   = allpairs[which( allpairs$Q < outliercut),7]

  # re-estimate probabilities excluding outliers
  print ("estimating p-values (2nd iteration)")
  pEstimates2 = estimateP (chromosomes,outname,numofbins=numofbins ,binrange=c(distancecutoff,maxinteractingdist),outliers=outliers,normPmeth = normPmeth)


  # score and filter interactions with new estimates
  print ("scoring interactions (2nd iteration)")
  allpairs   = scoreAndFilter(chromosomes = chromosomes,outname = outname,
                              mindist = distancecutoff, maxdist = maxinteractingdist,
                              averageDepth = pEstimates2$averageDepth,
                              spline = pEstimates2$spline,
                              N      = pEstimates2$N,corrMethod,
                              normPmeth = normPmeth)
  
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
    #if (file.exists(peaksizecount)) file.remove(peaksizecount)
    if (file.exists(pairsbedpe)) file.remove(pairsbedpe)
    if (file.exists(bedpefile)) file.remove(bedpefile)
    if (file.exists(bedfile)) file.remove(bedfile)
  }
}

##################################### print to logfile #####################################

logfile = paste(opt["outname"],".mango.log",sep="")
write(Sys.time(),file=logfile,append=TRUE)
write("Analyzed by Mango using the following parameters:",file=logfile,append=TRUE)
for (key in keys(args))
{
  write(paste( key, ":",opt[key]),file=logfile,append=TRUE)
}
write("",file=logfile,append=TRUE)
write("With the following results:",file=logfile,append=TRUE)
print (resultshash)
for (key in keys(resultshash))
{
  write(paste( key, ":",resultshash[[key]]),file=logfile,append=TRUE)
}

print("done")
