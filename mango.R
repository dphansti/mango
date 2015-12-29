# runs mango chia pet analysis pipeline

# Version info
Mangoversion = "1.0.5"

# Load Packages
suppressPackageStartupMessages(library("Rcpp"))
suppressPackageStartupMessages(library("hash"))
suppressPackageStartupMessages(library("mango"))
suppressPackageStartupMessages(library("optparse"))

print ("Starting mango ChIA PET analysis tool")
Sys.time()
set.seed(1)

##################################### read commandline paramters #####################################

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
  make_option(c("--bedtoolspath"),  default="NULL",help="full path to bedtools"),
  make_option(c("--macs2path"),  default="NULL",help="full path to macs2path"),
  make_option(c("--bowtiepath"),  default="NULL",help="full path to bowtiepath"),
  make_option(c("--verboseoutput"),  default="FALSE",help="if true output file will have more columns as well as a header row describing the columns"),
  
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
  
  make_option(c("--MACS_qvalue"),  default="0.05",help="MACS values"),
  make_option(c("--MACS_shiftsize"),  default="NULL",help="MACS shiftize.  NULL allows MACS to determine it"),
  make_option(c("--peakslop"),  default="500",help="Number of basespairs to extend peaks on both sides"),
  make_option(c("--peakinput"),  default="NULL",help="user supplied peaks file"),
  make_option(c("--blacklist"),  default="NULL",help="BED file of regions to remove from MACS peaks"),
  make_option(c("--gsize"),  default="hs",help="mappable genome size or effective genome size for MACS2"),
  
  #---------- STAGE 5 PARAMETERS ----------#
  
  make_option(c("--distcutrangemin"),  default="1000",help="range in which to look for the self-ligation distance"),
  make_option(c("--distcutrangemax"),  default="100000",help="range in which to look for the self-ligation distance"),
  make_option(c("--biascut"),  default="0.05",help="Self ligation bias cutoff"),
  make_option(c("--numofbins"),  default="50",help="number of bins for probability calculations"),
  make_option(c("--corrMethod"),  default="BH",help="multiple hypothesis tersting correction method"),
  make_option(c("--maxinteractingdist"),  default="1000000",help="maximum disance allowed for an interaction"),
  make_option(c("--FDR"),  default="0.05",help="FDR cutoff for interactions"),
  make_option(c("--extendreads"),  default="120",help="how many bp to extend reads towards peak"),
  make_option(c("--minPETS"),  default="2",help="minimum number of PETs required for an interaction (applied after FDR filtering)"),
  make_option(c("--reportallpairs"),  default="FALSE",help="Should all pairs be reported or just significant pairs"),
  make_option(c("--MHT"),  default="all",help="How should mutliple hypothsesis testing be done?  Correct for 'all' possible pairs of loci or only those 'found' with at least 1 PET")  
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

# check dependencies

# first look in PATH
progs = c("bedtools","macs2","bowtie")
Paths = DefinePaths(progs = progs)
bedtoolspath  = Paths[1]
macs2path     = Paths[2]
bowtiepath    = Paths[3]


# next, look in arguments
if (opt["bedtoolspath"] != "NULL")
{
  bedtoolspath = opt["bedtoolspath"]
}
if (opt["macs2path"] != "NULL")
{
  macs2path = opt["macs2path"]
}
if (opt["bowtiepath"] != "NULL")
{
  bowtiepath = opt["bowtiepath"]
}

# get software versions
bedtoolsversion = system(paste(bedtoolspath,"--version"),intern=TRUE)[1]
macs2version    = system(paste(macs2path,"--version 2>&1"),intern=TRUE)[1]
bowtieversion   = system(paste(bowtiepath,"--version"),intern=TRUE)[1]

# break if dependencies not found
Paths = c(bedtoolspath,macs2path,bowtiepath)
i = 0
pathsOK = T
for (p in Paths)
{
  i = i+ 1
  if (p == "notfound")
  {
    pathsOK = F
    print ("! Configuration Error !")
    print (paste("     Path to ",progs[i]," not in PATH or arguments",sep=""))
    print ("     Please add to PATH or arguments and try again")
    print ("")
  }
}
if (pathsOK == F)
{
  stop ("Exiting Mango.R pipeline.  Check system PATH")
}

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

logfile = paste(as.character(opt["outname"]),".mango.log",sep="")
if (file.exists(logfile) ==TRUE){file.remove(logfile)}
starttime = paste("Analysis start time:" , as.character(Sys.time()))
write(starttime,file=logfile,append=TRUE)

##################################### read in arguments #####################################

lines = c()
if (as.character(opt$argsfile) != "NULL")
{
  if (file.exists(as.character(opt$argsfile)) == FALSE)
  {
    stop (paste("Exiting Mango.R pipeline.  No file named ",opt$argsfile,sep=""))
  }
  lines = readLines(as.character(opt$argsfile))
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
  checkRequired(opt,c("fastq1","fastq2"))

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
  checkRequired(opt,c("bowtieref"))
  
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
  
  # split by chromosome and position
  print ("removing duplicate PETs")
  distancesplit = 10000000
  rmdupresults = removeDups(bedpefile,outname,distancesplit)
  
  # add info to results hash for log file
  resultshash[["duplicate PETs"]] = rmdupresults[1]
  resultshash[["nonduplicate PETs"]] = rmdupresults[2]
  resultshash[["interchromosomal PETs"]] = rmdupresults[3]
  resultshash[["intrachromosomal PETs"]] = rmdupresults[4]
  resultshash[["aligned PETs"]] = rmdupresults[5]
  
  # remove temporary files
  for (f in (5:length(rmdupresults)))
  {
    if (file.exists(rmdupresults[f])==TRUE){file.remove(rmdupresults[f])} 
  }
  
#   # split by chrom and sort bedpe
#   print ("sorting bedpe")
#   if (file.exists(bedpefilesort)){file.remove(bedpefilesort)}

#   
#   # this as the way to sort the whol file at once
#   #external_sort(bedpefile, bedpefilesort)
#   
#   # filter duplicates
#   print ("removing PCR duplicates")
#   if (file.exists(bedpefilesortrmdup)){file.remove(bedpefilesortrmdup)}
#   rmdupresults = removeDupBedpe(bedpefilesort,bedpefilesortrmdup,renamePets=TRUE);
#   resultshash[["duplicate PETs"]] = rmdupresults[1]
#   resultshash[["nonduplicate PETs"]] = rmdupresults[2]
#   resultshash[["interchromosomal PETs"]] = rmdupresults[3]
#   resultshash[["intrachromosomal PETs"]] = rmdupresults[4]
}

##################################### call peaks #####################################

if (4 %in% opt$stages)
{
  checkRequired(opt,c("bedtoolsgenome"))
  
  # gather arguments
  outname         = as.character(opt["outname"])
  MACS_qvalue     = as.character(opt["MACS_qvalue"])
  bedtoolsgenome  = as.character(opt["bedtoolsgenome"])
  peakslop        = as.character(opt["peakslop"])
  peakinput       = as.character(opt["peakinput"])
  MACS_shiftsize  = as.character(opt["MACS_shiftsize"])
  blacklist       = as.character(opt["blacklist"])
  gsize           = ascharacter(opt["gsize"])
  
  # filenames
  bedpefilesortrmdup = paste(outname ,".rmdup.bedpe",sep="")
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
    callpeaks(macs2path=macs2path,tagAlignfile,outname,qvalue=MACS_qvalue,
             bedtoolspath=bedtoolspath,bedtoolsgenome=bedtoolsgenome,
             peakslop=peakslop,MACS_shiftsize,gsize=gsize)
  }
  
  # extend and merge peaks according to peakslop
  print ("extending peaks")
  peakcounts = extendpeaks(peaksfile,peaksfileslop,bedtoolspath=bedtoolspath,
             bedtoolsgenome=bedtoolsgenome,peakslop=peakslop,blacklist=blacklist)
  resultshash[["peaks"]] = peakcounts[1]
  resultshash[["mergedpeaks"]] = peakcounts[2]
}


##################################### new group/score/filter pairs #####################################

if (5 %in% opt$stages)
{
  checkRequired(opt,c("outname","bedtoolsgenome"))
  
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
  chrominclude      = as.character(opt["chrominclude"])
  chromexclude      = as.character(opt["chromexclude"])
  reportallpairs    = as.character(opt["reportallpairs"])
  corrMethod = as.character(opt["corrMethod"])
  MHT    = as.character(opt["MHT"])
  extendreads = as.numeric(opt["extendreads"])
  verboseoutput = as.character(opt["verboseoutput"])

  # filenames
  tagAlignfile       = paste(outname,".tagAlign",sep="")
  tagAlignfileExt    = paste(outname ,".tagAlign.extended.bed",sep="")
  temppeakoverlap    = paste(outname ,".temppeakoverlap.bed",sep="")
  peaksfile          = paste(outname ,"_peaks.narrowPeak",sep="")
  peaksfileslop      = paste(outname ,"_peaks.slopPeak",sep="")
  peaksfileslopdepth = paste(outname ,"_peaks.slopPeak.depth",sep="")
  bedpefilesortrmdup = paste(outname ,".rmdup.bedpe",sep="")
  distancefile       = paste(outname ,".distance",sep="")
  distancecutpdf     = paste(outname ,".distance.pdf",sep="")
  modelspdf          = paste(outname ,".models.pdf",sep="")
  allpairsfile       = paste(outname ,".interactions.all.mango",sep="")
  fdrpairsfile       = paste(outname ,".interactions.fdr.mango",sep="")
  
  # counting reads per peak
  print ("counting reads per peak")
  if (file.exists(tagAlignfileExt) ==TRUE){file.remove(tagAlignfileExt)}
  if (file.exists(temppeakoverlap) ==TRUE){file.remove(temppeakoverlap)}
  DeterminePeakDepths(bedtools=bedtoolspath,bedtoolsgenome=bedtoolsgenome,extendreads=extendreads,tagAlignfile=tagAlignfile,
                  tagAlignfileExt=tagAlignfileExt,peaksfileslop=peaksfileslop,temppeakoverlap=temppeakoverlap)
  if (file.exists(tagAlignfileExt) ==TRUE){file.remove(tagAlignfileExt)}
  if (file.exists(temppeakoverlap) ==TRUE){file.remove(temppeakoverlap)}
  
  # build a file of just distances and same / dif
  print ("determining self-ligation distance")
  makeDistanceFile(bedpefilesortrmdup,distancefile,
                   distcutrangemin,
                   distcutrangemax)
  
  # calculate bias and cutoff
  distancecutoff = calcDistBias(distancefile,distancecutpdf=distancecutpdf,
                                range=c(distcutrangemin,distcutrangemax),
                                biascut= biascut)
  #distancecutoff = 20236
  print (paste("self-ligation cutoff =",distancecutoff))
    
  # group PETs into interactions
  print ("grouping PETs into interactions")
  chromosomes = groupPairs(bedpefilesortrmdup=bedpefilesortrmdup,
                           outname=outname,
                           peaksfile=peaksfileslop,
                           bedtoolspath = bedtoolspath,
                           bedtoolsgenome = bedtoolsgenome,
                           extendreads=extendreads,peaksfileslopdepth=peaksfileslopdepth,
                           verbose=FALSE)
  
  # filter out unwanted chromosomes
  originalchroms = chromosomes
  
  # get chromosomes from bedtools
  bedtoolsgenomeinfo = read.table(bedtoolsgenome,header=FALSE,sep="\t")
  chromosomes = bedtoolsgenomeinfo[,1]
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
  
  print ("modeling PETs based on peak depth and distance")
  #--------------- Gather IAB data ---------------#
  
  # gather all putative interactions
  putpairs = combineputativepairs(chromosomes,outname)
  
  # calculate interaction distances
  putpairs$distances = abs( (putpairs[,2] + putpairs[,3] ) / 2 - (putpairs[,5] + putpairs[,6] ) / 2  )
  
  # filter out putative interactions that don't fall into distance range
  putpairs = putpairs[which(putpairs$distances < maxinteractingdist & putpairs$distances > distancecutoff),]
  
  # calculate depths
  putpairs$depths = calcDepths(putpairs[,10:11],type="product")
  
  totalcombos = 0
  for (reps in (1:2))
  {
    #--------------- Distance Normalization ---------------#
    
    # determine borders to distance bins
    distanceborders = binmaker(putpairs$distances,binmethod="equalocc",numberbins=numofbins)
    
    # model IAB vs distance
    distance_IAB_model = model_chia(x=putpairs$distances,y=putpairs[,12],borders=distanceborders,yvals=TRUE)
    distance_IAB_model_file   = paste(outname ,".distance_IAB_model.",reps, ".text",sep="")
    write.table(distance_IAB_model,file=distance_IAB_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
    distance_IAB_spline =   smooth.spline(log10(distance_IAB_model[,1]),distance_IAB_model[,3],spar=.75)
    
    #--------------- Depth Normalization ---------------#
    
    # determine borders to depth bins
    depthborders = binmaker(putpairs$depths,binmethod="equalocc",numberbins=numofbins)
    
    # model IAB vs depth
    depth_IAB_model = model_chia(x=putpairs$depths,y=putpairs[,12],borders=depthborders,yvals=TRUE)
    depth_IAB_model_file   = paste(outname ,".depth_IAB_model.",reps, ".text",sep="")
    write.table(depth_IAB_model,file=depth_IAB_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
    depth_IAB_spline =   smooth.spline(log10(depth_IAB_model[,1]),depth_IAB_model[,3],spar=.75)
    
    # model Combos vs distance
    meanofx_dist  = rep(0,numofbins)
    sumofy_dist   = rep(0,numofbins)
    pvals_dist    = rep(0,numofbins)
    sumofx_dist   = rep(0,numofbins)
    countofx_dist = rep(0,numofbins)
    meanofx_depth  = rep(0,numofbins)
    sumofy_depth   = rep(0,numofbins)
    pvals_depth    = rep(0,numofbins)
    sumofx_depth   = rep(0,numofbins)
    countofx_depth = rep(0,numofbins)

    for (chrom in chromosomes)
    {
      # make combos
      combos = makecombos(chrom,outname,mindist=distancecutoff,maxdist=maxinteractingdist)
      
      # calculate distances
      combos$distance = abs( (combos[,2] + combos[,3] ) / 2 - (combos[,5] + combos[,6] ) / 2  )
      
      # calculate depths
      combos$depths = calcDepths(combos[,7:8],type="product")
      
      # model Combo vs distance
      distance_combo_model_chrom = model_chia(x=combos$distance,y=NA,borders=distanceborders,yvals=FALSE)
      
      sumofy_dist   = sumofy_dist   + distance_combo_model_chrom[,2]
      sumofx_dist   = sumofx_dist   + distance_combo_model_chrom[,4]
      countofx_dist = countofx_dist + distance_combo_model_chrom[,5]
       
      # model Combo vs depth
      depth_combo_model_chrom  = model_chia(x=combos$depths,y=NA,borders=depthborders,yvals=FALSE)
      
      sumofy_depth   = sumofy_depth   + depth_combo_model_chrom[,2]
      sumofx_depth   = sumofx_depth   + depth_combo_model_chrom[,4]
      countofx_depth = countofx_depth + depth_combo_model_chrom[,5]
    }
    
    # combine data from all chromosomes
    depth_combo_model    = cbind(sumofx_depth/countofx_depth, sumofy_depth, sumofy_depth/sum(sumofy_depth))
    distance_combo_model = cbind(sumofx_dist /countofx_dist,  sumofy_dist, sumofy_dist/sum(sumofy_dist))
    
    depth_combo_model_file   = paste(outname ,".depth_combo_model.",reps, ".text",sep="")
    write.table(depth_combo_model,file=depth_combo_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
    
    distance_combo_model_file   = paste(outname ,".distance_combo_model.",reps, ".text",sep="")
    write.table(distance_combo_model,file=distance_combo_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
    
    depth_combo_spline    =   smooth.spline(log10(depth_combo_model[,1]),depth_combo_model[,3],spar=.75)
    distance_combo_spline =   smooth.spline(log10(distance_combo_model[,1]),distance_combo_model[,3],spar=.75)

    if (reps == 2)
    {
      #--------------- Gather IAB data again ---------------#
      
      # gather all putative interactions
      putpairs = combineputativepairs(chromosomes,outname)
      
      # calculate interaction distances
      putpairs$distances = abs( (putpairs[,2] + putpairs[,3] ) / 2 - (putpairs[,5] + putpairs[,6] ) / 2  )
      
      # filter out putative interactions that don't fall into distance range
      putpairs = putpairs[which(putpairs$distances < maxinteractingdist & putpairs$distances > distancecutoff),]
      
      # calculate depths
      putpairs$depths = calcDepths(putpairs[,10:11],type="product")
    }

    #--------------- Score putative interactions ---------------#
    
#    putpairstemfile       = paste(outname ,".putpairs",sep="")
#    write.table(putpairs,file=putpairstemfile,quote = FALSE, sep = "\t",row.names = FALSE)
    
    # Assing the four probabilities
    putpairs$P_IAB_distance    = predict(distance_IAB_spline, log10(putpairs$distances))$y
    putpairs$P_combos_distance = predict(distance_combo_spline,log10(putpairs$distances))$y
    putpairs$P_IAB_depth       = predict(depth_IAB_spline,log10(putpairs$depths))$y
    putpairs$P_combos_depth    = predict(depth_combo_spline,log10(putpairs$depths))$y
    
#     # fix negative values
#     putpairs$P_IAB_distance[which(putpairs$P_IAB_distance <= 0)] = 
#       min(putpairs$P_IAB_distance[which(putpairs$P_IAB_distance > 0)])
#     putpairs$P_combos_distance[which(putpairs$P_combos_distance <= 0)] = 
#       min(putpairs$P_combos_distance[which(putpairs$P_combos_distance > 0)])
#     putpairs$P_IAB_depth[which(putpairs$P_IAB_depth <= 0)] = 
#       min(putpairs$P_IAB_depth[which(putpairs$P_IAB_depth > 0)])
#     putpairs$P_combos_depth[which(putpairs$P_combos_depth <= 0)] = 
#       min(putpairs$P_combos_depth[which(putpairs$P_combos_depth > 0)])

    # cap values to min and max
    putpairs$P_IAB_distance[which(putpairs$P_IAB_distance <= min(distance_IAB_model[,3]))] = min(distance_IAB_model[,3])
    putpairs$P_IAB_distance[which(putpairs$P_IAB_distance >= max(distance_IAB_model[,3]))] = max(distance_IAB_model[,3]) 

    putpairs$P_combos_distance[which(putpairs$P_combos_distance <= min(distance_combo_model[,3]))] = min(distance_combo_model[,3])
    putpairs$P_combos_distance[which(putpairs$P_combos_distance >= max(distance_combo_model[,3]))] = max(distance_combo_model[,3])

    putpairs$P_IAB_depth[which(putpairs$P_IAB_depth <=  min(depth_IAB_model[,3]))] =  min(depth_IAB_model[,3])  
    putpairs$P_IAB_depth[which(putpairs$P_IAB_depth >=  max(depth_IAB_model[,3]))] =  max(depth_IAB_model[,3])
                                                                                          
    putpairs$P_combos_depth[which(putpairs$P_combos_depth <= min(depth_combo_model[,3]))] =  min(depth_combo_model[,3])   
    putpairs$P_combos_depth[which(putpairs$P_combos_depth >= max(depth_combo_model[,3]))] =  max(depth_combo_model[,3])   

    # calculate the binomial probability
    totalcombos =  sum(sumofy_dist)
    putpairs$p_binom           = (putpairs$P_IAB_distance * putpairs$P_IAB_depth) / 
      (putpairs$P_combos_distance * putpairs$P_combos_depth * totalcombos)
   
    # calculate the total IABs
    totalIAB = sum(distance_IAB_model[,2])

    # calculate the final interaction P values
    putpairs$P = apply(cbind(putpairs$V12,rep(totalIAB,nrow(putpairs)),putpairs$p_binom),1,calcP)

    if (reps == 1)
    {
      putpairs[which(putpairs$P < 1/totalcombos ),]$V12 = 0
    }
  }
  
  #--------------- Correct for multiple hypothesis testing ---------------#

  print ("correcting for multiple hypothesis testing")
  n=nrow(putpairs)
  if (MHT == "all")
  {
    n=totalcombos
  }

  putpairs$Q = p.adjust(putpairs$P,method=corrMethod,n=n)

  #--------------- Organize  data ---------------#
  
  pairnames = paste("pair_",(1:nrow(putpairs)),sep="")
  putpairs = cbind(putpairs[,c(1,2,3,4,5,6)],pairnames,putpairs[,c(10,11,12,14,16,17,18,19,20,21,22)])
  names(putpairs) = c("chrom1","start1","end1","chrom2","start2","end2","name",
                      "peak1","peak2","PETs","distance",
                      "P_IAB_distance","P_combos_distance","P_IAB_depth","P_combos_depth",
                      "p_binom","P","Q")
  
  #--------------- Filter interactions ---------------#

  sig = putpairs[which(putpairs$Q < FDR & putpairs$PETs >= minPETS),]

  #--------------- Write outputs ---------------#
  
  print ("writing output files")
  # write results to output
  if (verboseoutput == TRUE)
  {
    if (reportallpairs == TRUE)
    {
      write.table(x=putpairs,file=allpairsfile,quote = FALSE, sep = "\t",row.names = FALSE)
    }
    write.table(x=sig,file=fdrpairsfile,quote = FALSE, sep = "\t",row.names = FALSE)
  }
  if (verboseoutput == FALSE)
  {
    if (reportallpairs == TRUE)
    {
      write.table(x=putpairs[,c("chrom1","start1","end1","chrom2","start2","end2","PETs","Q")],file=allpairsfile,quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
    }
    write.table(x=sig[,c("chrom1","start1","end1","chrom2","start2","end2","PETs","Q")],file=fdrpairsfile,quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
  } 
  resultshash[["putative interactions"]]    = nrow(putpairs)
  resultshash[["significant interactions"]] = nrow(sig)

  #--------------- Make plots ---------------#
  
  print ("plotting results")
  
  # plot models  
  pdf(modelspdf)
  par(mfrow=c(2,2))
  par(mgp=c(3,.3,0))
  plot(log10(distance_IAB_model[,1]),   distance_IAB_model[,3],pch=19,col="dodgerblue2",xlab="distance (bp)",ylab="IAB")
  lines(x=log10(distance_IAB_model[,1]),   predict(distance_IAB_spline,log10(distance_IAB_model[,1]))$y)   
  
  plot(log10(distance_combo_model[,1]), distance_combo_model[,3],  pch=19,col="firebrick2",xlab="distance (bp)",ylab="# combos")  
  lines(x=log10(distance_combo_model[,1]),   predict(distance_combo_spline,log10(distance_combo_model[,1]))$y)  
  
  plot(log10(depth_IAB_model[,1])   ,  depth_IAB_model[,3],    pch=19,col="dodgerblue2",xlab="depth (p1 * p2)",ylab="IAB") 
  lines(x=log10(depth_IAB_model[,1]),   predict(depth_IAB_spline,log10(depth_IAB_model[,1]))$y)
  
  plot(log10(depth_combo_model[,1]),   depth_combo_model[,3],  pch=19,col="firebrick2",xlab="depth (p1 * p2)",ylab="# combos") 
  lines(x=log10(depth_combo_model[,1]),   predict(depth_combo_spline,log10(depth_combo_model[,1]))$y)
  dev.off()
  
  #--------------- Delete temporary files ---------------#
  
  # clean up extra files
  print ("deleting temporary files")
  if (file.exists(distancefile)) file.remove(distancefile)
  for (chrom in originalchroms)
  {
    peaksizecount  = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
    pairsbedpe     = paste(outname,"." ,chrom, ".pairs.bedpe",sep="")
    bedpefile      = paste(outname,"." ,chrom, ".bedpe",sep="")
    bedfile        = paste(outname,"." ,chrom, ".bed",sep="")
    overlapfile    = paste(outname,"." ,chrom, ".bedNpeak",sep="")
    if (file.exists(peaksizecount)) file.remove(peaksizecount)
    if (file.exists(pairsbedpe)) file.remove(pairsbedpe)
    if (file.exists(bedpefile)) file.remove(bedpefile)
    if (file.exists(bedfile)) file.remove(bedfile)
    if (file.exists(overlapfile)) file.remove(overlapfile)
  }
} 

##################################### Make Log file #####################################

print ("writing to log file")

stoptime = paste("Analysis end time:" , as.character(Sys.time()))
write(stoptime,file=logfile,append=TRUE)
write("",file=logfile,append=TRUE)

write("Software Versions:",file=logfile,append=TRUE)
write(paste("mango version",":",Mangoversion),file=logfile,append=TRUE)
write(paste("bedtools version",":",bedtoolsversion),file=logfile,append=TRUE)
write(paste("macs2 version",":",macs2version),file=logfile,append=TRUE)
write(paste("bowtie version",":",bowtieversion),file=logfile,append=TRUE)

write("",file=logfile,append=TRUE)
write("Analyzed by Mango using the following parameters:",file=logfile,append=TRUE)

for (key in names(opt))
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










