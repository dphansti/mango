# runs mango chia pet analysis pipeline
usr = "doug"

##################################### paths to externals #####################################

#! Need to figure out how to get R to use the same environment as bash (i.e. same PATH)

fastqs = c("data/NH.K562_RAD21_K562_std_2.1_1.head.fastq",
           "data/NH.K562_RAD21_K562_std_2.1_2.head.fastq")

expname = "NH.K562_RAD21_K562_std_2.1"

if (usr == "alan")
{
  bowtiepath = "/home/aboyle/bowtie-1.0.1/bowtie"
  bowtieref  = "/home/aboyle/bowtie-1.0.1/indexes/hg19"
  outdir   = "/home/aboyle/mango2014test/"
}

if (usr == "doug")
{
  bowtiepath = "/Users/dougphanstiel/Tools/bowtie-1.0.0/bowtie"
  bowtieref  = "/Users/dougphanstiel/Tools/bowtie-1.0.0/indexes/hg19"
  bedtoolspath = "/Users/dougphanstiel/Tools/bedtools-2.17.0/bin/bedtools"
  bedtoolsgenome = "/Users/dougphanstiel/Tools/bedtools-2.17.0/genomes/human.hg19.genome"
  macs2path  = "/usr/local/bin/macs2"
  outdir     = "/Users/dougphanstiel/Desktop/mango2014test/"
  bigfastqs = c(paste(outdir,"NH.K562_RAD21_K562_std_2.1_1.fastq",sep=""),
                paste(outdir,"NH.K562_RAD21_K562_std_2.1_2.fastq",sep=""))
  fastqs = bigfastqs
}

if (usr == "aster")
{
  bowtiepath = "/Users/dougphanstiel/tools/bowtie-0.12.7/bowtie"
  bowtieref  = "/Users/dougphanstiel/tools/bowtie-0.12.7/indexes/hg19"
  outdir   = "/Volumes/HD3/projects/newmango/data/"
  bigfastqs = fastqs = c("/Volumes/HD3/projects/newmango/data/NH.K562_RAD21_K562_std_2.1_1.fastq",
                         "/Volumes/HD3/projects/newmango/data/NH.K562_RAD21_K562_std_2.1_2.fastq")
  #fastqs = bigfastqs
}

outname  = paste(outdir,expname,sep="")

##################################### initialization #####################################

library('mango')
print ("Starting mango ChIA PET analysis tool")
Sys.time()
set.seed(1)

##################################### parse fastqs #####################################

if (1 %in% stages)
{
  print ("finding linkers")
  parseFastq(fastqs[1],fastqs[2],basename = outname,minlength = 15,
           maxlength = 25, keepempty=TRUE)
}

###################################### align reads #####################################

if (2 %in% stages)
{
  print ("aligning reads")
  # filenames
  fastq1 = paste(outname ,"_1.same.fastq",sep="")
  fastq2 = paste(outname ,"_2.same.fastq",sep="")
  sam1 = paste(outname ,"_1.same.sam",sep="")
  sam2 = paste(outname ,"_2.same.sam",sep="")
  
  # align both ends of each PET
  alignBowtie(fastq=fastq1,output=sam1,bowtiepath=bowtiepath,bowtieref=bowtieref)
  alignBowtie(fastq=fastq2,output=sam2,bowtiepath=bowtiepath,bowtieref=bowtieref)
}

##################################### filter reads #####################################

if (3 %in% stages)
{
  # filenames
  bedpefile          = paste(outname ,".bedpe",sep="")
  bedpefilesort      = paste(outname ,".sort.bedpe",sep="")
  bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")
  
  # build bedpe
  print ("building bedpe")
  buildBedpe(sam1 =sam1, sam2 = sam2, bedpefile = bedpefile);
  
  # sort bedpe
  print ("sorting bedpe")
  external_sort(bedpefile, bedpefilesort)
  
  # filter duplicates
  print ("removing PCR duplicates")
  removeDupBedpe(bedpefilesort,bedpefilesortrmdup,renamePets=TRUE);
}

##################################### call peaks #####################################

if (4 %in% stages)
{
  # filenames
  bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")
  tagAlignfile  = paste(outname,".tagAlign",sep="")
  peaksfile     = paste(outname,"_peaks.narrowPeak",sep="")
  peaksfileslop = paste(outname,"_peaks.slopPeak",sep="")
  
  # reverse strands for peak calling
  buildTagAlign(bedpefilesortrmdup ,tagAlignfile )
  
  # call peaks 
  callpeaks(macs2path,tagAlignfile,outname,pvalue=MACS_pvalue,
            bedtoolspath=bedtoolspath,bedtoolsgenome=bedtoolsgenome,peakslop=peakslop)
  
  # extend and merge peaks according to peakslop
  extendpeaks(peaksfile,peaksfileslop,bedtoolspath=bedtoolspath,bedtoolsgenome=bedtoolsgenome,peakslop=peakslop)
}

##################################### group pairs #####################################

if (5 %in% stages)
{
  
  # filenames
  peaksfile     = paste(outname,"_peaks.narrowPeak",sep="")
  peaksfileslop = paste(outname,"_peaks.slopPeak",sep="")
  bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")
  distancefile       = paste(outname ,".distance",sep="")
  distancecutpdf     = paste(outname ,".distance.pdf",sep="")
  
  # build a file of just distances and same / dif
  makeDistanceFile(bedpefilesortrmdup,distancefile,distcutrange[1],distcutrange[2])
  
  # calculate bias and cutoff
  distancecutoff = calcDistBias(distancefile,distancecutpdf=distancecutpdf,range=distcutrange,biascut=biascut)
  
  # group PETs into interactions
  chromosomes = groupPairs(bedpefilesortrmdup=bedpefilesortrmdup,outname=outname,
                           peaksfile=peaksfileslop,verbose=FALSE,distancecutoff=distancecutoff)
  
  # score pairs with P value
  scorePairs(chromosomes,outname=outname,min_distance = distancecutoff,maxPval=maxPval,
             numofbins = numofbins,binrange=binrange,
             corrMethod=corrMethod,verbose=TRUE)
}

Sys.time()
print("done")



# ##################################### rewire reads #####################################
# 
# print ("rewiring PETs")
# 
# # make a pdf to print results
# pdf(paste(outname ,".rw.results.pdf",sep=""),height=10.5,width=10.5)
# par(mfrow=c(3,3))
# 
# rwrpetsfile  = paste(outname,".rwr.bedpe",sep="")
# obspetsfile  = paste(outname,".obs.bedpe",sep="")
# rwrreadsfile = paste(outname,".rwr.bed",sep="")
# obsreadsfile = paste(outname,".obs.bed",sep="")
# if (file.exists(rwrpetsfile)) file.remove(rwrpetsfile)
# if (file.exists(obspetsfile)) file.remove(obspetsfile)
# if (file.exists(rwrreadsfile)) file.remove(rwrreadsfile)
# if (file.exists(obsreadsfile)) file.remove(obsreadsfile)
# 
# # split into reads
# firstsample = TRUE
# counter = 0
# for (chrom in filestorewire)
# {
#   f = paste(outname,".",chrom,sep="")
#   print (f)
#   readsfile   = paste(f,".peakfilt.bed",sep="")
#   petsfile    = paste(f,".peakfilt.bedpe",sep="")
#   densities   = rewire(readsfile,petsfile,obspetsfile,rwrpetsfile,obsreadsfile,rwrreadsfile,counter=counter,smallreps=20,bigreps=10)
# 
#   # keep track of read numbers
#   counter = densities[[4]]
#   
#   # remove files for individual chromosomes
#   #if (file.exists(readsfile)) file.remove(readsfile)
#   #if (file.exists(petsfile)) file.remove(petsfile)
#   
#   # plot the distributions
#   par(ann=F)
#   plot(densities[[1]],col='red',xaxt='n')
#   axis(side=1,at=seq(1,8,by=1),labels = 10^seq(1,8,by=1),las=2)
#   points(densities[[2]],col='black',type='l')
#   mtext("distance",side=1,line=3.5,font=2)
#   mtext(paste("Chromosome", densities[[3]]),side=3,line=1,font=1.0,cex=2)
#   legend("topright",inset = .05,legend=c("observed","rewired"),col=c("red","black"),pch=15)
# }
# dev.off()
# 
# 
# 
# ##################################### call pairs #####################################
## make a read file
#print ("splitting PETs into reads")
#PETsToReads(bedpefilesortrmdup,bedfilesortrmdup);

##################################### determine min distance #####################################
##################################### Split and Filter #####################################

# # split bed and bedpe by chrom
# print ("splitting PETs and reads by chromosome")
# chromSplitResults = splitBedpe(bedpefilesortrmdup,outname);
# # filter bed for those that overlap peaks
# filterPeakOverlap(chroms=filestorewire,outname=outname,peaksfileslop=peaksfileslop,
#                   bedtoolspath=bedtoolspath,bedpe=FALSE,verbose=TRUE)
# 
# # filter bedpe for those that overlap peaks
# filterPeakOverlap(chroms=filestorewire,outname=outname,peaksfileslop=peaksfileslop,
#                   bedtoolspath=bedtoolspath,bedpe=TRUE,verbose=TRUE)


