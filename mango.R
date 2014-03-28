# runs mango chia pet analysis pipeline
#location = "work"

##################################### paths to externals #####################################

bowtiepath = "/Users/dougphanstiel/Tools/bowtie-1.0.0/bowtie"
bowtieref  = "/Users/dougphanstiel/Tools/bowtie-1.0.0/indexes/hg19"
#macs2path = 

fastqs = c("data/NH.K562_RAD21_K562_std_2.1_1.head.fastq",
           "data/NH.K562_RAD21_K562_std_2.1_2.head.fastq")
expname = "NH.K562_RAD21_K562_std_2.1.head"
outdir   = "/Users/dougphanstiel/Desktop/mango2014test/"
outname  = paste(outdir,expname,sep="")
linkers = c("GTTGGATAAG","GTTGGAATGT")
minlength = 15
maxlength = 25
keepempty=FALSE

# if (location == "home")
# {
#   bowtiepath = "/Users/dougphanstiel/Tools/bowtie-1.0.0/bowtie"
#   bowtieref  = "/Users/dougphanstiel/Tools/bowtie-1.0.0/indexes/hg19"
#   fastqs = c("/Users/dougphanstiel/Desktop/mango2014test/NH.K562_RAD21_K562_std_2.1_1.fastq",
#              "/Users/dougphanstiel/Desktop/mango2014test/NH.K562_RAD21_K562_std_2.1_2.fastq")
#   expname = "NH.K562_RAD21_K562_std_2.1"
#   outdir   = "/Users/dougphanstiel/Desktop/mango2014test/"
#   outname  = paste(outdir,expname,sep="")
#   linkers = c("GTTGGATAAG","GTTGGAATGT")
#   minlength = 15
#   maxlength = 25
#   keepempty=FALSE
# }
# if (location == "work")
# {
#   bowtiepath = "/Users/dougphanstiel/tools/bowtie-0.12.7/bowtie"
#   bowtieref  = "/Users/dougphanstiel/tools/bowtie-0.12.7/indexes/hg19"
#   macs2path = "/Users/dougphanstiel/tools/macs2/bin/macs2"
#   fastqs = c("/Volumes/HD3/projects/newmango/data/NH.K562_RAD21_K562_std_2.1_1.head.fastq",
#              "/Volumes/HD3/projects/newmango/data/NH.K562_RAD21_K562_std_2.1_2.head.fastq")
#   expname = "NH.K562_RAD21_K562_std_2.1"
#   outdir   = "/Volumes/HD3/projects/newmango/data/"
#   outname  = paste(outdir,expname,sep="")
#   linkers = c("GTTGGATAAG","GTTGGAATGT")
#   minlength = 15
#   maxlength = 25
#   keepempty=FALSE
# }

##################################### initialization #####################################

library('mango')
print ("Starting mango ChIA PET analysis tool")
Sys.time()
set.seed(1)

##################################### parse fastqs #####################################

print ("finding linkers")
parseFastq(fastqs[1],fastqs[2],basename = outname,minlength = 15,
           maxlength = 25, keepempty=TRUE)

###################################### align reads #####################################

print ("aligning reads")
# filenames
fastq1 = paste(outname ,"_1.same.fastq",sep="")
fastq2 = paste(outname ,"_2.same.fastq",sep="")
sam1 = paste(outname ,"_1.same.sam",sep="")
sam2 = paste(outname ,"_2.same.sam",sep="")

# align both ends of each PET
alignBowtie(fastq=fastq1,output=sam1,bowtiepath=bowtiepath,bowtieref=bowtieref)
alignBowtie(fastq=fastq2,output=sam2,bowtiepath=bowtiepath,bowtieref=bowtieref)

##################################### filter reads #####################################


# filenames
bedpefile          = paste(outname ,".bedpe",sep="")
bedpefilesort      = paste(outname ,".sort.bedpe",sep="")
bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")

# build bedpe
print ("building bedpe")
buildBedpe(sam1 =sam1, sam2 = sam2, bedpefile = bedpefile);

# sort bedpe (look into the C++ library STXXL to avoid using unix sort command)
print ("sorting bedpe")
command = paste("cat",bedpefile, " | sort > ",bedpefilesort)
system (command)

# filter duplicates
print ("removing PCR duplicates")
removeDupBedpe(bedpefilesort,bedpefilesortrmdup);

# separate by chrom (making reads and pets files)
print ("splitting PETs and reads by chromosome")
filestorewire = splitBedpe(bedpefilesortrmdup,outname);

##################################### rewire reads #####################################

print ("rewiring PETs")

# make a pdf to print results
pdf(paste(outname ,".rw.results.pdf",sep=""),height=10.5,width=10.5)
par(mfrow=c(3,3))

rwrpetsfile  = paste(outname,".rwr.bedpe",sep="")
obspetsfile  = paste(outname,".obs.bedpe",sep="")
if (file.exists(rwrpetsfile)) file.remove(rwrpetsfile)
if (file.exists(obspetsfile)) file.remove(obspetsfile)

# split into reads
for (f in filestorewire)
{
  readsfile   = paste(f,".bed",sep="")
  petsfile  = paste(f,".bedpe",sep="")
  densities   = rewire(readsfile,petsfile,obspetsfile,rwrpetsfile)
  
  par(ann=F)
  plot(densities[[1]],col='red',xaxt='n')
  axis(side=1,at=seq(1,8,by=1),labels = 10^seq(1,8,by=1),las=2)
  points(densities[[2]],col='black',type='l')
  mtext("distance",side=1,line=3.5,font=2)
  mtext(paste("Chromosome", densities[[3]]),side=3,line=1,font=1.0,cex=2)
  legend("topright",inset = .05,legend=c("observed","rewired"),col=c("red","black"),pch=15) 
}
dev.off()

##################################### call peaks #####################################

#callpeaks()





##################################### group pairs #####################################





##################################### call pairs #####################################

Sys.time()


