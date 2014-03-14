#' runs mango chia pet analysis pipeline
#'
#' This function find and trims linkers from fastq files and output 'same' and 'chim' files
#'  
#' @param arg1 first argument
#' @export
#' @examples
#' 
 
##################################### load libraries #####################################

library('logspline')

##################################### set var #####################################


set.seed(1)
bowtiepath = ""
bowtieref  = ""

##################################### parse fastqs #####################################

fastqs = c("/Volumes/HD3/projects/newmango/data/NH.K562_RAD21_K562_std_2.1_1.head.fastq",
           "/Volumes/HD3/projects/newmango/data/NH.K562_RAD21_K562_std_2.1_2.head.fastq")

fastqs = c("~/Desktop/mango2014test/NH.K562_RAD21_K562_std_2.1_1.head.fastq",
           "~/Desktop/mango2014test/NH.K562_RAD21_K562_std_2.1_2.head.fastq")

expname = "NH.K562_RAD21_K562_std_2.1"
outdir = NULL
outname="mango"
fastqsin = fastqs
fastqsout = NULL
linkers = c("GTTGGATAAG","GTTGGAATGT")
minlength = 15
maxlength = 25
keepempty=FALSE
linesatatime = 1000

parseFastq(fastqsin,outdir=NULL,outname=outname,linkers = linkers,
          minlegnth = 15,maxlength = 25,keepempty=TRUE,linesatatime = 1000)


##################################### align reads #####################################

fastq1 = "~/Desktop/mango2014test/mango_1.same.fastq"
fastq2 = "~/Desktop/mango2014test/mango_2.same.fastq"
sam1 = "~/Desktop/mango2014test/mango_1.same.sam"
sam2 = "~/Desktop/mango2014test/mango_2.same.sam"

bowtiepath = "/Users/dougphanstiel/Tools/bowtie-1.0.0/bowtie"
bowtieref  = "/Users/dougphanstiel/Tools/bowtie-1.0.0/indexes/hg19"

# align both ends of each PET
alignBowtie(fastq=fastq1,output=sam1,bowtiepath=bowtiepath,bowtieref=bowtieref)
alignBowtie(fastq=fastq2,output=sam2,bowtiepath=bowtiepath,bowtieref=bowtieref)


##################################### filter reads #####################################

samfilt1 = "~/Desktop/mango2014test/mango_1.same.filt.sam"
samfilt2 = "~/Desktop/mango2014test/mango_2.same.filt.sam"
bedpebase = "~/Desktop/mango2014test/mango"
path     = "~/Desktop/mango2014test/"

# filter for PETs that have both ends uniquely mapped
filterReads(sam1=sam1,sam2=sam2,samfilt1=samfilt1,samfilt1=samfilt2)

# build bedpes
buildBedpes(samfilt1=samfilt1,samfilt2=samfilt2,bedpebase)

# remove duplicates
removeDups(path=path,expname="mango")

# split into reads
petsToReads(path=path,expname="mango")

##################################### rewire reads #####################################

# split into reads
rewire(path=path,expname="mango")

##################################### call peaks #####################################


##################################### group pairs #####################################





##################################### call pairs #####################################




