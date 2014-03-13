#' runs mango chia pet analysis pipeline
#'
#' This function find and trims linkers from fastq files and output 'same' and 'chim' files
#'  
#' @param arg1 first argument
#' @export
#' @examples
#' 
 
##################################### load libraries #####################################



##################################### set var #####################################

bowtiepath = ""
bowtieref  = ""

##################################### parse fastqs #####################################

fastqs = c("/Volumes/HD3/projects/newmango/data/NH.K562_RAD21_K562_std_2.1_1.head.fastq",
           "/Volumes/HD3/projects/newmango/data/NH.K562_RAD21_K562_std_2.1_2.head.fastq")

outname = "NH.K562_RAD21_K562_std_2.1"
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

bowtiepath = ""
bowtieref  = ""

# align both ends of each PET
alignBowtie(fastq=fastq1,bowtiepath=bowtiepath,bowtieref=bowtieref)
alignBowtie(fastq=fastq2,bowtiepath=bowtiepath,bowtieref=bowtieref)


##################################### filter reads #####################################




##################################### shuffle reads #####################################





##################################### call peaks #####################################

> dir.create("results") # Note: all output data will be written to directory ✬results✬
> buildindex(basename="./results/tair10chr.fasta", reference="./data/tair10chr.fasta") # Build indexed reference genome
> targets <- read.


##################################### group pairs #####################################





##################################### call pairs #####################################




