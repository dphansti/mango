readParameters <-function(args=hash(),argsfile=NULL,verbose=FALSE)
{
  
  if (is.null(argsfile) == TRUE)
  {
    lines = c(
      
      "###### GENERAL PARAMETERS ######",
      
      "stages            = 1:5                   # stages of the pipeline to execute",
      "prefix            = mango                 # prefix for all output files",
      "outdir            = NULL",
      "chrominclude      = NULL                    # comma separated list of chromosomes to use (e.g. chr1,chr2,chr3,...).  Only these chromosomes will be processed",
      "chromexclude      = NULL                    # comma separated list of chromosomes to exclude (e.g. chrM,chrY)",
      
      "###### STAGE 1 PARAMETERS ######",
      
      "linkerA           = GTTGGATAAG              # linker sequences to look for",
      "linkerB           = GTTGGAATGT              # linker sequences to look for",
      "minlength         = 15                      # min length of reads after linker trimming",
      "maxlength         = 25                      # max length of reads after linker trimming",
      "keepempty         = FALSE                   # Should reads with no linker be kept",
      
      "###### STAGE 2 PARAMETERS ######",
      
      "shortreads           = TRUE                 # should bowtie alignments be done using paramter for very short reads (~20 bp)", 
      
      "###### STAGE 4 PARAMETERS ######",
      
      "MACS_pvalue       = 0.00001                  # MACS values",
      "peakslop          = 500                      # Number of basespairs to extend peaks on both sides",
      "peakinput         = NULL                     # name of user supplied peaks file",
      
      "###### STAGE 5 PARAMETERS ######",
      
      "distcutrangemin   = 1000                     # range in which to look for the self-ligation distance",
      "distcutrangemax   = 100000                   # range in which to look for the self-ligation distance",
      "biascut           = 0.05                     # Self ligation bias cutoff",
      "maxPval           = 0.01                     # P-value cutoff",
      "numofbins         = 30                       # number of bins for probability calculations",
      "corrMethod        = BY                       # multiple hypothesis tersting correction method",
      "maxinteractingdist= 10000000                 # maximum disance allowed for an interaction",
      "FDR               = 0.01                     # FDR cutoff for interactions",
      "minPETS           = 2                        # minimum number of PETs required for an interaction (applied after FDR filtering)",
      "reportallpairs    = FALSE                    # Should all pairs be reported or just significant pairs"
    )
  }
  
  if (is.null(argsfile) == FALSE)
  {
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
      args[[arginfo[1]]] = arginfo[2]
      if (verbose == TRUE)
      {
        print(paste(arginfo[1] , "        " , arginfo[2]))
      }
    }
  }
  return (args)
}

  
  
  