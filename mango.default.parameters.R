# Prefix to all files

###### GENERAL PARAMETERS ######

stages = (1:5)                            # stages of the pipeline to execute

###### STAGE 1 PARAMETERS ######

linkers   = c("GTTGGATAAG","GTTGGAATGT")  # linker sequences to look for
minlength = 15                            # min length of reads after linker trimming
maxlength = 25                            # max length of reads after linker trimming
keepempty = FALSE                         # Should reads with no linker be kept

###### STAGE 4 PARAMETERS ######

MACS_pvalue = .00001                      # MACS values
peakslop    = 500                         # Number of basespairs to extend peaks on both sides

###### STAGE 5 PARAMETERS ######

distcutrange   = c(1000,100000)           # range in which to look for the self-ligation distance
biascut = 0.05                            # Self ligation bias cutoff
maxPval = 0.01                            # P-value cutoff
numofbins = 50                            # number of bins for probability calculations
binrange=c(1000,250000000)                # range over which to construct bins for probability calculations
corrMethod = "BH"                         # multiple hypothesis tersting correction method
