defaultvariables <-function(file="mango.default.parameters")
{
  # Prefix to all files
  expname = "NH.K562_RAD21_K562_std_2.1"
  
  # linker sequences to look for
  linkers = c("GTTGGATAAG","GTTGGAATGT")
  
  # min and max length of reads after linker trimming
  minlength = 15
  maxlength = 25
  
  # Should reads with no linker be kept
  keepempty=FALSE
  
  # MACS values
  MACS_pvalue = .00001
  
  # Number of basespairs to extend peaks on both sides
  peakslop = 500
  
  # range in which to look for the self-ligation distance
  distcutrange   = c(1000,100000)
  
  # Self ligation bias cutoff
  biascut = 0.05
  
  # P-value cutoff
  maxPval = 0.01
  
  # binning for probability calcualtions
  numofbins = 50
  binrange=c(1000,250000000)
  
  # multiple hypothesis tersting correction method
  corrMethod = "BH"
  
  # stages
  stages = (1:5)
}