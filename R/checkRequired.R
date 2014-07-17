# Define a function that checks stage for neccesary arguments
checkRequired <-function(opt,listofargs)
{
  for (key in listofargs )
  {
    if ( as.character(opt[key]) == "NULL")
    {
      print (paste("missing required argument:",key))
      stop ("Exiting Mango.R pipeline.  Check required arguments")
      
    }
  }  
}

