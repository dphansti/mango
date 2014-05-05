# Define a function that checks stage for neccesary arguments
checkRequired <-function(args,listofargs)
{
  for (key in listofargs )
  {
    if (has.key(key,args) == FALSE)
    {
      print (paste("missing required argument:",key))
      break
    }
  }  
}
