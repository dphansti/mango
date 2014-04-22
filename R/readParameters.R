readParameters <-function(file="mango.default.parameters.R")
{
  lines = readLines(file)
  
  for (line in lines)
  {
    if (line == "")
    {
      print("2")
    }
    else
    {
    print("1")
    }
  }
}