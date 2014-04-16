#include <Rcpp.h>
using namespace Rcpp;


// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)
// For more on using Rcpp click the Help button on the editor toolbar

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <bitset>
#include <map>
#include "mergesort.h"
using namespace std;


// Define a function that joins vectors of strings
// [[Rcpp::export]]
std::string vector_join( const std::vector<std::string>& v, const std::string& token ){
    ostringstream result;
    for (std::vector<std::string>::const_iterator i = v.begin(); i != v.end(); i++){
        if (i != v.begin()) result << token;
        result << *i;
    }
    return result.str();
}

// Define a function that splits strings into vector
// [[Rcpp::export]]
std::vector<std::string> string_split( const std::string& s, const std::string& delimiter ){
    std::vector<std::string> result;
    std::string::size_type from = 0;
    std::string::size_type to = 0;
    
    while ( to != std::string::npos ){
        to = s.find( delimiter, from );
        if ( from < s.size() && from != to ){
            result.push_back( s.substr( from, to - from ) );
        }
        from = to + delimiter.size();
    }
    return result;
}

// [[Rcpp::export]]
std::string parseFastq(std::string fastq1, std::string fastq2,std::string basename,
              int minlength = 15,int maxlength = 25,
              bool keepempty = false, bool verbose = true,
              std::string linker1 = "GTTGGATAAG" , std::string linker2 = "GTTGGAATGT")
{

   // arguments
    ifstream file1(fastq1.c_str());
    ifstream file2(fastq2.c_str());
    ofstream same1 ( (basename + "_1.same.fastq").c_str() );
    ofstream same2 ( (basename + "_2.same.fastq").c_str() );
    ofstream chim1 ( (basename + "_1.chim.fastq").c_str() );
    ofstream chim2 ( (basename + "_2.chim.fastq").c_str() );
    
    // define variables
    std::string fqline1;
    std::string fqline2;
    int linecount = 0;
    int i = 0;
    std::vector<std::string> lines1;
    std::vector<std::string> lines2;
    
    while (getline(file1, fqline1))
    {
        // read lines and increment counters
        getline(file2, fqline2);
        i++;
        linecount++;

        // add lines to list
        lines1.push_back(fqline1);
        lines2.push_back(fqline2);
        
        // if list length is 4 perform operations, print, and clear lists
        if ( i == 4 )
        {
            
            // find the position of the linkers
            size_t r1l1found = lines1[1].find(linker1);
            size_t r1l2found = lines1[1].find(linker2);
            size_t r2l1found = lines2[1].find(linker1);
            size_t r2l2found = lines2[1].find(linker2);
            
            // determine the linker type (0 = none, 3 = both)
            // read 1
            int r1linker;
            if (r1l1found == -1 & r1l2found == -1)
            {
                r1linker = 0;
            }
            else if (r1l1found != -1 & r1l2found == -1)
            {
                r1linker = 1;
                lines1[1] =  lines1[1].substr (0,r1l1found);
                lines1[3] =  lines1[3].substr (0,r1l1found);
            }
            else if (r1l1found == -1 & r1l2found != -1)
            {
                r1linker = 2;
                lines1[1] =  lines1[1].substr (0,r1l2found);
                lines1[3] =  lines1[3].substr (0,r1l2found);
            }
            else if (r1l1found != -1 & r1l2found != -1)
            {
                r1linker = 3;
            }
            
            // read 2
            int r2linker;
            if (r2l1found == -1 & r2l2found == -1)
            {
                r2linker = 0;
            }
            else if (r2l1found != -1 & r2l2found == -1)
            {
                r2linker = 1;
                lines2[1] =  lines2[1].substr (0,r2l1found);
                lines2[3] =  lines2[3].substr (0,r2l1found);
            }
            else if (r2l1found == -1 & r2l2found != -1)
            {
                r2linker = 2;
                lines2[1] =  lines2[1].substr (0,r2l2found);
                lines2[3] =  lines2[3].substr (0,r2l2found);
            }
            else if (r2l1found != -1 & r2l2found != -1)
            {
                r2linker = 3;
            }
            
            // determine pairtype
            std::string pairtype = "unknown";
            if ((r1linker == 1 && r2linker == 1) || (r1linker == 2 && r2linker == 2))
            {
                    pairtype = "same";
            }
            if ((r1linker == 1 && r2linker == 2) || (r1linker == 2 && r2linker == 1))
            {
                pairtype = "chim";
            }
            if (r1linker == 3 || r2linker == 3)
            {
                pairtype = "ambi";
            }
            if (keepempty == true)
            {
                if ((r1linker == 0 && r2linker == 1) ||
                    (r1linker == 0 && r2linker == 2) ||
                    (r1linker == 1 && r2linker == 0) ||
                    (r1linker == 2 && r2linker == 0) ||
                    (r1linker == 0 && r2linker == 0))
                {
                    pairtype = "same";
                }
            }
            if (keepempty == false)
            {
                if (r1linker == 0 || r2linker == 0)
                {
                    pairtype = "ambi";
                }
            }
            
            // determine if they pass the size requirements and print to output
            if ((lines1[1].length() >= minlength ) &&  (lines1[1].length() <= maxlength ) &&
                (lines2[1].length() >= minlength ) &&  (lines2[1].length() <= maxlength ))
            {
                if (pairtype == "same")
                {
                    same1 << vector_join(lines1,"\n");
                    same2 << vector_join(lines2,"\n");
                    same1 << "\n";
                    same2 << "\n";

                }
                if (pairtype == "chim")
                {
                    chim1 << vector_join(lines1,"\n");
                    chim2 << vector_join(lines2,"\n");
                    chim1 << "\n";
                    chim2 << "\n";
                }
            }
            
            
            // reset lines
            i = 0;
            lines1.clear();
            lines2.clear();

        }
        
        // num % 2 computes the remainder when num is divided by 2
        if ( linecount % 1000000 == 0 )
        {
            cout << linecount;
            cout << "\n";
            
        }

    }
    

    
    // close streams
    same1.close();
    same2.close();
    chim1.close();
    chim2.close();
    file1.close();
    file2.close();
     
    return "linker parsing complete";
}


// Define a function that converts bitflag to 0s and 1s (in reverse order for some reason)
std::vector< int > get_bits( unsigned long x ) {
    std::string chars( std::bitset< sizeof(long) * CHAR_BIT >( x )
                      .to_string( char(0), char(1) ) );
    return std::vector< int >( chars.begin(), chars.end() );
}

// Define a function that returns the strand
std::string get_strand( std::vector< int > bitflag, int strandbit = 4 ) {
    std::string strand = "+";
    if (bitflag[63 - strandbit] == 0)
    {
        strand = "-";
    }
    return strand;
}

// Define a function that converts string to int
int StringToInt( std::string Text ) {
    int output;
    if ( ! (istringstream(Text) >> output) ) output = 0;
    return output;
}

// Define a function that converts string to int
std::string IntToString( int Number ) {
    std::string Result;          // string which will contain the result
    ostringstream convert;   // stream used for the conversion
    convert << Number;      // insert the textual representation of 'Number' in the characters in the stream
    Result = convert.str(); // set 'Result' to the contents of the stream
    return Result;
}


// Define a function that builds a bedpe file rom 2 sam file
// [[Rcpp::export]]
void buildBedpe(std::string sam1, std::string sam2,std::string bedpefile)
{
    
    // arguments
    ifstream file1(sam1.c_str());
    ifstream file2(sam2.c_str());
    ofstream bedpefilestream ( bedpefile.c_str() );

    // define variables
    std::string line1;
    std::string line2;
    int linecount = 0;
    
    while (getline(file1, line1))
    {
        // read lines and increment counter
        getline(file2, line2);
        linecount++;
        
        // split lines
        std::vector<std::string> e1 = string_split(line1,"\t");
        std::vector<std::string> e2 = string_split(line2,"\t");
        
        // get info for file 1
        std::string name1 = e1[0];
        name1 = string_split(name1,"_")[0];
        name1 = string_split(name1," ")[0];
        int bitflag1 = StringToInt(e1[1]);
        std::vector<int> bits1 = get_bits(bitflag1);
        std::string strand1 = get_strand(bits1);
        std::string sequence1 = e1[9];
        std::string chrom1 = e1[2];
        int start1 = StringToInt(e1[3]) -1;
        int stop1  = start1 + sequence1.length();
        
        // get info for file 2
        std::string name2 = e2[0];
        name2 = string_split(name2,"_")[0];
        name2 = string_split(name2," ")[0];
        int bitflag2 = StringToInt(e2[1]);
        std::vector<int> bits2 = get_bits(bitflag2);
        std::string strand2 = get_strand(bits2);
        std::string sequence2 = e2[9];
        std::string chrom2 = e2[2];
        int start2 = StringToInt(e2[3]) -1;
        int stop2  = start2 + sequence2.length();
        
        // check that read names match
        if (name1 != name2)
        {
            cout << "Error: read names of PET ends do not match";
            break;
        }
        
        // determine which read goes first
        bool reorder = false;
        if ((chrom1 == chrom2) & (start1 > start2) )
        {
            reorder = true;
        }
        if ((chrom1 != chrom2) & (chrom1 > chrom2) )
        {
            reorder = true;
        }
        
        // print out results
        if (reorder == false)
        {
            std::vector<std::string> outputvector;
            outputvector.push_back(chrom1);
            outputvector.push_back(IntToString(start1));
            outputvector.push_back(IntToString(stop1));
            outputvector.push_back(chrom2);
            outputvector.push_back(IntToString(start2));
            outputvector.push_back(IntToString(stop2));
            outputvector.push_back(name1);
            outputvector.push_back(".");
            outputvector.push_back(strand1);
            outputvector.push_back(strand2);
            std::string outputstring = vector_join(outputvector,"\t");
            bedpefilestream << outputstring;
            bedpefilestream << "\n";
        }
        
        if (reorder == true)
        {
            std::vector<std::string> outputvector;
            outputvector.push_back(chrom2);
            outputvector.push_back(IntToString(start2));
            outputvector.push_back(IntToString(stop2));
            outputvector.push_back(chrom1);
            outputvector.push_back(IntToString(start1));
            outputvector.push_back(IntToString(stop1));
            outputvector.push_back(name1);
            outputvector.push_back(".");
            outputvector.push_back(strand2);
            outputvector.push_back(strand1);
            std::string outputstring = vector_join(outputvector,"\t");
            bedpefilestream << outputstring;
            bedpefilestream << "\n";
        }
    }
}


// Define a function that removes duplicates from a bedpe file
// [[Rcpp::export]]
void removeDupBedpe(std::string infile,std::string outfile)
{
    // arguments
    ifstream file1(infile.c_str());
    ofstream bedpefilestream (outfile.c_str());
    
    // read in file line by line store currentline and last line
    std::string lastline;
    std::string currline;
    while (getline(file1, currline))
    {
        // split lines
        std::vector<std::string> lastEall = string_split(lastline,"\t");
        std::vector<std::string> currEall = string_split(currline,"\t");
        std::vector<std::string> lastE;
        std::vector<std::string> currE;
        
        // print first line
        if (lastEall.size() == 0)
        {
            bedpefilestream << currline;
            bedpefilestream << "\n";
            lastline = currline;
            continue;
        }
        
        // make shorter line for comparison
        for (int i=0; i<6; i++)
        {
            lastE.push_back(lastEall[i]);
            currE.push_back(currEall[i]);
        }
        std::string currlineshort = vector_join(currE,"_");
        std::string lastlineshort = vector_join(lastE,"_");
        
        // if duplicates continue
        if (lastlineshort == currlineshort)
        {
            lastline = currline;
            continue;
        }
        
        // print out non duplicates
        bedpefilestream << currline;
        bedpefilestream << "\n";
        
        // update last line
        lastline = currline;
    }
}






//// Define a function that puts pairs together
//void findPairs(std::string overlapfile, std::string petpairsfile,std::string peakpairsfile)
//{
//  
//}
//


// Define a class to keep track of peak information
class peak{
public:
  std::string name;
  std::string chrom;
  int start;
  int end;
  int intra;
};

// Define a class to keep track of pair information
class chiapair{
public:
  std::string pairname;
  std::string p1name;
  std::string p2name;
  std::string p1chrom;
  std::string p2chrom;
  int p1start;
  int p1end;
  int p2start;
  int p2end;
  int p1intra;
  int p2intra;
  int linking;
  int distance;
};



// Define a function that puts pairs together
// [[Rcpp::export]]
void findPairs(std::string overlapfile, std::string petpairsfile,std::string interactionfile)
{
  // (1) Read in overlap info
  
  // streams
  ifstream input(overlapfile.c_str());

  std::string line;
  //std::map<std::string, peak> peakinfodict;
  std::map<std::string, peak> peakinfodict;
  std::map<std::string, std::vector<std::string> > readpeakdict;
  
  // read in file line by line store currentline and last line
    while (getline(input, line))
    {
        // split lines
        std::vector<std::string> currEall = string_split(line,"\t");
        string readname = currEall[3];
        string peakname = currEall[9];
      
        // add info to peak dict
        if (peakinfodict.find(peakname) == peakinfodict.end())
        {
            peak p = *(new peak());
            peakinfodict.insert(std::pair<string,peak>(peakname,p));
            peakinfodict[peakname].name  = currEall[9];
            peakinfodict[peakname].chrom = currEall[6];
            peakinfodict[peakname].start = atoi( currEall[7].c_str() );
            peakinfodict[peakname].end   = atoi( currEall[8].c_str() );
            peakinfodict[peakname].intra = 0;
        }
        peakinfodict[peakname].intra++;
        
        // add info to readpeak dict
        if (readpeakdict.find(peakname) == readpeakdict.end())
        {
            std::vector<std::string> v = *(new std::vector<std::string>);
            readpeakdict.insert(std::pair<string,std::vector<std::string> > (readname, v));
        }
        readpeakdict[readname].push_back(peakname);
    }
  input.close();
  
  // (1) Go through PETs and make interactions
  ifstream inputpets(petpairsfile.c_str());

  std::string line2;
  std::map<std::string, chiapair> pairdict;
  
  // read in file line by line store currentline and last line
  while (getline(inputpets, line2))
  {
      // split lines
      std::vector<std::string> currEall = string_split(line2,"\t");
      string readname  = currEall[6];
      string r1 = readname + ".1";
      string r2 = readname + ".2";
      std::vector<std::string> p1s = readpeakdict[r1];
      std::vector<std::string> p2s = readpeakdict[r2];
      
      // interate through all combinations of peaks
      for (std::vector<std::string>::iterator p1 = p1s.begin() ; p1 != p1s.end() ; ++p1)
      {
        for (std::vector<std::string>::iterator p2 = p2s.begin() ; p2 != p2s.end() ; ++p2)
        {
          peak thep1 = peakinfodict[*p1];
          peak thep2 = peakinfodict[*p2];
          
          // swith the order of the peaks based on start position
          if (peakinfodict[*p1].start > peakinfodict[*p2].start   )
          {
            thep2 = peakinfodict[*p1];
            thep1 = peakinfodict[*p2];
          }

          // join the peak names for the name of the pair
          std::string pairname =  thep1.name + "_" + thep2.name;

          // add info to pair dict
          if (pairdict.find(pairname) == pairdict.end())
          {
            chiapair pa = *(new chiapair());
            pairdict.insert(std::pair<string,chiapair>(pairname,pa));
            pairdict[pairname].pairname  = pairname; 
            pairdict[pairname].p1name  = thep1.name; 
            pairdict[pairname].p2name  = thep2.name; 
            pairdict[pairname].p1chrom  = thep1.chrom; 
            pairdict[pairname].p2chrom  = thep2.chrom; 
            pairdict[pairname].p1start  = thep1.start; 
            pairdict[pairname].p1end  = thep1.end;
            pairdict[pairname].p2start  = thep2.start; 
            pairdict[pairname].p2end  = thep2.end; 
            pairdict[pairname].p1intra  = thep1.intra; 
            pairdict[pairname].p2intra  = thep2.intra; 
            pairdict[pairname].linking  = 0; 
            pairdict[pairname].distance  = thep2.start - thep1.end; 
          }

          //   STOP POINT 
           
           cout << IntToString(thep1.start) + " " + IntToString(thep1.end) + " " +
                   IntToString(thep2.start) + " " + IntToString(thep2.end);
          
         // cout << line2 + "   "  + readname + "  " + *p1 + "   " +  *p2;
          cout << "\n";
        }
      }
      
  } 
  inputpets.close();
  
  // print out info to pair file
  
  
  
  
  
}



// Define a function splits bed file by chromosome
// [[Rcpp::export]]
std::vector<std::string> splitBedbyChrom(std::string bedfile,std::string outnamebase)
{
    
    // streams
    ifstream file1 (bedfile.c_str());
    std::map<std::string, std::ofstream*> readoutput;
    
    std::string line;
    // read in file line by line store currentline and last line
    while (getline(file1, line))
    {
        // split lines
        std::vector<std::string> currEall = string_split(line,"\t");
        
        // add output string to dcit if neccesary
        std::string chrom = currEall[0];
        if ( (readoutput.find(chrom) == readoutput.end()) & (chrom != "*" )  ) {
            std::string outname = outnamebase + "." + chrom  + ".bed";
            readoutput[chrom] = new std::ofstream(outname.c_str());
        }
  
        // print reads
        if (chrom != "*")
        {
            *readoutput[chrom] << line;
            *readoutput[chrom] << "\n";
        }
    }
    
    
    // close reads files
    std::vector<std::string> chromosomes;
    for (std::map<std::string, std::ofstream*>::iterator i = readoutput.begin() ; i != readoutput.end() ; i ++ ) {
      i->second->close();
      chromosomes.push_back(i->first);
    }
    
    return (chromosomes);
}


// Define a function splits bedpe file into reads and PETs by chromosome
// [[Rcpp::export]]
std::vector<std::string> splitBedpe(std::string bedpein,std::string outnamebase, bool printreads = true , bool printpets = true)
{
    // keep track of output files
    std::vector<std::string> outputvector;
    
    // streams
    ifstream file1(bedpein.c_str());
    std::map<std::string, std::ofstream*> readoutput;
    std::map<std::string, std::ofstream*> petsoutput;
    
    std::string line;
    // read in file line by line store currentline and last line
    while (getline(file1, line))
    {
        // split lines
        std::vector<std::string> currEall = string_split(line,"\t");
        
        // make reads
        std::vector<std::string> read1vec;
        read1vec.push_back(currEall[0]);
        read1vec.push_back(currEall[1]);
        read1vec.push_back(currEall[2]);
        read1vec.push_back(currEall[6]);
        read1vec.push_back(currEall[7]);
        read1vec.push_back(currEall[8]);
        
        std::vector<std::string> read2vec;
        read2vec.push_back(currEall[3]);
        read2vec.push_back(currEall[4]);
        read2vec.push_back(currEall[5]);
        read2vec.push_back(currEall[6]);
        read2vec.push_back(currEall[7]);
        read2vec.push_back(currEall[9]);
        
        std::string read1 = vector_join(read1vec,"\t");
        std::string read2 = vector_join(read2vec,"\t");
        
        // get chromosome of each read
        std::string chrom1 = currEall[0];
        std::string chrom2 = currEall[3];
        
        if (printreads == true)
        {
          // print reads
          
          if ( (readoutput.find(chrom1) == readoutput.end()) & (chrom1 != "*" )  ) {  
              std::string outname = outnamebase + "." + chrom1 + ".bed";
              readoutput[chrom1] = new std::ofstream(outname.c_str());
          }
          
          if ( (readoutput.find(chrom2) == readoutput.end()) & (chrom2 != "*" ) ) {
              std::string outname = outnamebase + "." + chrom2 + ".bed";
              readoutput[chrom2] = new std::ofstream(outname.c_str());
          }

          // print reads
          if (chrom1 != "*")
          {
              *readoutput[chrom1] << read1;
              *readoutput[chrom1] << "\n";
          }
          if (chrom2 != "*")
          {
              *readoutput[chrom2] << read2;
              *readoutput[chrom2] << "\n";
          }
        }
        
        if ((chrom1 == "*") | (chrom2 == "*"))
        {
            continue;
        }
        
        if ( petsoutput.find(chrom1) == petsoutput.end() ) {
            
            std::string outname = outnamebase + "." + chrom1 + ".bedpe";
            petsoutput[chrom1] = new std::ofstream(outname.c_str());
            outputvector.push_back( chrom1);
            
        }
        *petsoutput[chrom1] << line;
        *petsoutput[chrom1] << "\n";
    }
    
    // close reads files
    for (std::map<std::string, std::ofstream*>::iterator i = readoutput.begin() ; i != readoutput.end() ; i ++ ) {
      i->second->close();
    }
    
    // close bedpe files
    for (std::map<std::string, std::ofstream*>::iterator i = petsoutput.begin() ; i != petsoutput.end() ; i ++ ) {
      i->second->close();
    }    

    return outputvector;
}

// [[Rcpp::export]]
void buildTagAlign(std::string bedpefile, std::string TagAlignfile) {
    
    // establish streams
    ifstream infile (bedpefile.c_str());
    ofstream outfile ( TagAlignfile.c_str() );

   int i = 0;
   std::string line;
   while (getline(infile, line))
    {
        // increment counters
        i++;
        
        // split line by tab
         std::vector<std::string> e = string_split(line,"\t");
        
        // reverse strands
        std::string newstrand1 = "-";
        if (e[8] == "-")
        {
          newstrand1 = "+";
        }
        e[8] = newstrand1;
        
        std::string newstrand2 = "-";
        if (e[9] == "-")
        {
          newstrand2 = "+";
        }
        e[9] = newstrand2;
        
        // print to output
        if (e[0] != "*")
        {
          std::vector<std::string> outputvector;
          outputvector.push_back(e[0]);
          outputvector.push_back(e[1]);
          outputvector.push_back(e[2]);
          outputvector.push_back(e[6]);
          outputvector.push_back(".");
          outputvector.push_back(e[8]);
          std::string outputstring = vector_join(outputvector,"\t");
          outfile << outputstring;
          outfile << "\n";
        }
        
        if (e[3] != "*")
        {
          std::vector<std::string> outputvector2;
          outputvector2.push_back(e[3]);
          outputvector2.push_back(e[4]);
          outputvector2.push_back(e[5]);
          outputvector2.push_back(e[6]);
          outputvector2.push_back(".");
          outputvector2.push_back(e[9]);
          std::string outputstring2 = vector_join(outputvector2,"\t");
          outfile << outputstring2;
          outfile << "\n";
        }
    }
    
    // close streams
    infile.close();
    outfile.close();
}

// Define a function to do an external sort
// [[Rcpp::export]]
void external_sort( std::string inputfile, std::string outputfile ){
    externalMergesort <string> externalMergeSorter(inputfile, outputfile, 1000000);
    return;

}





