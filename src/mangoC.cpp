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
std::vector< int > parseFastq(std::string fastq1, std::string fastq2,std::string basename,
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
    
    // keep track of PET types
    int samecount = 0;
    int chimcount = 0;
    int ambicount = 0;
    
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
            
            // add to counters
            if (pairtype == "same")
            {
              samecount++;
            }
            if (pairtype == "chim")
            {
              chimcount++;
            }
            if (pairtype == "ambi")
            {
              ambicount++;
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
     
    // report results
    std::vector< int > parsingresults;
    parsingresults.push_back(samecount);
    parsingresults.push_back(chimcount);
    parsingresults.push_back(ambicount);

    return parsingresults;
}

// Define a function that returns the strand
std::string get_strand( unsigned long x ) {
    std::string strand = "+";
    if ( x & 0x10 )
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

// Define a function that converts int to string 
std::string IntToString( int Number ) {
    std::string Result;          // string which will contain the result
    ostringstream convert;   // stream used for the conversion
    convert << Number;      // insert the textual representation of 'Number' in the characters in the stream
    Result = convert.str(); // set 'Result' to the contents of the stream
    return Result;
}


// Define a function that converts number to string
template <typename T>
std::string NumberToString ( T Number )
{
  stringstream ss;
	ss << Number;
	return ss.str();
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
        name1 = string_split(name1,"#")[0];
        int bitflag1 = StringToInt(e1[1]);
        std::string strand1 = get_strand(bitflag1);
        std::string sequence1 = e1[9];
        std::string chrom1 = e1[2];
        int start1 = StringToInt(e1[3]) -1;
        int stop1  = start1 + sequence1.length();
        
        // get info for file 2
        std::string name2 = e2[0];
        name2 = string_split(name2,"_")[0];
        name2 = string_split(name2," ")[0];
        name2 = string_split(name2,"#")[0];
        int bitflag2 = StringToInt(e2[1]);
        std::string strand2 = get_strand(bitflag2);
        std::string sequence2 = e2[9];
        std::string chrom2 = e2[2];
        int start2 = StringToInt(e2[3]) -1;
        int stop2  = start2 + sequence2.length();
        
        // skip double stars
        if ((chrom1 == "*") & (chrom2 == "*"))
        {
          continue;
        }
        
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
        if ((chrom1 != chrom2) & (chrom1 == "*") )
        {
            reorder = true;
        }
        if ((chrom1 != chrom2) & (chrom2 == "*") )
        {
            reorder = false;
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


//// Define a function that removes duplicates from a bedpe file
//// [[Rcpp::export]]
//std::vector< int > removeDupBedpe(std::string infile,std::string outfile , bool renamePets = true)
//{
//    // initialize petnumer
//    int petnumber  = 0;
//    int alllines   = 0;
//    int duplines   = 0;
//    int interchromosomal = 0;
//    int intrachromosomal = 0;
//    
//    // arguments
//    ifstream file1(infile.c_str());
//    ofstream bedpefilestream (outfile.c_str());
// 
//    // read in file line by line store currentline and last line
//    std::string lastline;
//    std::string currline;
//    while (getline(file1, currline))
//    {
//        // increment counter 
//        alllines++;
//        
//        // split lines
//        std::vector<std::string> lastEall = string_split(lastline,"\t");
//        std::vector<std::string> currEall = string_split(currline,"\t");
//        std::vector<std::string> lastE;
//        std::vector<std::string> currE;
//        
//        // print first line
//        if (lastEall.size() == 0)
//        {
//            if (renamePets == true)
//            {
//              petnumber ++;
//              currEall[6] = "obs_" +  IntToString(petnumber);
//            }
//          
//            bedpefilestream << vector_join(currEall,"\t");;
//            bedpefilestream << "\n";
//            lastline = currline;
//            continue;
//        }
//        
//        // make shorter line for comparison
//        for (int i=0; i<6; i++)
//        {
//            lastE.push_back(lastEall[i]);
//            currE.push_back(currEall[i]);
//        }
//        std::string currlineshort = vector_join(currE,"_");
//        std::string lastlineshort = vector_join(lastE,"_");
//        
//        // if duplicates continue
//        if (lastlineshort == currlineshort)
//        {
//            lastline = currline;
//            duplines++;
//            continue;
//        }
//        
//        // print out non duplicates
//        if (renamePets == true)
//        {
//          petnumber++;
//          currEall[6] = "obs_" +  IntToString(petnumber);
//        }
//        
//        bedpefilestream << vector_join(currEall,"\t");;
//        bedpefilestream << "\n";
//        
//        if ((currEall[0] == currEall[4]) & (currEall[0]  != "*") & (currEall[4]  != "*")  )
//        {
//          interchromosomal++;
//        }
//        if (currEall[0] != currEall[4]  & (currEall[0]  != "*") & (currEall[4]  != "*") )
//        {
//          intrachromosomal++;
//        }
//
//        // update last line
//        lastline = currline;
//    }
//    
//    // close files
//    file1.close();
//    bedpefilestream.close();
//    
//    
//    // report results
//    std::vector< int > rmdupresults;
//    rmdupresults.push_back(duplines);
//    rmdupresults.push_back(petnumber);
//    rmdupresults.push_back(interchromosomal);
//    rmdupresults.push_back(intrachromosomal);
//    return(rmdupresults);
//}

// Define a class to keep track of peak information
class peak{
public:
  std::string name;
  std::string chrom;
  int start;
  int end;
  int intra;
  int peakdepth;
  std::set<std::string> PETs;
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
  int p1depth;
  int p2depth;
  int p1intra;
  int p2intra;
  int linking;
  int distance;
};



// Define a function that puts pairs together
// [[Rcpp::export]]
void findPairs(std::string overlapfile, std::string petpairsfile,std::string interactionfile,std::string peakscount,std::string peaksfileslopdepth)
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
        std::vector<std::string> readnamestuff = string_split(readname,".");
        std::string readnamenonumeber = readnamestuff[0];
      
        // add info to peak dict
        if (peakinfodict.find(peakname) == peakinfodict.end())
        {
            peak p = *(new peak());
            peakinfodict.insert(std::pair<string,peak>(peakname,p));
            peakinfodict[peakname].name  = currEall[9];
            peakinfodict[peakname].chrom = currEall[6];
            peakinfodict[peakname].start = atoi( currEall[7].c_str() );
            peakinfodict[peakname].end   = atoi( currEall[8].c_str() );
            peakinfodict[peakname].peakdepth = 0;
            //peakinfodict[peakname].intra = 0;
        }
        peakinfodict[peakname].PETs.insert(readnamenonumeber);
        //peakinfodict[peakname].intra++;
        
        // add info to readpeak dict
        if (readpeakdict.find(peakname) == readpeakdict.end())
        {
            std::vector<std::string> v = *(new std::vector<std::string>);
            readpeakdict.insert(std::pair<string,std::vector<std::string> > (readname, v));
        }
        readpeakdict[readname].push_back(peakname);
    }
  input.close();
  
  // Add peak depth info
  ifstream peakfile(peaksfileslopdepth.c_str());
  while (getline(peakfile, line))
  {
    // split lines
    std::vector<std::string> currEall = string_split(line,"\t");
    std::string peakname  = currEall[3];
    int peakdepth = atoi( currEall[4].c_str() );
    
    if (peakinfodict.find(peakname) != peakinfodict.end())
    {
      peakinfodict[peakname].peakdepth = peakdepth;
    }
  }
  peakfile.close();
  
//  // count number of PETs in each peak
//  for (std::map<std::string, peak>::iterator chippeak = peakinfodict.begin() ; chippeak != peakinfodict.end() ; ++chippeak )
//  {
//    // cout << IntToString(chippeak->second.intra) + " \t " +  IntToString(chippeak->second.PETs.size()) + "\n";
//    chippeak->second.intra = chippeak->second.PETs.size();
//    
//    // clear the hash to save memory
//    chippeak->second.PETs.clear();
//  }
  
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
          std::string pairname =  thep1.name + ":" + thep2.name;

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
              pairdict[pairname].p1depth = thep1.peakdepth;
              pairdict[pairname].p2depth = thep2.peakdepth;           
              //pairdict[pairname].p1intra  = thep1.intra; 
              //pairdict[pairname].p2intra  = thep2.intra; 
              pairdict[pairname].linking  = 0; 
              pairdict[pairname].distance  = thep2.start - thep1.end; 
            }
            pairdict[pairname].linking++; 
        }
      }  
  } 
  inputpets.close();
  
  // print out info to pair file
  ofstream pairsfilestream (interactionfile.c_str());
  for (std::map<std::string, chiapair>::iterator cp = pairdict.begin() ; cp != pairdict.end() ; ++cp ) {
      
      pairsfilestream << cp->second.p1chrom;
      pairsfilestream << "\t";
      pairsfilestream << cp->second.p1start;
      pairsfilestream << "\t";
      pairsfilestream << cp->second.p1end;
      pairsfilestream << "\t";
      pairsfilestream << cp->second.p2chrom;
      pairsfilestream << "\t";
      pairsfilestream << cp->second.p2start;
      pairsfilestream << "\t";
      pairsfilestream << cp->second.p2end;
      pairsfilestream << "\t";      
      pairsfilestream << cp->second.pairname;
      pairsfilestream << "\t";
      pairsfilestream << cp->second.p1name;
      pairsfilestream << "\t";
      pairsfilestream << cp->second.p2name;
      pairsfilestream << "\t";
      pairsfilestream << cp->second.p1depth;
      pairsfilestream << "\t";
      pairsfilestream << cp->second.p2depth;
      pairsfilestream << "\t";    
//      pairsfilestream << cp->second.p1intra;
//      pairsfilestream << "\t";
//      pairsfilestream << cp->second.p2intra;
//      pairsfilestream << "\t";
      pairsfilestream << cp->second.linking;
      pairsfilestream << "\t";
      pairsfilestream << cp->second.distance;
      pairsfilestream << "\n";
    } 
    
  // close output stream
  pairsfilestream.close();
  
  // print out info to peak file
  ofstream peaksfilestream (peakscount.c_str());
  for (std::map<std::string, peak>::iterator chippeak = peakinfodict.begin() ; chippeak != peakinfodict.end() ; ++chippeak )
  {
      peaksfilestream << chippeak->second.chrom;
      peaksfilestream << "\t";
      peaksfilestream << chippeak->second.start;
      peaksfilestream << "\t";
      peaksfilestream << chippeak->second.end;
      peaksfilestream << "\t";
      peaksfilestream << chippeak->second.name;
      peaksfilestream << "\t";
      peaksfilestream << chippeak->second.peakdepth;
      peaksfilestream << "\t.\n";
//      peaksfilestream << chippeak->second.intra;
//      peaksfilestream << "\t.\n";
  }
  peaksfilestream.close();

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
void makeDistanceFile(std::string bedpefilesortrmdup,std::string distancefile,int mindist, int maxdist)
{
    // streams
    ifstream filein  (bedpefilesortrmdup.c_str());
    ofstream fileout (distancefile.c_str());

    
    // read in file line by line and make same dif calls
    std::string line;
    while (getline(filein, line))
    {
        // split lines
        std::vector<std::string> currEall = string_split(line,"\t");
        
        // skip unmapped and inter chrom
        if ((currEall[0] != currEall[3]) || (currEall[0] == "*")  ||  (currEall[3] == "*"))
        {
          continue;
        }
        
        // determine distance
        std:string distance = IntToString((StringToInt(currEall[5]) + StringToInt(currEall[4]))
        / 2 - (StringToInt(currEall[2]) + StringToInt(currEall[1])) / 2);
        
        // determine orientation
        std::string pairtype = "D";
        if (currEall[8] == currEall[9])
        {
          pairtype = "S";
        }
        
        
        if (StringToInt(distance) > mindist & StringToInt(distance) < maxdist)
        {
          fileout << distance + "\t" + pairtype + "\n";
        }
    
        // make reads
        std::vector<std::string> read1vec;
        read1vec.push_back(currEall[0]);
        read1vec.push_back(currEall[1]);
    }

    // close files
    filein.close();
    fileout.close();
}



// Define a function that joins file (normally files previously split by chromosome)
// [[Rcpp::export]]
void  joinchromfiles(std::vector<std::string> sortedchromfiles,std::string bedpefilesort)
{
  // open output stream
  ofstream fileout(bedpefilesort.c_str());
  
  // read in each input stream
  for (int i=0; i< sortedchromfiles.size() ; i++)
  {
    // open input stream
    ifstream filein(sortedchromfiles[i].c_str());
    
    std::string line;
    // read in file line by line store currentline and last line
    while (getline(filein, line))
    {
      fileout << line;      
      fileout << "\n";
    }
    filein.close();
  }
  fileout.close();
}

// Define a function the collects info from a peak / tagAlign overlap
// [[Rcpp::export]]
void DeterminePeakDepthsC(std::string temppeakoverlap,std::string peaksfileslopdepth)
{
  // streams
  ifstream fileIN(temppeakoverlap.c_str());
  ofstream fileOUT(peaksfileslopdepth.c_str());
  
  // make map of peaks
  std::map<std::string, int > peaksmap;
  
  // read in file line by line
  std::string line;
  while (getline(fileIN, line))
  {
    // split lines and determine bin
    std::vector<std::string> currEall = string_split(line,"\t");
    std::string peakchrom = currEall[6];
    std::string peakstart = currEall[7];
    std::string peakend   = currEall[8];
    std::string peakname  = currEall[9];
  
    // peak name
    std::string longname = peakchrom + "\t" + peakstart + "\t" + peakend + "\t" + peakname;
  
    // add peak to map
    if ( peaksmap.find(longname) == peaksmap.end() ) {  
      peaksmap[longname] = 0;
    }
    
    // increment count
    peaksmap[longname]++;
    
  }
  
  // close input file
  fileIN.close();
  
  // print out info
  for (std::map<std::string, int >::const_iterator longname = peaksmap.begin() ; longname != peaksmap.end() ; longname ++ ){
    fileOUT << longname->first + "\t" + NumberToString(peaksmap[longname->first]) + "\t.";
    fileOUT << "\n";
  }
  
  // close output file
  fileOUT.close();

}





// Define a function removes duplicates from a bedpe file
// [[Rcpp::export]]
std::vector< std::string > removeDups(std::string bedpein,std::string outnamebase,double distancesplit)
{
  // (1) split PETs by chrom and position 
  
  // keep track of output files
  std::vector<std::string> outputvectorPETs;
  
  // streams
  ifstream file1(bedpein.c_str());
  std::map<std::string, std::ofstream*> petsoutput;
  
  std::string line;
  while (getline(file1, line))
  {
    // split lines and determine bin
    std::vector<std::string> currEall = string_split(line,"\t");
    std::string chrom = currEall[0];
    double pos = StringToInt(currEall[1]);
    int bin = pos / distancesplit;
    std::string binstring = NumberToString(bin);
    
    // set output file name
    std::string outname = outnamebase + "." + chrom + "." + binstring + ".bedpe";
    
    // check if output file name exists (and make it fi neccesary)
    if ( petsoutput.find(outname) == petsoutput.end() ) {  
      petsoutput[outname] = new std::ofstream(outname.c_str());
      outputvectorPETs.push_back( outname);  
    }
    
    // print to output file
    *petsoutput[outname] << line;
    *petsoutput[outname] << "\n";
  }
  
  // close input stream
  file1.close();
  
  // close bedpe files streams
  for (std::map<std::string, std::ofstream*>::iterator i = petsoutput.begin() ; i != petsoutput.end() ; i ++ ) {
    i->second->close();
  }    
  
  // initialize counters
  int nondups  = 0;
  int alllines   = 0;
  int duplines   = 0;
  int interchromosomal = 0;
  int intrachromosomal = 0;
  
  // (2) read through each file and only print out non duplicates
  // open input stream
  std::string outputname = outnamebase + ".rmdup.bedpe";
  ofstream finaloutput(outputname.c_str());
  
  for (std::vector<std::string>::const_iterator i = outputvectorPETs.begin() ; i != outputvectorPETs.end() ; i ++ ) {
    
    // make new map
    std::map<std::string, int> PETmap;
    
    // open input stream
    ifstream file1(i->c_str());
    
    std::string line;
    while (getline(file1, line))
    {
      alllines++;
      
      // split lines and determine bin
      std::vector<std::string> currEall = string_split(line,"\t");
      std::string chrom1  = currEall[0];
      std::string pos1  = currEall[1];
      std::string chrom2 = currEall[3];
      std::string pos2 = currEall[4];
      std::string uniqcode = pos1 + "_" + chrom2 + "_" + pos2;

      // check if read has been seen before
      if ( PETmap.find(uniqcode) == PETmap.end() ) {  
        PETmap[uniqcode] = 0;
      }
      PETmap[uniqcode]++;
      
      if (PETmap[uniqcode] > 1)
      {
        duplines++;
      }
      
      // if it is the first instance print it to the output file
      if (PETmap[uniqcode] == 1)
      {
        nondups++;
        if (chrom1 == chrom2)
        {
          intrachromosomal++;
        }
        if (chrom1 != chrom2)
        {
          interchromosomal++;
        }
        finaloutput << line;
        finaloutput << "\n";
      }
    }    
    // close input stream
    file1.close();
    
  }
  // close input stream
  finaloutput.close();
  
  // report results
  std::vector< std::string > rmdupresults;
  rmdupresults.push_back(NumberToString(duplines));
  rmdupresults.push_back(NumberToString(nondups));
  rmdupresults.push_back(NumberToString(interchromosomal));
  rmdupresults.push_back(NumberToString(intrachromosomal));
  rmdupresults.push_back(NumberToString(alllines));
  
  for (std::vector<std::string>::const_iterator i = outputvectorPETs.begin() ; i != outputvectorPETs.end() ; i ++ ) { 
    rmdupresults.push_back(*i);
    }
  
  return(rmdupresults);
}



// Define a function splits bedpe file into reads and PETs by chromosome
// [[Rcpp::export]]
std::vector<std::string> splitBedpe(std::string bedpein,std::string outnamebase, bool printreads = true , bool printpets = true, bool skipstars=true,bool skipinter=true)
{
    // keep track of output files
    std::vector<std::string> outputvectorPETs;
    std::vector<std::string> outputvectorReads;
    
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
        read1vec.push_back(currEall[6] + ".1");
        read1vec.push_back(currEall[7]);
        read1vec.push_back(currEall[8]);
        
        std::vector<std::string> read2vec;
        read2vec.push_back(currEall[3]);
        read2vec.push_back(currEall[4]);
        read2vec.push_back(currEall[5]);
        read2vec.push_back(currEall[6] + ".2");
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
              outputvectorReads.push_back( chrom1);  
          }
          
          if ( (readoutput.find(chrom2) == readoutput.end()) & (chrom2 != "*" ) ) {
              std::string outname = outnamebase + "." + chrom2 + ".bed";
              readoutput[chrom2] = new std::ofstream(outname.c_str());
              outputvectorReads.push_back( chrom2);  
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
        
        if (skipstars == true)
        {
          if ((chrom1 == "*") | (chrom2 == "*"))
          {
              continue;
          }
        }
        
        if (skipinter == true)
        {
          if (chrom1 != chrom2)
          {
              continue;
          }
        }
        
        if (printpets == true)
        {
          if ( petsoutput.find(chrom1) == petsoutput.end() ) {  
              std::string outname = outnamebase + "." + chrom1 + ".bedpe";
              petsoutput[chrom1] = new std::ofstream(outname.c_str());
              outputvectorPETs.push_back( chrom1);  
          }
          *petsoutput[chrom1] << line;
          *petsoutput[chrom1] << "\n";
        }
    }
    
    // close reads files
    for (std::map<std::string, std::ofstream*>::iterator i = readoutput.begin() ; i != readoutput.end() ; i ++ ) {
      i->second->close();
    }
    
    // close bedpe files
    for (std::map<std::string, std::ofstream*>::iterator i = petsoutput.begin() ; i != petsoutput.end() ; i ++ ) {
      i->second->close();
    }    
    
    // combine the outputs
    std::vector<std::string> ReadAndPETchroms;
    ReadAndPETchroms.reserve( outputvectorReads.size() + outputvectorPETs.size() ); // preallocate memory
    ReadAndPETchroms.insert( ReadAndPETchroms.end(), outputvectorReads.begin(), outputvectorReads.end() );
    ReadAndPETchroms.insert( ReadAndPETchroms.end(), outputvectorPETs.begin(), outputvectorPETs.end() );


    return ReadAndPETchroms;
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
    externalMergesort <string> externalMergeSorter(inputfile, outputfile, 5000000);
    return;

}


// Define a function the filters out every second line from a file
// [[Rcpp::export]]
void everyotherline(std::string overlapin, std::string overlapout) {
    
    // establish streams
    ifstream infile  (overlapin.c_str());
    ofstream outfile (overlapout.c_str() );

   int i = 0;
   std::string line;
   while (getline(infile, line))
    {
        // increment counters
        i++;
        
        if (i == 1)
        {
          outfile << line;
          outfile << "\n";
        }
        if (i == 2)
        {
          i = 0;
        }
    }
    infile.close();
    outfile.close();
}



// Define a function that adds Q values and filters results
// [[Rcpp::export]]
void AddQvals(std::string interactionfile, std::string interactionfilefinal,std::vector<double> Q,double maxPval )
{
    // establish streams
    ifstream infile  (interactionfile.c_str());
    ofstream outfile (interactionfilefinal.c_str() );
  
    int i = -2;
    std::string line;
    while (getline(infile, line))
    {
      
      // increment counters
      i++;
      
      if (i == -1)
      {
        outfile << line + "\t" + "adjP";
        outfile << "\n";
        continue;
      }
      
      
      double Qvalue = Q[i];
        
      //if (Qvalue >= maxPval)
      //{
        outfile << line + "\t" + NumberToString(Qvalue);
        outfile << "\n";
      //}
    }
    infile.close();
    outfile.close();
  
}




