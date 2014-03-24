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
using namespace std;

// Define a function that joins vectors of strings
// [[Rcpp::export]]
string vector_join( const vector<string>& v, const string& token ){
    ostringstream result;
    for (typename vector<string>::const_iterator i = v.begin(); i != v.end(); i++){
        if (i != v.begin()) result << token;
        result << *i;
    }
    return result.str();
}

// Define a function that splits strings into vector
// [[Rcpp::export]]
vector<string> string_split( const string& s, const string& delimiter ){
    vector<string> result;
    string::size_type from = 0;
    string::size_type to = 0;
    
    while ( to != string::npos ){
        to = s.find( delimiter, from );
        if ( from < s.size() && from != to ){
            result.push_back( s.substr( from, to - from ) );
        }
        from = to + delimiter.size();
    }
    return result;
}

// [[Rcpp::export]]
string parseFastq(string fastq1, string fastq2,string basename,
              int minlength = 15,int maxlength = 25,
              bool keepempty = false, bool verbose = true,
              string linker1 = "GTTGGATAAG" , string linker2 = "GTTGGAATGT")
{

   // arguments
    ifstream file1(fastq1);
    ifstream file2(fastq2);
    ofstream same1 ( basename + "_1.same.fastq");
    ofstream same2 ( basename + "_2.same.fastq");
    ofstream chim1 ( basename + "_1.chim.fastq");
    ofstream chim2 ( basename + "_2.chim.fastq");
    
    // define variables
    string fqline1;
    string fqline2;
    int linecount = 0;
    int i = 0;
    vector<string> lines1;
    vector<string> lines2;
    
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
            string pairtype = "unknown";
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
vector< int > get_bits( unsigned long x ) {
    std::string chars( std::bitset< sizeof(long) * CHAR_BIT >( x )
                      .to_string( char(0), char(1) ) );
    return std::vector< int >( chars.begin(), chars.end() );
}

// Define a function that returns the strand
string get_strand( vector< int > bitflag, int strandbit = 4 ) {
    string strand = "+";
    if (bitflag[63 - strandbit] == 0)
    {
        strand = "-";
    }
    return strand;
}

// Define a function that converts string to int
int StringToInt( string Text ) {
    int output;
    if ( ! (istringstream(Text) >> output) ) output = 0;
    return output;
}

// Define a function that converts string to int
string IntToString( int Number ) {
    string Result;          // string which will contain the result
    ostringstream convert;   // stream used for the conversion
    convert << Number;      // insert the textual representation of 'Number' in the characters in the stream
    Result = convert.str(); // set 'Result' to the contents of the stream
    return Result;
}


// Define a function that builds a bedpe file rom 2 sam file
// [[Rcpp::export]]
void buildBedpe(string sam1, string sam2,string bedpefile)
{
    
    // arguments
    ifstream file1(sam1);
    ifstream file2(sam2);
    ofstream bedpefilestream ( bedpefile);

    // define variables
    string line1;
    string line2;
    int linecount = 0;
    
    while (getline(file1, line1))
    {
        // read lines and increment counter
        getline(file2, line2);
        linecount++;
        
        // split lines
        vector<string> e1 = string_split(line1,"\t");
        vector<string> e2 = string_split(line2,"\t");
        
        // get info for file 1
        string name1 = e1[0];
        name1 = string_split(name1,"_")[0];
        name1 = string_split(name1," ")[0];
        int bitflag1 = StringToInt(e1[1]);
        vector<int> bits1 = get_bits(bitflag1);
        string strand1 = get_strand(bits1);
        string sequence1 = e1[9];
        string chrom1 = e1[2];
        int start1 = StringToInt(e1[3]) -1;
        int stop1  = start1 + sequence1.length();
        
        // get info for file 2
        string name2 = e2[0];
        name2 = string_split(name2,"_")[0];
        name2 = string_split(name2," ")[0];
        int bitflag2 = StringToInt(e2[1]);
        vector<int> bits2 = get_bits(bitflag2);
        string strand2 = get_strand(bits2);
        string sequence2 = e2[9];
        string chrom2 = e2[2];
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
            vector<string> outputvector;
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
            string outputstring = vector_join(outputvector,"\t");
            bedpefilestream << outputstring;
            bedpefilestream << "\n";
        }
        
        if (reorder == true)
        {
            vector<string> outputvector;
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
            string outputstring = vector_join(outputvector,"\t");
            bedpefilestream << outputstring;
            bedpefilestream << "\n";
        }
    }
}


// Define a function that removes duplicates from a bedpe file
// [[Rcpp::export]]
void removeDupBedpe(string infile,string outfile)
{
    // arguments
    ifstream file1(infile);
    ofstream bedpefilestream (outfile);
    
    // read in file line by line store currentline and last line
    string lastline;
    string currline;
    while (getline(file1, currline))
    {
        // split lines
        vector<string> lastEall = string_split(lastline,"\t");
        vector<string> currEall = string_split(currline,"\t");
        vector<string> lastE;
        vector<string> currE;
        
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
        string currlineshort = vector_join(currE,"_");
        string lastlineshort = vector_join(lastE,"_");
        
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

// Define a function splits bedpe file into reads and PETs by chromosome
// [[Rcpp::export]]
vector<string> splitBedpe(string bedpein,string outnamebase)
{
    vector<string> outputvector;
    
    // streams
    ifstream file1(bedpein);
    std::map<std::string, std::ofstream*> readoutput;
    std::map<std::string, std::ofstream*> petsoutput;
    
    string line;
    // read in file line by line store currentline and last line
    while (getline(file1, line))
    {
        // split lines
        vector<string> currEall = string_split(line,"\t");
        
        // make reads
        vector<string> read1vec;
        read1vec.push_back(currEall[0]);
        read1vec.push_back(currEall[1]);
        read1vec.push_back(currEall[2]);
        read1vec.push_back(currEall[6]);
        read1vec.push_back(currEall[7]);
        read1vec.push_back(currEall[8]);
        
        vector<string> read2vec;
        read2vec.push_back(currEall[3]);
        read2vec.push_back(currEall[4]);
        read2vec.push_back(currEall[5]);
        read2vec.push_back(currEall[6]);
        read2vec.push_back(currEall[7]);
        read2vec.push_back(currEall[9]);
        
        string read1 = vector_join(read1vec,"\t");
        string read2 = vector_join(read2vec,"\t");
        
        // print reads
        string chrom1 = currEall[0];
        if ( (readoutput.find(chrom1) == readoutput.end()) & (chrom1 != "*" )  ) {
            
            string outname = outnamebase + "." + chrom1 + ".bed";
            readoutput[chrom1] = new std::ofstream(outname);
        }
        string chrom2 = currEall[3];
        if ( (readoutput.find(chrom2) == readoutput.end()) & (chrom2 != "*" ) ) {
            
            string outname = outnamebase + "." + chrom2 + ".bed";
            readoutput[chrom2] = new std::ofstream(outname);
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
        
        if ((chrom1 == "*") | (chrom2 == "*"))
        {
            continue;
        }
        
        if ( petsoutput.find(chrom1) == petsoutput.end() ) {
            
            string outname = outnamebase + "." + chrom1 + ".bedpe";
            petsoutput[chrom1] = new std::ofstream(outname);
            outputvector.push_back(outnamebase + "." + chrom1);
            
        }
        *petsoutput[chrom1] << line;
        *petsoutput[chrom1] << "\n";
    }
    
    // Close all of the files
    for(auto& pair : readoutput) {
        delete pair.second;
        pair.second = 0;
    }
    // Close all of the files
    for(auto& pair : petsoutput) {
        delete pair.second;
        pair.second = 0;
    }
    return outputvector;
}





