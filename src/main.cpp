#include "tinyxml.h"
#include <iostream>
#include <sstream>

#ifdef _WIN32
#include <conio.h>
#endif

#include <string>
#include "System.h"
using namespace std;

void usage();
int main(int argc,char *argv[])
{
  //
  // Instantiate the System collection. This holds all information 
  // about reaction systems and all molecular data.
  //
  mesmer::System System ;
  string inputfilename, outputfilename;
  if(argc<1)
  {
    usage();
    return -1;
  }

  inputfilename=argv[1];
  if(argc>2)
    outputfilename = argv[2];

  TiXmlDocument doc( inputfilename.c_str() );
  if ( !doc.LoadFile() )
  {
    cerr << "Could not load file " << inputfilename << endl;
    return -2;
  }
  TiXmlElement* root = doc.RootElement();
  if(root->ValueStr()!="me:mesmer")
  {
    cerr << inputfilename << " is not a MESMER file" << endl;
    return -2;
  }

  //
  // Parse input.
  // 
  if(!System.parse(root))
    return -2;
  clog << inputfilename << " successfully parsed. Now calculating..." <<endl;
  // 
  // Begin calculation.
  //
  System.calculate() ;

  // Save XML doc to a new file   
  if(outputfilename.empty())
  {
    string::size_type pos = inputfilename.rfind('.');
    if(pos!=string::npos)
      outputfilename = inputfilename.replace(pos, 1, "out.");
    else
      outputfilename = inputfilename + "out";
  }
  if(!doc.SaveFile(outputfilename))
  {
    cerr << "There was an error when writing " << outputfilename << endl;
    return -3;
  }

/*  #ifdef _DEBUG
  //CM keep window open
  cout << "Press any key to finish" <<endl;
  getch();
  #endif
*/

  return 0 ;
}

void usage()
{
cerr << "mesmer inputfilename [outputfilename]\n"
     << "The default outputfilename is the input name\
 with 'out' added before the extension.\n\
Any existing file will be overwritten." << endl;
}

