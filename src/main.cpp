#include <iostream>
#include <sstream>

#ifdef _WIN32
#include <conio.h>
#endif

#include "Persistence.h"
#include "System.h"

#include "XMLPersist.h"
using namespace std;

void usage();
int main(int argc,char *argv[])
{
  //
  // Instantiate the System collection. This holds all information 
  // about reaction systems and all molecular data.
  // It can be accessed via a namespace variable mesmer::pSys
  //
  mesmer::System System ;
  mesmer::pSys = &System;

  string inputfilename, outputfilename;
  if(argc<2)
  {
    usage();
    return -1;
  }

  inputfilename=argv[1];
  if(argc>2)
    outputfilename = argv[2];

  //This is where the type of IO is decided.
  //Opens the data file and checks that its root element is me:mesmer.
  mesmer::PersistPtr ppIOPtr = mesmer::XMLPersist::Create(inputfilename, "me:mesmer");

  //
  // Parse input.
  // 
  if(!ppIOPtr || !System.parse(ppIOPtr))
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
  if(!ppIOPtr->SaveFile(outputfilename))
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

