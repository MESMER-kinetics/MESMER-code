#include "tinyxml.h"
#include <iostream>
#include <sstream>

#ifdef _WIN32
#include <conio.h>
#endif

#include <string>
#include "System.h"
using namespace std;

int main(int argc,char *argv[])
{
  //
  // Instantiate the System collection. This holds all information 
  // about reaction systems and all molecular data.
  //
  mesmer::System System ;
  string filename;
  if(argc>1)
    filename=argv[1];
  else
    filename="test.xml";
  TiXmlDocument doc( filename.c_str() );

  if ( !doc.LoadFile() )
  {
    cerr << "Could not load file " << filename << endl;
    exit( 1 );
  }
  TiXmlElement* root = doc.RootElement();
  if(root->ValueStr()!="me:mesmer")
  {
    cerr << filename << " is not a MESMER file" << endl;
    exit(1);
  }

  //
  // Parse input.
  // 
  if(!System.parse(root))
    exit(1);

  // 
  // Begin calculation.
  //
  System.calculate() ;

  #ifdef _DEBUG
  //CM keep window open
  cout << "Press any key to finish" <<endl;
  getch();
  #endif
  exit(0);

  return 0 ;
}

