#include <iostream>
#include <sstream>

#ifdef _WIN32
#include <conio.h>
#endif

#include "System.h"

using namespace std ;
using namespace Constants ;
using namespace mesmer ;

void usage();
bool hasNoAlpha(string& kick_it);

int main(int argc,char *argv[])
{
  //
  // Instantiate the System collection. This holds all information
  // about reaction systems and all molecular data.
  //
  System _sys ;

  TimeCount events; unsigned int timeElapsed;

  //-------------------------------
  // process command line arguments
  string inputfilename, outputfilename;
  if(argc<2) {usage(); return -1;}
  inputfilename = argv[1];
  if (argc > 2) outputfilename = argv[2];

  //This is where the type of IO is decided.
  //Opens the data file and checks that its root element is me:mesmer.
  PersistPtr ppIOPtr = XMLPersist::XmlLoad(inputfilename, "me:mesmer");

  //------------
  // Parse input
    {
      stringstream errorMsg;
      string thisEvent = "Parse input xml file";
      string parseTimeStamp = events.setTimeStamp(thisEvent);
      errorMsg << thisEvent << " at " << parseTimeStamp << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
    }

  if(!ppIOPtr || !_sys.parse(ppIOPtr))
    return -2;

  //------------------
  // Begin calculation
  {
    stringstream errorMsg;
    string thisEvent = "Calculate EGME";
    errorMsg << inputfilename << " successfully parsed. " << thisEvent << " at " << events.setTimeStamp(thisEvent) << endl;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
  }

  _sys.calculate() ;

  //--------------------------------
  // Save XML document to a new file
  string saveTimeStamp(0);
  {
    stringstream errorMsg;
    string thisEvent = "Save XML document to a new file";
    saveTimeStamp = events.setTimeStamp(thisEvent, timeElapsed);
    errorMsg << thisEvent << " at " << saveTimeStamp  << " -- Time elapsed: " << timeElapsed << " seconds.\n";
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
  }

  if(outputfilename.empty())
  {
    vector<string> subStrings;
    string::iterator anchor = inputfilename.begin();
    string thisString;
    //populate the string vector
    while(1){
      if (anchor == inputfilename.end()){
        subStrings.push_back(thisString); // populating the vector 'subStrings' with the current string
        break; //break if reaches the end of inputfilename
      }
      else{
        if (*anchor == '.'){
          ++anchor; //go forward without populating the string 'thisString'
          subStrings.push_back(thisString); // populating the vector 'subStrings' with the current string
          thisString.clear(); //clearing the current string
        }
        else{
          thisString += *anchor;
          ++anchor;
        }
      }
    }

    //analyse the string vector
    outputfilename = subStrings[0] + '.'; // no need to check the first sub-string
    string fromTimeStamp;
    for (unsigned int i = 1; i < subStrings.size(); ++i){
      if (i == (subStrings.size() - 1)){//check the file extension
        if (!subStrings[i].compare("xml") && !subStrings[i].compare("XML")){
          stringstream errorMsg;
          errorMsg << "The file extension of the filename has to be XML\nCannot recognise: " << inputfilename;
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        }
        else{
          if (fromTimeStamp.size()){
            outputfilename += fromTimeStamp;
            outputfilename += '.';
          }
          outputfilename += saveTimeStamp;
          outputfilename += ".xml";
        }
      }
      else{ //not file extension (one of the 2nd, 3rd, 4th sub-strings), check if it is date format (no alphabets)
        if (hasNoAlpha(subStrings[i])){//assume anything has no alphabets as a date number
          fromTimeStamp = subStrings[i];
        }
        else{ //has alphabets
          outputfilename += subStrings[i];
          outputfilename += '.';
        }
      }
    }
  }
  if(!ppIOPtr->XmlSaveFile(outputfilename))
  {
    stringstream errorMsg;
    errorMsg << "There was an error when writing " << outputfilename;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
  }

/*  #ifdef _DEBUG
  //CM keep window open
  cout << "Press any key to finish" <<endl;
  getch();
  #endif
*/

  return 0 ;
}

bool hasNoAlpha(string& kick_it){
  string::iterator anchor = kick_it.begin();
  while(1){
    if (anchor == kick_it.end()){
      return true; //break if reaches the end of kick_it
    }
    else{
      if (isalpha(*anchor)) return false;
      ++anchor;
    }
  }
}

void usage()
{
  stringstream errorMsg;
  errorMsg  << "#----- mesmer inputfilename [outputfilename] -----#\n";
  errorMsg  << "The default outputfilename is the input name\n";
  errorMsg  << "with 'out' added before the extension.\n";
  errorMsg  << "  Any existing file will be overwritten.";
  obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
}

