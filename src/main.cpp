#include <iostream>
#include <sstream>

#ifdef _WIN32
#include <conio.h>
#endif

#include "System.h"

void usage();
bool hasNoAlpha(std::string& kick_it);

int main(int argc,char *argv[])
{
  //
  // Instantiate the System collection. This holds all information
  // about reaction systems and all molecular data.
  //
  mesmer::System _sys ;

  mesmer::TimeCount events; std::string thisEvent; unsigned int timeElapsed; 

  //-------------------------------
  // process command line arguments
  std::string inputfilename, outputfilename;
  if(argc<2) {usage(); return -1;}
  inputfilename = argv[1];
  if (argc > 2) outputfilename = argv[2];

  //This is where the type of IO is decided.
  //Opens the data file and checks that its root element is me:mesmer.
  mesmer::PersistPtr ppIOPtr = mesmer::XMLPersist::XmlCreate(inputfilename, "me:mesmer");

  //------------
  // Parse input
  thisEvent = "Parse input xml file";
  std::string parseTimeStamp = events.setTimeStamp(thisEvent);
  std::cout << thisEvent << " at " << parseTimeStamp << std::endl;
  if(!ppIOPtr || !_sys.parse(ppIOPtr))
    return -2;

  //------------------
  // Begin calculation
  std::clog << inputfilename << " successfully parsed. Now calculating..." << std::endl;

  thisEvent = "Calculate EGME";
  std::cout << thisEvent << " at " << events.setTimeStamp(thisEvent) << std::endl;
  _sys.calculate() ;

  //--------------------------------
  // Save XML document to a new file
  thisEvent = "Save XML document to a new file";
  std::string saveTimeStamp = events.setTimeStamp(thisEvent, timeElapsed);
  std::cout << thisEvent << " at " << saveTimeStamp  << " -- Time elapsed: " << timeElapsed << " seconds.\n";
  if(outputfilename.empty())
  {
    std::vector<std::string> subStrings;
    std::string::iterator anchor = inputfilename.begin();
    std::string thisString;
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
    std::string fromTimeStamp;
    for (unsigned int i = 1; i < subStrings.size(); ++i){
      if (i == (subStrings.size() - 1)){//check the file extension
        if (!subStrings[i].compare("xml") && !subStrings[i].compare("XML")){
          std::stringstream errorMsg;
          errorMsg << "The file extension of the filename has to be XML\nCannot recognise: " << inputfilename;
          mesmer::obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), mesmer::obError);
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
    std::stringstream errorMsg;
    errorMsg << "There was an error when writing " << outputfilename;
    mesmer::obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), mesmer::obError);
  }

/*  #ifdef _DEBUG
  //CM keep window open
  std::cout << "Press any key to finish" <<std::endl;
  getch();
  #endif
*/

  return 0 ;
}

bool hasNoAlpha(std::string& kick_it){
  std::string::iterator anchor = kick_it.begin();
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
  std::stringstream errorMsg;
  errorMsg  << "#----- mesmer inputfilename [outputfilename] -----#\n";
  errorMsg  << "The default outputfilename is the input name\n";
  errorMsg  << "with 'out' added before the extension.\n";
  errorMsg  << "  Any existing file will be overwritten.";
  mesmer::obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), mesmer::obInfo);
}

