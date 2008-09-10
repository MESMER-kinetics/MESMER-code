#include <iostream>
#include <fstream>
#include <sstream>
#include "System.h"

using namespace std ;
using namespace Constants ;
using namespace mesmer ;

void usage();
string version();
bool QACompare(string infilename);
string duplicateFileName(const string& inName, const string& suffix, const string& newTimeStamp = "");

int main(int argc,char *argv[])
{
  if(argc<2)
  {
    usage();
    return 0;
  }

  // process command line arguments
  string infilename, outfilename, testfilename, logfilename;
  bool nocalc=false, notimestamp=false, usecout=false, updatemols=true, qatest=false, nologging=false, changetestname=false;

  for(int iarg=1; iarg<argc;++iarg)
  {
    const char* p = argv[iarg];
    if(*p=='-')
    {
      switch(*++p)
      {
      case '?':
        usage();
        return 0;
      case 'g':
        nologging=true;
        cerr << "Logging is turned off: no info or test output" << endl;
        break;
      case 'l':
        updatemols=false;
        break;
      case 'n':
        notimestamp=true;
        break;
      case 'N':
        changetestname=true;
        break;
      case 'o': //output filename
        if(!*++p) //if no digit after -o use next arg if its not an option
        {
          if(++iarg<argc && *argv[iarg]!='-')
            p = argv[iarg];
          else //-o option has no filename
          {
            --iarg;
            usecout = true;
            break;
          }
        }
        outfilename = p;
        break;
      case 'p': //just parse the input file - no calculation
        nocalc=true;
        break;
      case 'q':
        qatest=true;
        break;
      case 'V':
        cout << version();
        return 0;
      case 'w': //error level for reporting
        if(!*++p) //if no digit after -w use next arg
          p=argv[++iarg];
        meErrorLog.SetOutputLevel((obMessageLevel)(*p-'0'));
        break;
      }
    }
    else{
      infilename = p;
      if (changetestname){
        string testSuffix(".test");
        string logSuffix(".log");
        testfilename = duplicateFileName(infilename, testSuffix);
        logfilename = duplicateFileName(infilename, logSuffix);
      }
      else{ // default test and log names
        testfilename = "mesmer.test";
        logfilename = "mesmer.log";
      }
    }
  }

  //-----------------------------------
  //Start the error handling system
  //
  ostream* plog=NULL;
  ofstream logstream(logfilename.c_str());
  if(!logstream)
    cerr << "Could not open " << logfilename << " for writing. Continuing without it." << endl;
  else
    plog = &logstream;

  //Write all messages to mesmer.log, if it was opened.
  meErrorLog.SetLogStream(plog);

  ofstream osout(testfilename.c_str());
  if(!osout)
    cerr << "Could not open " << testfilename << " for writing. Continuing without it." << endl;

  //Allow writing of messages to cerr, cwarn and cinfo. Send ctest to file mesmer.test
  //Original streams restored when this goes out of context
  OStreamRedirector osr(&meErrorLog, &osout, nologging);

  //-------------------------------

  //
  // Instantiate the System collection. This holds all information
  // about reaction systems and all molecular data.
  //
  System _sys ;

  TimeCount events; unsigned int timeElapsed;

  //This is where the type of IO is decided.
  //Opens the data file and checks that its root element is me:mesmer.
  PersistPtr ppIOPtr = XMLPersist::XmlLoad(infilename, "me:mesmer");
  if(!ppIOPtr)
    return -1;
  //------------
  // Parse input file

  {
    string thisEvent;
    if(infilename.empty())
      thisEvent = "Parsing xml from stdin...";
    else
      thisEvent = "Parsing input xml file... ";
    cerr << thisEvent; //usually output
    cinfo << infilename << endl;
    events.setTimeStamp(thisEvent);
  }

  if(!ppIOPtr || !_sys.parse(ppIOPtr))
  {
    cerr << "System parse failed." << endl;
    return -2;
  }

  /* Get a list of molecules which were not present in the data file
  but which were sucessfully recovered from the Library.
  Unless the -l option is used, insert them into the data tree and
  save to file at the end of main(). */
  vector<PersistPtr> libMols = _sys.getLibraryMols();
  if(updatemols && !libMols.empty())
  {
    for(size_t i=libMols.size();i;--i) //backwards so mols come out in order they went in
      ppIOPtr->XmlCopyElement(libMols[i-1]);
    //ppIOPtr->XmlSaveFile("E:/My Documents/MSVC/MESMER/examples/DiagramTest/testout.xml"); //***TEMPORARY****
  }

  meErrorLog.SetContext("");
  //------------------
  // Begin calculation
  {
    string thisEvent = "Calculate EGME";
    cinfo << "\nFile: \"" << infilename << "\" successfully parsed.\n" << thisEvent << endl;
    events.setTimeStamp(thisEvent);
  }

  if (nocalc) return 0;

  clog << "\nNow calculating..." << infilename << endl;

  switch (_sys.m_Flags.searchMethod){
    case 2:
      _sys.fitting() ;
      break;
    case 1:
      _sys.gridSearch();
      break;
    default:
      double chiSquare(1000.0);
      _sys.calculate(chiSquare) ;
  }

  meErrorLog.SetContext(__FUNCTION__);
  //--------------------------------
  // Save XML document to a new file
  string thisEvent = "Save XML document to a new file";
  string currentTimeStamp = events.setTimeStamp(thisEvent, timeElapsed);
  string saveTimeStamp = '.' + currentTimeStamp;
  cinfo << " -- Total time elapsed: " << timeElapsed << " seconds.\n" << endl;

  if(!usecout && outfilename.empty() && !infilename.empty())
  {
    string XMLsuffix(".xml");
    outfilename = duplicateFileName(infilename, XMLsuffix, saveTimeStamp);
  }

  if(!ppIOPtr->XmlSaveFile(outfilename))
    cerr << "There was an error when writing " << outfilename;
  else
  {
    meErrorLog.SetContext(""); //so no ***Error prefix
    if(!outfilename.empty())
      cerr << "System saved";
  }
  if(qatest)
  {
    osout.close();
    if(!QACompare(infilename))
    {
      cerr << "QA test failed" << endl;
      return -5;
    }
    else
      cerr << "QA test successful" << endl;
  }
  return 0 ;
}

string duplicateFileName(const string& inName, const string& suffix, const string& newTimeStamp){
  //construct duplicatedName from inName by adding or replacing a timestamp
  string duplicatedName = inName;
  string oldTimeStamp;
  for(;;) // remove extensions and timestamps
  {
    string::size_type pos = duplicatedName.rfind('.'); // search the string starting from the end
    //Break if no more dots or if part of current or parent directory aliases
    if(pos==string::npos || duplicatedName[pos+1]=='/' || duplicatedName[pos+1]=='\\')
      break;
    //Save last timestamp (that of the input file)
    if(oldTimeStamp.empty() && !duplicatedName.compare(pos, 3, ".20"))
      oldTimeStamp = duplicatedName.substr(pos, 15);
    duplicatedName.erase(pos);
  }
  if(!newTimeStamp.empty())
    duplicatedName += oldTimeStamp + newTimeStamp;
  duplicatedName += suffix;
  return duplicatedName;
}

void usage()
{
  cout << "Usage:\n"
    "In Windows command line:\n"
    "     mesmer infilename [options]\n"
    "In UNIX/Linux console:\n"
    "     ./mesmer.out infilename [options]\n"
    "If there is no infilename, stdin is used.\n"
    "Options:\n"
    " -o<outfilename>\n"
    "     If -o has no outfilename, stdout is used.\n"
    "     If there is no -o option, the output file name is constructed\n"
    "     from the input file name by adding or replacing a timstamp.\n"
    " -n  No timestamp in output file name.\n"
    " -N  Use the input filename for test and log, default is mesmer.test and mesmer.log.\n"
    " -p  Parse the input file only - no calculation.\n"
    " -w# Display only warning messages that are at least as severe as:\n"
    "       0 No messages\n"
    "       1 Errors\n"
    "       2 Warnings\n"
    "       3 Information\n"
    "       4 Audit messages\n"
    " -?  Display this help text\n"
    " -V  Output Mesmer version.\n\n"
    "For example:\n"
    "  mesmer HSO2.xml     will read HSO2.xml and write the output to\n"
    "                      something like HSO2.20071031_134701.xml\n"
    "                      Using this in turn as an input file will produce\n"
    "                      something like HSO2.20071031_134701.20071212_200423.xml\n"
    "  mesmer HSO2.xml -n  will read HSO2.xml and write the output back to it\n"
    "  mesmer -o -w0       will silently read from cin and write to cout\n"
    << endl;
  // "-q  Do a QA test.
}
string version()
{
  stringstream ss;
  ss << "Mesmer " << "v0.1" << " compiled: -- "  << __DATE__ << " -- " << __TIME__ << endl;
  return ss.str();
}

bool QACompare(string infilename)
{
  //Read the baseline file
  string::size_type pos = infilename.find_last_of("/\\:");
  infilename.erase(pos+1); //may erase everything if has no path
  infilename.append(TestSubFolder);
  ifstream QAfile(infilename.c_str());
  ifstream CurrentTest("mesmer.test");
  if(!QAfile || !CurrentTest)
  {  
    cerr << "Cannot open " << infilename << " or current mesmer.test" << endl;
    return false;
  }
  typedef std::istreambuf_iterator<char> istreambuf_iterator;
  return std::equal(
    istreambuf_iterator(QAfile),
    istreambuf_iterator(),
    istreambuf_iterator(CurrentTest)); 
  // return true even if there are extra characters at the end of outstream
  //but that doesn't matter here.
}

/*
Mesmer outputs:
                        Source                           Destination  
--------------------------------------------------------------------------------------
Main output       Functions in IPersist                Decided by -o option  
                  like XmlWriteElement                 Usually file,or cout for piping

Error messages    cerr <<...    or                     Console unless -w0
                  meErrorLog.ThrowError( , ,obError)   and mesmer.log

Warning messages  cwarn <<...   or                     Console unless -w0 or w1
                  meErrorLog.ThrowError( , ,obWarn)    and mesmer.log

Logging messages  cinfo <<...   or                     Console if -w2
                  meErrorLog.ThrowError( , ,obInfo)    and mesmer.log

QA test output    ctest <<...                          mesmer.test

Temporary debug   cout <<...                           File when redirected
                  Use for bulky                        from commandline with >
                  debug output                         (Otherwise console)

Temporary debug   clog <<...                           Console                    
                  Use for short debug output


Output to mesmer.log can be turned off with the -g option.
The files will be overwrtten each run.
A comparison with a previously stored version of mesmer.test can be made with the -q option.

*/
