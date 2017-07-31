/*  Copyright (C) 2009-2018 by
Struan H. Robertson, Chris Morley, Chi-Hsiu Liang,
David R. Glowacki, Michael J. Pilling and contributors
This file is a part of
Mesmer: Master Equation Solver for Multi-Energy well Reactions

Mesmer is free software: you can redistribute it and/or modify
it under the terms of the GNU Public License as published by
the Free Software Foundation, either version 2 of the License,
or (at your option) any later version.

Mesmer is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Public License for more details.

You should have received a copy of the GNU Public License
along with Mesmer.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "System.h"
#include "ParallelManager.h"
#include "TimeCounter.h"

using namespace std;
using namespace Constants;
using namespace mesmer;

void usage();
string version();
void banner(size_t nRanks);
bool QACompare(string infilename, bool NOptionUsed);
string duplicateFileName(const string& inName, const string& suffix, const string& newTimeStamp = "");
string replaceFilename(const string& inName, const string& newFilename);

int main(int argc, char **argv)
{
  //
  // MPI initiation.
  //

  ParallelManager parallelManager(argc, argv);

  //
  // The following invocation is required by the QD  library to fix a problem 
  // with the extended precision methods on x86 platforms.
  //
  unsigned int old_cw;
  fpu_fix_start(&old_cw);

  // Visual C++ for Windows has used 3 digit exponent when outputting
  // scientific format numbers even when it was not necessary. From
  // Visual Studio 2015 it uses the more standard 2-digit exponent
  // and the following code sets this for earlier Visual C++ compilers.
  // The QA baselines for Windows in Mesmer 4 and earlier are
  // consequently different from those for Mesmer 5 after Aug 2015.
  // To get original behaviour with older compilers define
  // USE_OLD_EXPONENT before compiling. 
#if !USE_OLD_EXPONENT && __MSC_VER && _MSC_VER<1900
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  if (argc < 2)
  {
    usage();
    parallelManager.finalize();
    return 0;
  }

  // process command line arguments
  string infilename, outfilename, testfilename, logfilename, punchfilename;
  vector<string> extraInfilenames;
  bool nocalc = false, usecout = false, updatemols = true, overwriteinput = false,
    qatest = false, nologging = false, changetestname = false;

  for (int iarg = 1; iarg < argc; ++iarg)
  {
    const char* p = argv[iarg];
    if (*p == '-')
    {
      switch (*++p)
      {
      case '?':
        usage();
        parallelManager.finalize();
        return 0;
      case 'g':
        nologging = true;
        cerr << "Logging is turned off: no info or test output" << endl;
        break;
      case 'l':
        updatemols = false;
        break;
      case 'n':
        overwriteinput = true;
        break;
      case 'N':
        changetestname = true;
        break;
      case 'o': //output filename
        if (!*++p) //if no digit after -o use next arg if its not an option
        {
          if (++iarg < argc && *argv[iarg] != '-')
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
        nocalc = true;
        break;
      case 'q':
        qatest = true;
        break;
      case 'T': //Table of plugins with descriptions
        cout << TopPlugin::List(TopPlugin::verbose) << endl;
        parallelManager.finalize();
        return 0;
      case 't': //Table of plugins
        cout << TopPlugin::List() << endl;
        parallelManager.finalize();
        return 0;
      case 'V':
        cout << version();
        parallelManager.finalize();
        return 0;
      case 'w': //error level for reporting
        if (!*++p) //if no digit after -w use next arg
          p = argv[++iarg];
        meErrorLog.SetOutputLevel((obMessageLevel)(*p - '0'));
        break;
      default:
        cerr << "The option -" << *p << " was not recognized" << endl;
        if (--argc < 2)
        {
          usage();
          parallelManager.finalize();
          return 0;
        }
      }
    }
    else {
      if (infilename.empty())
        infilename = p;
      else
        extraInfilenames.push_back(p);
    }
  }

  if (changetestname) {
    string testSuffix(".test");
    string logSuffix(".log");
    string punchSuffix(".punch");
    testfilename = duplicateFileName(infilename, testSuffix);
    logfilename = duplicateFileName(infilename, logSuffix);
    punchfilename = duplicateFileName(infilename, punchSuffix);
  }
  else { // default test and log names
    testfilename = "mesmer.test";
    logfilename = "mesmer.log";
    punchfilename = "mesmer.punch";
  }

  if (!usecout && outfilename.empty()) {
    string repname = "mesmer_out.xml";
    outfilename = overwriteinput ? infilename : replaceFilename(infilename, repname);
  }

  //-----------------------------------
  //Start the error handling system
  //
  ostream* plog = NULL;
  ofstream logstream(logfilename.c_str());
  if (!logstream)
    cerr << "Could not open " << logfilename << " for writing. Continuing without it." << endl;
  else
    plog = &logstream;

  //Write all messages to mesmer.log, if it was opened.
  meErrorLog.SetLogStream(plog);

  ofstream osout(testfilename.c_str());
  if (!osout)
    cerr << "Could not open " << testfilename << " for writing. Continuing without it." << endl;

  // Allow writing of messages to cerr, cwarn and cinfo. Send ctest to file mesmer.test
  // Original streams restored when this goes out of context.
  // SHR, 17/Jan/2017: The nologging feature is used here to stop slaves from writing
  // to output during parallel execution.
  OStreamRedirector osr(&meErrorLog, &osout, (nologging || parallelManager.rank()));

  banner(parallelManager.size());
  //if(parallelManager.size() >1)
  //  cerr << "Using " << parallelManager.size() << " processes" << endl;

  //-----------------------------------------------
  //Get the top level directory from an environment variable,
  //or if that fails, from two levels above the current directory
  //This contains librarymols.xml, defaults.xml.
  const char* pdir = getenv("MESMER_DIR");
  if (!pdir || *pdir == '\0') //env var absent or empty
    pdir = "../..";
  string MesmerDir(pdir);
  string::size_type pos = MesmerDir.find(';');
  if (pos != string::npos)
    MesmerDir.erase(pos); //Use the first directory in the env var

  cinfo << "defaults.xml, librarymols.xml are in " << MesmerDir << endl;
  //-------------------------------
  //
  // Instantiate the System collection. This holds all information
  // about reaction systems and all molecular data.
  //
  System _sys(MesmerDir + "/librarymols.xml", &parallelManager);
  if (!_sys.initialize())
    cerr << "Failed to create System object" << endl;

  _sys.m_Flags.punchFileName = punchfilename;

  TimeCount events; unsigned int timeElapsed;

  //This is where the type of IO is decided.
  //Opens the data file and checks that its root element is me:mesmer.
  PersistPtr ppIOPtr = XMLPersist::XmlLoad(infilename, MesmerDir + "/defaults.xml", "me:mesmer");
  if (!ppIOPtr) {
    parallelManager.finalize();
    return -1;
  }

  //Incorporate the additional input files into the main datafile
  vector<string>::reverse_iterator fileitr;
  for (fileitr = extraInfilenames.rbegin(); fileitr != extraInfilenames.rend(); ++fileitr)
  {
    if (!ppIOPtr->XmlInclude(*fileitr)) {
      parallelManager.finalize();
      return -1;
    }
  }

  //------------
  if (nocalc)
    _sys.assignMolTypes(ppIOPtr);

  // Parse input file
  {
    string thisEvent;
    if (infilename.empty())
      thisEvent = "\nParsing xml from stdin...";
    else
      thisEvent = "\nParsing input xml file...\n" + infilename + '\n';
    cinfo << thisEvent << endl; //usually output
    events.setTimeStamp(thisEvent);
  }
  int returnCode = 0;
  try {
    if (!ppIOPtr || !_sys.parse(ppIOPtr))
    {
      if (plog)
        plog->flush();
      cerr << "\nSystem parse failed." << endl;
      returnCode = -2;
      if (!nocalc) {
        parallelManager.finalize();
        return returnCode;
      }
    }
    _sys.WriteMetadata(infilename);
    if (plog)
      plog->flush();

    if (!nocalc)
    {
      //------------------
      // Begin calculation
      {
        string thisEvent = "Main Calculation begins\n";
        cinfo << "\nFile: \"" << infilename << "\" successfully parsed.\n" << thisEvent;
        events.setTimeStamp(thisEvent);
      }

      cinfo << "Now calculating..." << endl;

      _sys.executeCalculation();
    }
  }
  catch (std::runtime_error& e)
  {
    cinfo.flush();
    clog.flush();
    cerr << e.what() << endl;
    cerr.flush();
    parallelManager.finalize();
    exit(-1);
  }
  catch (std::logic_error&) {} // Outputs XML before terminating (for debugging)

  //--------------------------------
  // Save XML document
  string thisEvent = "Save XML document";
  string currentTimeStamp = events.setTimeStamp(thisEvent, timeElapsed);
  string saveTimeStamp = '.' + currentTimeStamp;
  cerr << " -- Total time elapsed: " << timeElapsed << " seconds.\n" << endl;

  //Any existing file with the same name as the one being written is renamed with a _prev suffix
  //Any old _prev file is not deleted unless the write of the new file is successful
  if (!usecout)
    rename(outfilename.c_str(), "temp");

  if (parallelManager.Master()) {
    if (!ppIOPtr->XmlSaveFile(outfilename))
    {
      cerr << "There was an error when writing " << outfilename;
      rename("temp", outfilename.c_str());
    }
    else
    {
      if (!outfilename.empty())
      {
        stringstream line;
        line << "System saved to " << outfilename << endl;
        line << "Total time elapsed: " << timeElapsed << " seconds." << endl;
        cinfo << line.str();

        if (!usecout)
        {
          string::size_type pos = outfilename.rfind('.');
          string prevname(outfilename);
          prevname = pos == string::npos ?
            prevname + "_prev" :
            prevname.replace(pos, 1, "_prev.");
          remove(prevname.c_str());
          rename("temp", prevname.c_str());
        }
      }
    }
  }

  if (qatest)
  {
    osout.close();
    if (!QACompare(infilename, changetestname))
    {
      cerr << "QA test *** FAILED ***" << endl;
      returnCode = -5;
    }
    else
      cerr << "QA test successful" << endl;
  }

  //
  // The following invocation is required by the QD  library to restore 
  // default floating point behaviour on x86 platforms.
  //
  fpu_fix_end(&old_cw);

  parallelManager.finalize();

  return returnCode;
}
/*
main return codes
0  ok
-1 runtime exception
-2 parse failed
-3 try to overwrite input file without -n option
-5 QA test failed
*/

string duplicateFileName(const string& inName, const string& suffix, const string& newTimeStamp) {
  //construct duplicatedName from inName by adding or replacing a timestamp
  string duplicatedName = inName;
  string oldTimeStamp;
  for (;;) // remove extensions and timestamps
  {
    string::size_type pos = duplicatedName.rfind('.'); // search the string starting from the end
    //Break if no more dots or if part of current or parent directory aliases
    if (pos == string::npos || duplicatedName[pos + 1] == '/' || duplicatedName[pos + 1] == '\\')
      break;
    //Save last timestamp (that of the input file)
    if (oldTimeStamp.empty() && !duplicatedName.compare(pos, 3, ".20"))
      oldTimeStamp = duplicatedName.substr(pos, 15);
    duplicatedName.erase(pos);
  }
  if (!newTimeStamp.empty())
    duplicatedName += oldTimeStamp + newTimeStamp;
  duplicatedName += suffix;
  return duplicatedName;
}

string replaceFilename(const string& inName, const string& newFilename) {
  string replacedName = inName;
  string::size_type posf = replacedName.rfind('/');
  string::size_type posr = replacedName.rfind('\\');
  if (posf != string::npos) {      // found '/'
    replacedName.erase(posf + 1);
    replacedName += newFilename;
  }
  else if (posr != string::npos) { // found '\\' 
    replacedName.erase(posr + 1);
    replacedName += newFilename;
  }
  else {                           // no '/' or '\\' found
    replacedName = newFilename;   // meaning the filename has no directory structure --> simply a filename
  }

  return replacedName;
}

void usage()
{
  cout << "Usage:\n"
    "In Windows command line:\n"
    "     mesmer infilename [options]\n"
    "In UNIX/Linux console:\n"
    "     ./mesmer.out infilename [options]\n"
    "If there is no infilename, stdin is used.\n"
    "All the chemistry and much of the program control is in the input file.\n\n"
    "Options:\n"
    " -o<outfilename>\n"
    "     If -o has no outfilename, stdout is used.\n"
    " -n  Allows the input file to be overwritten when there is no -o option,\n"
    "       otherwise the output file is \"mesmer_out.xml\"\n"
    "     Any file about to be overwritten is saved with suffix _prev to its name. \n"
    " -N  Use infilename for test and log, default is mesmer.test and mesmer.log.\n"
    " -p  Parse the input file only - no calculation.\n"
    " -q  Do a QA test. (Compares mesmer.test with a definitive version.)\n"
    " -w# Display only warning messages that are at least as severe as:\n"
    "       0 No messages\n"
    "       1 Errors\n"
    "       2 Warnings\n"
    "       3 Information\n"
    " -?  Display this help text\n"
    " -V  Output Mesmer version.\n"
    " -t  Output a table of available plugins.\n"
    " -T  As -t but with descriptions, if available.\n\n"
    << endl;
}
string version()
{
  stringstream ss;
  ss << "Mesmer v" << MESMER_VERSION << " compiled: -- " << __DATE__ << " -- " << __TIME__ << endl;
  return ss.str();
}

string platform()
{
  string OS;
#ifdef WIN32
  OS = "Windows";
#else
  OS = "Linux";
#endif
  return OS;
}

void banner(size_t nRanks)
{
  cinfo << endl;
  cinfo << "            /|      /|                        \\   " << endl;
  cinfo << "           / |     / | ________________________\\  " << endl;
  cinfo << "          /  |    /  |  _   _             _   _   " << endl;
  cinfo << "         /   |   /   | / \\ / \\ |\\     /| / \\ / \\  " << endl;
  cinfo << "        /    |  /    | \\_/ \\_  | \\   / | \\_/ |_/  " << endl;
  cinfo << "       /     | /     | /     \\ |  \\ /  | /   | \\  " << endl;
  cinfo << "    \\_/      |/      | \\_/ \\_/ |   V   | \\_/ |  \\ " << endl;
  cinfo << endl;
  cinfo << "     Mesmer: Master Equation Solver for Multi-Energy well Reactions" << endl;
  cinfo << "     " << version() << endl;
  cinfo << "     Copyright (C) 2009-2018 by" << endl;
  cinfo << "     Struan H. Robertson, Chris Morley, Chi-Hsiu Liang," << endl;
  cinfo << "     David R. Glowacki, Michael J. Pilling and contributors" << endl;
  cinfo << endl;
  cinfo << "     Mesmer is free software: you can redistribute it and/or modify" << endl;
  cinfo << "     it under the terms of the GNU Public License as published by" << endl;
  cinfo << "     the Free Software Foundation, either version 2 of the License," << endl;
  cinfo << "     or (at your option) any later version." << endl;
  cinfo << endl;
  cinfo << "     But WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
  cinfo << "     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
  cinfo << "     GNU Public License for more details." << endl;
  cinfo << endl;
  cinfo << "    You should have received a copy of the GNU Public License" << endl;
  cinfo << "    along with Mesmer.  If not, see <http://www.gnu.org/licenses/>." << endl;
  cinfo << endl;
  cinfo << "                     Platform:        " << platform() << endl;
  cinfo << "                     Execution Date:  " << date();
  cinfo << "                     Number of ranks: " << nRanks << endl;
  cinfo << endl;
}

bool QACompare(string infilename, bool NOptionUsed)
{
  string::size_type pos = infilename.find_last_of("/\\:");
  string filename = pos != string::npos ? infilename.substr(pos) : infilename; //name without path
  infilename.erase(pos + 1); //path; may erase everything if has no path
  pos = filename.find('.');
  filename.erase(pos + 1).append("test"); // e.g. "ipropyl_reservoir.test"
  if (!NOptionUsed)
    filename = "mesmer.test";

  //Access the baseline file
  infilename.append(TestFolder).append(filename);
  ifstream QAfile(infilename.c_str());
  ifstream CurrentTest(filename.c_str());
  if (!QAfile || !CurrentTest)
  {
    cerr << "Cannot open " << infilename << " or " << filename << endl;
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
Source                                                 Destination
------------------------------------------------------------------------------------------
Main XML output, including functions in                Decided by -o option, usually file,
IPersist like XmlWriteElement()                        or cout for piping

Error messages    cerr <<...    or                     Console (unless -w0)
meErrorLog.ThrowError( , ,obError)   and mesmer.log

Warning messages  cwarn <<...   or                     Console (unless -w0 or w1)
meErrorLog.ThrowError( , ,obWarn)    and mesmer.log

Logging messages  cinfo <<...   or                     Console if -w2
meErrorLog.ThrowError( , ,obInfo)    and mesmer.log

QA test output    ctest <<...                          mesmer.test

Temporary debug   cout <<...                           File when redirected from
Use for bulky debug output                             commandline with >
																											 (Otherwise console)

Temporary debug   clog <<...                           Console
Use for short debug output


Output to mesmer.log can be turned off with the -g option.
The files will be overwrtten each run.
A comparison with a previously stored version of mesmer.test can be made with the -q option.

*/
