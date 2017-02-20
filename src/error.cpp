/**********************************************************************
error.cpp - Handle error messages.
Details at end of error.cpp

Copyright (C) 2007 by Chris Morley
Based on original code from OpenBabel 
Copyright (C) 2002 by Stefan Kebekus, 2002-2006 by Geoffrey R. Hutchison
***********************************************************************/

#include <algorithm>
#include "error.h"

using namespace std;

namespace mesmer
{

  //Global output streams which by default go to clog
  //They may be redirected by an OStreamRedirector object 
  ostream cwarn(clog.rdbuf());
  ostream cinfo(clog.rdbuf());
  ostream ctest(clog.rdbuf());
	ostringstream stest("");

  //Global message handler
  MessageHandler meErrorLog;
  //Static member variable initialization
  std::vector<std::string> ErrorContext::stack;

  //////////////////////////////////////////////////////////////////
  MessageHandler::MessageHandler(ostream* logstream) :
  _outputLevel(obWarning), _logStream(logstream), _logging(true)
  {}

  ///////////////////////////////////////////////////////////////////
  void MessageHandler::ThrowError(const std::string &context, 
    const std::string &errorMsg,
    obMessageLevel level,
    errorQualifier qualifier) const
  {
    if (!_logging || errorMsg.empty())
      return;
    string txt;
    if(errorMsg.size()>1)
    { 
      if(!context.empty())
        txt = GetLevelText(level) + "In " + context;
      else if (!_defaultContext.empty())
        txt = "In " + _defaultContext + ": ";
    }
    txt += errorMsg;

    //if no new line at end, add one
    //if(txt[txt.size()-1]!='\n')
    //  txt += '\n';

    //ignore errors with a onceOnly qualifier if they have occurred before
    static vector<string> _messageList;
    if(qualifier==onceOnly)
    {
      if(find(_messageList.begin(), _messageList.end(), txt)==_messageList.end())
        _messageList.push_back(txt); //first time
      else
        return; //subsequent times - no error message
    }

    if(_logStream)
      *_logStream << txt << flush;

    if (level <= _outputLevel)
      clog << txt << flush; //write to console
  }
  ////////////////////////////////////////////////////////////////////
  ErrorContext::ErrorContext(const std::string& txt)
  {
    stack.push_back(txt);
    meErrorLog.SetContext(txt);
  }

  ErrorContext::~ErrorContext()
  {
    if(!stack.empty())
    {
      stack.pop_back();
      meErrorLog.SetContext(stack.empty() ? "" : stack.back());
    }
    else
      cerr << "Context Stack empty!" << endl;
  }

  ////////////////////////////////////////////////////////////////////
  int obLogBuf::sync()
  {
    //Being called after every << which is inconvenient.
    //Pass the string to the handler only when it ends in newline 
    string s = str();
    if(!s.empty() && s[s.size()-1]=='\n')
    {
      _handler->ThrowError("", str(), _messageLevel, _qualifier);
      str(std::string()); // clear the buffer
    }
    _qualifier = always;
    return 0;
  }

  ////////////////////////////////////////////////////////////////////
  OStreamRedirector::OStreamRedirector(MessageHandler* handler, std::ostream* Newctest, bool nologging) :
  _cerrbuf(obError, handler), _cwarnbuf(obWarning, handler), _cinfobuf(obInfo, handler)
  {
    //Save original buffers for output streams...
    _oldcerrbuf = std::cerr.rdbuf();
    _oldcwarnbuf = cwarn.rdbuf();
    _oldcinfobuf = cinfo.rdbuf();

    // and replace them with buffers which will send messages to the message handler
    std::cerr.rdbuf(&_cerrbuf);
    cwarn.rdbuf(&_cwarnbuf);
    cinfo.rdbuf(&_cinfobuf);

    //Optional redirection of ctest.
    //This is ordinary redirection, messages are not sent to throwError.
    if(Newctest)
    {
      _oldctestbuf = ctest.rdbuf();
      ctest.rdbuf(Newctest->rdbuf());
    }
    if(nologging)
    {  
      //turn off logging and test streams
      std::cerr.rdbuf(NULL);
      cwarn.rdbuf(NULL);
      cinfo.rdbuf(NULL);
      ctest.rdbuf(NULL);
    }
  }
  OStreamRedirector::~OStreamRedirector()
  {
    if(_oldcerrbuf)
      std::cerr.rdbuf(_oldcerrbuf);
    if(_oldcwarnbuf)
      cwarn.rdbuf(_oldcwarnbuf);
    if(_oldcinfobuf)
      cinfo.rdbuf(_oldcinfobuf);
    if(_oldctestbuf)
      ctest.rdbuf(_oldctestbuf);
  }

} // end namespace mesmer

/*
MessageHandler class controls the display of error messages. There is a
global instance meErrorLog, but other instances can be made if necessary.

All the error messages a handler receives are sent to an output stream that is
specified either when the MessageHandler object is constructed or by calling 
SetLogStream(). It is probably attached to a logfile.

Additionally, messages graded above a certain severity (see obMessageLevel) are
sent to cout (usually the console). The level is set by SetOutputLevel().

In the code, the most flexible way of sending an error message is to use
MessageHandler::ThrowError(string context, string msg, obMessageLevel level);
The context string gives locally significant information, such as the routine
where the call is made from (use __FUNCTION__) or the name of the structure
being worked on. If context is empty, a default context, set by calling
SetContext() is used instead. So a call:
meErrorLog.ThrowError(molID, "Molecule already defined", obError);
may appear in the log file and the console as:
Error in n-pentyl: Molecule already defined

An alternative method of sending error messages from code is to write to
cerr, or some extra similar output streams, cwarn and cinfo. The severities 
for these are obError, obWarning and obInfo respectively and the context is 
empty, so that a default context is used. Using these output stream simplifies
the formatting of the message and the output of numerical data. It is possible
to group several messages under one context header by not flushing the buffer,
i.e. don't use endl.

To use this to send messages, make an instance of OStreamRedirector, specifying
the message handler, e.g. OStreamRedirector osr(&meErrorLog). The streams will
be restored to their old destinations when the OStreamRedirector instance goes
out of scope. 

The reuse of cerr means that old code written using it will automatically use
the new system. This eases its introduction, since it is operational without
having to modify the message calling statements.

MessageHandler::ThrowError() now has an optional parameter which if set to onceOnly
will prevent the same error message being output more than once. Note that the
error message here includes the context, so that if several molecules contained
the same error, then there would be an error message for each, even if onceOnly had
been used. When using cerr, cwarn and cinfo with <<, onceOnly can be set by the
manipulator 'once', e.g.
  cerr << "Plotting error" << once << endl;
The whole of the text sent to error before endl is regarded as the error message
which is checked for duplication.

The optional parameter on the OStreamRedirector constructor is the stream to which
ctest should be redirected. This output is independent and not sent to 
mesmer.log or the console.

Recording of error messages can be switched off and on using 
MessageHandler::StopLogging() and StartLogging(). Logging is on by default
with a severity level of obInfo which outputs messages of severities obInfo,
obWarning and obError.
*/



//! \file oberror.cpp
//! \brief Handle error messages, warnings, notices, etc.
//!  Implements MessageHandler class.
