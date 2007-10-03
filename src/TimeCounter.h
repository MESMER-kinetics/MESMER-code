#ifndef GUARD_TimeCounter_h
#define GUARD_TimeCounter_h
//Time count header of Mesmer
//This header must be be p_thread safe and parallel computing safe
#include <iostream>
#include <time.h>
#include <vector>

//------------------------
//Usage:
//
//    Variable decalaration
//    \code
//    TimeCount events;
//    string thisEvent;
//    \endcode
//
//    Time stamping
//    \code
//    thisEvent = "Build Collison Matrix";
//    cout << thisEvent << " at " << events.setTimeStamp(thisEvent) << endl;
//    \endcode
//
//    Time stamp dumping
//    \code
//    cout << events << endl;
//    \endcode
//
//------------------------

using namespace std;

template<typename T>
string toString(T t)
{
  ostringstream s; s << t; return s.str();
}

namespace mesmer
{
  class EventObj
  {
    public:
      time_t timeStamp;
      string stampName;
      EventObj(const time_t& tt, const string& name):timeStamp(tt), stampName(name){}
  };

  class TimeCount
  {
    private:

      vector<EventObj> TimeMap;
      typedef vector<EventObj>::iterator TimeIter;

    public:

      string setTimeStamp(const string& timeStampName)
      {
        time_t tnow; time(&tnow);
        EventObj thisEvent(tnow, timeStampName);
        TimeMap.push_back(thisEvent);
        struct tm* timeinfo; timeinfo = localtime (&tnow);
        char buffer [20]; strftime (buffer,20,"%Y%m%d_%H%M%S", timeinfo);
        string myTime(buffer);
        if (TimeMap.size() > 1)
        {
          unsigned long seconds = TimeMap[TimeMap.size()-1].timeStamp - TimeMap[TimeMap.size()-2].timeStamp;
          myTime += " -- Time elapsed: "; myTime += toString(seconds); myTime += " seconds.";
        }
        return myTime;
      }

      string getTimeStamp(const string& timeStampName)
      {
        time_t tStamp;
        TimeIter curTimeIter = TimeMap.begin();
        while(curTimeIter != TimeMap.end()){
          if (curTimeIter->stampName == timeStampName)
          {
            tStamp = curTimeIter->timeStamp;
            break;
          }
        }
        if (curTimeIter == TimeMap.end()){
          string ErrTime = "xxxxxxxx_xxxxxx";
          return ErrTime;
        }

        struct tm* timeinfo; timeinfo = localtime (&tStamp);
        char buffer [20]; strftime (buffer,20,"%Y%m%d_%H%M%S", timeinfo);
        string myTime(buffer); return myTime;
      }

      friend ostream& operator <<(ostream &os, TimeCount& MTC);
  };

  ostream& operator<< (ostream &os, TimeCount& MTC) //function to dump all the time stamps
  {
    os << "\nTime stamps:\n";

    for (TimeCount::TimeIter curTimeIter = MTC.TimeMap.begin(); curTimeIter != MTC.TimeMap.end(); ++curTimeIter)
    {
      time_t tStamp; tStamp = curTimeIter->timeStamp;
      struct tm* timeinfo; timeinfo = gmtime(&tStamp);
      char buffer [20]; strftime (buffer,20,"%Y%m%d_%H%M%S", timeinfo); string myTime(buffer);
      os << myTime << " -- " << curTimeIter->stampName << endl;
    }
    return os;
  }


};

//---------------------------------------------------------------------------------------------------------
//The idea:
//
//    The format of a time stamp should look like "20071002_093401" and will be attached to the file.
//    The file name will be "<reaction_name>.<time_previous>.<time_current>.xml" to indicate
//    that this file is completed at time <time_current> and it was created from the calculation of the date
//    from a file completed at <time_previous>.
//
//    So, practically a Mesmer XML file will contain a reaction name, a previous time stamp, and a current
//    (completed) time stamp. If the file was created from a file without a time stamp, it will set
//    its <time_previous> to the initiation time of its calculation.
//
//    For example, a filename could look like:
//
//    pentyl_isomerization.20071002_093401.20071002_094825.xml
//
//    It could have initiated itself from from a file called something like the followings:
//
//    pentyl_isomerization.20071002_092341.20071002_093401.xml
//    pentyl_isomerization.20071002_093401.xml
//    pentyl_isomerization.xml
//
//    In this way, the files can be sorted easily.
//---------------------------------------------------------------------------------------------------------

#endif //GUARD_TimeCounter_h
