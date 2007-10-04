#ifndef GUARD_TimeCounter_h
#define GUARD_TimeCounter_h
//Time count header of Mesmer
//This header must be be p_thread safe and parallel computing safe
#include <sstream>
#include <iostream>
#include <time.h>
#include <vector>

//------------------------
//Usage:
//
//    Variable decalaration
//    \code
//    TimeCount events; unsigned int timeElapsed;
//    std::string thisEvent;
//    \endcode
//
//    Time stamping
//    \code
//    thisEvent = "Build Collison Matrix";
//    cout << thisEvent << " at " << events.setTimeStamp(thisEvent) << endl;
//    \endcode
//    
//    OR
//    
//    \code
//    thisEvent = "Build Collison Matrix";
//    std::cout << thisEvent << " at " << events.setTimeStamp(thisEvent, timeElapsed) << " -- Time elapsed: " << timeElapsed << " seconds.\n";
//    \endcode
//
//    Time stamp dumping
//    \code
//    cout << events << endl;
//    \endcode
//
//------------------------


template<typename T>
std::string toString(T t)
{
  std::ostringstream s; s << t; return s.str();
}

namespace mesmer
{
  class EventObj
  {
    public:
      time_t timeStamp;
      std::string stampName;
      EventObj(const time_t& tt, const std::string& name):timeStamp(tt), stampName(name){}
  };

  class TimeCount
  {
    private:

      std::vector<EventObj> TimeMap;
      typedef std::vector<EventObj>::iterator TimeIter;

    public:

      std::string setTimeStamp(const std::string& timeStampName);
      std::string TimeCount::setTimeStamp(const std::string& timeStampName, unsigned int &timeElapsed);
      std::string getTimeStamp(const std::string& timeStampName);
      std::string TimeCount::getTimeStamp(const std::string& timeStampName, unsigned int &timeElapsed);
      friend std::ostream& operator<<(std::ostream& os, TimeCount& MTC);
  };
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
