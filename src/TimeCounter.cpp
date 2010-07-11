//Time count function for Mesmer
#include "TimeCounter.h"

namespace mesmer
{
  std::string TimeCount::setTimeStamp(const std::string& timeStampName, unsigned int &timeElapsed)
  {
    std::time_t tnow; std::time(&tnow);
    EventObj thisEvent(tnow, timeStampName);
    TimeMap.push_back(thisEvent);
    struct std::tm* timeinfo; timeinfo = std::localtime (&tnow);
    char buffer [20]; std::strftime(buffer,20,"%Y%m%d_%H%M%S", timeinfo);
    std::string myTime(buffer);
    if (TimeMap.size() > 1) {
      timeElapsed = static_cast<unsigned int>(TimeMap[TimeMap.size()-1].timeStamp - TimeMap[TimeMap.size()-2].timeStamp) ;
    }
    else timeElapsed = 0;
    return myTime;
  }

  std::string TimeCount::setTimeStamp(const std::string& timeStampName)
  {
    std::time_t tnow; std::time(&tnow);
    EventObj thisEvent(tnow, timeStampName);
    TimeMap.push_back(thisEvent);
    struct std::tm* timeinfo; timeinfo = std::localtime (&tnow);
    char buffer [20]; std::strftime(buffer,20,"%Y%m%d_%H%M%S", timeinfo);
    std::string myTime(buffer);
    return myTime;
  }

  std::string TimeCount::getTimeStamp(const std::string& timeStampName)
  {
    std::time_t tStamp; bool foundTimeStamp = false;
    for (unsigned int i = 0; i < TimeMap.size(); ++i){
      if (TimeMap[i].stampName == timeStampName){
        tStamp = TimeMap[i].timeStamp;
        foundTimeStamp = true;
        break;
      }
    }
    if (!foundTimeStamp){
      std::string ErrTime = "xxxxxxxx_xxxxxx";
      return ErrTime;
    }
    
    struct std::tm* timeinfo; timeinfo = std::localtime(&tStamp);
    char buffer [20]; std::strftime(buffer,20,"%Y%m%d_%H%M%S", timeinfo);
    std::string myTime(buffer); return myTime;
  }

  std::string TimeCount::getTimeStamp(const std::string& timeStampName, unsigned int &timeElapsed)
  {
    std::time_t tStamp; bool foundTimeStamp = false;
    for (unsigned int i = 0; i < TimeMap.size(); ++i){
      if (TimeMap[i].stampName == timeStampName){
        tStamp = TimeMap[i].timeStamp;
        foundTimeStamp = true;
        if (i > 0) timeElapsed = static_cast<unsigned int>(TimeMap[i].timeStamp - TimeMap[i-1].timeStamp);
        else timeElapsed = 0;
        break;
      }
    }
    if (!foundTimeStamp){
      std::string ErrTime = "xxxxxxxx_xxxxxx";
      return ErrTime;
    }
    
    struct std::tm* timeinfo; timeinfo = std::localtime(&tStamp);
    char buffer [20]; std::strftime(buffer,20,"%Y%m%d_%H%M%S", timeinfo);
    std::string myTime(buffer); return myTime;
  }

  std::ostream& operator<<(std::ostream& os, TimeCount& MTC) //function to dump all the time stamps
  {
    os << "\nTime stamps:\n";

    for (TimeCount::TimeIter curTimeIter = MTC.TimeMap.begin(); curTimeIter != MTC.TimeMap.end(); ++curTimeIter)
    {
      std::time_t tStamp; tStamp = curTimeIter->timeStamp;
      struct std::tm* timeinfo; timeinfo = std::gmtime(&tStamp);
      char buffer [20]; std::strftime(buffer,20,"%Y%m%d_%H%M%S", timeinfo); std::string myTime(buffer);
      os << myTime << " -- " << curTimeIter->stampName << std::endl;
    }
    return os;
  }
};


