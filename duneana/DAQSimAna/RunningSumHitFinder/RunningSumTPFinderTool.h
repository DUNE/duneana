#ifndef RunningSumTPFinderTool_h
#define RunningSumTPFinderTool_h

#include <vector>
#include <iostream>

class RunningSumTPFinderTool {
 
 public:
  struct Hit
  {
  Hit(int _channel, int _startTime, int _charge, int _timeOverThreshold)
  : channel(_channel),
      startTime(_startTime),
      charge(_charge),
      timeOverThreshold(_timeOverThreshold)
    {}
    int channel;
    int startTime;
    int charge;
    int timeOverThreshold;
  };

  virtual ~RunningSumTPFinderTool() =default;

  virtual std::vector<RunningSumTPFinderTool::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& collection_samples) = 0;
 
};

#endif // include guard
