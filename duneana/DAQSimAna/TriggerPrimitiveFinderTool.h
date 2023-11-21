#ifndef TriggerPrimitiveFinderTool_h
#define TriggerPrimitiveFinderTool_h

#include <vector>
#include <iostream>

class TriggerPrimitiveFinderTool {
 
public:
    struct Hit
    {
    Hit(int _channel, int _startTime, int _peakCharge, int _SADC, int _peakTime,  int _timeOverThreshold)
            : channel(_channel),
              startTime(_startTime),
              peakCharge(_peakCharge),
	      SADC(_SADC),
	      peakTime(_peakTime),
              timeOverThreshold(_timeOverThreshold)
            {}
        int channel;
        int startTime;
        int peakCharge;
        int SADC;
        int peakTime; 
        int timeOverThreshold;
    };

    virtual ~TriggerPrimitiveFinderTool() =default;

    virtual std::vector<TriggerPrimitiveFinderTool::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& collection_samples) = 0;
 
};

#endif // include guard
