// A template for creating classes to find trigger primitives for
// testing DAQ algorithms

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "duneana/DAQSimAna/TriggerPrimitiveFinder/TriggerPrimitiveFinderTool.h"

class TriggerPrimitiveFinderTemplate : public TriggerPrimitiveFinderTool {
public:
    explicit TriggerPrimitiveFinderTemplate(fhicl::ParameterSet const & p);

    virtual std::vector<TriggerPrimitiveFinderTool::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& collection_samples);
    

private:

    // Example setting. This can be read from the input fcl file: see
    // the implementation of the constructor
    int m_threshold;
};


TriggerPrimitiveFinderTemplate::TriggerPrimitiveFinderTemplate(fhicl::ParameterSet const & p)
    // set m_threshold from the fcl file.
    // p.get() can take a second argument which is the default value
    : m_threshold(p.get<unsigned int>("Threshold"))
{
    // Do any other constructing you want here
}

// collection_samples is a vector of waveforms. Each waveform is
// itself a vector of all the ADC counts, in order, from a single
// collection wire. So, eg collection_samples[2] is a vector of the
// ADC samples in time order from the 2th collection wire
//
// (Trigger primitives are per-wire objects, so you'll mostly only act
// on one channel at a time, but you're given all of them here so you
// can do things like removing coherent noise)
//
// channel_numbers is a vector of the channel id of each channel in collection_samples, in the same order
std::vector<TriggerPrimitiveFinderTool::Hit>
TriggerPrimitiveFinderTemplate::findHits(const std::vector<unsigned int>& channel_numbers, 
                                   const std::vector<std::vector<short>>& collection_samples)
{
    // Make a return vector of TriggerPrimitiveFinderTool::Hit (see TriggerPrimitiveFinderTool.h)
    auto hits=std::vector<TriggerPrimitiveFinderTool::Hit>();

    // Construct Hits like:
    //
    // 1. TriggerPrimitiveFinderTool::Hit hit(channel, startTime, charge, timeOverThreshold)
    // or
    // 2. TriggerPrimitiveFinderTool::Hit hit;
    //    hit.channel=...;
    //    hit.startTime=...;
    //    hit.charge=...;
    //    hit.timeOverThreshold=...;
    //
    // Then add each one to the vector of hits

    return hits;
}

DEFINE_ART_CLASS_TOOL(TriggerPrimitiveFinderTemplate)
