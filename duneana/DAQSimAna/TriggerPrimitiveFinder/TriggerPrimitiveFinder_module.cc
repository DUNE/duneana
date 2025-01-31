////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveFinder
// Plugin Type: producer (art v2_10_03)
// File:        TriggerPrimitiveFinder_module.cc
//
// Generated at Tue Jun  5 06:03:24 2018 by Philip Rodrigues using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/WireReadout.h"
#include "lardata/ArtDataHelper/HitCreator.h"

#include "duneana/DAQSimAna/TriggerPrimitiveFinder/TriggerPrimitiveFinderTool.h"

#include <memory>

class TriggerPrimitiveFinder;


class TriggerPrimitiveFinder : public art::EDProducer {
public:
  explicit TriggerPrimitiveFinder(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TriggerPrimitiveFinder(TriggerPrimitiveFinder const &) = delete;
  TriggerPrimitiveFinder(TriggerPrimitiveFinder &&) = delete;
  TriggerPrimitiveFinder & operator = (TriggerPrimitiveFinder const &) = delete;
  TriggerPrimitiveFinder & operator = (TriggerPrimitiveFinder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:
    // The module name of the raw digits we're reading in
    std::string m_inputTag;
    // The actual Service that's doing the trigger primitive finding
    std::unique_ptr<TriggerPrimitiveFinderTool> m_finder;
};


TriggerPrimitiveFinder::TriggerPrimitiveFinder(fhicl::ParameterSet const & p)
    : EDProducer{p}, m_inputTag(p.get<std::string>("InputTag", "daq")),
      m_finder{art::make_tool<TriggerPrimitiveFinderTool>(p.get<fhicl::ParameterSet>("finder"))}
{
    produces<std::vector<recob::Hit>>();
    produces<art::Assns<raw::RawDigit, recob::Hit>>();
}

void TriggerPrimitiveFinder::produce(art::Event & e)
{
    auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(m_inputTag);
    auto& digits_in =*digits_handle;

    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    std::vector<std::vector<short>> collection_samples;
    std::vector<unsigned int> channel_numbers;
    std::map<raw::ChannelID_t, const raw::RawDigit*> chanToDigit;
    for(auto&& digit: digits_in){
        // Select just the collection channels for the primitive-finding algorithm
        const geo::SigType_t sigType = wireReadout.SignalType(digit.Channel());
        if(sigType==geo::kCollection){
            chanToDigit[digit.Channel()]=&digit;
            channel_numbers.push_back(digit.Channel());
            collection_samples.push_back(digit.ADCs());
        }
    }

    // Pass the full list of collection channels to the hit finding algorithm
    std::vector<TriggerPrimitiveFinderTool::Hit> hits=m_finder->findHits(channel_numbers, collection_samples);

    // Loop over the returned trigger primitives and turn them into recob::Hits
    recob::HitCollectionCreator hcol(e, false /* doWireAssns */, true /* doRawDigitAssns */);
    for(auto const& hit : hits){
        const raw::RawDigit* digit=chanToDigit[hit.channel];
        if(!digit){
            std::cout << "No digit with channel " << hit.channel << " found. Did you set the channel correctly?" << std::endl;
        }
        std::vector<geo::WireID> wids = wireReadout.ChannelToWire(hit.channel);
        geo::WireID wid = wids[0];

        recob::HitCreator lar_hit(*digit,                           //RAW DIGIT REFERENCE.
                              wid,                                  //WIRE ID.
                              hit.startTime,                        //START TICK.
                              hit.startTime+hit.timeOverThreshold,  //END TICK. 
                              hit.timeOverThreshold,                //RMS.
                              hit.startTime,                        //PEAK_TIME.
                              0,                                    //SIGMA_PEAK_TIME.
                              0,                                    //PEAK_AMPLITUDE.
                              0,                                    //SIGMA_PEAK_AMPLITUDE.
                              hit.charge,                           //HIT_INTEGRAL.
                              0,                                    //HIT_SIGMA_INTEGRAL.
                              hit.charge,                           //SUMMED CHARGE ROI. 
			      hit.charge,                           //SUMMED CHARGE HIT. TO BE FIXED 
                              0,                                    //MULTIPLICITY.
                              0,                                    //LOCAL_INDEX.
                              0,                                    //WIRE ID.
                              0                                     //DEGREES OF FREEDOM.
            );
        hcol.emplace_back(std::move(lar_hit), art::Ptr<raw::RawDigit>{digits_handle, 0});
    }
    hcol.put_into(e);
}

DEFINE_ART_MODULE(TriggerPrimitiveFinder)
