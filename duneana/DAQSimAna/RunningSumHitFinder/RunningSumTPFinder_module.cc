////////////////////////////////////////////////////////////////////////
// Class:       RunningSumTPFinder
// Plugin Type: producer (art v2_10_03)
// File:        RunningSumTPFinder_module.cc
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

#include "duneana/DAQSimAna/RunningSumHitFinder/RunningSumTPFinderTool.h"

#include <memory>

class RunningSumTPFinder;


class RunningSumTPFinder : public art::EDProducer {
public:
  explicit RunningSumTPFinder(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RunningSumTPFinder(RunningSumTPFinder const &) = delete;
  RunningSumTPFinder(RunningSumTPFinder &&) = delete;
  RunningSumTPFinder & operator = (RunningSumTPFinder const &) = delete;
  RunningSumTPFinder & operator = (RunningSumTPFinder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:
    // The module name of the raw digits we're reading in
    std::string m_inputTag;
    // The actual Service that's doing the trigger primitive finding
    std::unique_ptr<RunningSumTPFinderTool> m_finder1;
    std::unique_ptr<RunningSumTPFinderTool> m_finder2;
};


RunningSumTPFinder::RunningSumTPFinder(fhicl::ParameterSet const & p)
  : EDProducer{p}, m_inputTag(p.get<std::string>("InputTag", "daq")),
  m_finder1{art::make_tool<RunningSumTPFinderTool>(p.get<fhicl::ParameterSet>("finder1"))},
  m_finder2{art::make_tool<RunningSumTPFinderTool>(p.get<fhicl::ParameterSet>("finder2"))}
{
    produces<std::vector<recob::Hit>>();
    produces<art::Assns<raw::RawDigit, recob::Hit>>();
}

void RunningSumTPFinder::produce(art::Event & e)
{
    auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(m_inputTag);
    auto& digits_in =*digits_handle;

    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    std::vector<std::vector<short>>  induction_samples;
    std::vector<std::vector<short>> collection_samples;
    std::vector<unsigned int>  induction_channel_numbers;
    std::vector<unsigned int> collection_channel_numbers;
    std::map<raw::ChannelID_t, const raw::RawDigit*> indChanToDigit;
    std::map<raw::ChannelID_t, const raw::RawDigit*> colChanToDigit;
    for(auto&& digit: digits_in){
        // Select just the collection channels for the primitive-finding algorithm
        const geo::SigType_t sigType = wireReadout.SignalType(digit.Channel());
        if(sigType==geo::kInduction){
            indChanToDigit[digit.Channel()]=&digit;
            induction_channel_numbers.push_back(digit.Channel());
            induction_samples.push_back(digit.ADCs());
        }
        if(sigType==geo::kCollection){
            colChanToDigit[digit.Channel()]=&digit;
            collection_channel_numbers.push_back(digit.Channel());
            collection_samples.push_back(digit.ADCs());
        }
    }

    // Pass the full list of collection channels to the hit finding algorithm
    std::vector<RunningSumTPFinderTool::Hit> hits1=m_finder1->findHits( induction_channel_numbers,  induction_samples);
    std::vector<RunningSumTPFinderTool::Hit> hits2=m_finder2->findHits(collection_channel_numbers, collection_samples);

    // Loop over the returned trigger primitives and turn them into recob::Hits
    recob::HitCollectionCreator hcol(e, false /* doWireAssns */, true /* doRawDigitAssns */);
    for(auto const& hit : hits1){
        const raw::RawDigit* digit=indChanToDigit[hit.channel];
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
                              hit.charge,                           //SUMMED CHARGE. ROI
       			      hit.charge,                           //SUMMED CHARGE.HIT TO BE FIXED 
                              0,                                    //MULTIPLICITY.
                              0,                                    //LOCAL_INDEX.
                              0,                                    //WIRE ID.
                              0                                     //DEGREES OF FREEDOM.
            );
        hcol.emplace_back(std::move(lar_hit), art::Ptr<raw::RawDigit>{digits_handle, 0});
    }
    // Loop over the returned trigger primitives and turn them into recob::Hits
    for(auto const& hit : hits2){
        const raw::RawDigit* digit=colChanToDigit[hit.channel];
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
                              hit.charge,                           //SUMMED CHARGE. ROI
			      hit.charge,                           //SUMMED CHARGE.HIT TO BE FIXED 
                              0,                                    //MULTIPLICITY.
                              0,                                    //LOCAL_INDEX.
                              0,                                    //WIRE ID.
                              0                                     //DEGREES OF FREEDOM.
            );
        hcol.emplace_back(std::move(lar_hit), art::Ptr<raw::RawDigit>{digits_handle, 0});
    }
    hcol.put_into(e);
}

DEFINE_ART_MODULE(RunningSumTPFinder)
