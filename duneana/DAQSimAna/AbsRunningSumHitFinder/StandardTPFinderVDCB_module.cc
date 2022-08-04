////////////////////////////////////////////////////////////////////////
// Class:       StandardTPFinderVDCB
// File:        StandardTPFinderVDCB_module.cc
//
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

#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardataobj/RawData/RawDigit.h"

#include "duneana/DAQSimAna/TriggerPrimitiveFinderTool.h"

#include <memory>

class StandardTPFinderVDCB;


class StandardTPFinderVDCB : public art::EDProducer {
public:
  explicit StandardTPFinderVDCB(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StandardTPFinderVDCB(StandardTPFinderVDCB const &) = delete;
  StandardTPFinderVDCB(StandardTPFinderVDCB &&) = delete;
  StandardTPFinderVDCB & operator = (StandardTPFinderVDCB const &) = delete;
  StandardTPFinderVDCB & operator = (StandardTPFinderVDCB &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:
  // The module name of the raw digits we're reading in
  std::string m_inputTag;
  // The actual Service that's doing the trigger primitive finding
  std::unique_ptr<TriggerPrimitiveFinderTool> m_finderCol;
  std::unique_ptr<TriggerPrimitiveFinderTool> m_finderInd;
};


StandardTPFinderVDCB::StandardTPFinderVDCB(fhicl::ParameterSet const & p)
  : EDProducer{p}, 
  m_inputTag(p.get<std::string>("InputTag", "daq")),
  m_finderCol{art::make_tool<TriggerPrimitiveFinderTool>(p.get<fhicl::ParameterSet>("finderCol"))},
  m_finderInd{art::make_tool<TriggerPrimitiveFinderTool>(p.get<fhicl::ParameterSet>("finderInd"))}
{
    produces<std::vector<recob::Hit>>();
    produces<art::Assns<raw::RawDigit, recob::Hit>>();
}

void StandardTPFinderVDCB::produce(art::Event & e)
{
    art::ServiceHandle<geo::Geometry> geo;

    std::vector<std::vector<short>>  induction_samples;
    std::vector<std::vector<short>> collection_samples;
    std::vector<unsigned int>  induction_channel_numbers;
    std::vector<unsigned int>  collection_channel_numbers;
    std::map<raw::ChannelID_t, const raw::RawDigit*> indChanToDigit;
    std::map<raw::ChannelID_t, const raw::RawDigit*> colChanToDigit;

    std::vector<const raw::RawDigit*> digits_handle;
    e.getView("tpcrawdecoder:daq", digits_handle);

    //quick hacky way of getting the associations to work even though this gives wrong ADC values.
    auto const &digits = e.getValidHandle<std::vector<raw::RawDigit>>(m_inputTag);   
 
    for (size_t c = 0; c < digits_handle.size(); ++c){
      const raw::RawDigit* digit = digits_handle[c];
      
      if ( geo->SignalType( digit->Channel() ) == geo::kInduction ){
	indChanToDigit[digit->Channel()] = digit; 
	induction_channel_numbers.push_back( digit->Channel() );
	induction_samples.push_back(digit->ADCs());
      }//get induction samples 
      
      if ( geo->SignalType( digit->Channel() ) == geo::kCollection) {
	colChanToDigit[digit->Channel()] = digit; 
	collection_channel_numbers.push_back( digit->Channel() );
	collection_samples.push_back( digit->ADCs() );
      }//get collection samples 
    }//run over raw digits
    
    // Pass the full list of collection channels to the hit finding algorithm
    std::vector<TriggerPrimitiveFinderTool::Hit> hits_col=m_finderCol->findHits(collection_channel_numbers, collection_samples);
    std::vector<TriggerPrimitiveFinderTool::Hit> hits_ind=m_finderInd->findHits(induction_channel_numbers,  induction_samples);
    
    // Loop over the returned trigger primitives and turn them into recob::Hits
    recob::HitCollectionCreator hcol(e, false /* doWireAssns */, true /* doRawDigitAssns */);
    for(auto const& hit : hits_col){
      const raw::RawDigit* digit=colChanToDigit[hit.channel];
      
      std::vector<geo::WireID> wids = geo->ChannelToWire(hit.channel);
      geo::WireID wid = wids[0];
      
      recob::HitCreator lar_hit(*digit,                                //RAW DIGIT REFERENCE.
				wid,                                       //WIRE ID.
				hit.startTime,                             //START TICK.
				hit.startTime+hit.timeOverThreshold,       //END TICK. 
				hit.timeOverThreshold,                     //RMS.
				hit.startTime + hit.timeOverThreshold*0.5, //PEAK_TIME.
				0,                                         //SIGMA_PEAK_TIME.
				hit.peakCharge,                            //PEAK_AMPLITUDE.
				0,                                         //SIGMA_PEAK_AMPLITUDE.
				hit.SADC,                                  //HIT_INTEGRAL.
				0,                                         //HIT_SIGMA_INTEGRAL.
				hit.SADC,                                  //SUMMED CHARGE. 
				0,                                         //MULTIPLICITY.
				0,                                         //LOCAL_INDEX.
				0,                                         //WIRE ID. (?)
				0                                          //DEGREES OF FREEDOM.
				);
      hcol.emplace_back(std::move(lar_hit), art::Ptr<raw::RawDigit>{digits, 0});
    }
    // Loop over the returned trigger primitives and turn them into recob::Hits
    for(auto const& hit : hits_ind){
      const raw::RawDigit* digit=indChanToDigit[hit.channel];
      std::vector<geo::WireID> wids = geo->ChannelToWire(hit.channel);
      geo::WireID wid = wids[0];
      
      recob::HitCreator lar_hit(*digit,                                //RAW DIGIT REFERENCE.
				wid,                                       //WIRE ID.
				hit.startTime,                             //START TICK.
				hit.startTime+hit.timeOverThreshold,       //END TICK. 
				hit.timeOverThreshold,                     //RMS.
				hit.startTime + hit.timeOverThreshold*0.5, //PEAK_TIME.
				0,                                         //SIGMA_PEAK_TIME.
				hit.peakCharge,                            //PEAK_AMPLITUDE.
				0,                                         //SIGMA_PEAK_AMPLITUDE.
				hit.SADC,                                  //HIT_INTEGRAL.
				0,                                         //HIT_SIGMA_INTEGRAL.
				hit.SADC,                                  //SUMMED CHARGE. 
				0,                                         //MULTIPLICITY.
				1,                                         //LOCAL_INDEX.
				0,                                         //WIRE ID.
				0                                          //DEGREES OF FREEDOM.
				);
      hcol.emplace_back(std::move(lar_hit), art::Ptr<raw::RawDigit>{digits, 0});
    }
    hcol.put_into(e);
}

DEFINE_ART_MODULE(StandardTPFinderVDCB)
