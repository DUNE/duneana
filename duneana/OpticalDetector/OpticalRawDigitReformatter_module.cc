// Ben Jones, MIT, 2013
//
// This module finds periods of time-localized activity
// from the optical system, called Flashes.

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/OpticalDetectorData/OpticalRawDigit.h"
#include "lardataobj/OpticalDetectorData/OpticalTypes.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/BeamGateInfo.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

// C++ Includes
#include <cstring>

namespace opdet {

  class OpticalRawDigitReformatter : public art::EDProducer {
  public:
    // Standard constructor and destructor for an ART module.
    explicit OpticalRawDigitReformatter(const fhicl::ParameterSet&);

    // The producer routine, called once per event.
    void produce(art::Event&);

  private:
    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpDet collection
    std::string fGenModule;

    std::vector<std::string> CategoryLabels;
  };

}

namespace opdet {
  DEFINE_ART_MODULE(OpticalRawDigitReformatter)
}

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpticalRawDigitReformatter::OpticalRawDigitReformatter(const fhicl::ParameterSet& pset)
    : EDProducer{pset}
  {
    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");
    fGenModule = pset.get<std::string>("GenModule");

    CategoryLabels.push_back("Undefined");
    CategoryLabels.push_back("HighGain");
    CategoryLabels.push_back("LowGain");
    CategoryLabels.push_back("LogicPulse");
    CategoryLabels.push_back("FEMCosmicHighGain");
    CategoryLabels.push_back("FEMCosmicLowGain");
    CategoryLabels.push_back("FEMCosmicLogicPulse");
    CategoryLabels.push_back("FEMBeamHighGain");
    CategoryLabels.push_back("FEMBeamLowGain");
    CategoryLabels.push_back("FEMBeamLogicPulse");
    CategoryLabels.push_back("BeamPMTTrigger");
    CategoryLabels.push_back("CosmicPMTTrigger");

    // One for each category
    for (auto label : CategoryLabels)
      produces<std::vector<raw::OpDetWaveform>>(label);
  }

  //-----------------------------------------------------------------------
  void OpticalRawDigitReformatter::produce(art::Event& evt)
  {

    // These are the storage pointers we will put in the event, one per category
    std::vector<std::unique_ptr<std::vector<raw::OpDetWaveform>>> RawOpDetVecs;
    for (unsigned int i = 0; i < CategoryLabels.size(); i++) {
      std::unique_ptr<std::vector<raw::OpDetWaveform>> tmp(new std::vector<raw::OpDetWaveform>);
      RawOpDetVecs.push_back(std::move(tmp));
    }

    std::vector<const sim::BeamGateInfo*> beamGateArray;
    try {
      evt.getView(fGenModule, beamGateArray);
    }
    catch (art::Exception const& err) {
      if (err.categoryCode() != art::errors::ProductNotFound) throw;
    }

    // Read in the OpticalRawDigit collection from the event.
    art::Handle<std::vector<optdata::OpticalRawDigit>> ordHandle;
    evt.getByLabel(fInputModule, ordHandle);
    std::vector<optdata::OpticalRawDigit> const& ord_vec(*ordHandle);

    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    for (auto ord : ord_vec) {
      optdata::Channel_t channel = ord.ChannelNumber();
      optdata::TimeSlice_t timeSlice = ord.TimeSlice();
      optdata::Frame_t frame = ord.Frame();
      optdata::Optical_Category_t category = ord.Category();

      // Use the optical clock to conver timeSlice and frame
      // to an absolute time
      double timeStamp = clock_data.OpticalClock().Time(timeSlice, frame);

      RawOpDetVecs[category]->push_back(raw::OpDetWaveform(timeStamp, channel, ord));
    }

    // Store results into the event
    for (unsigned int i = 0; i < CategoryLabels.size(); i++) {
      // Only store collections which contain waveforms, assign the label
      if (RawOpDetVecs[i]->size() > 0) { evt.put(std::move(RawOpDetVecs[i]), CategoryLabels[i]); }
    }
  }

} // namespace opdet
