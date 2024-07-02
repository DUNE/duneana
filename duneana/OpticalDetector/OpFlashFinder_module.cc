// -*- mode: c++; c-basic-offset: 2; -*-
// Ben Jones, MIT, 2013
//
// This module finds periods of time-localized activity
// from the optical system, called Flashes, using OpHits as an input.
//
// Modified to make it more detector agnostic
// by Gleb Sinev, Duke, 2015
//

// LArSoft includes
#include "larana/OpticalDetector/OpFlashAlg.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes

// C++ Includes
#include <memory>
#include <string>

namespace opdet {

  class OpFlashFinder : public art::EDProducer {
  public:
    // Standard constructor and destructor for an ART module.
    explicit OpFlashFinder(const fhicl::ParameterSet&);

    // The producer routine, called once per event.
    void produce(art::Event&);

  private:
    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpHit collection

    Int_t fBinWidth;
    Float_t fFlashThreshold;
    Float_t fWidthTolerance;
    Double_t fTrigCoinc;
  };

}

namespace opdet {
  DEFINE_ART_MODULE(OpFlashFinder)
}

namespace opdet {

  //----------------------------------------------------------------------------
  // Constructor
  OpFlashFinder::OpFlashFinder(const fhicl::ParameterSet& pset) : EDProducer{pset}
  {

    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");

    fBinWidth = pset.get<int>("BinWidth");
    fFlashThreshold = pset.get<float>("FlashThreshold");
    fWidthTolerance = pset.get<float>("WidthTolerance");
    fTrigCoinc = pset.get<double>("TrigCoinc");

    produces<std::vector<recob::OpFlash>>();
    produces<art::Assns<recob::OpFlash, recob::OpHit>>();
  }

  //----------------------------------------------------------------------------
  void OpFlashFinder::produce(art::Event& evt)
  {
    // These are the storage pointers we will put in the event
    auto flashPtr = std::make_unique<std::vector<recob::OpFlash>>();
    auto assnPtr = std::make_unique<art::Assns<recob::OpFlash, recob::OpHit>>();

    // This will keep track of what flashes will assoc to what ophits
    // at the end of processing
    std::vector<std::vector<int>> assocList;

    auto const& geometry(*lar::providerFrom<geo::Geometry>());

    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    // Get OpHits from the event
    auto const opHitHandle = evt.getValidHandle<std::vector<recob::OpHit>>(fInputModule);

    RunFlashFinder(*opHitHandle,
                   *flashPtr,
                   assocList,
                   fBinWidth,
                   geometry,
                   fFlashThreshold,
                   fWidthTolerance,
                   clock_data,
                   fTrigCoinc);

    // Make the associations which we noted we need
    for (size_t i = 0; i != assocList.size(); ++i) {
      art::PtrVector<recob::OpHit> opHitPtrVector;
      for (size_t const hitIndex : assocList.at(i)) {
        opHitPtrVector.emplace_back(opHitHandle, hitIndex);
      }

      util::CreateAssn(evt, *flashPtr, opHitPtrVector, *(assnPtr.get()), i);
    }

    evt.put(std::move(flashPtr));
    evt.put(std::move(assnPtr));
  }

} // namespace opdet
