////////////////////////////////////////////////////////////////////////
//
//  Very simple module to filter events without all of the expected
//  products in the event.
//
//  Leigh Whitehead - leigh.howard.whitehead@cern.ch
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/OpFlash.h"

namespace mpf{
  class MissingProductFilter;
}

class mpf::MissingProductFilter : public art::EDFilter{
public:
  
  explicit MissingProductFilter(fhicl::ParameterSet const& pset);
  virtual ~MissingProductFilter();
  
  void beginJob() override;
  bool filter(art::Event& evt) override;
  void endJob() override;

private:

  std::vector<std::string> fModules;
};
  
//-----------------------------------------------------------------------
mpf::MissingProductFilter::MissingProductFilter(fhicl::ParameterSet const& pset):
  EDFilter(pset)
{
  fModules = pset.get<std::vector<std::string>>("Modules");
}

//-----------------------------------------------------------------------
mpf::MissingProductFilter::~MissingProductFilter(){}

//-----------------------------------------------------------------------
void mpf::MissingProductFilter::beginJob() {}

//-----------------------------------------------------------------------
bool mpf::MissingProductFilter::filter(art::Event& evt){

  bool prodMissing = false;

  // Can we do it all together?
  for(const std::string &s : fModules){

    try{
      if(s == "reco3d"){
        evt.getValidHandle<std::vector<recob::SpacePoint>>(s);
      }
      else if(s == "pandora"){
        evt.getValidHandle<std::vector<recob::PFParticle>>(s);
      }
      else if(s == "pandoraTrack"){
        evt.getValidHandle<std::vector<recob::Track>>(s);
      }
      else if(s == "pandoraShower"){
        evt.getValidHandle<std::vector<recob::Shower>>(s);
      }
      else if(s == "opflash" || s == "opflashInternal" || s == "opflashExternal"){
        evt.getValidHandle<std::vector<recob::OpFlash>>(s);
      }
      else{
        std::cerr << "Product " << s << " is not considered by the module and will be ignored." << std::endl;
      }
    }
    catch(...){
      prodMissing = true;
      std::cout << "Product " << s << " is missing. Event will be filtered out." << std::endl;
      break;
    }

  }

  if(!prodMissing){
    std::cout << "All products present, keeping this event" << std::endl;
  }
  return !prodMissing;

}

void mpf::MissingProductFilter::endJob() {}
 
DEFINE_ART_MODULE(mpf::MissingProductFilter)
