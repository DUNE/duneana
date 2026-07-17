////////////////////////////////////////////////////////////////////////
// Class:       RawDigitAna
// Plugin Type: analyzer
// File:        RawDigitAna_module.cc
//
// Generated at Wed Jun 12 by Sergio Manthey Corchado
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "TTree.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace rawdigit {
  class RawDigitAna;
}


class rawdigit::RawDigitAna : public art::EDAnalyzer {
public:
  explicit RawDigitAna(fhicl::ParameterSet const& p);

  RawDigitAna(RawDigitAna const&) = delete;
  RawDigitAna(RawDigitAna&&) = delete;
  RawDigitAna& operator=(RawDigitAna const&) = delete;
  RawDigitAna& operator=(RawDigitAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;
  void reset();

private:

  art::InputTag fModuleLabel;
  TTree* fTree;
  int fRun, fSubrun, fEvent;
  int fView, fChannel;
  std::vector<short> fADC;
  std::vector<int> fTPCs;

};


rawdigit::RawDigitAna::RawDigitAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fModuleLabel(p.get< art::InputTag >("ModuleLabel"))
{
}

void rawdigit::RawDigitAna::beginJob(){

  reset();
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("rawdigitTree","Tree with rawdigit info");
  fTree->Branch("run",     &fRun);
  fTree->Branch("subrun",  &fSubrun);
  fTree->Branch("event",   &fEvent);
  fTree->Branch("tpc",     &fTPCs);
  fTree->Branch("channel", &fChannel);
  fTree->Branch("view",    &fView);
  fTree->Branch("adc",     &fADC);

}  

void rawdigit::RawDigitAna::endJob(){

}
void rawdigit::RawDigitAna::analyze(art::Event const& e)
{
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.id().event();
  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  auto const& rawdigits = *(e.getValidHandle<std::vector<raw::RawDigit>>(fModuleLabel));
  for (raw::RawDigit const& rd : rawdigits) {
    fChannel = rd.Channel();
    fView = static_cast<int>(wireReadout.View(rd.Channel()));
    fTPCs.clear();
    for (auto const& wid : wireReadout.ChannelToWire(rd.Channel())) {
      fTPCs.push_back(wid.TPC);
    }

    fADC.assign(rd.ADCs().begin(), rd.ADCs().end());
    fTree->Fill();
  }  
}
void rawdigit::RawDigitAna::reset(){
  fChannel = -999;
  fView    = -999;
  fTPCs.clear();
  fADC.clear();
}

DEFINE_ART_MODULE(rawdigit::RawDigitAna)
