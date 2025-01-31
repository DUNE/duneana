////////////////////////////////////////////////////////////////////////
// Class:       WireAnaTree
// Plugin Type: analyzer (Unknown Unknown)
// File:        WireAnaTree_module.cc
//
// Generated at Fri Mar 29 14:02:45 2024 by Aaron Higuera Pichardo using cetskelgen
// from  version .
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
#include "lardataobj/RecoBase/Wire.h"


constexpr int kMaxNumberCh = 41471;
constexpr int kMaxTicks= 8000;

namespace wire {
  class WireAnaTree;
}


class wire::WireAnaTree : public art::EDAnalyzer {
public:
  explicit WireAnaTree(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireAnaTree(WireAnaTree const&) = delete;
  WireAnaTree(WireAnaTree&&) = delete;
  WireAnaTree& operator=(WireAnaTree const&) = delete;
  WireAnaTree& operator=(WireAnaTree&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;
  void reset();

private:

  art::InputTag fCalWireModuleLabel;
  TTree *fTree;
  int fRun;
  int fSubrun;
  int fEvent;
  int fW_plane[kMaxNumberCh];
  int fW_ch[kMaxNumberCh];
  float fW_signal[kMaxNumberCh][kMaxTicks];

};


wire::WireAnaTree::WireAnaTree(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fCalWireModuleLabel(p.get< art::InputTag >("CalWireModuleLabel"))
{
}

void wire::WireAnaTree::beginJob(){

  reset();
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("wireTree","Tree with wire info");
  fTree->Branch("event", &fEvent);
  fTree->Branch("subrun", &fSubrun);
  fTree->Branch("run", &fRun);
  fTree->Branch("w_plane", &fW_plane, Form("w_plane[%d]/I", kMaxNumberCh));
  fTree->Branch("w_ch", &fW_ch, Form("w_ch[%d]/I", kMaxNumberCh));
  fTree->Branch("w_signal", &fW_signal, Form("w_signal[%d][%d]/F", kMaxNumberCh, kMaxTicks));

}  

void wire::WireAnaTree::endJob(){

}
void wire::WireAnaTree::analyze(art::Event const& e)
{
  // Print a hello message.
  mf::LogInfo("WireAnaTree") << "Hello, WireAnaTree! ";
  // Implementation of required member function here.
  auto const& Wires = *(e.getValidHandle<std::vector<recob::Wire>>(fCalWireModuleLabel));
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.id().event();

  // Print general info about event and wires
  mf::LogInfo("WireAnaTree") << "Run: " << fRun << " Subrun: " << fSubrun << " Event: " << fEvent;
  mf::LogInfo("WireAnaTree") << "Number of wires: " << Wires.size();
  unsigned int idx =0;
  for (recob::Wire const& wire: Wires) {
    unsigned int jdx =0;
    if (idx >= kMaxNumberCh) break;
    for ( const float &signal : wire.Signal() ){
      if (jdx >= kMaxTicks) break;
      fW_signal[idx][jdx]= signal;
      ++ jdx;
    }
    fW_plane[idx] = wire.View();
    fW_ch[idx] = wire.Channel();
    ++ idx;
  }
  fTree->Fill();
}

void wire::WireAnaTree::reset(){

  for ( unsigned int i=0; i<kMaxNumberCh; ++i){
    fW_plane[i] = -999;   
    fW_ch[i] = -999;   
    for ( unsigned int j=0; j<kMaxTicks; ++j){
      fW_signal[i][j] = -999.0;
    }
  }
}
DEFINE_ART_MODULE(wire::WireAnaTree)
