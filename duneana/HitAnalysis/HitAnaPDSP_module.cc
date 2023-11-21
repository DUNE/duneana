////////////////////////////////////////////////////////////////////////
// Class:       HitAnaPDSP
// Plugin Type: analyzer (art v3_02_06)
// File:        HitAnaPDSP_module.cc
//
// Generated at Tue Aug 20 16:19:40 2019 by Tingjun Yang using cetskelgen
// from cetlib version v3_07_02.
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

#include "lardataobj/RecoBase/Hit.h"

#include "TTree.h"

namespace pdsp {
  class HitAnaPDSP;
}

class pdsp::HitAnaPDSP : public art::EDAnalyzer {
public:
  explicit HitAnaPDSP(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HitAnaPDSP(HitAnaPDSP const&) = delete;
  HitAnaPDSP(HitAnaPDSP&&) = delete;
  HitAnaPDSP& operator=(HitAnaPDSP const&) = delete;
  HitAnaPDSP& operator=(HitAnaPDSP&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  art::InputTag fHitModuleLabel;
  
  TTree *ftree;
  int run;
  int subrun;
  int event;
  std::vector<short> channel;
  std::vector<short> tpc;
  std::vector<short> plane;
  std::vector<short> wire;
  std::vector<double> charge;
  std::vector<double> peakt;
  std::vector<double> rms;
  std::vector<double> startt;
  std::vector<double> endt;
  std::vector<short> origin;
  std::vector<int> pdg;
};


pdsp::HitAnaPDSP::HitAnaPDSP(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  ,
  fHitModuleLabel(p.get< art::InputTag >("HitModuleLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdsp::HitAnaPDSP::analyze(art::Event const& e)
{
  run = e.run();
  subrun = e.subRun();
  event = e.id().event();
  channel.clear();
  tpc.clear();
  plane.clear();
  wire.clear();
  charge.clear();
  peakt.clear();
  rms.clear();
  startt.clear();
  endt.clear();
  origin.clear();
  pdg.clear();

   // Reconstruciton information
  auto hitListHandle = e.getHandle < std::vector < recob::Hit > >(fHitModuleLabel);
  if (!hitListHandle) {
    return;
  }

  for (auto const & hit : *hitListHandle){
    channel.push_back(hit.Channel());
    tpc.push_back(hit.WireID().TPC);
    plane.push_back(hit.WireID().Plane);
    wire.push_back(hit.WireID().Wire);
    charge.push_back(hit.Integral());
    peakt.push_back(hit.PeakTime());
    rms.push_back(hit.RMS());
    startt.push_back(hit.StartTick());
    endt.push_back(hit.EndTick());
  }

  if (!channel.empty()) ftree->Fill();

}

void pdsp::HitAnaPDSP::beginJob()
{
  art::ServiceHandle<art::TFileService> fileServiceHandle;
  ftree = fileServiceHandle->make<TTree>("ftree", "hit info");
  ftree->Branch("run", &run, "run/I");
  ftree->Branch("event", &event, "event/I");
  ftree->Branch("channel", &channel);
  ftree->Branch("tpc", &tpc);
  ftree->Branch("plane", &plane);
  ftree->Branch("wire", &wire);
  ftree->Branch("charge", &charge);
  ftree->Branch("peakt", &peakt);
  ftree->Branch("rms", &rms);
  ftree->Branch("startt", &startt);
  ftree->Branch("endt", &endt);
  ftree->Branch("origin", &origin);
  ftree->Branch("pdf", &pdg);
}

DEFINE_ART_MODULE(pdsp::HitAnaPDSP)
