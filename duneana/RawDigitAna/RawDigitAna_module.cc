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


constexpr int kMaxNumberCh = 41471;
constexpr int kMaxTicks= 8000;

namespace rawdigit {
  class RawDigitAna;
}


class rawdigit::RawDigitAna : public art::EDAnalyzer {
public:
  explicit RawDigitAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
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
  TTree *fTree;
  int fRun;
  int fSubrun;
  int fEvent;
  // int fW_plane[kMaxNumberCh];
  int fW_ch[kMaxNumberCh];
  float fW_signal[kMaxNumberCh][kMaxTicks];

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
  fTree->Branch("event", &fEvent);
  fTree->Branch("subrun", &fSubrun);
  fTree->Branch("run", &fRun);
  // fTree->Branch("w_plane", &fW_plane, Form("w_plane[%d]/I", kMaxNumberCh));
  fTree->Branch("w_ch", &fW_ch, Form("w_ch[%d]/I", kMaxNumberCh));
  fTree->Branch("w_signal", &fW_signal, Form("w_signal[%d][%d]/F", kMaxNumberCh, kMaxTicks));

}  

void rawdigit::RawDigitAna::endJob(){

}
void rawdigit::RawDigitAna::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  auto const& rawdigits = *(e.getValidHandle<std::vector<raw::RawDigit>>(fModuleLabel));
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.id().event();

  unsigned int idx = 0;
  for (raw::RawDigit const& rawdigit: rawdigits) {
    // std::cout << "Compression: " << rawdigit.Compression() << std::endl;
    if (idx >= kMaxNumberCh) {
      mf::LogWarning("RawDigitAna") << "Number of channels in rawdigit exceeds kMaxNumberCh. Only saving the first " << kMaxNumberCh << " channels.";
      break;
    }
    unsigned int jdx = 0;
    for ( int adc = 0; adc < int(rawdigit.NADC()); ++adc){
      if (jdx >= kMaxTicks) {
        mf::LogWarning("RawDigitAna") << "Number of ticks in rawdigit exceeds kMaxTicks. Only saving the first " << kMaxTicks << " ticks.";
        break;
      }
      fW_signal[idx][jdx]= rawdigit.ADC(adc);
      ++ jdx;
    }
    // fW_plane[idx] = rawdigit.View();
    fW_ch[idx] = rawdigit.Channel();
    ++ idx;
  }
  fTree->Fill();
}

void rawdigit::RawDigitAna::reset(){

  for ( unsigned int i=0; i<kMaxNumberCh; ++i){
    // fW_plane[i] = -999;   
    fW_ch[i] = -1;   
    for ( unsigned int j=0; j<kMaxTicks; ++j){
      fW_signal[i][j] = 0;
    }
  }
}

DEFINE_ART_MODULE(rawdigit::RawDigitAna)