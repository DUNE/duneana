////////////////////////////////////////////////////////////////////////
// Class:       GenRecoValidator
// Plugin Type: analyzer (Unknown Unknown)
// File:        GenRecoValidator_module.cc
//
// Generated at Mon Jan 29 08:43:39 2024 by Matteo Galli using cetskelgen
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
#include <TTree.h>

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimPhotons.h"

namespace dune {

  class GenRecoValidator;
}


class dune::GenRecoValidator : public art::EDAnalyzer {
public:
  explicit GenRecoValidator(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GenRecoValidator(GenRecoValidator const&) = delete;
  GenRecoValidator(GenRecoValidator&&) = delete;
  GenRecoValidator& operator=(GenRecoValidator const&) = delete;
  GenRecoValidator& operator=(GenRecoValidator&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  TTree *fTree;

  unsigned int fEventID;
  unsigned int Nstep;
  unsigned int Nmctruth;
  unsigned int Nhits;
  Int_t no_hits;
  unsigned int Nphotons;
  int nphotons;
  int no_photons;

  std::vector<float> edep;
  std::vector<float> num_photons;
  std::vector<float> num_electrons;
  std::vector<float> nuPDG_truth;
  std::vector<float> ccnc_truth;
  std::vector<float> mode_truth;
  std::vector<float> W_truth;           
  std::vector<float> X_truth;          
  std::vector<float> Y_truth;           
  std::vector<float> hitnuc_truth;      
  std::vector<float> theta_truth;       
  std::vector<float> enu_truth;         
  std::vector<float> nuvtxx_truth;      
  std::vector<float> nuvtxy_truth;      
  std::vector<float> nuvtxz_truth;      
  std::vector<float> time_truth;        
  std::vector<float> nu_dcosx_truth;    
  std::vector<float> nu_dcosy_truth;    
  std::vector<float> nu_dcosz_truth;    
  std::vector<float> zenith_truth;      
  std::vector<float> azimut_truth;      
  std::vector<float> beamangleX_truth;  
  std::vector<float> beamangleY_truth;  
  std::vector<float> beamangleZ_truth;  
  std::vector<float> lep_mom_truth;     
  std::vector<float> lep_dcosx_truth;   
  std::vector<float> lep_dcosy_truth;   
  std::vector<float> lep_dcosz_truth;   
  std::vector<float> nuWeight_truth;    
  std::vector<float> genie_no_primaries;
  std::vector<float> hit_channel;      
  std::vector<float> hit_tpc;     
  std::vector<float> hit_plane;        
  std::vector<float> hit_wire;         
  std::vector<float> hit_peakT;        
  std::vector<float> hit_charge;       
  std::vector<float> hit_ph;           
  std::vector<float> hit_startT;       
  std::vector<float> hit_endT;         
  std::vector<float> hit_width;        
  std::vector<float> hit_rms;          
  std::vector<float> hit_goodnessOfFit;
  std::vector<float> hit_multiplicity; 

  std::string fIonAndScintLabel;
  std::string fGenieGenModuleLabel;
  std::string fHitsModuleLabel;
  std::string fPhotonsLiteModuleLabel;
  bool fSaveIonAndScintInfo; 
  bool fSaveGenieGenInfo;
  bool fSaveHitsInfo;
  bool fSavePhotonsInfo;

};


dune::GenRecoValidator::GenRecoValidator(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fIonAndScintLabel(p.get<std::string>("IonAndScintLabel")),
  fGenieGenModuleLabel(p.get<std::string>("GenieGenModuleLabel")),
  fHitsModuleLabel(p.get<std::string>("HitsModuleLabel")),
  fPhotonsLiteModuleLabel(p.get<std::string>("PhotonsLiteModuleLabel")),
  fSaveIonAndScintInfo(p.get<bool>("SaveIonAndScintInfo")),
  fSaveGenieGenInfo(p.get<bool>("SaveGenieGenInfo")),
  fSaveHitsInfo(p.get<bool>("SaveHitsInfo")),
  fSavePhotonsInfo(p.get<bool>("SavePhotonsInfo"))

  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void dune::GenRecoValidator::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // Set the event ID
  fEventID = e.id().event();

  // Set all counters to 0 for the current event
  Nstep = 0;
  Nmctruth = 0;
  Nhits = 0;
  no_hits = 0;
  Nphotons = 0;
  nphotons = 0;
  no_photons = 0;

  // Get energy deposit
  if (fSaveIonAndScintInfo){
    art::ValidHandle<std::vector<sim::SimEnergyDeposit>> simHandle = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fIonAndScintLabel);
    std::vector<art::Ptr<sim::SimEnergyDeposit>> edeplist;
  
    if (simHandle.isValid())
      art::fill_ptr_vector(edeplist, simHandle);
  
    Nstep = edeplist.size();
    edep.clear();
    num_photons.clear();
    num_electrons.clear();
  
    edep.resize(Nstep);
    num_photons.resize(Nstep);
    num_electrons.resize(Nstep);
    
    for (size_t i=0; i<Nstep; i++)
    {
      edep[i]          = edeplist[i]->Energy();
      num_photons[i]   = edeplist[i]->NumPhotons();
      num_electrons[i] = edeplist[i]->NumElectrons();
    }
  } // end fSaveIonAndScintInfo

  // Get MC truth information
  if (fSaveGenieGenInfo){
    art::ValidHandle<std::vector<simb::MCTruth>> mctruthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fGenieGenModuleLabel);
    std::vector<art::Ptr<simb::MCTruth>> mclist;
  
    if (mctruthHandle.isValid())
      art::fill_ptr_vector(mclist, mctruthHandle);
  
    Nmctruth = mclist.size();
    nuPDG_truth.clear();
    ccnc_truth.clear();
    mode_truth.clear();
    W_truth.clear(); 
    X_truth.clear(); 
    Y_truth.clear(); 
    hitnuc_truth.clear(); 
    theta_truth.clear(); 
    enu_truth.clear(); 
    nuvtxx_truth.clear(); 
    nuvtxy_truth.clear(); 
    nuvtxz_truth.clear(); 
    time_truth.clear(); 
    nu_dcosx_truth.clear(); 
    nu_dcosy_truth.clear(); 
    nu_dcosz_truth.clear(); 
    zenith_truth.clear(); 
    azimut_truth.clear(); 
    beamangleX_truth.clear(); 
    beamangleY_truth.clear(); 
    beamangleZ_truth.clear(); 
    lep_mom_truth.clear(); 
    lep_dcosx_truth.clear(); 
    lep_dcosy_truth.clear(); 
    lep_dcosz_truth.clear(); 
    nuWeight_truth.clear(); 
    genie_no_primaries.clear();
  
  
    nuPDG_truth.resize(Nmctruth);
    ccnc_truth.resize(Nmctruth);
    mode_truth.resize(Nmctruth);
    W_truth.resize(Nmctruth); 
    X_truth.resize(Nmctruth); 
    Y_truth.resize(Nmctruth); 
    hitnuc_truth.resize(Nmctruth); 
    theta_truth.resize(Nmctruth); 
    enu_truth.resize(Nmctruth); 
    nuvtxx_truth.resize(Nmctruth); 
    nuvtxy_truth.resize(Nmctruth); 
    nuvtxz_truth.resize(Nmctruth); 
    time_truth.resize(Nmctruth); 
    nu_dcosx_truth.resize(Nmctruth); 
    nu_dcosy_truth.resize(Nmctruth); 
    nu_dcosz_truth.resize(Nmctruth); 
    zenith_truth.resize(Nmctruth); 
    azimut_truth.resize(Nmctruth); 
    beamangleX_truth.resize(Nmctruth); 
    beamangleY_truth.resize(Nmctruth); 
    beamangleZ_truth.resize(Nmctruth); 
    lep_mom_truth.resize(Nmctruth); 
    lep_dcosx_truth.resize(Nmctruth); 
    lep_dcosy_truth.resize(Nmctruth); 
    lep_dcosz_truth.resize(Nmctruth); 
    nuWeight_truth.resize(Nmctruth); 
    genie_no_primaries.resize(Nmctruth);
    
  
    for (size_t i=0; i<Nmctruth; i++)
    {
      nuPDG_truth[i]        = mclist[i]->GetNeutrino().Nu().PdgCode();
      ccnc_truth[i]         = mclist[i]->GetNeutrino().CCNC();
      mode_truth[i]         = mclist[i]->GetNeutrino().Mode();
      W_truth[i]            = mclist[i]->GetNeutrino().W();
      X_truth[i]            = mclist[i]->GetNeutrino().X();
      Y_truth[i]            = mclist[i]->GetNeutrino().Y();
      hitnuc_truth[i]       = mclist[i]->GetNeutrino().HitNuc();
  	  theta_truth[i]        = mclist[i]->GetNeutrino().Theta();
      enu_truth[i]          = mclist[i]->GetNeutrino().Nu().E();
      nuvtxx_truth[i]       = mclist[i]->GetNeutrino().Nu().Vx();
      nuvtxy_truth[i]       = mclist[i]->GetNeutrino().Nu().Vy();
      nuvtxz_truth[i]       = mclist[i]->GetNeutrino().Nu().Vz();
  	  time_truth[i]         = mclist[i]->GetNeutrino().Nu().T();
      nu_dcosx_truth[i]     = mclist[i]->GetNeutrino().Nu().Px()/mclist[i]->GetNeutrino().Nu().P();
      nu_dcosy_truth[i]     = mclist[i]->GetNeutrino().Nu().Py()/mclist[i]->GetNeutrino().Nu().P();
      nu_dcosz_truth[i]     = mclist[i]->GetNeutrino().Nu().Pz()/mclist[i]->GetNeutrino().Nu().P();
  	  zenith_truth[i]       = TMath::ACos(mclist[i]->GetNeutrino().Nu().Pz()/mclist[i]->GetNeutrino().Nu().P());
  	  azimut_truth[i]       = TMath::ATan((mclist[i]->GetNeutrino().Nu().Py()/mclist[i]->GetNeutrino().Nu().P())/(mclist[i]->GetNeutrino().Nu().Pz()/mclist[i]->GetNeutrino().Nu().P()));
  	  beamangleX_truth[i]   = TMath::ACos(mclist[i]->GetNeutrino().Nu().Px()/mclist[i]->GetNeutrino().Nu().P());
  	  beamangleY_truth[i]   = TMath::ACos(mclist[i]->GetNeutrino().Nu().Py()/mclist[i]->GetNeutrino().Nu().P());
  	  beamangleZ_truth[i]   = TMath::ACos(mclist[i]->GetNeutrino().Nu().Pz()/mclist[i]->GetNeutrino().Nu().P());
      lep_mom_truth[i]      = mclist[i]->GetNeutrino().Lepton().P();
      lep_dcosx_truth[i]    = mclist[i]->GetNeutrino().Lepton().Px()/mclist[i]->GetNeutrino().Lepton().P();
      lep_dcosy_truth[i]    = mclist[i]->GetNeutrino().Lepton().Py()/mclist[i]->GetNeutrino().Lepton().P();
      lep_dcosz_truth[i]    = mclist[i]->GetNeutrino().Lepton().Pz()/mclist[i]->GetNeutrino().Lepton().P();
      genie_no_primaries[i] = mclist[i]->NParticles();
      
      auto gt = e.getValidHandle< std::vector<simb::GTruth> >("generator");
  
      if ( gt ){
        auto gtruth = (*gt)[0];
        nuWeight_truth[i]  = gtruth.fweight;
      }
    }
  } // end fSaveGenieGenInfo

  // Get Hits information
  if (fSaveHitsInfo){
    art::ValidHandle<std::vector<recob::Hit>> hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);
    std::vector<art::Ptr<recob::Hit>> hitlist;
  
    if (hitHandle.isValid())
      art::fill_ptr_vector(hitlist, hitHandle);
  
    Nhits = hitlist.size();
    hit_channel.clear();
    hit_tpc.clear();
    hit_plane.clear();
    hit_wire.clear();
    hit_peakT.clear();
    hit_charge.clear();
    hit_ph.clear();
    hit_startT.clear();
    hit_endT.clear();
    hit_width.clear();
    hit_rms.clear();
    hit_goodnessOfFit.clear();
    hit_multiplicity.clear();
  
    hit_channel.resize(Nhits);
    hit_tpc.resize(Nhits);
    hit_plane.resize(Nhits);
    hit_wire.resize(Nhits);
    hit_peakT.resize(Nhits);
    hit_charge.resize(Nhits);
    hit_ph.resize(Nhits);
    hit_startT.resize(Nhits);
    hit_endT.resize(Nhits);
    hit_width.resize(Nhits);
    hit_rms.resize(Nhits);
    hit_goodnessOfFit.resize(Nhits);
    hit_multiplicity.resize(Nhits);
  
    no_hits = (int) Nhits;
  
    for (size_t i=0; i<Nhits; i++)
    {
      hit_channel[i]       = hitlist[i]->Channel();
      hit_tpc[i]           = hitlist[i]->WireID().TPC;
      hit_plane[i]         = hitlist[i]->WireID().Plane;
      hit_wire[i]          = hitlist[i]->WireID().Wire;
      hit_peakT[i]         = hitlist[i]->PeakTime();
      hit_charge[i]        = hitlist[i]->Integral();
      hit_ph[i]            = hitlist[i]->PeakAmplitude();
      hit_startT[i]        = hitlist[i]->StartTick();
      hit_endT[i]          = hitlist[i]->EndTick();
      hit_width[i]         = hitlist[i]->EndTick() - hitlist[i]->StartTick();
      hit_rms[i]           = hitlist[i]->RMS();
      hit_goodnessOfFit[i] = hitlist[i]->GoodnessOfFit();
      hit_multiplicity[i]  = hitlist[i]->Multiplicity();
    }
  } // end fSaveHitsInfo

  // Get Photons Lite information
  if (fSavePhotonsInfo){
    art::ValidHandle<std::vector<sim::SimPhotonsLite>> pliteHandle = e.getValidHandle<std::vector<sim::SimPhotonsLite>>(fPhotonsLiteModuleLabel);
    std::vector<art::Ptr<sim::SimPhotonsLite>> plitelist;
  
    if (pliteHandle.isValid())
      art::fill_ptr_vector(plitelist, pliteHandle);
  
    Nphotons = plitelist.size();
  
    for(size_t i = 0; i < Nphotons; i++){
      auto spl = pliteHandle->at(i);
      for(auto const& it : spl.DetectedPhotons){
        nphotons += it.second;
      }
    }
  
    no_photons = nphotons;
  } // end SavePhotonsInfo

  // Fill Tree
  fTree->Fill();
}

void dune::GenRecoValidator::beginJob()
{
  // Get TFileService to create the output TTree for us
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Output TTree");
  
  // Add branches to TTree
  fTree->Branch("eventID",      &fEventID);
  fTree->Branch("nstep",        &Nstep);
  fTree->Branch("e_dep",         &edep);
  fTree->Branch("num_photons",   &num_photons);
  fTree->Branch("num_electrons", &num_electrons);
  fTree->Branch("nuPDG_truth",    &nuPDG_truth);
  fTree->Branch("ccnc_truth",   &ccnc_truth);
  fTree->Branch("mode_truth",   &mode_truth);
  fTree->Branch("W_truth",            &W_truth           );
  fTree->Branch("X_truth",            &X_truth           );
  fTree->Branch("Y_truth",            &Y_truth           );
  fTree->Branch("hitnuc_truth",       &hitnuc_truth      );
  fTree->Branch("theta_truth",        &theta_truth       );
  fTree->Branch("enu_truth",          &enu_truth         );
  fTree->Branch("nuvtxx_truth",       &nuvtxx_truth      );
  fTree->Branch("nuvtxy_truth",       &nuvtxy_truth      );
  fTree->Branch("nuvtxz_truth",       &nuvtxz_truth      );
  fTree->Branch("time_truth",         &time_truth        );
  fTree->Branch("nu_dcosx_truth",     &nu_dcosx_truth    );
  fTree->Branch("nu_dcosy_truth",     &nu_dcosy_truth    );
  fTree->Branch("nu_dcosz_truth",     &nu_dcosz_truth    );
  fTree->Branch("zenith_truth",       &zenith_truth      );
  fTree->Branch("azimut_truth",       &azimut_truth      );
  fTree->Branch("beamangleX_truth",   &beamangleX_truth  );
  fTree->Branch("beamangleY_truth",   &beamangleY_truth  );
  fTree->Branch("beamangleZ_truth",   &beamangleZ_truth  );
  fTree->Branch("lep_mom_truth",      &lep_mom_truth     );
  fTree->Branch("lep_dcosx_truth",    &lep_dcosx_truth   );
  fTree->Branch("lep_dcosy_truth",    &lep_dcosy_truth   );
  fTree->Branch("lep_dcosz_truth",    &lep_dcosz_truth   );
  fTree->Branch("nuWeight_truth",     &nuWeight_truth    );
  fTree->Branch("genie_no_primaries", &genie_no_primaries);
  fTree->Branch("no_hits",           &no_hits           );
  fTree->Branch("hit_channel",       &hit_channel       );
  fTree->Branch("hit_tpc",           &hit_tpc           );
  fTree->Branch("hit_plane",         &hit_plane         );
  fTree->Branch("hit_wire",          &hit_wire          );
  fTree->Branch("hit_peakT",         &hit_peakT         );
  fTree->Branch("hit_charge",        &hit_charge        );
  fTree->Branch("hit_ph",            &hit_ph            );
  fTree->Branch("hit_startT",        &hit_startT        );
  fTree->Branch("hit_endT",          &hit_endT          );
  fTree->Branch("hit_width",         &hit_width         );
  fTree->Branch("hit_rms",           &hit_rms           );
  fTree->Branch("hit_goodnessOfFit", &hit_goodnessOfFit );
  fTree->Branch("hit_multiplicity",  &hit_multiplicity  );
  fTree->Branch("no_photons",        &no_photons        );
}

void dune::GenRecoValidator::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(dune::GenRecoValidator)
