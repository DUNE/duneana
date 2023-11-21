////////////////////////////////////////////////////////////////////////
// Class:       TrackingEfficiency
// Module Type: analyzer
// File:        TrackingEfficiency_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
//
//  Using AnaTree as base, have written a module to compare reconstruction
//  algorithm efficiencies.
//
//  Thomas Karl Warburton
//  k.warburton@sheffield.ac.uk
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TLorentzVector.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TNtuple.h"
#include "TFile.h"

//standard library includes
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>
#include <limits> // std::numeric_limits<>

namespace TrackingEfficiency {
  class TrackingEfficiency;
}

class TrackingEfficiency::TrackingEfficiency : public art::EDAnalyzer {
public:
  explicit TrackingEfficiency(fhicl::ParameterSet const & p);
  virtual ~TrackingEfficiency();

  void analyze(art::Event const & e) override;

  void beginRun(art::Run const& run) override;
  void beginJob() override;
  void endJob() override;
  void endRun(art::Run const&) override;
  //void reconfigure(fhicl::ParameterSet const & p) ;
  
  // ----------- Declare my structs ----------
  struct EfficHists {
    //EfficHists();
    // EfficHists(const std::string& subdir);
    TH1D *Length;
    TH1D *Energy;
    TH1D *Energy_Depos;
    TH1D *Theta;
    TH1D *Theta_XZ;
    TH1D *Theta_YZ;
    TH1D *Phi;
    TH2D *Phi_v_Theta;
    TH1D *Eta_XY;
    TH1D *Eta_ZY;
  } MCAllEffic, MatchedAllEffic, MCProtonEffic, MatchedProtonEffic;
  struct TEfficiencies {
    //TEfficiencies();
    //TEfficiencies(const std::string& subdir);
    TEfficiency *Effic_Length;
    TEfficiency *Effic_Energy;
    TEfficiency *Effic_Energy_Depos;
    TEfficiency *Effic_Theta;
    TEfficiency *Effic_Theta_XZ;
    TEfficiency *Effic_Theta_YZ;
    TEfficiency *Effic_Phi;
    TEfficiency *Effic_Phi_v_Theta;
    TEfficiency *Effic_Eta_XY;
    TEfficiency *Effic_Eta_ZY;
  } AllEfficiencies, ProtonEfficiencies;
  // ----------- Declare my structs ----------

private:
  
  void ResetVars();
  void MCTruthInformation ( detinfo::DetectorClocksData const& clockData,
                            const simb::MCParticle *particle, double &Energy, double &EnergyDeposited, double &TPCLength,
			    double &Theta_XZ, double &Theta_YZ, double &Eta_XY, double &Eta_ZY, double &Theta, double &Phi,
			    int &MCPdgCode, int &MCTrackId );
  void HistoFiller ( struct TrackingEfficiency::TrackingEfficiency::EfficHists *HistPtr, 
		     double Energy, double EnergyDeposited, double TPCLength, 
		     double Theta_XZ, double Theta_YZ, double Eta_XY, double Eta_ZY, double Theta, double Phi );
  void Make_Efficiencies ( struct TrackingEfficiency::TrackingEfficiency::TEfficiencies *AllHistPtr, bool matched, 
			   double Energy, double EnergyDeposited, double TPCLength, 
			   double Theta_XZ, double Theta_YZ, double Eta_XY, double Eta_ZY, double Theta, double Phi );
  void EfficPlot_1D ( TH1D *AllHist, TH1D *MatchHist, TEfficiency *EfficHist );
  void EfficPlot_2D ( TH2D *AllHist, TH2D *MatchHist, TEfficiency *EfficHist );

  TTree* fTree;

  double MCTruthT0, MCTruthTickT0, MCTruthTrackID;
  double PhotonCounterT0, PhotonCounterTickT0, PhotonCounterID;
  double TrackTheta_XZ, TrackTheta_YZ, TrackEta_XY, TrackEta_ZY, TrackTheta, TrackPhi;
  double TrackLength; 

  double MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi;
  double MCTPCLength, MCEnergy, MCEnergyDeposited;
  int MCPdgCode, MCTrackId;
  int MatchedTrackID;
  
  bool MCMatched, AllChargedTrack;
  // Handles
  art::ServiceHandle<geo::Geometry> geom;  
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  double XDriftVelocity;
  double WindowSize;
  // Parameter List
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fMCTruthT0ModuleLabel;
  std::string fPhotonT0ModuleLabel;

 };
// ********************************** Begin Run *******************************************************
void TrackingEfficiency::TrackingEfficiency::beginRun(art::Run const& run) {

}
// *********************************** Begin Job ********************************************************
void TrackingEfficiency::TrackingEfficiency::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("TrackingTree","analysis tree");
  fTree->Branch("MCTPCLength"      ,&MCTPCLength      ,"MCTPCLength/D"      );
  fTree->Branch("MCEnergy"         ,&MCEnergy         ,"MCEnergy/D"         );
  fTree->Branch("MCEnergyDeposited",&MCEnergyDeposited,"MCEnergyDeposited/D");
  fTree->Branch("MCTheta"          ,&MCTheta          ,"MCTheta/D"          );
  fTree->Branch("MCTheta_XZ"       ,&MCTheta_XZ       ,"MCTheta_XZ/D"       );
  fTree->Branch("MCTheta_YZ"       ,&MCTheta_YZ       ,"MCTheta/D"          );
  fTree->Branch("MCPhi"            ,&MCPhi            ,"MCPhi/D"            );
  fTree->Branch("MCEta_XY"         ,&MCEta_XY         ,"MCEta_XY/D"         );
  fTree->Branch("MCEta_ZY"         ,&MCEta_ZY         ,"MCEta_ZY/D"         );
  fTree->Branch("MCMatched"        ,&MCMatched        ,"MCMatched/B"        );
  fTree->Branch("MCPdgCode"        ,&MCPdgCode        ,"MCPdgCode/I"        );
  fTree->Branch("MCTrackId"        ,&MCTrackId        ,"MCTrackId/I"        );
  fTree->Branch("MatchedTrackID"   ,&MatchedTrackID   ,"MatchedTrackID/I"   );


  
  MCAllEffic.Length       = tfs->make<TH1D>("Length_MCAll"      ,"All charged Monte Carlo particle track lengths;Length (cm);Efficiency", 50, 0, 500 );
  MCAllEffic.Energy       = tfs->make<TH1D>("Energy_MCAll"      ,"All charged Monte Carlo particle energies;Energy (GeV);Efficiency"    , 50, 0, 20  );
  MCAllEffic.Energy_Depos = tfs->make<TH1D>("Energy_Depos_MCAll","All charged Monte Carlo particle energy deposited;Energy Deposited (GeV);Efficiency",50, 0, 3.5  );
  MCAllEffic.Theta        = tfs->make<TH1D>("Theta_MCAll"       ,"All charged Monte Carlo particle theta;Theta (rad);Efficiency",    50, 0, TMath::Pi() );
  MCAllEffic.Theta_XZ     = tfs->make<TH1D>("Theta_XZ_MCAll"    ,"All charged Monte Carlo particle thetaXZ;ThetaxZ (rad);Efficiency",50, -TMath::Pi(), TMath::Pi() ); 
  MCAllEffic.Theta_YZ     = tfs->make<TH1D>("Theta_YZ_MCAll"    ,"All charged Monte Carlo particle thetaYZ;ThetaYZ (rad);Efficiency",50, -TMath::Pi(), 0 ); 
  MCAllEffic.Phi          = tfs->make<TH1D>("Phi_MCAll"         ,"All charged Monte Carlo particle phi;Phi (rad);Efficiency",        50, -TMath::Pi(), 0 ); 
  MCAllEffic.Eta_XY       = tfs->make<TH1D>("Eta_XY_MCAll"      ,"All charged Monte Carlo particle etaXZ;EtaxY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  MCAllEffic.Eta_ZY       = tfs->make<TH1D>("Eta_ZY_MCAll"      ,"All charged Monte Carlo particle etaYZ;EtaZY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  MCAllEffic.Phi_v_Theta  = tfs->make<TH2D>("Phi_v_Theta_MCAll" ,"All charged Monte Carlo phi versus theta;Theta (rad);Phi (rad)",   50, 0, TMath::Pi(), 50, -TMath::Pi(), 0);
  
  MatchedAllEffic.Length       = tfs->make<TH1D>("Length_MatchAll"      ,"Matched all charged Monte Carlo particle track lengths;Length (cm);Efficiency",50, 0, 500 );
  MatchedAllEffic.Energy       = tfs->make<TH1D>("Energy_MatchAll"      ,"Matched all charged Monte Carlo particle energies;Energy (GeV);Efficiency"    ,50, 0, 20  );
  MatchedAllEffic.Energy_Depos = tfs->make<TH1D>("Energy_Depos_MatchAll","Matched all charged Monte Carlo particle energy deposited;Energy Deposited (GeV);Efficiency",50, 0, 3.5 );
  MatchedAllEffic.Theta        = tfs->make<TH1D>("Theta_MatchAll"       ,"Matched all charged Monte Carlo particle theta;Theta (rad);Efficiency",    50, 0, TMath::Pi() );
  MatchedAllEffic.Theta_XZ     = tfs->make<TH1D>("Theta_XZ_MatchAll"    ,"Matched all charged Monte Carlo particle thetaXZ;ThetaxZ (rad);Efficiency",50, -TMath::Pi(), TMath::Pi() ); 
  MatchedAllEffic.Theta_YZ     = tfs->make<TH1D>("Theta_YZ_MatchAll"    ,"Matched all charged Monte Carlo particle thetaYZ;ThetaYZ (rad);Efficiency",50, -TMath::Pi(), 0 ); 
  MatchedAllEffic.Phi          = tfs->make<TH1D>("Phi_MatchAll"         ,"Matched all charged Monte Carlo particle phi;Phi (rad);Efficiency",        50, -TMath::Pi(), 0 ); 
  MatchedAllEffic.Eta_XY       = tfs->make<TH1D>("Eta_XY_MatchAll"      ,"Matched all charged Monte Carlo particle etaXZ;EtaxY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  MatchedAllEffic.Eta_ZY       = tfs->make<TH1D>("Eta_ZY_MatchAll"      ,"Matched all charged Monte Carlo particle etaYZ;EtaZY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  MatchedAllEffic.Phi_v_Theta  = tfs->make<TH2D>("Phi_v_Theta_MatchAll" ,"Matched all charged Monte Carlo phi versus theta;Theta (rad);Phi (rad)",   50, 0, TMath::Pi(), 50, -TMath::Pi(), 0);

  AllEfficiencies.Effic_Length       = tfs->make<TEfficiency>("Effic_Length_All"      ,"Efficiency vs Particle Length for all charged particles;Length (cm);Efficiency",   50, 0, 500 );
  AllEfficiencies.Effic_Energy       = tfs->make<TEfficiency>("Effic_Energy_All"      ,"Efficiency vs Particle Energy for all charged particles;Energy (GeV);Efficiency",  50, 0, 20  );
  AllEfficiencies.Effic_Energy_Depos = tfs->make<TEfficiency>("Effic_Energy_Depos_All","Efficiency vs Particle Energy Deposited for all charged particles;Energy Deposited (GeV);Efficiency",50, 0, 3.5 );
  AllEfficiencies.Effic_Theta        = tfs->make<TEfficiency>("Effic_Theta_All"       ,"Efficiency vs Particle Theta for all charged particles;Theta (rad);Efficiency",    50, 0, TMath::Pi()  );
  AllEfficiencies.Effic_Theta_XZ     = tfs->make<TEfficiency>("Effic_Theta_XZ_All"    ,"Efficiency vs Particle ThetaXZ for all charged particles;ThetaxZ (rad);Efficiency",50, -TMath::Pi(), TMath::Pi() );
  AllEfficiencies.Effic_Theta_YZ     = tfs->make<TEfficiency>("Effic_Theta_YZ_All"    ,"Efficiency vs Particle ThetaYZ for all charged particles;ThetaYZ (rad);Efficiency",50, -TMath::Pi(), 0 );
  AllEfficiencies.Effic_Phi          = tfs->make<TEfficiency>("Effic_Phi_All"         ,"Efficiency vs Particle Phi for all charged particles;Phi (rad);Efficiency",        50, -TMath::Pi(), 0 );
  AllEfficiencies.Effic_Eta_XY       = tfs->make<TEfficiency>("Effic_Eta_XY_All"      ,"Efficiency vs Particle EtaXZ for all charged particles;EtaxY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  AllEfficiencies.Effic_Eta_ZY       = tfs->make<TEfficiency>("Effic_Eta_ZY_All"      ,"Efficiency vs Particle EtaYZ for all charged particles;EtaZY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  AllEfficiencies.Effic_Phi_v_Theta  = tfs->make<TEfficiency>("Effic_Phi_v_Theta_All" ,"Efficiency of Phi versus theta for all charged particles;Theta (rad);Phi (rad)",   50, 0, TMath::Pi(), 50, -TMath::Pi(), 0);
   
  // --------------------- Now for Protons ---------------------------
  MCProtonEffic.Length       = tfs->make<TH1D>("Length_MCProton"      ,"Proton charged Monte Carlo particle track lengths;Length (cm);Efficiency", 50, 0, 500 );
  MCProtonEffic.Energy       = tfs->make<TH1D>("Energy_MCProton"      ,"Proton charged Monte Carlo particle energies;Energy (GeV);Efficiency"    , 50, 0, 20  );
  MCProtonEffic.Energy_Depos = tfs->make<TH1D>("Energy_Depos_MCProton","Proton charged Monte Carlo particle energy deposited;Energy Deposited (GeV);Efficiency",50, 0, 3.5 );
  MCProtonEffic.Theta        = tfs->make<TH1D>("Theta_MCProton"       ,"Proton charged Monte Carlo particle theta;Theta (rad);Efficiency",    50, 0, TMath::Pi() );
  MCProtonEffic.Theta_XZ     = tfs->make<TH1D>("Theta_XZ_MCProton"    ,"Proton charged Monte Carlo particle thetaXZ;ThetaxZ (rad);Efficiency",50, -TMath::Pi(), TMath::Pi() ); 
  MCProtonEffic.Theta_YZ     = tfs->make<TH1D>("Theta_YZ_MCProton"    ,"Proton charged Monte Carlo particle thetaYZ;ThetaYZ (rad);Efficiency",50, -TMath::Pi(), 0 ); 
  MCProtonEffic.Phi          = tfs->make<TH1D>("Phi_MCProton"         ,"Proton charged Monte Carlo particle phi;Phi (rad);Efficiency",        50, -TMath::Pi(), 0 ); 
  MCProtonEffic.Eta_XY       = tfs->make<TH1D>("Eta_XY_MCProton"      ,"Proton charged Monte Carlo particle etaXZ;EtaxY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  MCProtonEffic.Eta_ZY       = tfs->make<TH1D>("Eta_ZY_MCProton"      ,"Proton charged Monte Carlo particle etaYZ;EtaZY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  MCProtonEffic.Phi_v_Theta  = tfs->make<TH2D>("Phi_v_Theta_MCProton" ,"Proton charged Monte Carlo phi versus theta;Theta (rad);Phi (rad)",   50, 0, TMath::Pi(), 50, -TMath::Pi(), 0);
  
  MatchedProtonEffic.Length       = tfs->make<TH1D>("Length_MatchProton"      ,"Matched all charged Monte Carlo particle track lengths;Length (cm);Efficiency",50, 0, 500 );
  MatchedProtonEffic.Energy       = tfs->make<TH1D>("Energy_MatchProton"      ,"Matched all charged Monte Carlo particle energies;Energy (GeV);Efficiency"    ,50, 0, 20  );
  MatchedProtonEffic.Energy_Depos = tfs->make<TH1D>("Energy_Depos_MatchProton","Matched all charged Monte Carlo particle energy deposited;Energy Deposited (GeV);Efficiency",50, 0, 3.5 );
  MatchedProtonEffic.Theta        = tfs->make<TH1D>("Theta_MatchProton"       ,"Matched all charged Monte Carlo particle theta;Theta (rad);Efficiency",    50, 0, TMath::Pi() );
  MatchedProtonEffic.Theta_XZ     = tfs->make<TH1D>("Theta_XZ_MatchProton"    ,"Matched all charged Monte Carlo particle thetaXZ;ThetaxZ (rad);Efficiency",50, -TMath::Pi(), TMath::Pi() ); 
  MatchedProtonEffic.Theta_YZ     = tfs->make<TH1D>("Theta_YZ_MatchProton"    ,"Matched all charged Monte Carlo particle thetaYZ;ThetaYZ (rad);Efficiency",50, -TMath::Pi(), 0 ); 
  MatchedProtonEffic.Phi          = tfs->make<TH1D>("Phi_MatchProton"         ,"Matched all charged Monte Carlo particle phi;Phi (rad);Efficiency",        50, -TMath::Pi(), 0 ); 
  MatchedProtonEffic.Eta_XY       = tfs->make<TH1D>("Eta_XY_MatchProton"      ,"Matched all charged Monte Carlo particle etaXZ;EtaxY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  MatchedProtonEffic.Eta_ZY       = tfs->make<TH1D>("Eta_ZY_MatchProton"      ,"Matched all charged Monte Carlo particle etaYZ;EtaZY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  MatchedProtonEffic.Phi_v_Theta  = tfs->make<TH2D>("Phi_v_Theta_MatchProton" ,"Matched all charged Monte Carlo phi versus theta;Theta (rad);Phi (rad)",   50, 0, TMath::Pi(), 50, -TMath::Pi(), 0);

  ProtonEfficiencies.Effic_Length       = tfs->make<TEfficiency>("Effic_Length_Proton"      ,"Efficiency vs Particle Length for all charged particles;Length (cm);Efficiency",   50, 0, 500 );
  ProtonEfficiencies.Effic_Energy       = tfs->make<TEfficiency>("Effic_Energy_Proton"      ,"Efficiency vs Particle Energy for all charged particles;Energy (GeV);Efficiency",  50, 0, 20  );
  ProtonEfficiencies.Effic_Energy_Depos = tfs->make<TEfficiency>("Effic_Energy_Depos_Proton","Efficiency vs Particle Energy Deposited for all charged particles;Energy Deposited (GeV);Efficiency",50, 0, 3.5 );
  ProtonEfficiencies.Effic_Theta        = tfs->make<TEfficiency>("Effic_Theta_Proton"       ,"Efficiency vs Particle Theta for all charged particles;Theta (rad);Efficiency",    50, 0, TMath::Pi() );
  ProtonEfficiencies.Effic_Theta_XZ     = tfs->make<TEfficiency>("Effic_Theta_XZ_Proton"    ,"Efficiency vs Particle ThetaXZ for all charged particles;ThetaxZ (rad);Efficiency",50, -TMath::Pi(), TMath::Pi() );
  ProtonEfficiencies.Effic_Theta_YZ     = tfs->make<TEfficiency>("Effic_Theta_YZ_Proton"    ,"Efficiency vs Particle ThetaYZ for all charged particles;ThetaYZ (rad);Efficiency",50, -TMath::Pi(), 0 );
  ProtonEfficiencies.Effic_Phi          = tfs->make<TEfficiency>("Effic_Phi_Proton"         ,"Efficiency vs Particle Phi for all charged particles;Phi (rad);Efficiency",        50, -TMath::Pi(), 0 );
  ProtonEfficiencies.Effic_Eta_XY       = tfs->make<TEfficiency>("Effic_Eta_XY_Proton"      ,"Efficiency vs Particle EtaXZ for all charged particles;EtaxY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  ProtonEfficiencies.Effic_Eta_ZY       = tfs->make<TEfficiency>("Effic_Eta_ZY_Proton"      ,"Efficiency vs Particle EtaYZ for all charged particles;EtaZY (rad);Efficiency",    50, -TMath::Pi(), TMath::Pi() );
  ProtonEfficiencies.Effic_Phi_v_Theta  = tfs->make<TEfficiency>("Effic_Phi_v_Theta_Proton" ,"Efficiency of Phi versus theta for all charged particles;Theta (rad);Phi (rad)",   50, 0, TMath::Pi(), 50, -TMath::Pi(), 0);
  
}
// ************************************ End Job *********************************************************
void TrackingEfficiency::TrackingEfficiency::endJob() {
}
// ************************************ End Run *********************************************************
void TrackingEfficiency::TrackingEfficiency::endRun(art::Run const&) {
}

// ********************************** pset param *******************************************************
TrackingEfficiency::TrackingEfficiency::TrackingEfficiency(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fHitsModuleLabel         ( pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel        ( pset.get< std::string >("TrackModuleLabel"))
  , fMCTruthT0ModuleLabel    ( pset.get< std::string >("MCTruthT0ModuleLabel"))
  , fPhotonT0ModuleLabel     ( pset.get< std::string >("PhotonT0ModuleLabel"))
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
  XDriftVelocity      = detProp.DriftVelocity()*1e-3; //cm/ns
  WindowSize          = detProp.NumberTimeSamples() * clockData.TPCClock().TickPeriod() * 1e3;
}
// ******************************************************************************************************
TrackingEfficiency::TrackingEfficiency::~TrackingEfficiency()
{
  // Clean up dynamic memory and other resources here.

}
// ************************************ Analyse *********************************************************
void TrackingEfficiency::TrackingEfficiency::analyze(art::Event const & evt)
{
  //std::cout << "\n\n************* New Event / Module running *************\n\n" << std::endl;
  // Implementation of required member function here. 
  std::vector<art::Ptr<recob::Track> > tracklist;
  auto trackListHandle = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  if (trackListHandle)
    art::fill_ptr_vector(tracklist, trackListHandle);
  
  const sim::ParticleList& plist = pi_serv->ParticleList();
  
  auto trackh = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  
  auto trackvh = evt.getHandle< std::vector< art::PtrVector < recob::Track > > >(fTrackModuleLabel);
  
  int NPart = 0;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
  if ( trackListHandle.isValid() ) { // Check that trackListHandle is valid.....
    art::FindManyP<recob::Hit>        fmht  (trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::T0>           fmt0  (trackListHandle, evt, fMCTruthT0ModuleLabel);
    art::FindMany<anab::T0>           fmphot(trackListHandle, evt, fPhotonT0ModuleLabel);
    int ntracks_reco=tracklist.size();      
  
    // **************************************************************
    // Want to loop through Monte Carlo Truth Particles. 
    // **************************************************************
    for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
      simb::MCParticle *particle = ipar->second;
      
      MCMatched       = false; // MCParticle not matched by defualt    
      AllChargedTrack = false; // MCParticle not charged by default
      MatchedTrackID  = -1;    // No reco track has this ID

      // ----- Check that have a charged particle -----
      if ( fabs(particle->PdgCode() ) == 13 )  AllChargedTrack = true;
      else if ( particle->PdgCode() == 2212 ) AllChargedTrack = true;
      /*
	else if ( fabs(particle->PdgCode() == 11) ) AllChargedTrack = true;  // Electron
	else if ( fabs(particle->PdgCode() == 211) ) AllChargedTrack = true; // Pions
	else if ( fabs(particle->PdgCode() == 321) ) AllChargedTrack = true; // Kaons
      //*/
      if (!AllChargedTrack) continue;
    
      // ---- Get MCTruth Information and check that MCParticle goes in TPC ---
      MCTruthInformation ( clockData, particle, MCEnergy, MCEnergyDeposited, MCTPCLength,
			   MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi,
			   MCPdgCode, MCTrackId );
      if ( MCTPCLength < 1. ) continue;
      ++NPart;
      //std::cout << "Looking at MCParticle " << MCTrackId << ", with PdgCode " << MCPdgCode << ", MCTPCLength " << MCTPCLength << std::endl;
      
      // ----- Fill the histograms for all the Monte Carlo Particles -----
      HistoFiller ( &MCAllEffic, MCEnergy, MCEnergyDeposited, MCTPCLength, 
		    MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi );
      if ( MCPdgCode == 2212 ) {
	HistoFiller ( &MCProtonEffic, MCEnergy, MCEnergyDeposited, MCTPCLength, 
		      MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi );
      } // Proton
          
      for(int Track=0; Track < ntracks_reco; ++Track){
	if ( MCMatched ) continue;
	// *****************************************************************************
	// T0 stuff - So can correct X Start/End positions and identify MCParticle!!
	// *****************************************************************************
	double TickT0 = 0;
	if ( fmt0.isValid() ) {
	  std::vector<const anab::T0*> T0s = fmt0.at(Track);
	  for (size_t t0size =0; t0size < T0s.size(); t0size++) {
	    MCTruthT0      = T0s[t0size]->Time();
	    MCTruthTickT0  = MCTruthT0 / sampling_rate(clockData);
	    MCTruthTrackID = T0s[t0size]->TriggerBits();
	  } // T0 size
	  TickT0 = MCTruthTickT0;
	} // T0 valid
	if ( fmphot.isValid() ) {
	  std::vector<const anab::T0*> PhotT0 = fmphot.at(Track);
	  for (size_t T0it=0; T0it<PhotT0.size(); ++T0it) {
	    PhotonCounterT0     = PhotT0[T0it]->Time();
	    PhotonCounterTickT0 = PhotonCounterT0 / sampling_rate(clockData);
	    PhotonCounterID     = PhotT0[T0it]->TriggerBits();
	  }
	  //TickT0 = PhotonCounterTickT0;
	} 
	// ************** END T0 stuff ***************
	
	// ***** Check whether MCParticle made this track *****
	if ( MCTrackId != MCTruthTrackID ) continue;
	      
	// ******************************************************************************************
	// Correct X and get track length etc now that we have matched a Track with an MCParticle!!
	// ******************************************************************************************
	// ----- Hit Stuff ------
	std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(Track);
	double Hit_Size = allHits.size();
	// ---- Correct X positions and get track positions.
        recob::Track::Point_t trackStart, trackEnd;
        std::tie(trackStart, trackEnd) = tracklist[Track]->Extent(); 
	trackStart.SetX( trackStart.X() - detProp.ConvertTicksToX( TickT0, allHits[Hit_Size-1]->WireID().Plane, allHits[Hit_Size-1]->WireID().TPC, allHits[Hit_Size-1]->WireID().Cryostat )); // Correct X, last entry is first 'hit'
        trackEnd.SetX( trackEnd.X() - detProp.ConvertTicksToX( TickT0, allHits[0]->WireID().Plane, allHits[0]->WireID().TPC, allHits[0]->WireID().Cryostat)); // Correct X, first entry is last 'hit'
	// ---- Get lengths and angles.
	art::Ptr<recob::Track> ptrack(trackh, Track);
	const recob::Track& track = *ptrack;
	TrackLength = track.Length();
	int ntraj = track.NumberTrajectoryPoints();
	if(ntraj > 0) {
	  auto dir = track.VertexDirection();
	  TrackTheta_XZ = std::atan2(dir.X(), dir.Z());
	  TrackTheta_YZ = std::atan2(dir.Y(), dir.Z());
	  TrackEta_XY   = std::atan2(dir.X(), dir.Y());
	  TrackEta_ZY   = std::atan2(dir.Z(), dir.Y());
	  TrackTheta    = dir.Theta();
	  TrackPhi      = dir.Phi();
	}
	// *********** Corrected X and track information **************
      
	// ***************************** Check how well reconstrcuted ******************************
	//std::cout << "Track " << tracklist[Track]->ID() << " had reco Length " << TrackLength << ", MCLength " << MCTPCLength << " ---> " << TrackLength / MCTPCLength << std::endl;
	if ( TrackLength / MCTPCLength > 0.5 &&
	     TrackLength / MCTPCLength < 1.5 &&
	     MCMatched != true && // If first time this MCParticle passed criteria
	     1 ) {

	  MCMatched = true; // Stop MCParticle passing criteria with multiple tracks
	  MatchedTrackID = tracklist[Track]->ID(); // Which trackid matched this MCParticle?
	  //std::cout << "Particle matched so fill numerator and denominator\n" << std::endl;
	  // --- Track is well reconstructed, so fill the all charged histograms and TEfficiencies!!! ----
	  HistoFiller ( &MatchedAllEffic, MCEnergy, MCEnergyDeposited, MCTPCLength, 
			MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi );
	  Make_Efficiencies ( &AllEfficiencies, true, 
			      MCEnergy, MCEnergyDeposited, MCTPCLength, 
			      MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi);
	  	
	  // --- If is a proton fill that efficiency plot too ----
	  if ( MCPdgCode == 2212 ) {
	    HistoFiller ( &MatchedProtonEffic, MCEnergy, MCEnergyDeposited, MCTPCLength, 
			  MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi );
	    Make_Efficiencies ( &ProtonEfficiencies, true, 
				MCEnergy, MCEnergyDeposited, MCTPCLength, 
				MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi);
	  } // Proton
	
	} // If well matched track
	// ***************************** Check how well reconstrcuted ******************************
      
      } // Track Loop
      if (!MCMatched) { // If charged MCParticle doesn't get matched want to fill just denominator.
	//std::cout << "MCParticle not matched so only fill denominator...\n" << std::endl;
	Make_Efficiencies ( &AllEfficiencies, false, 
			    MCEnergy, MCEnergyDeposited, MCTPCLength, 
			    MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi);
	// --- If a proton, also want to fill that one too! ---
	if ( MCPdgCode == 2212 ) {
	  Make_Efficiencies ( &ProtonEfficiencies, false, 
			      MCEnergy, MCEnergyDeposited, MCTPCLength, 
			      MCTheta_XZ, MCTheta_YZ, MCEta_XY, MCEta_ZY, MCTheta, MCPhi);
	} // Proton
      } // If not matched
      // ******** Fill Tree for each MCParticle **********
      fTree->Fill();
    } // Loop over MCParticles
  } // if trackListHandle.isValid()
  //std::cout << "\nThis event had " << NPart << " particles." << std::endl;
}
// *********************************** Monte Carlo Truth Extraction ********************************************************
void TrackingEfficiency::TrackingEfficiency::MCTruthInformation (detinfo::DetectorClocksData const& clockData,
                                                                 const simb::MCParticle *particle, double &Energy, double &EnergyDeposited, double &TPCLength,
								  double &Theta_XZ, double &Theta_YZ, double &Eta_XY, double &Eta_ZY, double &Theta, double &Phi,
								  int &MCPdgCode, int &MCTrackId ) {
  int numberTrajectoryPoints = particle->NumberTrajectoryPoints(); // Looking at each MC hit
  //double TPCLengthHits[numberTrajectoryPoints];
  std::vector<double> TPCLengthHits(numberTrajectoryPoints, 0);
  bool BeenInVolume = false;
  int FirstHit=0, LastHit=0;
  TPCLength = 0;

  MCPdgCode       = particle->PdgCode();
  MCTrackId       = particle->TrackId();
  
  for(unsigned int MCHit=0; MCHit < TPCLengthHits.size(); ++MCHit) {
    if (MCHit!=0) TPCLengthHits[MCHit] = pow ( pow( (particle->Vx(MCHit-1)-particle->Vx(MCHit)),2)
					       + pow( (particle->Vy(MCHit-1)-particle->Vy(MCHit)),2)
					       + pow( (particle->Vz(MCHit-1)-particle->Vz(MCHit)),2)
					       , 0.5 );
    // --- Check if hit is in TPC...
    auto const position = geo::vect::toPoint(particle->Position(MCHit).Vect());
    geo::TPCID tpcid = geom->FindTPCAtPosition(position);
    if (tpcid.isValid) { 
      // -- Check if hit is within drift window...
      geo::TPCGeo      const& tpc  = geom->TPC(tpcid);
      double XPlanePosition      = tpc.Plane(0).GetCenter().X();
      double DriftTimeCorrection = fabs( position.X() - XPlanePosition ) / XDriftVelocity;
      double TimeAtPlane         = particle->T() + DriftTimeCorrection;
      if ( TimeAtPlane < trigger_offset(clockData)
           || TimeAtPlane > trigger_offset(clockData) + WindowSize
	   ) continue;
      // -- Good hit in TPC
      LastHit = MCHit;
      if ( !BeenInVolume ) {
	BeenInVolume = true;
	FirstHit = MCHit;
      }
    } // TPC.valid	    
    /*
    if (MCTrackId == 54 ) {
      std::cout << MCHit << " of " << numberTrajectoryPoints << " x,y,z " << tmpPosition[0] << " " << tmpPosition[1] << " " << tmpPosition[2] << " " << TPCLengthHits[MCHit] << ", t " << particle->T() << ". In volume? " << BeenInVolume << " first/last hit " << FirstHit << " " << LastHit << std::endl;
    }
    //*/
  } // MCTrajPoints
  
  Energy          = particle->E(FirstHit);
  EnergyDeposited = particle->E(FirstHit) - particle->E(LastHit);
  for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) {
    TPCLength += TPCLengthHits[Hit];
    //if (MCTrackId == 54 ) std::cout << FirstHit << " " << LastHit << ", " << Hit << ". dx " << TPCLengthHits[Hit] << " Len " << TPCLength << std::endl;
  }
    TLorentzVector& momentumStart  = (TLorentzVector&)particle->Momentum(FirstHit);   
  TVector3 mcstartmom = particle->Momentum(FirstHit).Vect();
  Theta_XZ = std::atan2(momentumStart.Px(), momentumStart.Pz());
  Theta_YZ = std::atan2(momentumStart.Py(), momentumStart.Pz());
  Eta_XY   = std::atan2(momentumStart.Px(), momentumStart.Py());
  Eta_ZY   = std::atan2(momentumStart.Pz(), momentumStart.Py());
  Theta    = mcstartmom.Theta();
  Phi      = mcstartmom.Phi();
  
} // MCTruthInformation

// ******************************** Efficiency Filler ****************************************************
void TrackingEfficiency::TrackingEfficiency::HistoFiller ( struct TrackingEfficiency::TrackingEfficiency::EfficHists *HistPtr, 
							   double Energy, double EnergyDeposited, double TPCLength, 
							   double Theta_XZ, double Theta_YZ, double Eta_XY, double Eta_ZY, double Theta, double Phi ) {
  HistPtr->Length       -> Fill( TPCLength       );
  HistPtr->Energy       -> Fill( Energy          );
  HistPtr->Energy_Depos -> Fill( EnergyDeposited );
  HistPtr->Theta        -> Fill( Theta           );
  HistPtr->Theta_XZ     -> Fill( Theta_XZ        );
  HistPtr->Theta_YZ     -> Fill( Theta_YZ        );
  HistPtr->Phi          -> Fill( Phi             );
  HistPtr->Phi_v_Theta  -> Fill( Theta, Phi      );
  HistPtr->Eta_XY       -> Fill( Eta_XY          );
  HistPtr->Eta_ZY       -> Fill( Eta_ZY          );
} // HistoFiller

// ******************************** Make Efficiencies ****************************************************
void TrackingEfficiency::TrackingEfficiency::Make_Efficiencies ( struct TrackingEfficiency::TrackingEfficiency::TEfficiencies *AllHistPtr, bool matched, 
								 double Energy, double EnergyDeposited, double TPCLength, 
								 double Theta_XZ, double Theta_YZ, double Eta_XY, double Eta_ZY, double Theta, double Phi ) {
  AllHistPtr->Effic_Length      -> Fill ( matched, TPCLength       );
  AllHistPtr->Effic_Energy      -> Fill ( matched, Energy          );
  AllHistPtr->Effic_Energy_Depos-> Fill ( matched, EnergyDeposited );
  AllHistPtr->Effic_Theta       -> Fill ( matched, Theta           );
  AllHistPtr->Effic_Theta_XZ    -> Fill ( matched, Theta_XZ        );
  AllHistPtr->Effic_Theta_YZ    -> Fill ( matched, Theta_YZ        );
  AllHistPtr->Effic_Phi         -> Fill ( matched, Phi             );
  AllHistPtr->Effic_Eta_XY      -> Fill ( matched, Eta_XY          );
  AllHistPtr->Effic_Eta_ZY      -> Fill ( matched, Eta_ZY          );

  AllHistPtr->Effic_Phi_v_Theta -> Fill ( matched, Theta, Phi      );
}
// ******************************** Reset Variables *****************************************************
void TrackingEfficiency::TrackingEfficiency::ResetVars() {
}
// ******************************** Define Module *****************************************************
DEFINE_ART_MODULE(TrackingEfficiency::TrackingEfficiency)
