////////////////////////////////////////////////////////////////////////
// Class:      ROIAna 
// Plugin Type: analyzer (art v3_00_00)
// File:        ROIAna_module.cc
// Written by Tejin Cai & Matthew Man
// Reach out for questions/issues/bugs
////////////////////////////////////////////////////////////////////////
#ifndef WIREANA_H
#define WIREANA_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
//#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "lardataobj/RawData/RDTimeStamp.h"


#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"

#include "art_root_io/TFileService.h"

// ROOT includes
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"

//Others
#define DEFAULT_VALUE -99999

namespace roiana {

  class ROIAna;

}

class roiana::ROIAna : public art::EDAnalyzer {
public:
  explicit ROIAna(fhicl::ParameterSet const& pset);
  ROIAna(ROIAna const&) = delete;
  ROIAna(ROIAna&&) = delete;
  ROIAna& operator=(ROIAna const&) = delete;
  ROIAna& operator=(ROIAna&&) = delete;
  virtual ~ROIAna() noexcept {};

  /////////////////////////////////////////////
  // Required functions.
  void analyze(art::Event const& evt) override;

  /////////////////////////////////////////////
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  /////////////////////////////////////////////
  std::map<int,bool> ProcessROIWide(std::vector<art::Ptr<raw::RawDigit>>& rawList);

  bool ispeak(float x) { return (x != 0); }
  bool isZero(float x) { return x == 0 ; }
private:
  ////////////////////////////////////////////
  // Labels
  const art::InputTag fWireProducerLabel; 
  const art::InputTag fRawProducerLabel;
  const art::InputTag fSimChannelLabel;
  const art::InputTag fSimulationProducerLabel;
  ////////////////////////////////////////////
  // Log Control
  int fLogLevel;


  /////////////////////////////////////////////
  // Geometry Options && Tool options
  int fNPlanes;
  int fNChanPerApa;
  int fNTicksPerWire;

  std::vector<std::string> fLabels, fInteraction;

  int fROI_Peak;
  int fROI_CH;

  std::map<raw::ChannelID_t, std::pair<art::Ptr<raw::RawDigit>, art::Ptr<sim::SimChannel>>> 
    ch_w_sc;
  std::map<int, std::string> 
    trkid_to_label_map;

  /////////////////////////////////////////////
  // Tree Variables

  int run;
  int subrun;
  int event;
  int MC;
  int flag;

  // --- MC ROI Tree Variables
  int PDGROITree;
  int TrackIDROITree;
  std::string Generator;
  float MCParticleEnergy;
  float TEnergyDepositedU;
  float TChargeDepositedU;
  float TEnergyDepositedV;
  float TChargeDepositedV;
  float TEnergyDepositedX;
  float TChargeDepositedX;
  float TEnergyDepositedROIU;
  float TChargeDepositedROIU;
  float TEnergyDepositedROIV;
  float TChargeDepositedROIV;
  float TEnergyDepositedROIX;
  float TChargeDepositedROIX;
  std::map<int,float> TrackIDEnergyMapU;
  std::map<int,float> TrackIDChargeMapU;
  std::map<int,float> TrackIDEnergyMapV;
  std::map<int,float> TrackIDChargeMapV;
  std::map<int,float> TrackIDEnergyMapX;
  std::map<int,float> TrackIDChargeMapX;
  std::map<int,float> TrackIDEnergyMapROIU;
  std::map<int,float> TrackIDChargeMapROIU;
  std::map<int,float> TrackIDEnergyMapROIV;
  std::map<int,float> TrackIDChargeMapROIV;
  std::map<int,float> TrackIDEnergyMapROIX;
  std::map<int,float> TrackIDChargeMapROIX;

  // --- MC Truth Tree Variables
  std::vector<int> TPart;
  std::vector<std::map<int,simb::MCParticle>> Parts = {};
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  // --- MC Interaction Tree Variables
  std::string Interaction;
  int PDG;
  float Energy;
  //std::vector<int> DaughterPDG;
  std::vector<double> Momentum,StartVertex,EndVertex;
  //std::vector<double> DaughterE,DaughterPx,DaughterPy,DaughterPz,DaughterStartVx,DaughterStartVy,DaughterStartVz,DaughterEndVx,DaughterEndVy,DaughterEndVz;
  
  /////////////////////////////////////////////
  // Private Functions
  template<class T>
  void SortWirePtrByChannel( std::vector<art::Ptr<T>> &vec, bool increasing );

  
  /////////////////////////////////////////////
  // ROI Metrics
  art::ServiceHandle<cheat::ParticleInventoryService> PIS;
  void ROIMetrics( std::map<int,bool> );

  /////////////////////////////////////////////
  // ROI Filter
  void ROIFilter( std::vector<art::Ptr<raw::RawDigit>>& rawList, std::map<int,bool> ret );

  /////////////////////////////////////////////
  // ROI Efficiencies
  void ROIEfficiencies( std::map<int,bool>, int , std::map<int,sim::IDE>);

  // Helper functions
  std::string PrintInColor ( std::string InputString, std::string MyString, int MyColor );
  void FillMyMaps    ( std::map< int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand );
  void FillMCInteractionTree(
    std::map< int, simb::MCParticle> &MCParticleList,
    std::vector<std::string> ProcessList,
    int fLogLevel);
  bool InMyMap(
    int TrID,
    std::map< int,
    simb::MCParticle> ParMap);
  int GetColor( std::string ColorName );

  /////////////////////////////////////////////
  // Declare output data
  //TTree *fTree;
  std::string fTreeName;
  TTree* fROITree;
  TTree* fMCTruthTree;
  TTree* fInteractionTree;

  TH1F* TrueEnergyDeposited;
  TH1F* TrueEnergyDepositedInROI; //Filled in TagROITruth
  TH1F* TrueEnergyDepositedRatio;
  TH1F* TrueChargeDeposited;
  TH1F* TrueChargeDepositedInROI; //Filled in TagROITruth
  TH1F* TrueChargeDepositedRatio;
  TH1F* DataReductionRate;
  TH1F* MarleySignalSensitivity;
  float fECMin; //minimum energy and charge to begin accumulation
  float fHistEnergyMax;
  float fHistChargeMax;
};





#endif
