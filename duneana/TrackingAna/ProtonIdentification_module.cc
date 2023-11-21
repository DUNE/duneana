////////////////////////////////////////////////////////////////////////
// Class:       ProtonIdentification
// Module Type: analyzer
// File:        ProtonIdentification_module.cc
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
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

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

const int kMaxHits = 250;

namespace ProtonIdentification {
  class ProtonIdentification;
}

class ProtonIdentification::ProtonIdentification : public art::EDAnalyzer {
public:
  explicit ProtonIdentification(fhicl::ParameterSet const & p);
  virtual ~ProtonIdentification();

  void analyze(art::Event const & e) override;

  void beginRun(art::Run const& run) override;
  void beginJob() override;
  void endJob() override;
  void endRun(art::Run const&) override;
  //void reconfigure(fhicl::ParameterSet const & p) ;

private:
  
  void ResetVars();
  void MCTruthInformation ( detinfo::DetectorClocksData const& clockData, const simb::MCParticle *particle );

  void  TrackBoundaries( TVector3 larStart, TVector3 larEnd );
  double CalcDist( double X1, double Y1, double Z1, double X2, double Y2, double Z2 );
  double CalcPIDA( std::vector<const anab::Calorimetry*> calos );
  double PIDAFunc( double dedx, double resrng );

  // The T0 variables
  double MCTruthT0, MCTruthTickT0, MCTruthTrackID;
  double PhotonCounterT0, PhotonCounterTickT0, PhotonCounterID;
  // Some angular variables
  double TrackTheta_XZ, TrackTheta_YZ, TrackEta_XY, TrackEta_ZY, TrackTheta, TrackPhi;
  double MCTheta_XZ   , MCTheta_YZ   , MCEta_XY   , MCEta_ZY   , MCTheta   , MCPhi   ;

  // The TPC boundaries
  // MinX, MaxX, MinY, MaxY, MinZ, MaxZ
  double boundaries[6]; // For use in finding detector boundaries.  
  
  // Some global variables for what I have
  int TrueMuons=0, TrueElectrons=0, TrueProtons=0, TrueGammas=0, TrueOthers=0;
  int AllBad = 0, MuonBad = 0, ProtonBad = 0, GammaBad = 0, ElectronBad=0, OtherBad = 0;
  int AllAll = 0, MuonAll = 0, ProtonAll = 0, GammaAll = 0, ElectronAll=0, OtherAll = 0;
  int AllRange = 0, MuonRange = 0, ProtonRange = 0, GammaRange = 0, ElectronRange=0, OtherRange = 0;

  // My exported TTree
  TTree *fRecoTree;
  // MC variables
  int    MCPdgCode, MCTrackId, MatchedTrackID, MCContainment, MCIsPrimary;
  double MCTPCLength, MCEnergy, MCEnergyDeposited, MCPIDA;
  bool   MCStartInTPC, MCEndInTPC, MCCorOrient;
  double MCStartX, MCStartY, MCStartZ, MCEndX, MCEndY, MCEndZ;
  // Reco variables
  double PIDA, PIDA_Plane0, PIDA_Plane1, PIDA_Plane2, StartFromEdge, EndFromEdge, AvdEdx, TrackLength;
  double CorrectedStartX, CorrectedStartY, CorrectedStartZ, CorrectedEndX, CorrectedEndY, CorrectedEndZ;
  bool   RecoStartInTPC, RecoEndInTPC;
  int    RecoContainment, TrackHits[3];
  int    CaloPlane0, CaloPlane1, CaloPlane2;
  double ResRPlane0[kMaxHits], ResRPlane1[kMaxHits], ResRPlane2[kMaxHits];
  double dEdxPlane0[kMaxHits], dEdxPlane1[kMaxHits], dEdxPlane2[kMaxHits];


  // A soley truth tree...
  TTree* fTrueTree;
  int True_Particles, True_PdgCode[1000], True_Contained[1000], True_ID[1000], True_Primary[1000];
  double True_Length[1000], True_StartX[1000], True_StartY[1000], True_StartZ[1000], True_EndX[1000], True_EndY[1000], True_EndZ[1000];

  // Handles
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  
  // Some variables
  double XDriftVelocity;
  double WindowSize;

  // Parameter List
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fMCTruthT0ModuleLabel;
  std::string fPhotonT0ModuleLabel;
  std::string fCounterT0ModuleLabel;
  std::string fCalorimetryModuleLabel;
  bool fUsePhotons;
  double PIDApower, fBoundaryEdge;
  int  Verbose;
};
// ********************************** Begin Run *******************************************************
void ProtonIdentification::ProtonIdentification::beginRun(art::Run const& run) {

}
// *********************************** Begin Job ********************************************************
void ProtonIdentification::ProtonIdentification::beginJob()
{
  // --------- Working out detector dimensions ----------------
  boundaries[0] = boundaries[2] = boundaries[4] = DBL_MAX;
  boundaries[1] = boundaries[3] = boundaries[5] = -DBL_MAX;

  auto const* geom = lar::providerFrom<geo::Geometry>();
  for (geo::TPCGeo const& TPC: geom->Iterate<geo::TPCGeo>()) {
    // get center in world coordinates
    auto const center = TPC.GetCenter();
    //double tpcDim[3] = {TPC.ActiveHalfWidth(), TPC.ActiveHalfHeight(), 0.5*TPC.ActiveLength() };
    double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
    if( center.X() - tpcDim[0] < boundaries[0] ) boundaries[0] = center.X() - tpcDim[0];
    if( center.X() + tpcDim[0] > boundaries[1] ) boundaries[1] = center.X() + tpcDim[0];
    if( center.Y() - tpcDim[1] < boundaries[2] ) boundaries[2] = center.Y() - tpcDim[1];
    if( center.Y() + tpcDim[1] > boundaries[3] ) boundaries[3] = center.Y() + tpcDim[1];
    if( center.Z() - tpcDim[2] < boundaries[4] ) boundaries[4] = center.Z() - tpcDim[2];
    if( center.Z() + tpcDim[2] > boundaries[5] ) boundaries[5] = center.Z() + tpcDim[2];
    //std::cout << boundaries[0] << " " << boundaries[1] << " " << boundaries[2] << " " << boundaries[3] << " " <<boundaries[4] << " " << boundaries[5] << std::endl;
  } // for all TPC
  //std::cout << boundaries[0] << " " << boundaries[1] << " " << boundaries[2] << " " << boundaries[3] << " " <<boundaries[4] << " " << boundaries[5] << std::endl;
  
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fRecoTree = tfs->make<TTree>("ReconstructedTree","analysis tree");  
  // The Truth information
  fRecoTree->Branch("MCPdgCode"        ,&MCPdgCode        ,"MCPdgCode/I"        );
  fRecoTree->Branch("MCTPCLength"      ,&MCTPCLength      ,"MCTPCLength/D"      );
  fRecoTree->Branch("MCTrackId"        ,&MCTrackId        ,"MCTrackId/I"        );
  fRecoTree->Branch("MCEnergy"         ,&MCEnergy         ,"MCEnergy/D"         );
  fRecoTree->Branch("MCEnergyDeposited",&MCEnergyDeposited,"MCEnergyDeposited/D");
  fRecoTree->Branch("MCStartInTPC"     ,&MCStartInTPC     ,"MCStartInTPC/B"     );
  fRecoTree->Branch("MCEndInTPC"       ,&MCEndInTPC       ,"MCEndInTPC/B"       );
  fRecoTree->Branch("MCCorOrient"      ,&MCCorOrient      ,"MCCorOrient/B"      );
  fRecoTree->Branch("MCContainment"    ,&MCContainment    ,"MCContainment/I"    );
  fRecoTree->Branch("MCIsPrimary"      ,&MCIsPrimary      ,"MCIsPrimary/I"      );
  fRecoTree->Branch("MCPIDA"           ,&MCPIDA           ,"MCPIDA/D"           );
  fRecoTree->Branch("MCStartX"         ,&MCStartX         ,"MCStartX/D"         );
  fRecoTree->Branch("MCStartY"         ,&MCStartY         ,"MCStartY/D"         );
  fRecoTree->Branch("MCStartZ"         ,&MCStartZ         ,"MCStartZ/D"         );
  fRecoTree->Branch("MCEndX"           ,&MCEndX           ,"MCEndX/D"           );
  fRecoTree->Branch("MCEndY"           ,&MCEndY           ,"MCEndY/D"           );
  fRecoTree->Branch("MCEndZ"           ,&MCEndZ           ,"MCEndZ/D"           );
  // Now for the reco stuff...
  fRecoTree->Branch("TrackLength"      ,&TrackLength      ,"TrackLength/D"      );
  fRecoTree->Branch("TrackHits"        ,&TrackHits        ,"TrackHits[3]/I"     );
  fRecoTree->Branch("MatchedTrackID"   ,&MatchedTrackID   ,"MatchedTrackID/I"   );
  fRecoTree->Branch("CorrectedStartX"  ,&CorrectedStartX  ,"CorrectedStartX/D"  );
  fRecoTree->Branch("CorrectedStartY"  ,&CorrectedStartY  ,"CorrectedStartY/D"  );
  fRecoTree->Branch("CorrectedStartZ"  ,&CorrectedStartZ  ,"CorrectedStartZ/D"  );
  fRecoTree->Branch("CorrectedEndX"    ,&CorrectedEndX    ,"CorrectedEndX/D"    );
  fRecoTree->Branch("CorrectedEndY"    ,&CorrectedEndY    ,"CorrectedEndY/D"    );
  fRecoTree->Branch("CorrectedEndZ"    ,&CorrectedEndZ    ,"CorrectedEndZ/D"    );
  fRecoTree->Branch("PIDA"             ,&PIDA             ,"PIDA/D"             );
  fRecoTree->Branch("PIDA_Plane0"      ,&PIDA_Plane0      ,"PIDA_Plane0/D"      );
  fRecoTree->Branch("PIDA_Plane1"      ,&PIDA_Plane1      ,"PIDA_Plane1/D"      );
  fRecoTree->Branch("PIDA_Plane2"      ,&PIDA_Plane2      ,"PIDA_Plane2/D"      );
  fRecoTree->Branch("RecoStartInTPC"   ,&RecoStartInTPC   ,"RecoStartInTPC/B"   );
  fRecoTree->Branch("RecoEndInTPC"     ,&RecoEndInTPC     ,"RecoEndInTPC/B"     );
  fRecoTree->Branch("StartFromEdge"    ,&StartFromEdge    ,"StartFromEdge/D"    );
  fRecoTree->Branch("EndFromEdge"      ,&EndFromEdge      ,"EndFromEdge/D"      );
  fRecoTree->Branch("RecoContainment"  ,&RecoContainment  ,"RecoContainment/I"  );
  // Calorimetry objects
  fRecoTree->Branch("CaloPlane0", &CaloPlane0, "CaloPlane0/I"            );
  fRecoTree->Branch("ResRPlane0", &ResRPlane0, "ResRPlane0[CaloPlane0]/D");
  fRecoTree->Branch("dEdxPlane0", &dEdxPlane0, "dEdxPlane0[CaloPlane0]/D");
  fRecoTree->Branch("CaloPlane1", &CaloPlane1, "CaloPlane1/I"            );
  fRecoTree->Branch("ResRPlane1", &ResRPlane1, "ResRPlane1[CaloPlane1]/D");
  fRecoTree->Branch("dEdxPlane1", &dEdxPlane1, "dEdxPlane1[CaloPlane1]/D");
  fRecoTree->Branch("CaloPlane2", &CaloPlane2, "CaloPlane2/I"            );
  fRecoTree->Branch("ResRPlane2", &ResRPlane2, "ResRPlane2[CaloPlane2]/D");
  fRecoTree->Branch("dEdxPlane2", &dEdxPlane2, "dEdxPlane2[CaloPlane2]/D");
  
  // ==== My truth tree
  fTrueTree = tfs->make<TTree>("TruthTree","truth_tree");
  fTrueTree->Branch("True_Particles", &True_Particles, "True_Particles/I" );
  fTrueTree->Branch("True_PdgCode"  , &True_PdgCode  , "True_PdgCode[True_Particles]/I"   );
  fTrueTree->Branch("True_ID"       , &True_ID       , "True_ID[True_Particles]/I"        );
  fTrueTree->Branch("True_Contained", &True_Contained, "True_Contained[True_Particles]/I" );
  fTrueTree->Branch("True_Primary"  , &True_Primary  , "True_Primary[True_Particles]/I"   );
  fTrueTree->Branch("True_Length"   , &True_Length   , "True_Length[True_Particles]/D"    );
  fTrueTree->Branch("True_StartX"   , &True_StartX   , "True_StartX[True_Particles]/D"    );
  fTrueTree->Branch("True_StartY"   , &True_StartY   , "True_StartY[True_Particles]/D"    );
  fTrueTree->Branch("True_StartZ"   , &True_StartZ   , "True_StartZ[True_Particles]/D"    );
  fTrueTree->Branch("True_EndX"     , &True_EndX     , "True_EndX[True_Particles]/D"      );
  fTrueTree->Branch("True_EndY"     , &True_EndY     , "True_EndY[True_Particles]/D"      );
  fTrueTree->Branch("True_EndZ"     , &True_EndZ     , "True_EndZ[True_Particles]/D"      );
}
// ************************************ End Job *********************************************************
void ProtonIdentification::ProtonIdentification::endJob() {
  std::cout << "\nFinished all the events for module " << fTrackModuleLabel << "..."
	    << "\nMonte Carlo's in the TPC had a total of " << TrueMuons << " Muons, " << TrueProtons << " Protons, " << TrueElectrons << " Electrons, " << TrueGammas << " Gammas, and " << TrueOthers << " others." 
	    << "\nHad a total of " << AllAll << " tracks. Comprised of Muons " << MuonAll << ", Protons " << ProtonAll << ", Gammas " << GammaAll << ", Electrons " << ElectronBad << ", Others " << OtherAll  
	    << "\n" << AllBad << " of which had high PIDA's, Muons " << MuonBad << ", Protons " << ProtonBad << ", Gammas " << GammaBad << ", Electrons " << ElectronBad << ", Others " << OtherBad
	    << "\n" << AllRange << " tracks were within PIDA range, Muons " << MuonRange << ", Protons " << ProtonRange << ", Gammas " << GammaRange << ", ElectronRange " << ElectronRange << ", Others " << OtherRange 
	    << std::endl;
}
// ************************************ End Run *********************************************************
void ProtonIdentification::ProtonIdentification::endRun(art::Run const&) {
}

// ********************************** pset param *******************************************************
ProtonIdentification::ProtonIdentification::ProtonIdentification(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fHitsModuleLabel         ( pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel        ( pset.get< std::string >("TrackModuleLabel"))
  , fMCTruthT0ModuleLabel    ( pset.get< std::string >("MCTruthT0ModuleLabel"))
  , fPhotonT0ModuleLabel     ( pset.get< std::string >("PhotonT0ModuleLabel"))
  , fCounterT0ModuleLabel    ( pset.get< std::string >("CounterT0ModuleLabel"))
  , fCalorimetryModuleLabel  ( pset.get< std::string >("CalorimetryModuleLabel"))
  , fUsePhotons              ( pset.get< bool        >("UsePhotons"))
  , PIDApower                ( pset.get< double      >("PIDApower"))
  , fBoundaryEdge            ( pset.get< double      >("BoundaryEdge"))
  , Verbose                  ( pset.get< int         >("Verbose"))
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
  XDriftVelocity      = detProp.DriftVelocity()*1e-3; //cm/ns
  WindowSize          = detProp.NumberTimeSamples() * clockData.TPCClock().TickPeriod() * 1e3;
}
// ******************************************************************************************************
ProtonIdentification::ProtonIdentification::~ProtonIdentification()
{
  // Clean up dynamic memory and other resources here.

}
// ************************************ Analyse *********************************************************
void ProtonIdentification::ProtonIdentification::analyze(art::Event const & evt)
{
  //std::cout << "\n\n************* New Event / Module running *************\n\n" << std::endl;
  // Implementation of required member function here. 
  std::vector<art::Ptr<recob::Track> > tracklist;
  auto trackListHandle = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  if (trackListHandle)
    art::fill_ptr_vector(tracklist, trackListHandle);

  auto trackh = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  
  std::vector<art::Ptr<raw::ExternalTrigger> > triglist;
  auto trigListHandle = evt.getHandle< std::vector<raw::ExternalTrigger> >(fCounterT0ModuleLabel);
  if (trigListHandle)
    art::fill_ptr_vector(triglist, trigListHandle);  
  
  const sim::ParticleList& plist = pi_serv->ParticleList();
  True_Particles = 0;
  for (int qq=0; qq<1000; ++qq) {
    True_PdgCode  [ qq ] = True_Contained[ qq ] = True_Length   [ qq ] = True_ID[qq] = True_Primary[qq] = 0;
    True_StartX   [ qq ] = True_StartY   [ qq ] = True_StartZ   [ qq ] = 0;
    True_EndX     [ qq ] = True_EndY     [ qq ] = True_EndZ     [ qq ] = 0;
  }
// Quickly get total number of MC Particles...
  for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
    // Some variables by each particle
    bool InTPC = false, StIn = false, EnIn = false;
    double ThisStX=0, ThisStY=0, ThisStZ=0, ThisEnX=0, ThisEnY=0, ThisEnZ=0;    
    double ThisLen=0;
    int    ThisPDG=0, ThisCon=-1, ThisID=-1, ThisPri=-1;
    // Get my particle and loop through points
    simb::MCParticle *particle = ipar->second;
    ThisPDG = particle->PdgCode();
    ThisID  = particle->TrackId();
    if (particle->Process() == "primary") ThisPri=1;
    else ThisPri=0;
    for ( unsigned int a=0; a<particle->NumberTrajectoryPoints(); ++a ) {
      // Get Positions and if in TPC
      auto const position = geo::vect::toPoint(particle->Position(a).Vect());
      geo::TPCID tpcid = geom->FindTPCAtPosition(position);
      if (tpcid.isValid && a == 0) 
	StIn = true;
      if (tpcid.isValid && a == (particle->NumberTrajectoryPoints() -1) )
	EnIn = true;
      // If not in TPC yet...
      if (!InTPC) {
	if (tpcid.isValid) {
	  InTPC = true;
          ThisStX = position.X();
          ThisStY = position.Y();
          ThisStZ = position.Z();
	  // Increment some total output variables....
	  if (abs(particle->PdgCode()) == 13 ) ++TrueMuons;
	  else if (abs(particle->PdgCode()) == 11 ) ++TrueElectrons;
	  else if (abs(particle->PdgCode()) == 22 ) ++TrueGammas;
	  else if (abs(particle->PdgCode()) == 2212 ) ++TrueProtons;
	  else ++TrueOthers;
	} // First time in TPC
      } else {
	// If already been in TPC
	if (tpcid.isValid) {
	  // And still in the TPC!
          ThisEnX = position.X();
          ThisEnY = position.Y();
          ThisEnZ = position.Z();
	} // Still in TPC
      } // Already in TPC
    } // Traj points
    if(!StIn && !EnIn ) ThisCon = 0;  // through track
    if(!StIn &&  EnIn ) ThisCon = 1;  // entering track
    if( StIn && !EnIn ) ThisCon = 2;  // escaping track
    if( StIn &&  EnIn ) ThisCon = 3;  // contained track
    ThisLen = CalcDist( ThisStX, ThisStY, ThisStZ, ThisEnX, ThisEnY, ThisEnZ );
    /*
    std::cout << "Finished particle info...InTPC " << InTPC << " PDG " << ThisPDG 
	      << ", start " << ThisStX << ", " << ThisStY << ", " << ThisStZ << "....End " << ThisEnX << ", " << ThisEnY << ", " << ThisEnZ
	      << "....StIn " << StIn << ", EnIn " << EnIn << ", ThisCon " << ThisCon << ", Len " << ThisLen
	      << std::endl;
    */
    if (InTPC && True_Particles < 1000) {
      True_PdgCode  [ True_Particles ] = ThisPDG;
      True_ID       [ True_Particles ] = ThisID;
      True_Contained[ True_Particles ] = ThisCon;
      True_Primary  [ True_Particles ] = ThisPri;
      True_Length   [ True_Particles ] = ThisLen;
      True_StartX   [ True_Particles ] = ThisStX;
      True_StartY   [ True_Particles ] = ThisStY;
      True_StartZ   [ True_Particles ] = ThisStZ;
      True_EndX     [ True_Particles ] = ThisEnX;
      True_EndY     [ True_Particles ] = ThisEnY;
      True_EndZ     [ True_Particles ] = ThisEnZ;
      ++True_Particles;
    }
  } // Particle list
  // Now I have gone through particle list fill the tree.
  fTrueTree -> Fill();
  
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
  if ( trackListHandle.isValid() ) { // Check that trackListHandle is valid.....
    art::FindManyP<recob::Hit>        fmht   (trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::T0>           fmt0   (trackListHandle, evt, fMCTruthT0ModuleLabel);
    art::FindMany<anab::T0>           fmphot (trackListHandle, evt, fPhotonT0ModuleLabel);
    art::FindMany<anab::Calorimetry>  fmcal  (trackListHandle, evt, fCalorimetryModuleLabel);
    int ntracks_reco=tracklist.size();      
    
    for(int Track=0; Track < ntracks_reco; ++Track){
      ResetVars(); // Reset variables.
      std::cout << "\n***** Looking at new track " << Track << " of " << ntracks_reco << ", event " << evt.event() << ", Module " << fTrackModuleLabel << std::endl;

      // Load the new track info, and the basic track properties.
      art::Ptr<recob::Track> ptrack(trackh, Track);
      const recob::Track& track = *ptrack;
      TrackLength = track.Length();
      unsigned int NumTraj = track.NumberTrajectoryPoints();
      // ---- Get lengths and angles.
      if (NumTraj < 2) continue;
      else {
	TVector3 dir = track.VertexDirection<TVector3>();
	TrackTheta_XZ = std::atan2(dir.X(), dir.Z());
	TrackTheta_YZ = std::atan2(dir.Y(), dir.Z());
	TrackEta_XY   = std::atan2(dir.X(), dir.Y());
	TrackEta_ZY   = std::atan2(dir.Z(), dir.Y());
	TrackTheta    = dir.Theta();
	TrackPhi      = dir.Phi();
      }
      // ----- Hit Stuff ------
      std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(Track);
      int NumHits = allHits.size();
      for (int hit=0; hit<NumHits; ++hit) {
	// The recob::Hit -> View() is backwards...
	if ( allHits[hit]->View() == 2 ) ++TrackHits[0];
	if ( allHits[hit]->View() == 1 ) ++TrackHits[1];
	if ( allHits[hit]->View() == 0 ) ++TrackHits[2];
      }
      std::cout << "There were " << TrackHits[0] << ", " << TrackHits[1] << ", " << TrackHits[2] << " one each plane. " << std::endl;
      
      // *****************************************************************************
      // T0 stuff - So can correct X Start/End positions and identify MCParticle!!
      // *****************************************************************************
      double TickT0 = -1.;
      if ( fmt0.isValid() ) {
	std::vector<const anab::T0*> T0s = fmt0.at(Track);
	for (size_t t0size =0; t0size < T0s.size(); t0size++) {
	  MCTruthT0      = T0s[t0size]->Time();
	  MCTruthTickT0  = MCTruthT0 / sampling_rate(clockData);
	  MCTruthTrackID = T0s[t0size]->TriggerBits();
	} // T0 size
	TickT0 = MCTruthTickT0;
      } // T0 valid
      // ========== Rework out the MCTruth TrackID...Want to see if negative TrackID...
      std::cout << "\n\nWorking out this TrackID" << std::endl;
      std::map<int,double> trkide;
      for(size_t h = 0; h < allHits.size(); ++h){
	art::Ptr<recob::Hit> hit = allHits[h];
	std::vector<sim::IDE> ides;
	std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
	for(size_t e = 0; e < TrackIDs.size(); ++e){
	  trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
	}
      }
      // double maxe = -1, tote = 0; // tote unused
      double maxe = -1;
      int TrackID = -1;
      for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
	// tote += ii->second; // unused
	if ((ii->second)>maxe){
	  maxe = ii->second;
	  TrackID = ii->first;
	  std::cout << "Highest E track was " << TrackID << ", deposited " << maxe << std::endl;
	}
      }
      std::cout << "Old Best track was " << MCTruthTrackID << ", New Best track is " << TrackID << std::endl;
      if ( fabs(MCTruthTrackID) != fabs(TrackID) ) std::cout << "!!!!!I HAVE A TOTALLY DIFFERENT TRACK!!!!!!\n" << std::endl;
      else if ( MCTruthTrackID != TrackID )        std::cout << "The TrackIDs are different signs.\n" << std::endl;
      else                                         std::cout << "The TrackIDs are the same.\n" << std::endl;
      MCTruthTrackID = TrackID;
      std::cout << "The MCTruthTrackID is now " << MCTruthTrackID << std::endl;
      // ========== Rework out the MCTruth TrackID...Want to see if negative TrackID...

      // If using photon detectors...
      if ( fmphot.isValid() && fUsePhotons) {
	std::vector<const anab::T0*> PhotT0 = fmphot.at(Track);
	for (size_t T0it=0; T0it<PhotT0.size(); ++T0it) {
	  PhotonCounterT0     = PhotT0[T0it]->Time();
	  PhotonCounterTickT0 = PhotonCounterT0 / sampling_rate(clockData);
	  PhotonCounterID     = PhotT0[T0it]->TriggerBits();
	}
	TickT0 = PhotonCounterTickT0;
      } 
      // ************** END T0 stuff ***************
      
      if (TickT0 == -1.) continue;
      double XCorFac = detProp.ConvertTicksToX( TickT0, 0, 0, 0 );
      std::cout << "The TickT0 is " << TickT0 << ", giving an x correction to each hit of around " << XCorFac << std::endl;
      // ******************************************************************************************
      // Correct X and get track length etc now that we have matched a Track with an MCParticle!!
      // ******************************************************************************************

      // ---- Correct X positions!
      std::vector < TVector3 > CorrectedLocations;
      CorrectedLocations.clear();
      for ( unsigned int point=0; point < NumTraj; ++point ) {
	const TVector3 ThisLoc = track.LocationAtPoint<TVector3>(point);
	TVector3 CorrectLoc = ThisLoc;
	CorrectLoc[0] = CorrectLoc[0] - detProp.ConvertTicksToX( TickT0, allHits[NumTraj-(1+point)]->WireID().Plane, allHits[NumTraj-(1+point)]->WireID().TPC, allHits[NumTraj-(1+point)]->WireID().Cryostat );
	CorrectedLocations.push_back(CorrectLoc);
      }
      CorrectedStartX = CorrectedLocations[0][0];
      CorrectedStartY = CorrectedLocations[0][1];
      CorrectedStartZ = CorrectedLocations[0][2];
      CorrectedEndX   = CorrectedLocations[NumTraj-1][0];
      CorrectedEndY   = CorrectedLocations[NumTraj-1][1];
      CorrectedEndZ   = CorrectedLocations[NumTraj-1][2];
      
      // **************************************************************
      // Determine what kind of particle actually caused the track.....
      
      int ii = 0;
      for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
	simb::MCParticle *MyParticle = ipar->second;
	++ii;
	if ( MyParticle->TrackId() != fabs(MCTruthTrackID) )
	  continue;
	MatchedTrackID = Track;
	// ---- Get MCTruth Information and check that MCParticle goes in TPC ---
        MCTruthInformation ( clockData, MyParticle );

	// Work out if the track is back to front...
	// MC Start -> Track
	double St_St = CalcDist( MCStartX, MCStartY, MCStartZ, CorrectedStartX, CorrectedStartY, CorrectedStartZ );
	double St_En = CalcDist( MCStartX, MCStartY, MCStartZ, CorrectedEndX  , CorrectedEndY  , CorrectedEndZ   );
	// MC End   -> Track
	double En_St = CalcDist( MCEndX  , MCEndY  , MCEndZ  , CorrectedStartX, CorrectedStartY, CorrectedStartZ );
	double En_En = CalcDist( MCEndX  , MCEndY  , MCEndZ  , CorrectedEndX  , CorrectedEndY  , CorrectedEndZ   );

	// If backwards...
	if ( (St_En < St_St) && (En_St < En_En) ) {
	  std::cout << "I think that this track is backwards." << std::endl;
	} else {
	  MCCorOrient = true;
	}
	if (Verbose) {
	  std::cout <<"MC Start   ("<< MCStartX        << ", " << MCStartY        << ", " << MCStartZ        << ")\n"
		    << "Reco Start ("<< CorrectedStartX << ", " << CorrectedStartY << ", " << CorrectedStartZ << ")\n"
		    << "MC End     ("<< MCEndX          << ", " << MCEndY          << ", " << MCEndZ          << ")\n"		    
		    << "Reco End   ("<< CorrectedEndX   << ", " << CorrectedEndY   << ", " << CorrectedEndZ   << ")\n"
		    << "Dist of MC Start to Track Start/End is " << St_St << ", " << St_En << "\n"
		    << "Dist of MC End   to Track Start/End is " << En_St << ", " << En_En << "\n";
	  if (!MCCorOrient)
	    std::cout << "The track is the wrong way around...Will need to correct later.\n";
	}
	break;
      } // Check what particle caused this track....
      // ****************************************************************
      
      ++AllAll;
      if (fabs(MCPdgCode) == 13  ) ++MuonAll;
      else if (fabs(MCPdgCode) == 11 ) ++ElectronAll;
      else if (MCPdgCode == 2212 ) ++ProtonAll;
      else if (MCPdgCode == 22   ) ++GammaAll;
      else ++OtherAll;
            
      // Want to select only tracks which stop in the detector.
      TrackBoundaries( CorrectedLocations[0], CorrectedLocations[NumTraj-1] );
      //if ( RecoContainment == 0 || RecoContainment == 2 ) continue;
     
      // ---- Make a vector of calorimetry objects...
      std::vector<const anab::Calorimetry*> calos = fmcal.at(Track);
      PIDA = CalcPIDA ( calos );
     
      // ************************************************************
      // Now see what values of PIDA each particle type has
      // ************************************************************
      int QuickThresh = 25;
      if (PIDA_Plane2 > QuickThresh ) {
	std::cout << "\nTrack " << Track << " has a PIDA value of " << PIDA << ", " << NumTraj << " traj points and " << NumHits << " hits, PdGCode " << MCPdgCode << " , MCTrackID " << MCTrackId << std::endl;
	++AllBad;
	if (fabs(MCPdgCode) == 13  ) ++MuonBad;
	else if (fabs(MCPdgCode) == 11) ++ElectronBad;
	else if (MCPdgCode == 2212 ) ++ProtonBad;
	else if (MCPdgCode == 22   ) ++GammaBad;
	else ++OtherBad;
      }
      if (PIDA_Plane2 > 15 && PIDA_Plane2 < 25 && (RecoContainment == 1 || RecoContainment == 3) ) {
	++AllRange;
	if (fabs(MCPdgCode) == 13  ) ++MuonRange;
	else if (fabs(MCPdgCode) == 11) ++ElectronRange;
	else if (MCPdgCode == 2212 ) ++ProtonRange;
	else if (MCPdgCode == 22   ) ++GammaRange;
	  else ++OtherRange;
      }
      // ******** Fill Tree for each MCParticle **********
      fRecoTree->Fill();

      // Write some output so I can check how particles do when subject to my macro cuts....
      std::cout << "\n==== When subject to macro cuts ===" << std::endl;
      // What kind of particle is it?
      if (fabs(MCPdgCode) == 13)        std::cout << "This is a muon." << std::endl;
      else if (fabs(MCPdgCode) == 2212) std::cout << "This is a proton." << std::endl;
      else {
	std::cout << "This is a " << MCPdgCode << ", not interested, so continuing." << std::endl;
	continue;
      }
      // MC cont cuts
      if (MCContainment == 0 || MCContainment == 2 ) std::cout << "Would get cut by MC cut on stopping particles." << std::endl;
      else if (MCContainment == 1 ) std::cout << "Would get cut by MC cut on contained particles." << std::endl;
      else if (MCContainment == 3 ) std::cout << "This is a fully contained track.." << std::endl;
      // Reco cont cuts
      if (RecoContainment == 0 || RecoContainment == 2 ) std::cout << "Would get cut by Reco cut on stopping particles." << std::endl;
      else if (RecoContainment == 1 ) std::cout << "Would get cut by Reco cut on contained particles." << std::endl;
      else if (RecoContainment == 3 ) std::cout << "This is a fully contained track." << std::endl;
      // Unreasonably high PIDA value.
      if (PIDA > 25)         std::cout << "Unreasonably high PIDA of " << PIDA << ". ";
      else                   std::cout << "Has a Reasonabe PIDA of " << PIDA << ". ";
      if (PIDA_Plane2 > 25 ) std::cout << "Unreasonably high PIDA_Plane2 of " << PIDA_Plane2 << std::endl;
      else                   std::cout << "Has a Reasonabe PIDA_Plane2 of " << PIDA_Plane2 << std::endl;
      if (PIDA > 14 && PIDA < 18) std::cout << "Would lie in the proton PIDA area" << std::endl;
      if (PIDA > 5  && PIDA < 9 ) std::cout << "Would lie in the muon PIDA area" << std::endl;
      // Minimum number of coll plane hits
      if (CaloPlane2 < 5)       std::cout << "Has less than 5 collection plane hits." << std::endl;
      else if (CaloPlane2 < 10) std::cout << "Has less than 10 collection plane hits." << std::endl;
      else                        std::cout << "Has more than 10 collection plane hits " << CaloPlane2 << std::endl;
      // Correct orientation
      if (!MCCorOrient ) std::cout << "The track was the wrong way around!" << std::endl;
      else               std::cout << "The track was the right way around!" << std::endl;
      // Missed end points...
      if ( fabs( CorrectedEndX - MCEndX ) > 2.5 ) std::cout << "Missed the end point of the particle in X" << std::endl;
      if ( fabs( CorrectedEndY - MCEndY ) > 2.5 ) std::cout << "Missed the end point of the particle in Y" << std::endl;
      if ( fabs( CorrectedEndZ - MCEndZ ) > 2.5 ) std::cout << "Missed the end point of the particle in Z" << std::endl;
    } // Loop over Tracks
  } // if trackListHandle.isValid()
} // Analyse

// ******************************** Calc PIDA  ****************************************************
double ProtonIdentification::ProtonIdentification::CalcPIDA ( std::vector<const anab::Calorimetry*> calos ) {
  double PIDA  = 0, dEdxSum = 0;
  int UsedHits = 0, TotHits = 0;

  // I need to decide which way to go through the hits...
  double SumdEdx_St = 0., SumdEdx_En = 0., ResRng_St = 0, ResRng_En = 0;
  unsigned int AvHits = 5;
  if (calos[2]->dEdx().size()) {
    ResRng_St  = calos[2]->ResidualRange()[0];
    ResRng_En  = calos[2]->ResidualRange()[calos[2]->dEdx().size()-1];
  }
  // If don't have 2*AvHits (10) hits on the collection plane then use mid point + 1 hit
  if ( calos[2]->dEdx().size() < (2*AvHits) ) {
    AvHits = 1 + (0.5 * calos[2]->dEdx().size());
  }
  for ( unsigned int PlHit=0; PlHit < calos[2]->dEdx().size(); ++PlHit ) { // loop through hits on the collection plane
    if ( PlHit <= AvHits ) {
      SumdEdx_St += calos[2]->dEdx()[PlHit];
    }
    if ( calos[2]->dEdx().size() - PlHit <= AvHits ) {
      SumdEdx_En += calos[2]->dEdx()[PlHit];
    }
    //std::cout << "Looking at hit " << PlHit << " of " << (int)calos[2]->dEdx().size() << "...SumdEdx_St = " << SumdEdx_St << ", and SumdEdx_En = " << SumdEdx_En << std::endl;
  }
  double AvdEdx_St = SumdEdx_St / AvHits;
  double AvdEdx_En = SumdEdx_En / AvHits;
  // The dEdx at the start of the track should be less than that at the end...
  bool LowResSt = false;
  if ( ResRng_St < ResRng_En )
    LowResSt = true;
  bool LowdEdxSt = false;
  if ( AvdEdx_St < AvdEdx_En )
    LowdEdxSt = true;
  
  if ( LowResSt && (!MCCorOrient) ) {
    std::cout << "Track backwards, but calorimetry correct so setting everything to true...." << std::endl;
    MCCorOrient = true;
    double TmpX, TmpY, TmpZ;
    TmpX            = CorrectedStartX; TmpY            = CorrectedStartY; TmpZ            = CorrectedStartZ;
    CorrectedStartX = CorrectedEndX  ; CorrectedStartY = CorrectedEndY  ; CorrectedStartZ = CorrectedEndZ  ;
    CorrectedEndX   = TmpX           ; CorrectedEndY   = TmpY           ; CorrectedEndZ   = TmpZ           ;
    std::cout << "Reco Start ("<< CorrectedStartX   << ", " << CorrectedStartY   << ", " << CorrectedStartZ   << ")\n"
	      << "Reco End   ("<< CorrectedEndX     << ", " << CorrectedEndY     << ", " << CorrectedEndZ     << ")\n";
  }
  std::cout << "AvdEdx_St is " << AvdEdx_St << ", and AvdEdx_En is " << AvdEdx_En << "....ResRng_St is " << ResRng_St << ", and ResRng_En is " << ResRng_En 
	    << " ====>>> LowResSt " << LowResSt << ", LowdEdxSt " << LowdEdxSt << ", and TruthOrient? " << MCCorOrient << std::endl;

  if ( !MCCorOrient ) {
    std::cout << "This track is wrong?" << std::endl;
  }

  // *********** How do I deal with backwards tracks....What do I do if only one is backwards?.....
  //  If the dEdx is the wrong way around then without truth I would assume that the track is backwards.
  //  This means that I should use whether the MC is correct as a later cut.
  //  So if MCCorOrient == false then get PIDA using start of 'track'
  for ( int Plane=0; Plane < (int)calos.size(); ++Plane ) { // Loop through planes
    double PlanePIDA=0; int PlaneHits=0;
    for ( int PlaneHit=0; PlaneHit < (int)calos[Plane]->dEdx().size(); ++PlaneHit ) { // loop through hits on each plane
      double ThisdEdx = calos[Plane]->dEdx()[PlaneHit];
      double ThisResR = calos[Plane]->ResidualRange()[PlaneHit];
      // Increment TotHits and dEdx sum
      dEdxSum += ThisdEdx;
      ++TotHits;
      // Fill the calorimetry objects
      if (Plane==0 && CaloPlane0<kMaxHits) {
	dEdxPlane0[CaloPlane0] = ThisdEdx;
	ResRPlane0[CaloPlane0] = ThisResR;
	++CaloPlane0;
      } else if (Plane==1 && CaloPlane1<kMaxHits) {
	dEdxPlane1[CaloPlane1] = ThisdEdx;
	ResRPlane1[CaloPlane1] = ThisResR;
	++CaloPlane1;
      } else if (Plane==2 && CaloPlane2<kMaxHits) {
	dEdxPlane2[CaloPlane2] = ThisdEdx;
	ResRPlane2[CaloPlane2] = ThisResR;
	++CaloPlane2;
      }
      // Output some calorimetry information
      /*
      double ThisRang = calos[Plane]->Range();
      double ThisKinE = calos[Plane]->KineticEnergy();
      TVector3 ThisPos = calos[Plane]->XYZ()[PlaneHit];
      std::cout << "Looking at hit " << PlaneHit << " of " << (int)calos[Plane]->dEdx().size() << "...dEdx is " << ThisdEdx << ", ResR is " 
		<< ThisResR << ", this Range is " << ThisRang << ", TrackLength is " << TrackLength << ", KinE is " << ThisKinE 
		<< ". This is pos is (" << ThisPos[0] << ", " << ThisPos[1] << ", " << ThisPos[2] << ")." 
		<< std::endl;
      //*/
      
      // ==== If MCCorOrient == true
      // Work out PIDA if ResRange < 30 cm
      if ( ThisResR < 30 ) { // Only want PIDA for last 30 cm
	PlanePIDA += PIDAFunc( ThisdEdx, ThisResR );
	++PlaneHits;
      } // If ResRange < 30 cm
      // ===== This is where I need to do things...
    } // Loop over hits on each plane
      // Increment whole track PIDA.
    PIDA     += PlanePIDA;
    UsedHits += PlaneHits;
    // Work out PIDA for this plane
    PlanePIDA = PlanePIDA/PlaneHits;
    if (Plane == 0 ) PIDA_Plane0 = PlanePIDA;
    else if (Plane == 1 ) PIDA_Plane1 = PlanePIDA;
    else if (Plane == 2 ) PIDA_Plane2 = PlanePIDA;
  } // Loop over planes
  
  if ( UsedHits ) // If had any hits, work out PIDA and calculate
    PIDA = PIDA / UsedHits;
  AvdEdx = dEdxSum / TotHits;
  return PIDA;
} // CalcPIDA

// ******************************** Track Boundaries ****************************************************
void ProtonIdentification::ProtonIdentification::TrackBoundaries ( TVector3 larStart, TVector3 larEnd ) {
  
  if (MCCorOrient == false) {
    TVector3 larTemp = larStart;
    larStart = larEnd;
    larEnd   = larTemp;
    std::cout << "Swapping positions for reco containment. " << std::endl;
  }
  std::cout << "The track positions are (" << larStart[0] << ", " << larStart[1] << ", " << larStart[2] << ") ("<< larEnd[0] << ", " << larEnd[1] << ", " << larEnd[2] << ") and boundaries are "
	    << boundaries[0] << " -> " << boundaries[1] << ", " << boundaries[2] << " -> " << boundaries[3] << ", " <<boundaries[4] << " -> " << boundaries[5] << std::endl;
 
  double XStartDiff = std::min( fabs(boundaries[0]-larStart[0]), fabs(boundaries[1]-larStart[0]) );
  double YStartDiff = std::min( fabs(boundaries[2]-larStart[1]), fabs(boundaries[3]-larStart[1]) );
  double ZStartDiff = std::min( fabs(boundaries[4]-larStart[2]), fabs(boundaries[5]-larStart[2]) );
  if ( larStart[0] < boundaries[0] || larStart[0] > boundaries[1] ) XStartDiff = -XStartDiff;
  if ( larStart[1] < boundaries[2] || larStart[1] > boundaries[3] ) YStartDiff = -YStartDiff;
  if ( larStart[2] < boundaries[4] || larStart[2] > boundaries[5] ) ZStartDiff = -ZStartDiff;
  RecoStartInTPC = false;
  if ( XStartDiff > fBoundaryEdge && XStartDiff > 0 &&
       YStartDiff > fBoundaryEdge && YStartDiff > 0 && 
       ZStartDiff > fBoundaryEdge && ZStartDiff > 0
       ) {
    RecoStartInTPC = true;
  }
  StartFromEdge = std::min( XStartDiff, std::min(YStartDiff,ZStartDiff) );
  
  double XEndDiff = std::min( fabs(boundaries[0]-larEnd[0]), fabs(boundaries[1]-larEnd[0]) );
  double YEndDiff = std::min( fabs(boundaries[2]-larEnd[1]), fabs(boundaries[3]-larEnd[1]) );
  double ZEndDiff = std::min( fabs(boundaries[4]-larEnd[2]), fabs(boundaries[5]-larEnd[2]) );
  if ( larEnd[0] < boundaries[0] || larEnd[0] > boundaries[1] ) XEndDiff = -XEndDiff;
  if ( larEnd[1] < boundaries[2] || larEnd[1] > boundaries[3] ) YEndDiff = -YEndDiff;
  if ( larEnd[2] < boundaries[4] || larEnd[2] > boundaries[5] ) ZEndDiff = -ZEndDiff;
  RecoEndInTPC = false;
  if ( XEndDiff > fBoundaryEdge && XEndDiff > 0 && 
       YEndDiff > fBoundaryEdge && YEndDiff > 0 && 
       ZEndDiff > fBoundaryEdge && ZEndDiff > 0
       ) {
    RecoEndInTPC = true;
  }  
  EndFromEdge = std::min( XEndDiff, std::min(YEndDiff,ZEndDiff) );
  
  // *** What to do if the track is reconstructed backwards though....
  // *** Do I want to look into that here or elsewhere?
  
  if(!RecoStartInTPC && !RecoEndInTPC ) RecoContainment = 0;  // Through track
  if(!RecoStartInTPC &&  RecoEndInTPC ) RecoContainment = 1;  // Entering track
  if( RecoStartInTPC && !RecoEndInTPC ) RecoContainment = 2;  // Escaping track
  if( RecoStartInTPC &&  RecoEndInTPC ) RecoContainment = 3;  // Contained track
  // If the timing is totally off...
  if( XStartDiff < 0 || XEndDiff < 0 ) {
    std::cout << "This particle was outside active vol...." << std::endl;
    if ( std::max( fabs(larStart[1]-MCStartY), fabs(larStart[2]-MCStartZ ) ) < 5 &&
	 std::max( fabs(larEnd[1]  -MCEndY)  , fabs(larEnd[2]  -MCEndZ   ) ) < 5
	 ) {
      std::cout << "But the Y and Z difference are small, so changing to MCContainment" << std::endl;
      RecoContainment = MCContainment;
    }
  }

  std::cout << "RecoStartInTPC? " << RecoStartInTPC << ", RecoEndInTPC? " << RecoEndInTPC
	    << ", StartFromEdge " << StartFromEdge << ", EndFromEdge " << EndFromEdge
	    << ", RecoContainment " << RecoContainment << ", MCContainment " << MCContainment	    
	    << std::endl;
} // TrackBoundaries
// *********************************** Monte Carlo Truth Extraction ********************************************************
void ProtonIdentification::ProtonIdentification::MCTruthInformation (detinfo::DetectorClocksData const& clockData,
                                                                     const simb::MCParticle *particle ) {
  unsigned int numberTrajectoryPoints = particle->NumberTrajectoryPoints(); // Looking at each MC hit
  //double TPCLengthHits[numberTrajectoryPoints];
  //double TPCEnDepos   [numberTrajectoryPoints];
  std::vector<double> TPCLengthHits(numberTrajectoryPoints, 0);
  std::vector<double> TPCEnDepos(numberTrajectoryPoints, 0);
  bool BeenInVolume = false;
  int FirstHit=0, LastHit=0;
    
  for(unsigned int MCHit=0; MCHit <  TPCLengthHits.size(); ++MCHit) {
    if (MCHit!=0) {
      TPCLengthHits[MCHit] = pow ( pow( (particle->Vx(MCHit-1)-particle->Vx(MCHit)),2)
				   + pow( (particle->Vy(MCHit-1)-particle->Vy(MCHit)),2)
				   + pow( (particle->Vz(MCHit-1)-particle->Vz(MCHit)),2)
				   , 0.5 );
      TPCEnDepos[MCHit] = 1000 * (particle->E(MCHit-1) - particle->E(MCHit));
    }
    // --- Check if hit is in TPC...
    auto const position = geo::vect::toPoint(particle->Position(MCHit).Vect());
    geo::TPCID tpcid = geom->FindTPCAtPosition(position);
    if (tpcid.isValid) { 
      if (MCHit == 0 ) MCStartInTPC = true;
      if (MCHit == numberTrajectoryPoints-1 ) MCEndInTPC = true;
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
      MCEndX  = particle->Vx(MCHit) ; MCEndY = particle->Vy(MCHit) ; MCEndZ = particle->Vz(MCHit);
      if ( !BeenInVolume ) {
	BeenInVolume = true;
	FirstHit = MCHit;
	MCStartX = particle->Vx(MCHit) ; MCStartY = particle->Vy(MCHit) ; MCStartZ = particle->Vz(MCHit);
      }
    } // TPC.valid
  } // MCTrajPoints
  // What is the true containment?
  if(!MCStartInTPC && !MCEndInTPC ) MCContainment = 0;  // through track
  if(!MCStartInTPC &&  MCEndInTPC ) MCContainment = 1;  // entering track
  if( MCStartInTPC && !MCEndInTPC ) MCContainment = 2;  // escaping track
  if( MCStartInTPC &&  MCEndInTPC ) MCContainment = 3;  // contained track
  std::cout << "This track is from a " << particle->PdgCode() <<", StartInTPC? " << MCStartInTPC << ", EndInTPC? " << MCEndInTPC << ", MCContainment " << MCContainment << std::endl;

  // Work out the energy deposited etc.
  MCEnergy          = particle->E(FirstHit);
  MCEnergyDeposited = particle->E(FirstHit) - particle->E(LastHit);
  for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) {
    MCTPCLength += TPCLengthHits[Hit];
  }
  // ********* Need to work out MC PIDA *********
  // I think I want to do this using SimChannels..
  double SumPIDA = 0;
  int    MCHits  = 0;
  double NewDist = 0;
  for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) {
    NewDist += TPCLengthHits[Hit];
    double ThdEdx = TPCEnDepos[Hit] / TPCLengthHits[Hit];
    double ThResR = MCTPCLength - NewDist;
    // If no dEdx then can't add to PIDA value
    if (TPCEnDepos[Hit] == 0) continue;
    double ThPIDA = PIDAFunc( ThdEdx, ThResR );
    if ( ThResR < 30 ) { 
      SumPIDA += ThPIDA;
      ++MCHits;
    }
    if (Verbose > 1) {
      std::cout << "Looking at Hit " << Hit << " to " << LastHit << ". Total length was " << MCTPCLength << ", now gone " << NewDist << ". "
		<< "dEdx here was " << ThdEdx << ", ResRange was " << ThResR << ", so PIDA = " << ThPIDA << " = " << ThdEdx << " * " << ThResR << "^"<<PIDApower<< ".\n"
		<< "The PIDA sum is now " << SumPIDA << ", and I have used " << MCHits << " hits."
		<< std::endl;
    }
  }
  MCPIDA = SumPIDA / MCHits;
  std::cout << "The MCPIDA is " << MCPIDA << std::endl;

  // Set the MC parameters for the track
  MCPdgCode       = particle->PdgCode();
  MCTrackId       = particle->TrackId();
  if (particle->Process() == "primary") {
    MCIsPrimary = 1;
    if (MCTruthTrackID < 0)
      MCIsPrimary = -1;
  } else 
    MCIsPrimary = 0;
  std::cout << "What kind of primary identifier? " << MCIsPrimary << std::endl;
  TLorentzVector& momentumStart  = (TLorentzVector&)particle->Momentum(FirstHit);   
  TVector3 mcstartmom = particle->Momentum(FirstHit).Vect();
  MCTheta_XZ = std::atan2(momentumStart.Px(), momentumStart.Pz());
  MCTheta_YZ = std::atan2(momentumStart.Py(), momentumStart.Pz());
  MCEta_XY   = std::atan2(momentumStart.Px(), momentumStart.Py());
  MCEta_ZY   = std::atan2(momentumStart.Pz(), momentumStart.Py());
  MCTheta    = mcstartmom.Theta();
  MCPhi      = mcstartmom.Phi();
  
  return;
} // MCTruthInformation
// ******************************** Calc the dists ******************************************************
double ProtonIdentification::ProtonIdentification::CalcDist( double X1, double Y1, double Z1, double X2, double Y2, double Z2 ) {
  double Sq = TMath::Power( X1-X2 , 2 ) + TMath::Power( Y1-Y2 , 2 ) + TMath::Power( Z1-Z2 , 2 );
  double Va = TMath::Power( Sq, 0.5 );
  return Va;
} // Calc the dists
// ******************************** PIDA func ******************************************************
double ProtonIdentification::ProtonIdentification::PIDAFunc( double dedx, double resrng ) {
  double Va = dedx * pow( resrng, PIDApower );
  return Va;
}
// ******************************** Reset Variables *****************************************************
void ProtonIdentification::ProtonIdentification::ResetVars() {
  // MC info
  MCTPCLength  = MCEnergy   = MCEnergyDeposited = MCPIDA = 0;
  MCPdgCode    = MCTrackId  = MatchedTrackID = MCContainment = 0;
  MCStartX = MCStartY = MCStartZ = MCEndX = MCEndY = MCEndZ = 0;
  // Reco info  
  PIDA = PIDA_Plane0 = PIDA_Plane1   = PIDA_Plane2   = StartFromEdge   = EndFromEdge  = TrackLength    = 0;
  CorrectedStartX = CorrectedStartY = CorrectedStartZ = CorrectedEndX = CorrectedEndY = CorrectedEndZ = 0;
  TrackHits[0]    = TrackHits[1]    = TrackHits[2] = RecoContainment = 0;
  // Booleans
  MCStartInTPC = MCEndInTPC = RecoStartInTPC = RecoEndInTPC = false;
  MCCorOrient  = false;
  // Calorimetry
  CaloPlane0 = CaloPlane1 = CaloPlane2 = 0;
  for (int ii=0; ii<kMaxHits; ++ii) { 
    ResRPlane0[ii] = ResRPlane1[ii] = ResRPlane2[ii] = 0;
    dEdxPlane0[ii] = dEdxPlane1[ii] = dEdxPlane2[ii] = 0;
  }
}
// ******************************** Define Module *****************************************************
DEFINE_ART_MODULE(ProtonIdentification::ProtonIdentification)
