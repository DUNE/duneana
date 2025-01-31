////////////////////////////////////////////////////////////////////////
// Class:       NeutronDecayN2Ana
// Module Type: analyzer
// File:        NeutronDecayN2Ana_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
//
//  Having a look at FD reconstruction quantities
//
//  Thomas Karl Warburton
//  k.warburton@sheffield.ac.uk
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TTree.h"

//standard library includes
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>
#include <limits> // std::numeric_limits<>

struct IDEYLess {
  bool operator()(const sim::IDE& first, const sim::IDE& second) {
    return first.y < second.y;
  }
};

#define MaxPart   100
#define MaxParent 25
#define MaxPrim   10
#define HolderVal 9999

namespace NeutronDecayN2Ana {
  class NeutronDecayN2Ana;
}

class NeutronDecayN2Ana::NeutronDecayN2Ana : public art::EDAnalyzer {
public:
  explicit NeutronDecayN2Ana(fhicl::ParameterSet const & p);
  virtual ~NeutronDecayN2Ana();

  void analyze(art::Event const & e) override;

  void beginRun(art::Run const& run) override;
  void beginJob() override;
  void endJob() override;
  void endRun(art::Run const&) override;
  //void reconfigure(fhicl::ParameterSet const & p) ;
  
private:
  // ------ My functions ------
  void ResetVars();
  void FillVars( std::vector<int> &TrackIDVec, int &numParts, float EDep[MaxPart], float DaughtEDep[MaxPart], float DecayEDep[MaxPart], float Start[MaxPart][4], float End[MaxPart][4],
		 int nParents[MaxPart], int Parent[MaxPart][MaxParent], int ParTrID[MaxPart][MaxParent], int PDG[MaxPart], int TrID[MaxPart], int Cont[MaxPart], int FromDecay[MaxPart],
		 int ThisID, unsigned int ThisTDC, sim::IDE ThisIDE, const simb::MCParticle& MCPart, bool Decay, bool OrigParticle, bool &Written );
  bool  UnAssignLoop( int nPart, float St[MaxPart][4], sim::IDE thisIDE, float NearE[MaxPart] );
  float CalcDist( float X1, float Y1, float Z1, float X2, float Y2, float Z2 );
  bool  IsInTPC( float X1, float Y1, float Z1 , float Bound[6]);
  // Handles
  art::ServiceHandle<geo::Geometry> geom;

  // Parameter List
  int   Verbosity;
  int   NearEDeps;
  float ActiveBounds[6]; // Cryostat boundaries ( neg x, pos x, neg y, pos y, neg z, pos z )

  std::map<int, const simb::MCParticle*> truthmap; // A map of the truth particles.
  std::vector<int> AllTrackIDs; // A vector of all of my stored TrackIDs
  bool BadEvent;
  int NEvent;

  TTree* fDecayTree;
  int Run;
  int Event;
  float DistEdge[3][2];
  float TotalEDep;
  float EDepNearEdge2;
  float EDepNearEdge5;
  float EDepNearEdge10;

  float TopX, TopY, TopZ, BotX, BotY, BotZ;

  // Primary particles
  int   nPrim, PrimPDG[MaxPrim];
  float PrimEn[MaxPrim], PrimMom[MaxPrim];
  // NumParts
  int nMuon, nPion, nPi0, nKaon, nElec, nProt;
  // PdgCode
  int MuonPDG[MaxPart], PionPDG[MaxPart], Pi0PDG[MaxPart], KaonPDG[MaxPart], ElecPDG[MaxPart], ProtPDG[MaxPart];
  // TrackID
  int MuonTrID[MaxPart], PionTrID[MaxPart], Pi0TrID[MaxPart], KaonTrID[MaxPart], ElecTrID[MaxPart], ProtTrID[MaxPart];
  // TrackID
  int MuonCont[MaxPart], PionCont[MaxPart], Pi0Cont[MaxPart], KaonCont[MaxPart], ElecCont[MaxPart], ProtCont[MaxPart];
  // From a decay?
  int MuonFromDecay[MaxPart], PionFromDecay[MaxPart], Pi0FromDecay[MaxPart], KaonFromDecay[MaxPart], ElecFromDecay[MaxPart], ProtFromDecay[MaxPart];
  // NumParents
  int nParentMuon[MaxPart], nParentPion[MaxPart], nParentPi0[MaxPart], nParentKaon[MaxPart], nParentElec[MaxPart], nParentProt[MaxPart];
  // Muon Parent PDGs
  int MuonParents[MaxPart][MaxParent], PionParents[MaxPart][MaxParent], Pi0Parents[MaxPart][MaxParent], KaonParents[MaxPart][MaxParent], ElecParents[MaxPart][MaxParent], ProtParents[MaxPart][MaxParent];
  // Muon Parent IDs
  int MuonParTrID[MaxPart][MaxParent], PionParTrID[MaxPart][MaxParent], Pi0ParTrID[MaxPart][MaxParent], KaonParTrID[MaxPart][MaxParent], ElecParTrID[MaxPart][MaxParent], ProtParTrID[MaxPart][MaxParent];
  // EDep
  float MuonEDep[MaxPart], PionEDep[MaxPart], Pi0EDep[MaxPart], KaonEDep[MaxPart], ElecEDep[MaxPart], ProtEDep[MaxPart];
  // Daughter EDep
  float MuonDaughtersEDep[MaxPart], PionDaughtersEDep[MaxPart], Pi0DaughtersEDep[MaxPart], KaonDaughtersEDep[MaxPart], ElecDaughtersEDep[MaxPart], ProtDaughtersEDep[MaxPart];
  // Decay EDep
  float MuonDecayEDep[MaxPart], PionDecayEDep[MaxPart], Pi0DecayEDep[MaxPart], KaonDecayEDep[MaxPart], ElecDecayEDep[MaxPart], ProtDecayEDep[MaxPart];
  // Near EDep
  float MuonNearEDep[MaxPart] , PionNearEDep[MaxPart] , Pi0NearEDep[MaxPart] , KaonNearEDep[MaxPart] , ElecNearEDep[MaxPart], ProtNearEDep[MaxPart];
  // Start
  float MuonStart[MaxPart][4], PionStart[MaxPart][4], Pi0Start[MaxPart][4], KaonStart[MaxPart][4], ElecStart[MaxPart][4], ProtStart[MaxPart][4];
  // End
  float MuonEnd[MaxPart][4]  , PionEnd[MaxPart][4]  , Pi0End[MaxPart][4]  , KaonEnd[MaxPart][4]  , ElecEnd[MaxPart][4], ProtEnd[MaxPart][4];
};
// ******************************** Reset Variables *****************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::ResetVars() {
  BadEvent = false;
  Run = Event = 0;
  TotalEDep = EDepNearEdge2 = EDepNearEdge5 = EDepNearEdge10 = 0;
  // DistEdge
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      DistEdge[i][j]=HolderVal;
  //Primary particles
  nPrim = 0;
  for (int pr=0; pr<MaxPrim; ++pr) {
    PrimPDG[pr] = PrimEn[pr] = PrimMom[pr] = 0;
  }
  // Each particle type
  nMuon = nPion = nPi0 = nKaon = nElec = nProt = 0;
  for (int i=0; i<MaxPart; ++i) {
    MuonPDG[i]     = PionPDG[i]     = Pi0PDG[i]     = KaonPDG[i]     = ElecPDG[i]     = ProtPDG[i]     = 0;
    MuonTrID[i]    = PionTrID[i]    = Pi0TrID[i]    = KaonTrID[i]    = ElecTrID[i]    = ProtTrID[i]    = 0;
    MuonCont[i]    = PionCont[i]    = Pi0Cont[i]    = KaonCont[i]    = ElecCont[i]    = ProtCont[i]    = -1;
    MuonEDep[i]    = PionEDep[i]    = Pi0EDep[i]    = KaonEDep[i]    = ElecEDep[i]    = ProtEDep[i]    = 0;
    nParentMuon[i] = nParentPion[i] = nParentPi0[i] = nParentKaon[i] = nParentElec[i] = nParentProt[i] = 0;
    for (int j=0; j<MaxParent; ++j) {
      MuonParents[i][j] = PionParents[i][j] = Pi0Parents[i][j] = KaonParents[i][j] = ElecParents[i][j] = ProtParents[i][j] = 0;
      MuonParTrID[i][j] = PionParTrID[i][j] = Pi0ParTrID[i][j] = KaonParTrID[i][j] = ElecParTrID[i][j] = ProtParTrID[i][j] = 0;
    }
    MuonDaughtersEDep[i] = PionDaughtersEDep[i] = Pi0DaughtersEDep[i] = KaonDaughtersEDep[i] = ElecDaughtersEDep[i] = ProtDaughtersEDep[i] = 0;
    MuonDecayEDep[i]     = PionDecayEDep[i]     = Pi0DecayEDep[i]     = KaonDecayEDep[i]     = ElecDecayEDep[i]     = ProtDecayEDep[i]     = 0;
    MuonNearEDep[i]      = PionNearEDep[i]      = Pi0NearEDep[i]      = KaonNearEDep[i]      = ElecNearEDep[i]      = ProtNearEDep[i]      = 0;
    for (int j=0; j<4; ++j) {
      MuonStart[i][j] = PionStart[i][j] = Pi0Start[i][j] = KaonStart[i][j] = ElecStart[i][j] = ProtStart[i][j] = HolderVal;
      MuonEnd[i][j]   = PionEnd[i][j]   = Pi0End[i][j]   = KaonEnd[i][j]   = ElecEnd[i][j]   = ProtEnd[i][j]   = HolderVal;
    }
  }
}
// ********************************** Begin Run *******************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::beginRun(art::Run const& run) {
}
// *********************************** Begin Job ********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::beginJob()
{
  NEvent = 0;
  // Build my Cryostat boundaries array...Taken from Tyler Alion in Geometry Core. Should still return the same values for uBoone.
  ActiveBounds[0] = ActiveBounds[2] = ActiveBounds[4] = DBL_MAX;
  ActiveBounds[1] = ActiveBounds[3] = ActiveBounds[5] = -DBL_MAX;

  TopX = TopY = TopZ = -DBL_MAX;
  BotX = BotY = BotZ = DBL_MAX;

  // ----- FixMe: Assume single cryostats ------
  auto const* geom = lar::providerFrom<geo::Geometry>();
  for (geo::TPCGeo const& TPC: geom->Iterate<geo::TPCGeo>()) {
    // get center in world coordinates
    auto const center = TPC.GetCenter();
    //double tpcDim[3] = {TPC.ActiveHalfWidth(), TPC.ActiveHalfHeight(), 0.5*TPC.ActiveLength() };
    double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() }; // Gives same result as what Matt has
    if( center.X() - tpcDim[0] < ActiveBounds[0] ) ActiveBounds[0] = center.X() - tpcDim[0];
    if( center.X() + tpcDim[0] > ActiveBounds[1] ) ActiveBounds[1] = center.X() + tpcDim[0];
    if( center.Y() - tpcDim[1] < ActiveBounds[2] ) ActiveBounds[2] = center.Y() - tpcDim[1];
    if( center.Y() + tpcDim[1] > ActiveBounds[3] ) ActiveBounds[3] = center.Y() + tpcDim[1];
    if( center.Z() - tpcDim[2] < ActiveBounds[4] ) ActiveBounds[4] = center.Z() - tpcDim[2];
    if( center.Z() + tpcDim[2] > ActiveBounds[5] ) ActiveBounds[5] = center.Z() + tpcDim[2];
  } // for all TPC

  // Going from what I see in the event displays...
  ActiveBounds[0] = -723;
  ActiveBounds[1] = 723;
  ActiveBounds[2] = -600;
  ActiveBounds[3] = 600;
  ActiveBounds[4] = -1;
  ActiveBounds[5] = 5809;
  std::cout << "Active Boundaries: "
	    << "\n\tx: " << ActiveBounds[0] << " to " << ActiveBounds[1]
	    << "\n\ty: " << ActiveBounds[2] << " to " << ActiveBounds[3]
	    << "\n\tz: " << ActiveBounds[4] << " to " << ActiveBounds[5]
	    << std::endl;
  /*
  std::cout << "\n\n******Total World*******\n\n" << std::endl;
  std::cout << "The total mass is " << geom->TotalMass() << std::endl;

  std::cout << "\n\n******vol TPC Active Inner*******\n\n" << std::endl;
  std::cout << "The total mass is " << geom->TotalMass("volTPCActiveInner") << std::endl;

  std::cout << "\n\n******vol TPC Inner*******\n\n" << std::endl;
  std::cout << "The total mass is " << geom->TotalMass("volTPCInner") << std::endl;

  std::cout << "\n\n******vol TPC Active Outer*******\n\n" << std::endl;
  std::cout << "The total mass is " << geom->TotalMass("volTPCActiveOuter") << std::endl;

  std::cout << "\n\n******vol TPC Outer*******\n\n" << std::endl;
  std::cout << "The total mass is " << geom->TotalMass("volTPCOuter") << std::endl;

  std::cout << "\n\n******vol Cryostat*******\n\n" << std::endl;
  std::cout << "The total mass is " << geom->TotalMass("volCryostat") << std::endl;

  std::cout << "\n\n*************\n\n" << std::endl;
  */
  /*
  double minx_, maxx_, miny_, maxy_, minz_, maxz_;
  minx_ = miny_ = minz_ = DBL_MAX;
  maxx_ = maxy_ = maxz_ = -DBL_MAX;
  for (unsigned int c=0; c<geom->Ncryostats(); c++) {
    const geo::CryostatGeo& cryostat=geom->Cryostat(c);
    for (unsigned int t=0; t<cryostat.NTPC(); t++) {
      geo::TPCID id;
      id.Cryostat=c;
      id.TPC=t;
      id.isValid=true;
      const geo::TPCGeo& tpc=cryostat.TPC(id);
      //std::cout << t << "\t" << (tpc.Length()/2) << ", " << (tpc.ActiveLength()/2) << std::endl;
      if (tpc.MinX()<minx_) minx_=tpc.MinX();
      if (tpc.MaxX()>maxx_) maxx_=tpc.MaxX();
      if (tpc.MinY()<miny_) miny_=tpc.MinY();
      if (tpc.MaxY()>maxy_) maxy_=tpc.MaxY();
      if (tpc.MinZ()<minz_) minz_=tpc.MinZ();
      if (tpc.MaxZ()>maxz_) maxz_=tpc.MaxZ();  
    }
  }
  std::cout << "Using Matt's method " << minx_ << ", " << maxx_ << "\t" <<  miny_ << ", " << maxy_ << "\t" <<  minz_ << ", " << maxz_ << std::endl;
  */

  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fDecayTree = tfs->make<TTree>("ReconstructedTree","analysis tree");

  fDecayTree->Branch("Run"           , &Run           , "Run/I"            );
  fDecayTree->Branch("Event"         , &Event         , "Event/I"          );
  fDecayTree->Branch("DistEdge"      , &DistEdge      , "DistEdge[3][2]/F" );
  fDecayTree->Branch("TotalEDep"     , &TotalEDep     , "TotalEDep/F"      );
  fDecayTree->Branch("EDepNearEdge2" , &EDepNearEdge2 , "EDepNearEdge2/F"   );
  fDecayTree->Branch("EDepNearEdge5" , &EDepNearEdge5 , "EDepNearEdge5/F"   );
  fDecayTree->Branch("EDepNearEdge10", &EDepNearEdge10, "EDepNearEdge10/F"  );

  // Primary particles
  fDecayTree->Branch("nPrim"    , &nPrim    , "nPrim/I"           );
  fDecayTree->Branch("PrimPDG"  , &PrimPDG  , "PrimPDG[nPrim]/I"  );
  fDecayTree->Branch("PrimEn"   , &PrimEn   , "PrimEn[nPrim]/F"   );
  fDecayTree->Branch("PrimMom"  , &PrimMom  , "PrimMom[nPrim]/F"  );
  // Muon
  fDecayTree->Branch("nMuon"            ,&nMuon            ,"nMuon/I"                   );
  fDecayTree->Branch("nParentMuon"      ,&nParentMuon      ,"nParentMuon[nMuon]/I"      );
  fDecayTree->Branch("MuonParents"      ,&MuonParents      ,"MuonParents[nMuon][25]/I" );
  fDecayTree->Branch("MuonParTrID"      ,&MuonParTrID      ,"MuonParTrID[nMuon][25]/I" );
  fDecayTree->Branch("MuonPDG"          ,&MuonPDG          ,"MuonPDG[nMuon]/I"          );
  fDecayTree->Branch("MuonTrID"         ,&MuonTrID         ,"MuonTrID[nMuon]/I"         );
  fDecayTree->Branch("MuonCont"         ,&MuonCont         ,"MuonCont[nMuon]/I"         );
  fDecayTree->Branch("MuonFromDecay"    ,&MuonFromDecay    ,"MuonFromDecay[nMuon]/I"    );
  fDecayTree->Branch("MuonEDep"         ,&MuonEDep         ,"MuonEDep[nMuon]/F"         );
  fDecayTree->Branch("MuonStart"        ,&MuonStart        ,"MuonStart[nMuon][4]/F"     );
  fDecayTree->Branch("MuonEnd"          ,&MuonEnd          ,"MuonEnd[nMuon][4]/F"       );
  fDecayTree->Branch("MuonDaughtersEDep",&MuonDaughtersEDep,"MuonDaughtersEDep[nMuon]/F");
  fDecayTree->Branch("MuonDecayEDep"    ,&MuonDecayEDep    ,"MuonDecayEDep[nMuon]/F"    );
  fDecayTree->Branch("MuonNearEDep"     ,&MuonNearEDep     ,"MuonNearEDep[nMuon]/F"     );
  // Pion
  fDecayTree->Branch("nPion"            ,&nPion            ,"nPion/I"                   );
  fDecayTree->Branch("nParentPion"      ,&nParentPion      ,"nParentPion[nPion]/I"      );
  fDecayTree->Branch("PionParents"      ,&PionParents      ,"PionParents[nPion][25]/I" );
  fDecayTree->Branch("PionParTrID"      ,&PionParTrID      ,"PionParTrID[nPion][25]/I" );
  fDecayTree->Branch("PionPDG"          ,&PionPDG          ,"PionPDG[nPion]/I"          );
  fDecayTree->Branch("PionTrID"         ,&PionTrID         ,"PionTrID[nPion]/I"         );
  fDecayTree->Branch("PionCont"         ,&PionCont         ,"PionCont[nPion]/I"         );
  fDecayTree->Branch("PionFromDecay"    ,&PionFromDecay    ,"PionFromDecay[nPion]/I"    );
  fDecayTree->Branch("PionEDep"         ,&PionEDep         ,"PionEDep[nPion]/F"         );
  fDecayTree->Branch("PionStart"        ,&PionStart        ,"PionStart[nPion][4]/F"     );
  fDecayTree->Branch("PionEnd"          ,&PionEnd          ,"PionEnd[nPion][4]/F"       );
  fDecayTree->Branch("PionDaughtersEDep",&PionDaughtersEDep,"PionDaughtersEDep[nPion]/F");
  fDecayTree->Branch("PionDecayEDep"    ,&PionDecayEDep    ,"PionDecayEDep[nPion]/F"    );
  fDecayTree->Branch("PionNearEDep"     ,&PionNearEDep     ,"PionNearEDep[nPion]/F"     );
  // Pi0
  fDecayTree->Branch("nPi0"            ,&nPi0            ,"nPi0/I"                  );
  fDecayTree->Branch("nParentPi0"      ,&nParentPi0      ,"nParentPi0[nPi0]/I"      );
  fDecayTree->Branch("Pi0Parents"      ,&Pi0Parents      ,"Pi0Parents[nPi0][25]/I" );
  fDecayTree->Branch("Pi0ParTrID"      ,&Pi0ParTrID      ,"Pi0ParTrID[nPi0][25]/I" );
  fDecayTree->Branch("Pi0PDG"          ,&Pi0PDG          ,"Pi0PDG[nPi0]/I"          );
  fDecayTree->Branch("Pi0TrID"         ,&Pi0TrID         ,"Pi0TrID[nPi0]/I"         );
  fDecayTree->Branch("Pi0Cont"         ,&Pi0Cont         ,"Pi0Cont[nPi0]/I"         );
  fDecayTree->Branch("Pi0FromDecay"    ,&Pi0FromDecay    ,"Pi0FromDecay[nPi0]/I"    );
  fDecayTree->Branch("Pi0EDep"         ,&Pi0EDep         ,"Pi0EDep[nPi0]/F"         );
  fDecayTree->Branch("Pi0Start"        ,&Pi0Start        ,"Pi0Start[nPi0][4]/F"     );
  fDecayTree->Branch("Pi0End"          ,&Pi0End          ,"Pi0End[nPi0][4]/F"       );
  fDecayTree->Branch("Pi0DaughtersEDep",&Pi0DaughtersEDep,"Pi0DaughtersEDep[nPi0]/F");
  fDecayTree->Branch("Pi0DecayEDep"    ,&Pi0DecayEDep    ,"Pi0DecayEDep[nPi0]/F"    );
  fDecayTree->Branch("Pi0NearEDep"     ,&Pi0NearEDep     ,"Pi0NearEDep[nPi0]/F"     );
  // Kaon
  fDecayTree->Branch("nKaon"            ,&nKaon            ,"nKaon/I"                   );
  fDecayTree->Branch("nParentKaon"      ,&nParentKaon      ,"nParentKaon[nKaon]/I"      );
  fDecayTree->Branch("KaonParents"      ,&KaonParents      ,"KaonParents[nKaon][25]/I" );
  fDecayTree->Branch("KaonParTrID"      ,&KaonParTrID      ,"KaonParTrID[nKaon][25]/I" );
  fDecayTree->Branch("KaonPDG"          ,&KaonPDG          ,"KaonPDG[nKaon]/I"          );
  fDecayTree->Branch("KaonTrID"         ,&KaonTrID         ,"KaonTrID[nKaon]/I"         );
  fDecayTree->Branch("KaonCont"         ,&KaonCont         ,"KaonCont[nKaon]/I"         );
  fDecayTree->Branch("KaonFromDecay"    ,&KaonFromDecay    ,"KaonFromDecay[nKaon]/I"    );
  fDecayTree->Branch("KaonEDep"         ,&KaonEDep         ,"KaonEDep[nKaon]/F"         );
  fDecayTree->Branch("KaonStart"        ,&KaonStart        ,"KaonStart[nKaon][4]/F"     );
  fDecayTree->Branch("KaonEnd"          ,&KaonEnd          ,"KaonEnd[nKaon][4]/F"       );
  fDecayTree->Branch("KaonDaughtersEDep",&KaonDaughtersEDep,"KaonDaughtersEDep[nKaon]/F");
  fDecayTree->Branch("KaonDecayEDep"    ,&KaonDecayEDep    ,"KaonDecayEDep[nKaon]/F"    );
  fDecayTree->Branch("KaonNearEDep"     ,&KaonNearEDep     ,"KaonNearEDep[nKaon]/F"     );
  // Electron
  fDecayTree->Branch("nElec"            ,&nElec            ,"nElec/I"                   );
  fDecayTree->Branch("nParentElec"      ,&nParentElec      ,"nParentElec[nElec]/I"      );
  fDecayTree->Branch("ElecParents"      ,&ElecParents      ,"ElecParents[nElec][25]/I" );
  fDecayTree->Branch("ElecParTrID"      ,&ElecParTrID      ,"ElecParTrID[nElec][25]/I" );
  fDecayTree->Branch("ElecFromDecay"    ,&ElecFromDecay    ,"ElecFromDecay[nElec]/I"    );
  fDecayTree->Branch("ElecPDG"          ,&ElecPDG          ,"ElecPDG[nElec]/I"          );
  fDecayTree->Branch("ElecTrID"         ,&ElecTrID         ,"ElecTrID[nElec]/I"         );
  fDecayTree->Branch("ElecCont"         ,&ElecCont         ,"ElecCont[nElec]/I"         );
  fDecayTree->Branch("ElecEDep"         ,&ElecEDep         ,"ElecEDep[nElec]/F"         );
  fDecayTree->Branch("ElecStart"        ,&ElecStart        ,"ElecStart[nElec][4]/F"     );
  fDecayTree->Branch("ElecEnd"          ,&ElecEnd          ,"ElecEnd[nElec][4]/F"       );
  fDecayTree->Branch("ElecDaughtersEDep",&ElecDaughtersEDep,"ElecDaughtersEDep[nElec]/F");
  fDecayTree->Branch("ElecDecayEDep"    ,&ElecDecayEDep    ,"ElecDecayEDep[nElec]/F"    );
  fDecayTree->Branch("ElecNearEDep"     ,&ElecNearEDep     ,"ElecNearEDep[nElec]/F"     );
  // Proton
  fDecayTree->Branch("nProt"            ,&nProt            ,"nProt/I"                   );
  fDecayTree->Branch("nParentProt"      ,&nParentProt      ,"nParentProt[nProt]/I"      );
  fDecayTree->Branch("ProtParents"      ,&ProtParents      ,"ProtParents[nProt][25]/I" );
  fDecayTree->Branch("ProtParTrID"      ,&ProtParTrID      ,"ProtParTrID[nProt][25]/I" );
  fDecayTree->Branch("ProtFromDecay"    ,&ProtFromDecay    ,"ProtFromDecay[nProt]/I"    );
  fDecayTree->Branch("ProtPDG"          ,&ProtPDG          ,"ProtPDG[nProt]/I"          );
  fDecayTree->Branch("ProtTrID"         ,&ProtTrID         ,"ProtTrID[nProt]/I"         );
  fDecayTree->Branch("ProtCont"         ,&ProtCont         ,"ProtCont[nProt]/I"         );
  fDecayTree->Branch("ProtEDep"         ,&ProtEDep         ,"ProtEDep[nProt]/F"         );
  fDecayTree->Branch("ProtStart"        ,&ProtStart        ,"ProtStart[nProt][4]/F"     );
  fDecayTree->Branch("ProtEnd"          ,&ProtEnd          ,"ProtEnd[nProt][4]/F"       );
  fDecayTree->Branch("ProtDaughtersEDep",&ProtDaughtersEDep,"ProtDaughtersEDep[nProt]/F");
  fDecayTree->Branch("ProtDecayEDep"    ,&ProtDecayEDep    ,"ProtDecayEDep[nProt]/F"    );
  fDecayTree->Branch("ProtNearEDep"     ,&ProtNearEDep     ,"ProtNearEDep[nProt]/F"     );
}
// ************************************ End Job *********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::endJob() {
  std::cout << "\nAfter all of that Top = " << TopX << ", " << TopY << ", " << TopZ << ". Bot = " << BotX << ", " << BotY << ", " << BotZ << std::endl;
}
// ************************************ End Run *********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::endRun(art::Run const&) {
}

// ********************************** pset param *******************************************************
NeutronDecayN2Ana::NeutronDecayN2Ana::NeutronDecayN2Ana(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , Verbosity( pset.get< int >( "Verbosity" ))
  , NearEDeps( pset.get< int >( "NearEDeps" ))
{

}
// ******************************************************************************************************
NeutronDecayN2Ana::NeutronDecayN2Ana::~NeutronDecayN2Ana()
{
  // Clean up dynamic memory and other resources here.

}
// ************************************ Analyse *********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::analyze(art::Event const & evt) {
  ++NEvent;
  ResetVars();
  Run   = evt.run();
  Event = evt.event();
  
  //if ( Event != 39 ) return;

  if (Verbosity)
    std::cout << "\n\n************* New Event / Module running - Run " << Run << ", Event " << Event << ", NEvent " << NEvent << " *************\n\n" << std::endl;

  // Any providers I need.
  auto const* geo = lar::providerFrom<geo::Geometry>();
  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  /*
  // Implementation of required member function here. 
  std::vector<art::Ptr<recob::Track> > tracklist;
  auto trackListHandle = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  if (trackListHandle)
    art::fill_ptr_vector(tracklist, trackListHandle);
  */
  // Make a map of MCParticles which I can access later.
  auto truth = evt.getHandle<std::vector<simb::MCParticle> >("largeant");
  truthmap.clear();
  for (size_t i=0; i<truth->size(); i++) {
    truthmap[truth->at(i).TrackId()]=&((*truth)[i]);
    if (truth->at(i).Mother() == 0) {
      PrimPDG[nPrim] = truth->at(i).PdgCode();
      PrimEn [nPrim] = truth->at(i).E(0);
      PrimMom[nPrim] = truth->at(i).P(0);
      ++nPrim;
      if (Verbosity) {
	std::cout << "The primary particle in this event was a " << PrimPDG[nPrim-1]<<" with Energy " << PrimEn[nPrim-1] 
		  << ", momentum " << PrimMom[nPrim-1]<< ", process " << truth->at(i).Process() << std::endl; 
      }
    }
  }
  if (Verbosity)
    std::cout << "There were " << nPrim << " primary particles." << std::endl;

  // Get a vector of sim channels.
  auto simchannels = evt.getHandle<std::vector<sim::SimChannel> >("largeant");
  
  // Make vectors to hold all of my particle TrackIDs
  std::vector<int> MuonVec, PionVec, Pi0Vec, KaonVec, ElecVec, ProtVec;
  AllTrackIDs.clear();

  // Want to loop through the simchannels, and the ides.
  int chanIt = 0;
  
  std::vector< sim::IDE > UnAssignVec;
  UnAssignVec.clear();

  for (auto const& simchannel:*simchannels) {
    // ------ Only want to look at collection plane hits ------
    //std::cout << "Looking at a new SimChannel, it was on channel " << simchannel.Channel() << ", which is plane " << geo->SignalType(simchannel.Channel()) << std::endl;
    if (wireReadout.SignalType(simchannel.Channel()) != geo::kCollection) continue;
    int tdcideIt = 0;
    // ------ Loop through all the IDEs for this channel ------
    for (auto const& tdcide:simchannel.TDCIDEMap()) {
      unsigned int tdc = tdcide.first;
      auto const& idevec=tdcide.second;
      if (Verbosity>1)
	std::cout << "TDC IDE " << tdcideIt << " of " << simchannel.TDCIDEMap().size() << ", has tdc " << tdc << ", and idevec of size " << idevec.size() << std::endl;
      int ideIt = 0;
      // ------ Look at each individual IDE ------
      for (auto const& ide:idevec) {
	int ideTrackID = ide.trackID;
	float ideEnergy = ide.energy;
	/*
	  if (ideTrackID == 1 && Verbosity == 2) {
	  float ideNumEl  = ide.numElectrons;
	  std::cout << "      IDE " << ideIt << " of " << idevec.size() << " has TrackId " << ideTrackID << ", energy " << ideEnergy << ", NumEl " << ideNumEl << ", tdc " << tdc << ". ";
	  std::cout << "Position " << idePos[0] << ", " << idePos[1] << ", " << idePos[2] << std::endl;
	  }
	//*/
	
	// ***** Find out what particle the ide is due to...
	const simb::MCParticle& Origpart=*( truthmap[ abs(ideTrackID) ] );
	int OrigPdgCode = Origpart.PdgCode();

        geo::TPCID tpcid=geo->FindTPCAtPosition(geo::Point_t{ide.x, ide.y, ide.z});
	if (!(geo->HasTPC(tpcid)) ) {
	  if (Verbosity)
	    std::cout << "Outside the Active volume I found at the top!" << std::endl;
	  continue;
	}
	  
	// ------ I want to work the closest IDE to an edge of the active volume ------
	// If I am writing out the distance to the edge of the active volume
	if ( DistEdge[0][0] > ( ide.x - ActiveBounds[0]) ) DistEdge[0][0] =  ide.x - ActiveBounds[0];
	if ( DistEdge[0][1] > (-ide.x + ActiveBounds[1]) ) DistEdge[0][1] = -ide.x + ActiveBounds[1];
	if ( DistEdge[1][0] > ( ide.y - ActiveBounds[2]) ) DistEdge[1][0] =  ide.y - ActiveBounds[2];
	if ( DistEdge[1][1] > (-ide.y + ActiveBounds[3]) ) DistEdge[1][1] = -ide.y + ActiveBounds[3];
	if ( DistEdge[2][0] > ( ide.z - ActiveBounds[4]) ) DistEdge[2][0] =  ide.z - ActiveBounds[4];
	if ( DistEdge[2][1] > (-ide.z + ActiveBounds[5]) ) DistEdge[2][1] = -ide.z + ActiveBounds[5];
	
	// Work out the EDeps within a distance to the detector edge.
	float XEDep = std::min( ide.x - ActiveBounds[0], -ide.x + ActiveBounds[1] );
	float YEDep = std::min( ide.y - ActiveBounds[2], -ide.y + ActiveBounds[3] );
	float ZEDep = std::min( ide.z - ActiveBounds[4], -ide.z + ActiveBounds[5] );
	float MEDep = std::min( XEDep, std::min( YEDep, ZEDep ) );
	// If the deposition was outside my active volume continue...
	if ( MEDep < 0 ) {
	  //std::cout << "Deposition outside active vol ("<<ide.x<<", " << ide.y<< ", "<< ide.z<<")." << std::endl;
	  continue;
	}
	if ( MEDep < 2 ) {
	  //std::cout << "Edep at " << ide.x << ", " << ide.y << " " << ide.z << "..." << MEDep << ", from a " << OrigPdgCode 
	  //	    << ", Energy = " << ide.energy << " ==>> " << EDepNearEdge2+ide.energy << std::endl;
	  EDepNearEdge2  += ide.energy;
	}
	if ( MEDep < 5 ) {
	  EDepNearEdge5  += ide.energy;
	}
	if ( MEDep < 10 ) {
	  EDepNearEdge10 += ide.energy;
	}
	
	if (ide.x > TopX) {
	  TopX = ide.x;
	  //std::cout << "Changed TopX to " << TopX << std::endl;
	}
	if (ide.y > TopY) {
	  TopY = ide.y;
	  //std::cout << "Changed TopY to " << TopY << std::endl;
	}
	if (ide.z > TopZ) {
	  TopZ = ide.z;
	  //std::cout << "Changed TopZ to " << TopZ << std::endl;
	}
	if (ide.x < BotX) BotX = ide.x;
	if (ide.y < BotY) BotY = ide.y;
	if (ide.z < BotZ) BotZ = ide.z;

	// ------ Add the energy deposition from this IDE to the sum of IDEs
	TotalEDep += ideEnergy;

	//std::cout << "IDE is at (" << ide.x << ", " << ide.y << ", " << ide.z << "). MinX is " << XEDep << ", MinY is " << YEDep << ", MinZ is " << ZEDep << "===> MinDist is " << MEDep
	//	  << "...is this less than DistToEdge? " << IsClose << " ====> EDepNearEdge is now " << EDepNearEdge
	//	  << std::endl;
	
	// If I am writing out the most +- depositions in X, Y, Z
	/*
	if ( DistEdge[0][0] > ide.x ) DistEdge[0][0] = ide.x;
	if ( DistEdge[0][1] < ide.x ) DistEdge[0][1] = ide.x;
	if ( DistEdge[1][0] > ide.y ) DistEdge[1][0] = ide.y;
	if ( DistEdge[1][1] < ide.y ) DistEdge[1][1] = ide.y;
	if ( DistEdge[2][0] > ide.z ) DistEdge[2][0] = ide.z;
	if ( DistEdge[2][1] < ide.z ) DistEdge[2][1] = ide.z;
	*/
	// ------ Now to work out which particles in particular I want to save more information about... ------
	// ----------------- I want to write out the IDE to the relevant parent particle type -----------------
	// ------- This means looping back through the IDE's parents until I find an interesting TrackID ------
	bool isDecay = false;
	bool WrittenOut = false;
	bool OrigPart = true;
	if (Verbosity>1)
	  std::cout << "\nLooking at IDE " << ideIt << ", ideTrackID is " << ideTrackID << ", it was due to a " 
		    << truthmap[ abs(ideTrackID) ]->PdgCode() << ", process " << truthmap[ abs(ideTrackID) ]->Process() 
		    << std::endl;
	while ( ideTrackID != 0 && !WrittenOut ) {
	  const simb::MCParticle& part=*( truthmap[ abs(ideTrackID) ] );
	  int PdgCode=part.PdgCode();
	  if ( PdgCode != OrigPdgCode || ideTrackID < 0 )
	    OrigPart = false;
	  // ========== Muons ==========
	  if      ( (PdgCode == -13  || PdgCode == 13) ) 
	    FillVars( MuonVec, nMuon, MuonEDep, MuonDaughtersEDep, MuonDecayEDep, MuonStart, MuonEnd, nParentMuon, MuonParents, MuonParTrID, MuonPDG, MuonTrID, MuonCont, MuonFromDecay,
		      ideTrackID, tdc, ide, part, OrigPart, isDecay, WrittenOut  );
	  // ========== Pions ==========
	  else if ( (PdgCode == -211 || PdgCode == 211) && part.Process() != "pi+Inelastic" && part.Process() != "pi-Inelastic" ) 
	    FillVars( PionVec, nPion, PionEDep, PionDaughtersEDep, PionDecayEDep, PionStart, PionEnd, nParentPion, PionParents, PionParTrID, PionPDG, PionTrID, PionCont, PionFromDecay,
		      ideTrackID, tdc, ide, part, OrigPart, isDecay, WrittenOut);
	  // ========== Pi0s  ==========
	  else if ( PdgCode == 111 )
	    FillVars( Pi0Vec , nPi0 , Pi0EDep , Pi0DaughtersEDep , Pi0DecayEDep , Pi0Start , Pi0End , nParentPi0 , Pi0Parents , Pi0ParTrID , Pi0PDG , Pi0TrID , Pi0Cont , Pi0FromDecay ,
		      ideTrackID, tdc, ide, part, OrigPart, isDecay, WrittenOut );
	  // ========== Kaons ===========
	  else if ( (PdgCode == 321 || PdgCode == -321) && part.Process() != "kaon+Inelastic" && part.Process() != "kaon-Inelastic" ) 
	    FillVars( KaonVec, nKaon, KaonEDep, KaonDaughtersEDep, KaonDecayEDep, KaonStart, KaonEnd, nParentKaon, KaonParents, KaonParTrID, KaonPDG, KaonTrID, KaonCont, KaonFromDecay, 
		      ideTrackID, tdc, ide, part, OrigPart, isDecay, WrittenOut );
	  // ========== Elecs ===========
	  else if ( (PdgCode == -11 || PdgCode == 11) ) {
	    // Electrons can shower straight away, so I want to treat all deposits as if it was from the electron.
	    OrigPart = true;
	    FillVars( ElecVec, nElec, ElecEDep, ElecDaughtersEDep, ElecDecayEDep, ElecStart, ElecEnd, nParentElec, ElecParents, ElecParTrID, ElecPDG, ElecTrID, ElecCont, ElecFromDecay,
		      ideTrackID, tdc, ide, part, OrigPart, isDecay, WrittenOut );
	  // ========== Prots ===========
	  } else if ( PdgCode == 2212 )
	    FillVars( ProtVec, nProt, ProtEDep, ProtDaughtersEDep, ProtDecayEDep, ProtStart, ProtEnd, nParentProt, ProtParents, ProtParTrID, ProtPDG, ProtTrID, ProtCont, ProtFromDecay,
		      ideTrackID, tdc, ide, part, OrigPart, isDecay, WrittenOut );
	  // ========== If still not one of my intersting particles I need to find this particles parent ==========
	  else {
	    ideTrackID = part.Mother();
	    if (ideTrackID==0) {
	      if (Verbosity > 1) {
		int ThrowIDE = ide.trackID;
		std::cout << "None of the particles in this chain are interesting, so skipping. It was due to TrackID " << ThrowIDE << ", PdGCode " << truthmap[ abs(ThrowIDE) ]->PdgCode() 
			  << ", energy " << ide.energy << ", and pos ("<<ide.x<<", "<<ide.y<<", "<<ide.z<<").\n"
			  << "  The ide parentage was PdG (TrackID) " << truthmap[ abs(ThrowIDE) ]->PdgCode() << " (" << ThrowIDE << ") ";
		while ( ThrowIDE != 0 ) {
		  ThrowIDE = truthmap[ abs(ThrowIDE) ]->Mother();
		  if ( ThrowIDE == 0) 
		    std::cout << ". end." << std::endl;
		  else
		    std::cout << " <-- " << truthmap[ abs(ThrowIDE) ]->PdgCode() << " (" << ThrowIDE << ") ";
		}
	      }
	      UnAssignVec.push_back( ide );
	      break;
	    }
	    PdgCode = truthmap[ abs(ideTrackID) ]->PdgCode();
	    // === Work out if a decay I care about ===
	    if ( part.PdgCode() != PdgCode ) isDecay = false; // Only reset to false if not a daughter of the same particle type eg not K+ -> K+
	    if ( ideTrackID > 0 && part.Process() == "Decay" && 
		 ( PdgCode == -13 || PdgCode == 13   || PdgCode == -211 || PdgCode == 211 || // muon or pion
		   PdgCode == 321 || PdgCode == -321 || PdgCode == -11  || PdgCode == 11  || // kaon or elec
		   PdgCode == 111 || // pi0
		   PdgCode == 2212   // proton
		   )
 		 )
	      { // If decay I care about.
		//std::cout << "This particle was from a decay!" << std::endl;
		isDecay = true;
	      }
	    if (Verbosity > 1) {
	      std::cout << "Not something interesting so moving backwards. The parent of that particle had trackID " << ideTrackID 
			<< ", was from a " << PdgCode << ", from a decay? " << isDecay << std::endl;
	    }
	  } // If not one of chosen particles.
	} // While loop
	if (BadEvent) {
	  if (Verbosity) {
 	    std::cout << "Had too many of one particle type, voiding this event now...nKaon " << nKaon << ", nElec " << nElec << ", nMuon " << nMuon 
		      << ", nPion " << nPion << ", nPi0 " << nPi0 << ", nProt " << nProt << "....Up to now there was " << TotalEDep << " energy deposited." << std::endl;
	  }
	  return;
	}
	ideIt++;
      } // Each IDE ( ide:idevec )
      ++tdcideIt;
    } // IDE vector for SimChannel ( tdcide:simcahnnel.TPCIDEMap() )
    ++chanIt;
  } // Loop through simchannels
  
  // ------------------------------------- Now loop through all of my unassigned IDEs -------------------------------------
  std::cout << "\nLooked through all of my IDEs, my UnAssignVec has size " << UnAssignVec.size() << std::endl; 
  int StillUnassign = 0;
  for (auto const& ide:UnAssignVec) {
    if (Verbosity > 1) {
      std::cout << "Going through my uninteresting particles...this IDE was due to TrackID " << ide.trackID << ", energy " << ide.energy 
		<< ", and pos ("<<ide.x<<", "<<ide.y<<", "<<ide.z<<"). " << std::endl;
      int ideTrId = ide.trackID;
      std::cout << " The ide parentage was PdG (TrackID) " << truthmap[ abs(ideTrId) ]->PdgCode() << " (" << ideTrId << ") ";
      while ( ideTrId != 0 ) {
	ideTrId = truthmap[ abs(ideTrId) ]->Mother();
	if ( ideTrId == 0) 
	  std::cout << ". end." << std::endl;
	else
	  std::cout << " <-- " << truthmap[ abs(ideTrId) ]->PdgCode() << " (" << ideTrId << ") ";
      }
    }
    bool FoundDep = false;
    // ==== Loop through kaons    
    FoundDep = UnAssignLoop( nKaon, KaonStart, ide, KaonNearEDep ); 
    if (FoundDep) continue;
    // ==== Loop through elecs    
    FoundDep = UnAssignLoop( nElec, ElecStart, ide, ElecNearEDep ); 
    if (FoundDep) continue;
    // ==== Loop through muons    
    FoundDep = UnAssignLoop( nMuon, MuonStart, ide, MuonNearEDep ); 
    if (FoundDep) continue;
    // ==== Loop through pions    
    FoundDep = UnAssignLoop( nPion, PionStart, ide, PionNearEDep ); 
    if (FoundDep) continue;
    // ==== Loop through pi0s    
    FoundDep = UnAssignLoop( nPi0, Pi0Start, ide, Pi0NearEDep ); 
    if (FoundDep) continue;
    // ==== Loop through protons
    FoundDep = UnAssignLoop( nProt, ProtStart, ide, ProtNearEDep ); 
    if (FoundDep) continue;
    // ==== If still not assigned this IDE
    if (Verbosity > 1) std::cout << "!!!This IDE was nowhere near anything...." << std::endl;
    ++StillUnassign;
  }
  std::cout << "There were still " << StillUnassign << " unassigned IDEs." << std::endl;
  // ------------------------------------- Now loop through all of my unassigned IDEs -------------------------------------

  // ------------------------------------- Output some information about my particles -------------------------------------
  if (Verbosity) {
    if (nKaon) {
      std::cout << "\nThere are kaons in this event!!" << std::endl;
      for (int KL=0; KL<nKaon; ++KL) {
	std::cout << "Kaon " << KL << " had properties: PDG " << KaonPDG[KL] << ", TrackID " << KaonTrID[KL] << ", EDep " << KaonEDep[KL] 
		  << ", DaughtEDep " << KaonDaughtersEDep[KL] << ", DecayEDep " << KaonDecayEDep[KL] << ", NearEDep " << KaonNearEDep[KL] << std::endl;
      }
    }
    if (nElec) {
      std::cout << "There are electrons in this event!!" << std::endl;
      for (int EL=0; EL<nElec; ++EL) {
	std::cout << "Elec " << EL << " had properties: PDG " << ElecPDG[EL] << ", TrackID " << ElecTrID[EL] << ", EDep " << ElecEDep[EL] 
		  << ", DaughtEDep " << ElecDaughtersEDep[EL] << ", DecayEDep " << ElecDecayEDep[EL] << ", NearEDep " << ElecNearEDep[EL] 
		  << "\nStart pos " << ElecStart[EL][0] << ", " << ElecStart[EL][1] << ", " << ElecStart[EL][2]
		  << ". + End pos " << ElecEnd  [EL][0] << ", " << ElecEnd  [EL][1] << ", " << ElecEnd  [EL][2]
		  << "\nThe parents of this electron were: ";
	for (int zz=0; zz<nParentElec[EL]; ++zz) {
	  std::cout << ElecParents[EL][zz] << " ("<<ElecParTrID[EL][zz]<<"), ";
	}
	std::cout << "."<< std::endl;
      }
    }
    if (nMuon) {
      std::cout << "There are " << nMuon << " muons in this event!!" << std::endl;
    }
    if (nPion) {
      std::cout << "There are " << nPion << " pions in this event!!" << std::endl;
    }
    if (nPi0) {
      std::cout << "There are " << nPi0  << " pi0's in this event!!" << std::endl;
    }
    if (nProt) {
      std::cout << "There are " << nProt << " protons in this event!!" << std::endl;
    }

    std::cout << "\nThere was " << EDepNearEdge2 << ", (" << EDepNearEdge5 << "), [" << EDepNearEdge10 << "] MeV EDep within 2, (5), [10] cm of walls" << std::endl;
  }
  // ------------------------------------- Output some information about my particles -------------------------------------
  
  // ------ Fill the Tree ------
  if (nKaon) {
    if (Verbosity)
      std::cout << "There were Kaons in the detector so writing out this event." << std::endl;
    fDecayTree->Fill();
  }
  return;
} // Analyse
// ******************************** Fill variables *****************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::FillVars( std::vector<int> &TrackIDVec, int &numParts, float EDep[MaxPart], float DaughtEDep[MaxPart], 
						     float DecayEDep[MaxPart], float Start[MaxPart][4], float End[MaxPart][4],
						     int nParents[MaxPart], int Parent[MaxPart][MaxParent], int ParTrID[MaxPart][MaxParent],
						     int PDG[MaxPart], int TrID[MaxPart], int Cont[MaxPart], int FromDecay[MaxPart],
						     int ThisID, unsigned int ThisTDC, sim::IDE ThisIDE, const simb::MCParticle& MCPart, 
						     bool OrigParticle, bool Decay, bool &Written ) {
  if (numParts+1 > MaxPart-1) {
    BadEvent = true;
    Written=true;
    return;
  }

  std::vector<int>::iterator it=std::find( TrackIDVec.begin(), TrackIDVec.end(), abs(ThisID) );
  int partNum = 0;
  // ------- Figure out which partNum I want to use --------
  if ( it==TrackIDVec.end() ) {
    TrackIDVec.push_back ( abs(ThisID) ); // Push back this particle type IDVec
    AllTrackIDs.push_back( abs(ThisID) ); // Push back all particle type IDVec
    PDG[numParts]  = MCPart.PdgCode();
    TrID[numParts] = abs( MCPart.TrackId() );
    if ( MCPart.Process() == "Decay" ) FromDecay[numParts] = 1;
    else FromDecay[numParts] = 0;
    // Work out if contained within TPC...
    bool StartIn = IsInTPC( MCPart.Vx(0) , MCPart.Vy(0) , MCPart.Vz(0) , ActiveBounds );
    bool EndIn   = IsInTPC( MCPart.EndX(), MCPart.EndY(), MCPart.EndZ(), ActiveBounds );
    if(!StartIn && !EndIn ) Cont[numParts] = 0;  // through track
    if(!StartIn &&  EndIn ) Cont[numParts] = 1;  // entering track
    if( StartIn && !EndIn ) Cont[numParts] = 2;  // escaping track
    if( StartIn &&  EndIn ) Cont[numParts] = 3;  // contained track
    if (Verbosity)
      std::cout << "\nPushing back a new ideTrackID " << abs(ThisID) << ", it was from a " << MCPart.PdgCode() << ", process " << MCPart.Process() 
		<< ", PartDecay? " << FromDecay[numParts] << ", deposition from a decay? " << Decay << "\n"
		<< "Starting location is " << MCPart.Vx(0) << ", " << MCPart.Vy(0) << ", " << MCPart.Vz(0)  << "==> InTPC? " << StartIn
		<< ". Ending location is " << MCPart.EndX()<< ", " << MCPart.EndY()<< ", " << MCPart.EndZ() << "==> InTPC? " << EndIn
		<< " =====>>> Cont? " << Cont[numParts]
		<< std::endl;
    // ---- Work out the particles ancestry ----
    Parent[numParts][0] = MCPart.Mother();
    int NumParent = 0;
    int ParentID  = Parent[numParts][0];
    while ( ParentID > 0 && NumParent < MaxParent) {
      int pdg    = truthmap[ ParentID ]->PdgCode();
      if (Verbosity==2)
	std::cout << "The particles parent was TrackID " << ParentID << ", PdgCode " << pdg << ", it was from process " << truthmap[ ParentID ]->Process() << std::endl;
      // If this particle had same parent as previously, then put this TrackID in my array instead
      if ( pdg == Parent[numParts][NumParent-1] ) {
	if (Verbosity==2)
	  std::cout << "Overwriting what was in ParentID as same particle." << std::endl;
	ParTrID[numParts][NumParent-1] = ParentID;
      } else {
	Parent [numParts][NumParent] = pdg;
	ParTrID[numParts][NumParent] = ParentID;
	++NumParent;
      }
      ParentID = truthmap[ ParentID ]->Mother();
    }
    
    if (Verbosity)
      std::cout << "There were a total of " << NumParent << " parents of this particle. The parentage was: " << PDG[numParts] << " ("<< TrID[numParts] <<") ";
    for (int aa=0; aa<NumParent; ++aa)
      if (Verbosity)
	std::cout << " <- " << Parent[numParts][aa] << " ("<<ParTrID[numParts][aa] << ") ";
    if (Verbosity) std::cout << "." << std::endl;
    
    nParents[numParts] = NumParent;
    partNum = numParts;
    ++numParts;
  } else {  
    partNum = it - TrackIDVec.begin();
  }
  
  // ----------- If TrackID is positive add it to the EDep from this track, if negative then it has to be from a daughter --------------
  if ( Verbosity > 1 ) {
    std::cout << "OrigParticle " << OrigParticle << ", TrackID " << ThisID << ", " << ThisIDE.trackID 
	      << ", PDgCode " << PDG[numParts-1] << ", " << MCPart.PdgCode() << ", decay " << Decay 
	      << std::endl;
  }
  if ( OrigParticle ) {
    EDep[partNum]       += ThisIDE.energy;
    
    bool RepSt = false;
    if (Start[partNum][1] == HolderVal)
      RepSt = true;
    else {
      float St_Ar = CalcDist( MCPart.Vx(0), MCPart.Vy(0), MCPart.Vz(0), Start[partNum][0], Start[partNum][1], Start[partNum][2] );
      float St_ID = CalcDist( MCPart.Vx(0), MCPart.Vy(0), MCPart.Vz(0), ThisIDE.x        , ThisIDE.y        , ThisIDE.z         );
      if ( St_ID < St_Ar )
	RepSt = true;
    }
    if (RepSt) {
      Start[partNum][0] = ThisIDE.x;
      Start[partNum][1] = ThisIDE.y;
      Start[partNum][2] = ThisIDE.z;
      Start[partNum][3] = ThisTDC;
    }
    // ------- Do I want to replace the end array?
    bool RepEn = false;
    if (End[partNum][1] == HolderVal)
      RepEn = true;
    else {
      float En_Ar = CalcDist( MCPart.EndX(), MCPart.EndY(), MCPart.EndZ(), End[partNum][0], End[partNum][1], End[partNum][2] );
      float En_ID = CalcDist( MCPart.EndX(), MCPart.EndY(), MCPart.EndZ(), ThisIDE.x      , ThisIDE.y      , ThisIDE.z       );
      if ( En_ID < En_Ar )
	RepEn = true;
    }
    if (RepEn) {
      End[partNum][0] = ThisIDE.x;
      End[partNum][1] = ThisIDE.y;
      End[partNum][2] = ThisIDE.z;
      End[partNum][3] = ThisTDC;
    }
    // --------- Work out the start / end of this track ----------
  } else if (Decay) {
    DecayEDep[partNum]  += ThisIDE.energy;
  } else {
    DaughtEDep[partNum] += ThisIDE.energy;
  }
  
  Written=true;
} // FillVars
// ******************************** UnAssignLoop *****************************************************
bool NeutronDecayN2Ana::NeutronDecayN2Ana::UnAssignLoop( int nPart, float St[MaxPart][4], sim::IDE thisIDE, float NearE[MaxPart] ) {
  for (int nP=0; nP < nPart; ++nP) {
    float ThisDist = CalcDist( St[nP][0], St[nP][1], St[nP][2], thisIDE.x, thisIDE.y, thisIDE.z );
    if (Verbosity > 1) {
      std::cout << "  Particle " << nP << " had initial position ("<<St[nP][0]<<", "<<St[nP][1]<<", "<<St[nP][2] <<")...How close was this to the IDE in question? " << ThisDist << std::endl;
    }
    if (ThisDist < NearEDeps) {
      NearE[nP] += thisIDE.energy;
      return true;
    }
  }
  return false;
} // UnAssignLoop
// ******************************** Define Module *****************************************************
float NeutronDecayN2Ana::NeutronDecayN2Ana::CalcDist( float X1, float Y1, float Z1, float X2, float Y2, float Z2 ) {
  float Sq = TMath::Power( X1-X2 , 2 ) + TMath::Power( Y1-Y2 , 2 ) + TMath::Power( Z1-Z2 , 2 );
  float Va = TMath::Power( Sq, 0.5 );
  return Va;
}
// ******************************** Define Module *****************************************************
bool  NeutronDecayN2Ana::NeutronDecayN2Ana::IsInTPC( float X1, float Y1, float Z1, float Bound[6] ) {
  bool InX = ( X1 > Bound[0] && X1 < Bound[1] );
  bool InY = ( Y1 > Bound[2] && Y1 < Bound[3] );
  bool InZ = ( Z1 > Bound[4] && Z1 < Bound[5] );
  bool Within = (InX && InY && InZ );
  /*
    std::cout << "Is ("<<X1<<", "<<Y1<<", "<<Z1<<") within the TPC? InX " << InX << ", InY " << InY << ", InZ " << InZ << ", within = " << Within
	    << "  X "<<Bound[0]<<" to "<<Bound[1]
	    << ", X "<<Bound[2]<<" to "<<Bound[3]
	    << ", X "<<Bound[4]<<" to "<<Bound[5]
	    <<std::endl;
  */
  return Within;
}
// ******************************** Define Module *****************************************************
DEFINE_ART_MODULE(NeutronDecayN2Ana::NeutronDecayN2Ana)
