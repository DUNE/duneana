////////////////////////////////////////////////////////////////////////
// Class:       NucleonDecayFilter
// Module Type: analyzer
// File:        NucleonDecayFilter_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
//
//  Making a filter module to reduce the number of events which nucleon
//   decay analyses have to look at.
//  For example, if there are no Kaon energy deposits, then shouldn't 
//   reconstruct a Kaon in the event!
//
//  Thomas Karl Warburton
//  k.warburton@sheffield.ac.uk
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/MCCheater/BackTracker.h"
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

#define MaxPart   100
#define MaxParent 100
#define HolderVal 9999

namespace filt{
  class NucleonDecayFilter : public art::EDFilter {
  public:
    explicit NucleonDecayFilter(fhicl::ParameterSet const & pset);
    virtual ~NucleonDecayFilter() {};
    virtual bool filter(art::Event& e);
    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob() ;

  private:
    // My functions
    void FillVars( std::vector<int> &TrackIDVec, int &numParts, int ThisID );

    // Global variables
    double ActiveBounds[6]; // Cryostat boundaries ( neg x, pos x, neg y, pos y, neg z, pos z )
    std::map<int, const simb::MCParticle*> truthmap; // A map of the truth particles.

    // Handles
    art::ServiceHandle<geo::Geometry> geom;
    geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

    // FHiCL parameters
    bool fWantMuon;
    bool fWantPion;
    bool fWantPi0;
    bool fWantKaon;
    bool fWantElec;

    // Tree variables
    TTree *MyOutTree;
    int Run;
    int SubRun;
    int Event;
    int nMuon, nPion, nPi0, nKaon, nElec;
     
  };
  // *************************************************************************************
  NucleonDecayFilter::NucleonDecayFilter::NucleonDecayFilter(fhicl::ParameterSet const & pset) : EDFilter{pset} {
    this->reconfigure(pset);
  }
  // ************************************ Reconfigure ************************************
  void NucleonDecayFilter::reconfigure(fhicl::ParameterSet const& pset) {
    fWantMuon = pset.get<bool>("WantMuon");
    fWantPion = pset.get<bool>("WantPion");
    fWantPi0  = pset.get<bool>("WantPi0" );
    fWantKaon = pset.get<bool>("WantKaon");
    fWantElec = pset.get<bool>("WantElec");
    std::cout << "Your filter is selecting only events with: "
	      << "Muons " << fWantMuon << ", "
	      << "Pions " << fWantPion << ", "
	      << "Pi0s "  << fWantPi0  << ", "
	      << "Kaons " << fWantKaon << ", "
	      << "Elecs " << fWantElec << "."
	      << std::endl;
  }
  // ************************************ Filter Func ************************************
  bool NucleonDecayFilter::filter(art::Event & e) {
    Run    = e.run();
    SubRun = e.subRun();
    Event  = e.event();
    nMuon  = nPion = nPi0 = nKaon = nElec = 0;

    // Any providers I need.
    auto const* geo = lar::providerFrom<geo::Geometry>();
    
    // Make a map of MCParticles which I can access later.
    auto truth = e.getHandle<std::vector<simb::MCParticle> >("largeant");
    truthmap.clear();
    for (size_t i=0; i<truth->size(); i++)
      truthmap[truth->at(i).TrackId()]=&((*truth)[i]);
    //std::cout << "The primary muon in this event had an energy of " << truthmap[1]->E() << std::endl;

    // Get a vector of sim channels.
    auto simchannels = e.getHandle<std::vector<sim::SimChannel> >("largeant");

    // Make vectors to hold all of my particle TrackIDs
    std::vector<int> MuonVec, PionVec, Pi0Vec, KaonVec, ElecVec;
    
    // ------ Now loop through all of my hits ------
    for (auto const& simchannel:*simchannels) {
    // ------ Only want to look at collection plane hits ------
      if (wireReadout.SignalType(simchannel.Channel()) != geo::kCollection) continue;
      // ------ Loop through all the IDEs for this channel ------
      for (auto const& tdcide:simchannel.TDCIDEMap()) {
	auto const& idevec=tdcide.second;
	// ------ Look at each individual IDE ------
	for (auto const& ide:idevec) {
	  int ideTrackID = ide.trackID;
	  // ------ Is the hit in a TPC? ------
          geo::TPCID tpcid=geo->FindTPCAtPosition(geo::Point_t{ide.x, ide.y, ide.z});
	  if (!(geo->HasTPC(tpcid)) ) {
	    //std::cout << "Outside the Active volume I found at the top!" << std::endl;
	    continue;
	  }
	  // ------ Now to work out which particles in particular I want to save more information about... ------
	  // ----------------- I want to write out the IDE to the relevant parent particle type -----------------
	  // ------- This means looping back through the IDE's parents until I find an interesting TrackID ------
	  bool WrittenOut = false;
	  while ( ideTrackID != 0 && !WrittenOut ) {
	    const simb::MCParticle& part=*( truthmap[ abs(ideTrackID) ] );
	    int PdgCode=part.PdgCode();
	    WrittenOut = true; // If one of my particles don't want to go further back.
	    // ========== Muons ==========
	    if      ( (PdgCode == -13  || PdgCode == 13)  ) 
	      FillVars( MuonVec, nMuon, ideTrackID );
	    // ========== Pions ==========
	    else if ( (PdgCode == -211 || PdgCode == 211) && part.Process() != "pi+Inelastic" && part.Process() != "pi-Inelastic" ) 
	      FillVars( PionVec, nPion, ideTrackID );
	    // ========== Pi0s  ==========
	    else if ( PdgCode == 111  )
	      FillVars( Pi0Vec , nPi0 , ideTrackID );
	    // ========== Kaons ===========
	    else if ( (PdgCode == 321 || PdgCode == -321) && part.Process() != "kaon+Inelastic" && part.Process() != "kaon-Inelastic" ) 
	      FillVars( KaonVec, nKaon, ideTrackID );
	    // ========== Elecs ===========
	    else if ( (PdgCode == -11 || PdgCode == 11)  )
	      FillVars( ElecVec, nElec, ideTrackID );
	    // ========== If still not one of my intersting particles I need to find this particles parent ==========
	    else {
	      ideTrackID = part.Mother();
	      PdgCode = truthmap[ abs(ideTrackID) ]->PdgCode();
	      WrittenOut = false; // If not one of my particles change back to false.
	    } // If not one of chosen particles.
	  } // While loop
	} // Each IDE ( ide:idevec )
      } // IDE vector for SimChannel ( tdcide:simcahnnel.TPCIDEMap() )
    } // Loop through simchannels
    /*
    std::cout << "Looking at Run " << Run << " SubRun " << SubRun << " Event " << Event
	      << ". There were " << nMuon << " muons, " << nPion << " pions, " << nPi0 << " pi0s, " << nKaon << " Kaons, " << nElec << " Electrons."
	      << std::endl;
    */
    // Fill Tree and decide whether to keep this event...
    MyOutTree->Fill();
    bool EventPasses = true;
    if ( (fWantMuon && nMuon == 0) || nMuon > MaxPart) EventPasses = false;
    if ( (fWantPion && nPion == 0) || nPi0  > MaxPart) EventPasses = false;
    if ( (fWantPi0  && nPi0  == 0) || nPion > MaxPart) EventPasses = false;
    if ( (fWantKaon && nKaon == 0) || nKaon > MaxPart) EventPasses = false;
    if ( (fWantElec && nElec == 0) || nElec > MaxPart) EventPasses = false;
    return EventPasses;
  }
  // ************************************* Begin Job *************************************
  void NucleonDecayFilter::beginJob() {
    ActiveBounds[0] = ActiveBounds[2] = ActiveBounds[4] = DBL_MAX;
    ActiveBounds[1] = ActiveBounds[3] = ActiveBounds[5] = -DBL_MAX;
    // ----- FixMe: Assume single cryostats ------
    auto const* geom = lar::providerFrom<geo::Geometry>();
    for (geo::TPCGeo const& TPC: geom->Iterate<geo::TPCGeo>()) {
      // get center in world coordinates
      auto const center = TPC.GetCenter();
      double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
      if( center.X() - tpcDim[0] < ActiveBounds[0] ) ActiveBounds[0] = center.X() - tpcDim[0];
      if( center.X() + tpcDim[0] > ActiveBounds[1] ) ActiveBounds[1] = center.X() + tpcDim[0];
      if( center.Y() - tpcDim[1] < ActiveBounds[2] ) ActiveBounds[2] = center.Y() - tpcDim[1];
      if( center.Y() + tpcDim[1] > ActiveBounds[3] ) ActiveBounds[3] = center.Y() + tpcDim[1];
      if( center.Z() - tpcDim[2] < ActiveBounds[4] ) ActiveBounds[4] = center.Z() - tpcDim[2];
      if( center.Z() + tpcDim[2] > ActiveBounds[5] ) ActiveBounds[5] = center.Z() + tpcDim[2];
    } // for all TPC
    std::cout << "Active Boundaries: "
	      << "\n\tx: " << ActiveBounds[0] << " to " << ActiveBounds[1]
	      << "\n\ty: " << ActiveBounds[2] << " to " << ActiveBounds[3]
	      << "\n\tz: " << ActiveBounds[4] << " to " << ActiveBounds[5]
	      << std::endl;
  
    art::ServiceHandle<art::TFileService> tfs;
    MyOutTree = tfs->make<TTree>("FilterTree","FilterTree");
    MyOutTree->Branch( "Run",   &Run   , "Run/I"    );
    MyOutTree->Branch( "SubRun",&SubRun, "SubRun/I" );
    MyOutTree->Branch( "Event", &Event , "Event/I"  );
    MyOutTree->Branch( "nMuon", &nMuon , "nMuon/I"  );
    MyOutTree->Branch( "nPion", &nPion , "nPion/I"  );
    MyOutTree->Branch( "nPi0" , &nPi0  , "nPi0/I"   );
    MyOutTree->Branch( "nKaon", &nKaon , "nKaon/I"  );
    MyOutTree->Branch( "nElec", &nElec , "nElec/I"  );
  }
  // *********************************** Fill Variable ***********************************
  void NucleonDecayFilter::FillVars( std::vector<int> &TrackIDVec, int &numParts, int ThisID ) {
    std::vector<int>::iterator it=std::find( TrackIDVec.begin(), TrackIDVec.end(), abs(ThisID) );
    if ( it==TrackIDVec.end() ) {
      TrackIDVec.push_back ( abs(ThisID) ); // Push back this particle type IDVec
      ++numParts;
    }
  }
  // *********************************** Define Module ***********************************
  DEFINE_ART_MODULE(NucleonDecayFilter)
}
