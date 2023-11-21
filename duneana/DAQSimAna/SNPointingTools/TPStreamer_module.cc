////////////////////////////////////////////////////////////////////////
// Class:       TPStreamer
// Module Type: analyzer
// File:        TPStreamer_module.cc
// Author:      K. Wawrowska
// Date:        October 2023
//
// Run custom hit finder and dump the three-view hit information 
// to a txt file with labels corresponding to the MC producer
////////////////////////////////////////////////////////////////////////


#include <fstream>

// Framework includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


enum PType{ kUnknown=0, kGen, kBgd}; //0 == noise hit, 1 == signal hit, 2 == radiological/background hit

class TPStreamer : public art::EDAnalyzer {

public:

  explicit TPStreamer(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  TPStreamer(TPStreamer const &) = delete;
  TPStreamer(TPStreamer &&) = delete;
  TPStreamer & operator = (TPStreamer const &) = delete;
  TPStreamer & operator = (TPStreamer &&) = delete;

  void analyze(art::Event const & evt) override;

private:

  //--- custom functions 
  void ResetVariables(); 
  PType WhichParType(int TrkID); // Particle type from track ID 
  void FillMyMaps( std::map< int, simb::MCParticle> &MyMap,
			art::FindManyP<simb::MCParticle> Assn, 
			art::Handle< std::vector<simb::MCTruth> > Handle,  
			std::map<int, int>* indexMap=nullptr);

  //--- producer labels
  std::string m_GenLabel; //generator used for signal events 
  std::string m_GeantLabel; //g4
  std::string m_HitLabel; // which hit finder to use 

  //--- time window around a hit for tagging (useful if filter algorithms elongate hit TOTs)
  raw::TDCtick_t m_HitWindow;  // in ticks 
  bool m_AbsRS; //needed for separate hit tagging to account for modified hit profile

  //--- files where we dump our collection and induction TPs
  std::string m_outputFilename;
  std::ofstream m_outputFile; 

  std::map< int, simb::MCParticle> GenParts;
  std::map< int, simb::MCParticle> BgdParts;

  //Mapping from track ID to particle type, for use in WhichParType() 
  std::map<int, PType> trkIDToPType; 
  std::vector<int> Hit_True_MainTrID;

  // --- Declare our services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  
}; // class TPStreamer

//......................................................
TPStreamer::TPStreamer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p), 
  m_GenLabel(      p.get<std::string>("GenLabel")),
  m_GeantLabel(    p.get<std::string>("GEANT4Label")),
  m_HitLabel(      p.get<std::string>("HitLabel")),
  m_HitWindow(     p.get<raw::TDCtick_t>("HitWindow", 20)), 
  m_AbsRS(         p.get<bool>("AbsRSHits", false)),
  m_outputFilename(p.get<std::string>("OutputFile")),
  m_outputFile(m_outputFilename)
{
  
} // TPStreamer constructor

void TPStreamer::ResetVariables()
{
  GenParts.clear();
  BgdParts.clear();
  trkIDToPType.clear(); 
  Hit_True_MainTrID.clear();

} // ResetVariables

void TPStreamer::FillMyMaps(std::map<int, simb::MCParticle> &MyMap,
				    art::FindManyP<simb::MCParticle> Assn, 
				    art::Handle< std::vector<simb::MCTruth> > Handle, 
				    std::map<int,int>* indexMap)
{
  for ( size_t L1=0; L1 < Handle->size(); ++L1 ) {
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      const simb::MCParticle ThisParticle = (*Assn.at(L1).at(L2));
      MyMap[ThisParticle.TrackId()] = ThisParticle;
      if(indexMap) indexMap->insert({ThisParticle.TrackId(), L1});
    }
  }
  return;

} // FillMyMaps

PType TPStreamer::WhichParType (int TrkID)
{
  PType ThisPType= kUnknown;
  auto const& it = trkIDToPType.find(TrkID);
  if(it!=trkIDToPType.end()){   ThisPType=it->second;  }
  return ThisPType; 
} // WhichParType

//......................................................
void TPStreamer::analyze(art::Event const & evt)
{

  ResetVariables(); 

  // --- First want to associate all the g4 tracks to their respective generators
  //--  Store the particles and their corresponding g4 track IDs in maps

  auto GenHandles = evt.getMany<std::vector<simb::MCTruth>>();
  for (auto const& handle: GenHandles){

    //Get the gen module labels for signal + radiologicals
    const std::string GenModuleLabel = handle.provenance()->moduleLabel();

    //Get a map between G4 Track IDs and signal MC parts. 
    if (GenModuleLabel == m_GenLabel){
      
      auto GenTrue = evt.getHandle< std::vector<simb::MCTruth> >(m_GenLabel);
      if (GenTrue){	
        art::FindManyP<simb::MCParticle> GenAssn( GenTrue, evt, m_GeantLabel); 
	      FillMyMaps( GenParts, GenAssn, GenTrue); 
      }
    }
    //Get a map between G4 Track IDs and bgd MC parts.
    else{
      auto BgdTrue = evt.getHandle< std::vector<simb::MCTruth> >(GenModuleLabel);   
      if (BgdTrue){                                                                                                           
	      art::FindManyP<simb::MCParticle> BgdAssn(BgdTrue, evt, m_GeantLabel);      

        //Create a temporary map for the specific background source 
        std::map<int, simb::MCParticle> tempBgdMap; 
        FillMyMaps(tempBgdMap, BgdAssn, BgdTrue);

        //Merge the temporary map with the full backgrounds map
        BgdParts.insert(tempBgdMap.begin(), tempBgdMap.end());                                                         
      }
    }
  }


  // Get any G4 tracks that weren't matched to signal or radio gen
  // and add them to signal map if they were most likely daughter particles (e.g photons from EM shower)
  
  art::ValidHandle<std::vector <simb::MCParticle> > mcParticles = evt.getValidHandle<std::vector <simb::MCParticle> >(m_GeantLabel);
  if (mcParticles.isValid()){
    std::map< int, simb::MCParticle> DaughterParts;

    for (unsigned int i = 0; i < mcParticles->size(); ++i){
      const simb::MCParticle trueParticle = mcParticles->at(i);
      // if (trueParticle.Mother() != 0){ 	DaughterParts[trueParticle.TrackId()] = trueParticle;      }
     
      // Check if the TrackId() is not in GenParts or BgdParts
      if (GenParts.find(trueParticle.TrackId()) == GenParts.end() && BgdParts.find(trueParticle.TrackId()) == BgdParts.end()){
        DaughterParts[trueParticle.TrackId()] = trueParticle;
      }
    }
    //Add daughter particles to signal map
    GenParts.insert(DaughterParts.begin(), DaughterParts.end()); 
  }
 
  // ---



  //Map each particle map to its corresponding enum tag 
  std::map<PType, std::map< int, simb::MCParticle >&> PTypeToMap{
    {kGen,  GenParts },
    {kBgd,  BgdParts }
  };
  
  //run over the particle assn map
  for (auto const& it : PTypeToMap){

    //particle tag e.g. kGen
    const PType p = it.first;
    // gen-g4 mapping e.g. GenParts
    auto const& m = it.second;

    //run over each row in e.g. GenParts
    for (auto const& it2 : m){
      //add a row to the trkIDToPType map consisting of [particle trk ID, kGen] 
      //trkIDToPType is a 2xn matrix
      trkIDToPType.insert( std::make_pair(it2.first, p));
    }
  }

  // --- Lift out the reco hits:
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(m_HitLabel);
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);


  
  for(size_t hit = 0; hit < reco_hits->size(); ++hit) {
    recob::Hit const& ThisHit = reco_hits->at(hit);   // current hit 

    //time window in which hits get mapped to IDEs
    raw::TDCtick_t WindowStart = ThisHit.StartTick() - m_HitWindow;
    raw::TDCtick_t WindowEnd = ThisHit.EndTick() + m_HitWindow;

    //Use separate window for AbsRS induction hit to account for the hits being elongated/merged 
    //PeakT typically corresponds to inflexion point -- i.e. hit EndT for unfiltered hit. 
    if ((ThisHit.View() != 2) && (m_AbsRS==true)){ WindowEnd = ThisHit.PeakTime() +m_HitWindow;  }
    
    //--- Ionization drift electrons (IDEs) associated with current hit
    std::vector<sim::TrackIDE> ThisHitIDE = bt_serv->ChannelToTrackIDEs(clockData, 
									ThisHit.Channel(), 
								        WindowStart,
									WindowEnd);

    
    //---Get the G4 track associated to the IDEs 
    double TopEFrac = -DBL_MAX;
    Hit_True_MainTrID.push_back(-1);     //in the case of noise hit when there's no track associated with the hit 

    if (ThisHitIDE.size()){
      for (size_t ideL=0; ideL < ThisHitIDE.size(); ++ideL){
        if (ThisHitIDE[ideL].energyFrac > TopEFrac){
          TopEFrac = ThisHitIDE[ideL].energyFrac;
          Hit_True_MainTrID.at(hit) = std::abs( ThisHitIDE[ideL].trackID );
        }
      }
    }
    PType ThisPType = WhichParType( Hit_True_MainTrID.at(hit));
  
    //dump the hits to a file with one hit per line in the following format:
    //event number, plane, start time, end time, peak time, TOT, channel, SADC, Peak ADC, MC producer

    // to match DAQ (just for simplicity), format needs to be
    // ([tp.time_start, tp.time_over_threshold, tp.time_peak, tp.channel, tp.adc_integral, tp.adc_peak, tp.detid, tp.flag, tp.version]) + truth

    // I want numbers to be printed out always as integers, not floats
    m_outputFile << std::fixed << std::setprecision(0) << std::setw(0) << std::setfill('0'); // TODO alignment
    m_outputFile
     << ThisHit.StartTick() << ' '
     << ThisHit.EndTick() - ThisHit.StartTick() << ' ' // TOT
		 << ThisHit.PeakTime() << ' ' 
     << ThisHit.Channel() << ' '
		 << ThisHit.SummedADC() << ' '  // integral
     << ThisHit.PeakAmplitude() << ' '
     << (int) ThisHit.Channel()/2560 << ' ' // detid
     << 1 << ' ' // type
     << 1 << ' ' // algorithm, currently kTPCDefault
     << 1 << ' ' // version
     << 0 << ' ' // flag
     // now MC trith
     << ThisPType << ' ' 
     << evt.event() << ' '
     << ThisHit.View() << ' '
     << 0 << ' ' // energy, to be implemented in BT
     << 0 << ' ' // nElectrons, to be implemented in BT
     << 0 << ' ' // trackId, to be implemented in BT
     << 0 << ' ' // trueX, to be implemented in BT
     << 0 << ' ' // trueY, to be implemented in BT
     << 0 << ' ' // trueZ, to be implemented in BT
     <<  std::endl; 

  } // Loop over reco_hits.
} // Analyze TPStreamer.


//......................................................
DEFINE_ART_MODULE(TPStreamer)
