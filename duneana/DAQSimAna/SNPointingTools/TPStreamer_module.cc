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

// adding all of the different backgrounds to have a complete set, can group them in the analysis
enum PType{ 
  kUnknown=0,
  kMarley=1,
  kAr39GenInLAr=2,
  kKr85GenInLAr=3,
  kAr42GenInLAr=4,
  kK42From42ArGenInLAr=5,
  kRn222ChainRn222GenInLAr=6,
  kRn222ChainPo218GenInLAr=7,
  kRn222ChainPb214GenInLAr=8,
  kRn222ChainBi214GenInLAr=9,
  kRn222ChainPb210GenInLAr=10,
  kK40GenInCPA=11,
  kU238ChainGenInCPA=12,
  kK42From42ArGenInCPA=13,
  kRn222ChainPo218GenInCPA=14,
  kRn222ChainPb214GenInCPA=15,
  kRn222ChainBi214GenInCPA=16,
  kRn222ChainPb210GenInCPA=17,
  kRn222ChainFromBi210GenInCPA=18,
  kCo60GenInAPA=19,
  kU238ChainGenInAPA=20,
  kRn222ChainGenInPDS=21,
  kNeutronGenInRock=22
};

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
  std::string m_MarleyLabel; //generator used for signal events 
  // add all background variables at once please
  std::string m_Ar39GenInLArLabel;
  std::string m_Kr85GenInLArLabel;
  std::string m_Ar42GenInLArLabel;
  std::string m_K42From42ArGenInLArLabel;
  std::string m_Rn222ChainRn222GenInLArLabel;
  std::string m_Rn222ChainPo218GenInLArLabel;
  std::string m_Rn222ChainPb214GenInLArLabel;
  std::string m_Rn222ChainBi214GenInLArLabel;
  std::string m_Rn222ChainPb210GenInLArLabel;
  std::string m_K40GenInCPALabel;
  std::string m_U238ChainGenInCPALabel;
  std::string m_K42From42ArGenInCPALabel;
  std::string m_Rn222ChainPo218GenInCPALabel;
  std::string m_Rn222ChainPb214GenInCPALabel;
  std::string m_Rn222ChainBi214GenInCPALabel;
  std::string m_Rn222ChainPb210GenInCPALabel;
  std::string m_Rn222ChainFromBi210GenInCPALabel;
  std::string m_Co60GenInAPALabel;
  std::string m_U238ChainGenInAPALabel;
  std::string m_Rn222ChainGenInPDSLabel;
  std::string m_NeutronGenInRockLabel;

  std::string m_GeantLabel; //g4
  std::string m_HitLabel; // which hit finder to use 

  //--- time window around a hit for tagging (useful if filter algorithms elongate hit TOTs)
  raw::TDCtick_t m_HitWindow;  // in ticks 
  bool m_AbsRS; //needed for separate hit tagging to account for modified hit profile

  //--- files where we dump our collection and induction TPs
  std::string m_outputFilename;
  std::ofstream m_outputFile; 

  std::map< int, simb::MCParticle> MarleyPMap;
  std::map< int, simb::MCParticle> BgdPMap;

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
  // gen labels to retrieve MCtruth
  m_MarleyLabel(                      p.get<std::string>("GenLabels.MarleyLabel")),
  m_Ar39GenInLArLabel(                p.get<std::string>("GenLabels.Ar39GenInLarLabel")),
  m_Kr85GenInLArLabel(                p.get<std::string>("GenLabels.Kr85GenInLarLabel")),
  m_Ar42GenInLArLabel(                p.get<std::string>("GenLabels.Ar42GenInLarLabel")),
  m_K42From42ArGenInLArLabel(         p.get<std::string>("GenLabels.K42From42ArGenInLarLabel")),
  m_Rn222ChainRn222GenInLArLabel(     p.get<std::string>("GenLabels.Rn222ChainRn222GenInLarLabel")),
  m_Rn222ChainPo218GenInLArLabel(     p.get<std::string>("GenLabels.Rn222ChainPo218GenInLarLabel")),
  m_Rn222ChainPb214GenInLArLabel(     p.get<std::string>("GenLabels.Rn222ChainPb214GenInLarLabel")),
  m_Rn222ChainBi214GenInLArLabel(     p.get<std::string>("GenLabels.Rn222ChainBi214GenInLarLabel")),
  m_Rn222ChainPb210GenInLArLabel(     p.get<std::string>("GenLabels.Rn222ChainPb210GenInLarLabel")),
  m_K40GenInCPALabel(                 p.get<std::string>("GenLabels.K40GenInCPALabel")),
  m_U238ChainGenInCPALabel(           p.get<std::string>("GenLabels.U238ChainGenInCPALabel")),
  m_K42From42ArGenInCPALabel(         p.get<std::string>("GenLabels.K42From42ArGenInCPALabel")),
  m_Rn222ChainPo218GenInCPALabel(     p.get<std::string>("GenLabels.Rn222ChainPo218GenInCPALabel")),
  m_Rn222ChainPb214GenInCPALabel(     p.get<std::string>("GenLabels.Rn222ChainPb214GenInCPALabel")),
  m_Rn222ChainBi214GenInCPALabel(     p.get<std::string>("GenLabels.Rn222ChainBi214GenInCPALabel")),
  m_Rn222ChainPb210GenInCPALabel(     p.get<std::string>("GenLabels.Rn222ChainPb210GenInCPALabel")),
  m_Rn222ChainFromBi210GenInCPALabel( p.get<std::string>("GenLabels.Rn222ChainFromBi210GenInCPALabel")),
  m_Co60GenInAPALabel(                p.get<std::string>("GenLabels.Co60GenInAPALabel")),
  m_U238ChainGenInAPALabel(           p.get<std::string>("GenLabels.U238ChainGenInAPALabel")),
  m_Rn222ChainGenInPDSLabel(          p.get<std::string>("GenLabels.Rn222ChainGenInPDSLabel")),
  m_NeutronGenInRockLabel(            p.get<std::string>("GenLabels.NeutronGenInRockLabel")),

  m_GeantLabel(    p.get<std::string>("GEANT4Label")),
  m_HitLabel(      p.get<std::string>("HitLabel")),
  m_HitWindow(     p.get<raw::TDCtick_t>("HitWindow", 20)), 
  m_AbsRS(         p.get<bool>("AbsRSHits", false)),
  m_outputFilename(p.get<std::string>("OutputFile")), // TODO add here reading of threshold to compose filename
  m_outputFile(m_outputFilename)
{
  
} // TPStreamer constructor

void TPStreamer::ResetVariables()
{
  MarleyPMap.clear();
  BgdPMap.clear();
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
  // create a quick output textfile called "trkIDToPType.txt" to check the mapping
  std::ofstream trkIDToPTypeFile("trkIDToPType.txt", std::ios_base::app);
  trkIDToPTypeFile << TrkID << " " << it->second << std::endl;
  trkIDToPTypeFile.close();
  // TODO remove this

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
    if (GenModuleLabel == m_MarleyLabel){
      
      auto GenTrue = evt.getHandle< std::vector<simb::MCTruth> >(m_MarleyLabel);
      if (GenTrue){	
        art::FindManyP<simb::MCParticle> GenAssn( GenTrue, evt, m_GeantLabel); 
	      FillMyMaps( MarleyPMap, GenAssn, GenTrue); 
      }
    }
    //Get a map between G4 Track IDs and bgd MC parts.
    else{
      // in this way, you only check that there is a handle to not classify noise as background
      // but the proper way is to have all the handles of different brackgrounds
      auto BgdTrue = evt.getHandle< std::vector<simb::MCTruth> >(GenModuleLabel);   
      if (BgdTrue){                                                                                                           
	      art::FindManyP<simb::MCParticle> BgdAssn(BgdTrue, evt, m_GeantLabel);      

        //Create a temporary map for the specific background source 
        std::map<int, simb::MCParticle> tempBgdMap; 
        FillMyMaps(tempBgdMap, BgdAssn, BgdTrue);

        //Merge the temporary map with the full backgrounds map
        BgdPMap.insert(tempBgdMap.begin(), tempBgdMap.end());                                                         
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
     
      // Check if the TrackId() is not in MarleyPMap or BgdPMap
      if (MarleyPMap.find(trueParticle.TrackId()) == MarleyPMap.end() && BgdPMap.find(trueParticle.TrackId()) == BgdPMap.end()){
        DaughterParts[trueParticle.TrackId()] = trueParticle;
      }
    }
    //Add daughter particles to signal map
    MarleyPMap.insert(DaughterParts.begin(), DaughterParts.end()); 
  }
 
  // ---



  //Map each particle map to its corresponding enum tag 
  std::map<PType, std::map< int, simb::MCParticle >&> PTypeToMap{
    {kMarley,  MarleyPMap },
    {kAr39GenInLAr, Ar39GenInLArPMap },
    {kKr85GenInLAr, Kr85GenInLArPMap },
    {kAr42GenInLAr, Ar42GenInLArPMap },
    {kK42From42ArGenInLAr, K42From42ArGenInLArPMap },
    {kRn222ChainRn222GenInLAr, Rn222ChainRn222GenInLArPMap },
    {kRn222ChainPo218GenInLAr, Rn222ChainPo218GenInLArPMap },
    {kRn222ChainPb214GenInLAr, Rn222ChainPb214GenInLArPMap },
    {kRn222ChainBi214GenInLAr, Rn222ChainBi214GenInLArPMap },
    {kRn222ChainPb210GenInLAr, Rn222ChainPb210GenInLArPMap },
    {kK40GenInCPA, K40GenInCPAPMap },
    {kU238ChainGenInCPA, U238ChainGenInCPAPMap },
    {kK42From42ArGenInCPA, K42From42ArGenInCPAPMap },
    {kRn222ChainPo218GenInCPA, Rn222ChainPo218GenInCPAPMap },
    {kRn222ChainPb214GenInCPA, Rn222ChainPb214GenInCPAPMap },
    {kRn222ChainBi214GenInCPA, Rn222ChainBi214GenInCPAPMap },
    {kRn222ChainPb210GenInCPA, Rn222ChainPb210GenInCPAPMap },
    {kRn222ChainFromBi210GenInCPA, Rn222ChainFromBi210GenInCPAPMap },
    {kCo60GenInAPA, Co60GenInAPAPMap },
    {kU238ChainGenInAPA, U238ChainGenInAPAPMap },
    {kRn222ChainGenInPDS, Rn222ChainGenInPDSPMap },
    {kNeutronGenInRock, NeutronGenInRockPMap }
  };
  
  //run over the particle assn map
  for (auto const& it : PTypeToMap){

    //particle tag e.g. kGen
    const PType p = it.first;
    // gen-g4 mapping e.g. MarleyPMap
    auto const& m = it.second;

    //run over each row in e.g. MarleyPMap
    for (auto const& it2 : m){
      //add a row to the trkIDToPType map consisting of [particle trk ID, kGen] 
      //trkIDToPType is now a 22xn matrix
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
