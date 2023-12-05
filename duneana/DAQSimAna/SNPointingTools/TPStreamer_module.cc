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

  std::map< int, simb::MCParticle> MarleyMap;
  //  I am not satisfied with this, but for the moment let's not change too much the structure
  // an update will be transparent to the user
  std::map< int, simb::MCParticle> Ar39GenInLArPMap;
  std::map< int, simb::MCParticle> Kr85GenInLArPMap;
  std::map< int, simb::MCParticle> Ar42GenInLArPMap;
  std::map< int, simb::MCParticle> K42From42ArGenInLArPMap;
  std::map< int, simb::MCParticle> Rn222ChainRn222GenInLArPMap;
  std::map< int, simb::MCParticle> Rn222ChainPo218GenInLArPMap;
  std::map< int, simb::MCParticle> Rn222ChainPb214GenInLArPMap;
  std::map< int, simb::MCParticle> Rn222ChainBi214GenInLArPMap;
  std::map< int, simb::MCParticle> Rn222ChainPb210GenInLArPMap;
  std::map< int, simb::MCParticle> K40GenInCPAPMap;
  std::map< int, simb::MCParticle> U238ChainGenInCPAPMap;
  std::map< int, simb::MCParticle> K42From42ArGenInCPAPMap;
  std::map< int, simb::MCParticle> Rn222ChainPo218GenInCPAPMap;
  std::map< int, simb::MCParticle> Rn222ChainPb214GenInCPAPMap;
  std::map< int, simb::MCParticle> Rn222ChainBi214GenInCPAPMap;
  std::map< int, simb::MCParticle> Rn222ChainPb210GenInCPAPMap;
  std::map< int, simb::MCParticle> Rn222ChainFromBi210GenInCPAPMap;
  std::map< int, simb::MCParticle> Co60GenInAPAPMap;
  std::map< int, simb::MCParticle> U238ChainGenInAPAPMap;
  std::map< int, simb::MCParticle> Rn222ChainGenInPDSPMap;
  std::map< int, simb::MCParticle> NeutronGenInRockPMap;
  
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
 // print to file the read labels
 std::ofstream labelsFile("labels.txt", std::ios_base::app);
  labelsFile << m_MarleyLabel << std::endl;
  labelsFile << m_Ar39GenInLArLabel << std::endl;
  labelsFile << m_Kr85GenInLArLabel << std::endl; 
  labelsFile << m_Ar42GenInLArLabel << std::endl;
  labelsFile << m_K42From42ArGenInLArLabel << std::endl;
  labelsFile << m_Rn222ChainRn222GenInLArLabel << std::endl;
  labelsFile << m_Rn222ChainPo218GenInLArLabel << std::endl;
  labelsFile << m_Rn222ChainPb214GenInLArLabel << std::endl;
  labelsFile << m_Rn222ChainBi214GenInLArLabel << std::endl;
  labelsFile << m_Rn222ChainPb210GenInLArLabel << std::endl;
  labelsFile << m_K40GenInCPALabel << std::endl;
  labelsFile << m_U238ChainGenInCPALabel << std::endl;
  labelsFile << m_K42From42ArGenInCPALabel << std::endl;
  labelsFile << m_Rn222ChainPo218GenInCPALabel << std::endl;
  labelsFile << m_Rn222ChainPb214GenInCPALabel << std::endl;
  labelsFile << m_Rn222ChainBi214GenInCPALabel << std::endl;
  labelsFile << m_Rn222ChainPb210GenInCPALabel << std::endl;
  labelsFile << m_Rn222ChainFromBi210GenInCPALabel << std::endl;
  labelsFile << m_Co60GenInAPALabel << std::endl;
  labelsFile << m_U238ChainGenInAPALabel << std::endl;
  labelsFile << m_Rn222ChainGenInPDSLabel << std::endl;
  labelsFile << m_NeutronGenInRockLabel << std::endl;
  labelsFile.close();
  // TODO remove this

} // TPStreamer constructor

void TPStreamer::ResetVariables()
{
  MarleyMap.clear();
  Ar39GenInLArPMap.clear();
  Kr85GenInLArPMap.clear();
  Ar42GenInLArPMap.clear();
  K42From42ArGenInLArPMap.clear();
  Rn222ChainRn222GenInLArPMap.clear();
  Rn222ChainPo218GenInLArPMap.clear();
  Rn222ChainPb214GenInLArPMap.clear();
  Rn222ChainBi214GenInLArPMap.clear();
  Rn222ChainPb210GenInLArPMap.clear();
  K40GenInCPAPMap.clear();
  U238ChainGenInCPAPMap.clear();
  K42From42ArGenInCPAPMap.clear();
  Rn222ChainPo218GenInCPAPMap.clear();
  Rn222ChainPb214GenInCPAPMap.clear();
  Rn222ChainBi214GenInCPAPMap.clear();
  Rn222ChainPb210GenInCPAPMap.clear();
  Rn222ChainFromBi210GenInCPAPMap.clear();
  Co60GenInAPAPMap.clear();
  U238ChainGenInAPAPMap.clear();
  Rn222ChainGenInPDSPMap.clear();
  NeutronGenInRockPMap.clear();
  trkIDToPType.clear(); 
  Hit_True_MainTrID.clear();

} // ResetVariables

void TPStreamer::FillMyMaps(std::map<int, simb::MCParticle> &MyMap,
				    art::FindManyP<simb::MCParticle> Assn, 
				    art::Handle< std::vector<simb::MCTruth> > Handle, 
				    std::map<int,int>* indexMap)
{
  std::ofstream AssnFile("Assn.txt", std::ios_base::app);
  for ( size_t L1=0; L1 < Handle->size(); ++L1 ) {
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      // print to file the size of this, with some details
      AssnFile << Assn.at(L1).size() << " " << L1 << " " << L2 << " " << Assn.at(L1).at(L2)->TrackId() << std::endl;
      // print module label in the association
      // AssnFile << Assn.provenance()->moduleLabel() << std::endl;
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

  if(it!=trkIDToPType.end()){   ThisPType=it->second; /*std::cout << "This PType in whichptype is " << ThisPType << std::endl;*/  }
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
    const std::string thisGenLabel = handle.provenance()->moduleLabel();
    // print out the module label in a file
    std::ofstream genLabelsFile("genLabels.txt", std::ios_base::app);
    genLabelsFile << thisGenLabel << std::endl;
    genLabelsFile.close();
    // TODO remove this

    //Get the particle assn for this gen module label
    // art::FindManyP<simb::MCParticle> Association(handle, evt, m_GeantLabel);
    
    std::ofstream MapSizes("MapSizes.txt", std::ios_base::app);

    if (thisGenLabel == m_MarleyLabel) {
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_MarleyLabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel);
        FillMyMaps( MarleyMap, Assn, thisHandle); 
      }
    }
    if (thisGenLabel == m_Ar39GenInLArLabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Ar39GenInLArLabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel);
        FillMyMaps( Ar39GenInLArPMap, Assn, thisHandle); 
        std::cout << "Filling map of ar39" << std::endl;
        MapSizes << "Size of Ar39 Map" << Ar39GenInLArPMap.size() << std::endl;
        std::cout << "Size of Ar39 Map" << Ar39GenInLArPMap.size() << std::endl;
      }
    }
    if (thisGenLabel == m_Kr85GenInLArLabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Kr85GenInLArLabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel);
        FillMyMaps( Kr85GenInLArPMap, Assn, thisHandle); 
        MapSizes << "Size of Kr85 Map" << Kr85GenInLArPMap.size() << std::endl;
        std::cout << "Size of Kr85 Map" << Kr85GenInLArPMap.size() << std::endl;
      }
    }
    if (thisGenLabel == m_Ar42GenInLArLabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Ar42GenInLArLabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel);
        FillMyMaps( Ar42GenInLArPMap, Assn, thisHandle); 
        MapSizes << "Size of Ar42 Map" << Ar42GenInLArPMap.size() << std::endl;
        std::cout << "Size of Ar42 Map" << Ar42GenInLArPMap.size() << std::endl;
      }
    }
    if (thisGenLabel == m_K42From42ArGenInLArLabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_K42From42ArGenInLArLabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel);
        FillMyMaps( K42From42ArGenInLArPMap, Assn, thisHandle); 
        MapSizes << "Size of K42From42Ar Map" << K42From42ArGenInLArPMap.size() << std::endl;
        std::cout << "Size of K42From42Ar Map" << K42From42ArGenInLArPMap.size() << std::endl;
      }
    }
    if (thisGenLabel == m_Rn222ChainRn222GenInLArLabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainRn222GenInLArLabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel);
        FillMyMaps( Rn222ChainRn222GenInLArPMap, Assn, thisHandle); 
        MapSizes << "Size of Rn222ChainRn222 Map" << Rn222ChainRn222GenInLArPMap.size() << std::endl;
        std::cout << "Size of Rn222ChainRn222 Map" << Rn222ChainRn222GenInLArPMap.size() << std::endl;
      }
    }
    if (thisGenLabel == m_Rn222ChainPo218GenInLArLabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainPo218GenInLArLabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel);
        FillMyMaps( Rn222ChainPo218GenInLArPMap, Assn, thisHandle); 
      }
    }
    if (thisGenLabel == m_Rn222ChainPb214GenInLArLabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainPb214GenInLArLabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel);
        FillMyMaps( Rn222ChainPb214GenInLArPMap, Assn, thisHandle); 
      }
    }
    if (thisGenLabel == m_Rn222ChainBi214GenInLArLabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainBi214GenInLArLabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel);
        FillMyMaps( Rn222ChainBi214GenInLArPMap, Assn, thisHandle); 
      }
    }
    if (thisGenLabel == m_Rn222ChainPb210GenInLArLabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainPb210GenInLArLabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel); 
        FillMyMaps( Rn222ChainPb210GenInLArPMap, Assn, thisHandle); 
      }
    }
    if (thisGenLabel == m_K40GenInCPALabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_K40GenInCPALabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel); 
        FillMyMaps( K40GenInCPAPMap, Assn, thisHandle); 
      }
    }
    if (thisGenLabel == m_U238ChainGenInCPALabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_U238ChainGenInCPALabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel); 
        FillMyMaps( U238ChainGenInCPAPMap, Assn, thisHandle); 
      }
    }
    if (thisGenLabel == m_K42From42ArGenInCPALabel){
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_K42From42ArGenInCPALabel);
      if (thisHandle){	
        art::FindManyP<simb::MCParticle> Assn( thisHandle, evt, m_GeantLabel); 
        FillMyMaps( K42From42ArGenInCPAPMap, Assn, thisHandle); 
      }
    }

    if (thisGenLabel == m_Rn222ChainPo218GenInCPALabel) {
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainPo218GenInCPALabel);
      if (thisHandle) {
        art::FindManyP<simb::MCParticle> Assn(thisHandle, evt, m_GeantLabel);
        FillMyMaps(Rn222ChainPo218GenInCPAPMap, Assn, thisHandle);
      }
    }
    if (thisGenLabel == m_Rn222ChainPb214GenInCPALabel) {
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainPb214GenInCPALabel);
      if (thisHandle) {
        art::FindManyP<simb::MCParticle> Assn(thisHandle, evt, m_GeantLabel);
        FillMyMaps(Rn222ChainPb214GenInCPAPMap, Assn, thisHandle);
      }
    }
    if (thisGenLabel == m_Rn222ChainBi214GenInCPALabel) {
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainBi214GenInCPALabel);
      if (thisHandle) {
        art::FindManyP<simb::MCParticle> Assn(thisHandle, evt, m_GeantLabel);
        FillMyMaps(Rn222ChainBi214GenInCPAPMap, Assn, thisHandle);
      }
    }
    if (thisGenLabel == m_Rn222ChainPb210GenInCPALabel) {
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainPb210GenInCPALabel);
      if (thisHandle) {
        art::FindManyP<simb::MCParticle> Assn(thisHandle, evt, m_GeantLabel);
        FillMyMaps(Rn222ChainPb210GenInCPAPMap, Assn, thisHandle);
      }
    }
    if (thisGenLabel == m_Rn222ChainFromBi210GenInCPALabel) {
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainFromBi210GenInCPALabel);
      if (thisHandle) {
        art::FindManyP<simb::MCParticle> Assn(thisHandle, evt, m_GeantLabel);
        FillMyMaps(Rn222ChainFromBi210GenInCPAPMap, Assn, thisHandle);
      }
    }
    if (thisGenLabel == m_Co60GenInAPALabel) {
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Co60GenInAPALabel);
      if (thisHandle) {
        art::FindManyP<simb::MCParticle> Assn(thisHandle, evt, m_GeantLabel);
        FillMyMaps(Co60GenInAPAPMap, Assn, thisHandle);
      }
    }
    if (thisGenLabel == m_U238ChainGenInAPALabel) {
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_U238ChainGenInAPALabel);
      if (thisHandle) {
        art::FindManyP<simb::MCParticle> Assn(thisHandle, evt, m_GeantLabel);
        FillMyMaps(U238ChainGenInAPAPMap, Assn, thisHandle);
      }
    }
    if (thisGenLabel == m_Rn222ChainGenInPDSLabel) {
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_Rn222ChainGenInPDSLabel);
      if (thisHandle) {
        art::FindManyP<simb::MCParticle> Assn(thisHandle, evt, m_GeantLabel);
        FillMyMaps(Rn222ChainGenInPDSPMap, Assn, thisHandle);
        std::cout << "Filling map of Rn222ChainGenInPDSPMap" << std::endl;
      }
    }
    if (thisGenLabel == m_NeutronGenInRockLabel) {
      auto thisHandle = evt.getHandle< std::vector<simb::MCTruth> >(m_NeutronGenInRockLabel);
      if (thisHandle) {
        art::FindManyP<simb::MCParticle> Assn(thisHandle, evt, m_GeantLabel);
        FillMyMaps(NeutronGenInRockPMap, Assn, thisHandle);
      }
    }
  }


  // Get any G4 tracks that weren't matched to signal or radio gen
  // and add them to signal map if they were most likely daughter particles (e.g photons from EM shower)
  
  art::ValidHandle<std::vector <simb::MCParticle> > mcParticles = evt.getValidHandle<std::vector <simb::MCParticle> >(m_GeantLabel);
  if (mcParticles.isValid()){
    std::map< int, simb::MCParticle> DaughterParts;

    for (unsigned int i = 0; i < mcParticles->size(); ++i) {
      const simb::MCParticle trueParticle = mcParticles->at(i);
      // if (trueParticle.Mother() != 0){ 	DaughterParts[trueParticle.TrackId()] = trueParticle;      } // why this was commented?

      // Check if the TrackId() is not in any map, then add it to the daughter map
      // TODO maybe make this a bit more efficient?
      std::vector<std::map<int, simb::MCParticle>> AllMaps = {
        MarleyMap, Ar39GenInLArPMap, Kr85GenInLArPMap, Ar42GenInLArPMap, K42From42ArGenInLArPMap, 
        Rn222ChainRn222GenInLArPMap, Rn222ChainPo218GenInLArPMap, Rn222ChainPb214GenInLArPMap, Rn222ChainBi214GenInLArPMap, Rn222ChainPb210GenInLArPMap, K40GenInCPAPMap, U238ChainGenInCPAPMap, 
        K42From42ArGenInCPAPMap, Rn222ChainPo218GenInCPAPMap, Rn222ChainPb214GenInCPAPMap, Rn222ChainBi214GenInCPAPMap, Rn222ChainPb210GenInCPAPMap, Rn222ChainFromBi210GenInCPAPMap, 
        Co60GenInAPAPMap, U238ChainGenInAPAPMap, Rn222ChainGenInPDSPMap, NeutronGenInRockPMap};
      
      bool found = false;
      for (auto const& it : AllMaps) {
        if (it.find(trueParticle.TrackId()) != it.end()) {
          found = true;
          break;
        }
      }
      if (!found) DaughterParts[trueParticle.TrackId()] = trueParticle;
    }
    //Add daughter particles to signal map. 
    // TODO check if daughter particles are supposed to be only for marley or also bgd
    MarleyMap.insert(DaughterParts.begin(), DaughterParts.end()); 
  }
 
  // ---



  //Map each particle map to its corresponding enum tag, maybe can avoid
  std::map<PType, std::map< int, simb::MCParticle >&> PTypeToMap{
    {kMarley,  MarleyMap },
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
    std::cout << "PType is " << p << std::endl;
    // gen-g4 mapping e.g. MarleyMap
    auto const& m = it.second;
    std::cout << "Size of map is " << m.size() << std::endl;

    //run over each row in e.g. MarleyMap
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
    
    // add multiple catches? TODO
    std::vector<const sim::IDE*> ThisSimIDE = bt_serv->HitToSimIDEs_Ps(clockData, ThisHit.Channel(), WindowStart, WindowEnd);
    // there should be more conditions here...

    double trueX = 0, trueY = 0, trueZ = 0, trueEnergy = 0, nElectrons = 0, trackId = 0;
    
    for(unsigned int i = 0; i < ThisSimIDE.size(); i++)
    {
      if(ThisSimIDE.at(i)->trackID==Hit_True_MainTrID.at(hit))
      {
        trueX = ThisSimIDE.at(i)->x;
        trueY = ThisSimIDE.at(i)->y;
        trueZ = ThisSimIDE.at(i)->z;
        trueEnergy = ThisSimIDE.at(i)->energy;
        nElectrons = ThisSimIDE.at(i)->numElectrons;
        break;
      }
    }

    // I want numbers to be printed out always as integers, not floats
    m_outputFile << std::fixed << std::setprecision(0) << std::setw(0) << std::setfill('0'); // TODO alignment
    // printing TPs
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
     << trueX << ' ' // trueX, to be implemented in BT
     << trueY << ' ' // trueY, to be implemented in BT
     << trueZ << ' ' // trueZ, to be implemented in BT
     << trueEnergy << ' ' // energy, to be implemented in BT
     << nElectrons << ' ' // nElectrons, to be implemented in BT OR SOMETHING ELSE
     << 0 << ' ' // trackId, to be implemented in BT OR SOMETHING ELSE
     <<  std::endl; 

  } // Loop over reco_hits.
} // Analyze TPStreamer.


//......................................................
DEFINE_ART_MODULE(TPStreamer)
