////////////////////////////////////////////////////////////////////////////////////
// Class:       SolarNuAna                                                        //
// Module Type: analyzer                                                          //
// File:        SolarNuAna_module.cc                                              //
//                                                                                //  
// Written by Sergio Manthey Corchado with guidence of Daniel Pershey             // 
// developed from Michael Baird's DAQSimAna_module                                //    
////////////////////////////////////////////////////////////////////////////////////

// C++ includes
#ifndef SolarNuAna_h
#define SolarNuAna_h 

// C++ includes
// ROOT includes
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom.h"
// #include "DisplacementVector3D.h"

// Framework includes (not all might be necessary)
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

class SolarNuAna : public art::EDAnalyzer {

public:
  // --- Standard constructor and destructor for an ART module.
  explicit SolarNuAna(fhicl::ParameterSet const & p);
  SolarNuAna(SolarNuAna const &)               = delete;
  SolarNuAna(SolarNuAna &&)                    = delete;
  SolarNuAna & operator = (SolarNuAna const &) = delete;
  SolarNuAna & operator = (SolarNuAna &&)      = delete;
  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;

private:
  // --- Some of our own functions.
  void              ResetVariables();
  long unsigned int WhichParType  ( int TrID );
  bool              InMyMap       ( int TrID, std::map< int, simb::MCParticle> ParMap );
  void              FillMyMaps    ( std::map< int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand );
  void              CalcAdjHits   ( std::vector< recob::Hit > MyVec,std::vector< std::vector<recob::Hit> >& Clusters,TH1I* MyHist, TH1F* MyADCIntHist, bool HeavDebug );
  void              PrintInColor  ( std::string MyString, int MyColor );
  int               GetColor      ( std::string MyString );
  std::string       str           ( int MyInt );
  std::string       str           ( float MyFloat );
  std::string       str           ( double MyDouble );
  
  // --- Our fcl parameter labels for the modules that made the data products
  std::string fRawDigitLabel,fHitLabel,fOpHitLabel,fOpDetWaveformLabel,fOpFlashLabel,fGEANTLabel; 

  // --- Input settings imported from the fcl
  std::string sPositionReco;
  std::vector<std::string> fLabels;
  int fMaxDetSizeY, fClusterMatchMinNHit, fGoalInd0MatchTime, fGoalInd1MatchTime;
  float fClusterMatchTime,fAdjClusterTime,fAdjClusterRad,fAdjOpFlashRad,fAdjOpFlashTime,fAdjOpFlashMaxPECut;
  bool fTestNewClReco, fDebug;
  
  // --- Our TTrees, and its associated variables.
  TTree* fSolarNuAnaTree;
  TTree* fMCTruthTree;
  std::string MGenLabel;
  int Run,SubRun,Event,Flag,MNHit,MGen,MInd0TPC,MInd1TPC,MInd0NHits,MInd1NHits,MMainID;
  float TNuQSqr,TNuE,TNuP,TNuX,TNuY,TNuZ,avX,avY,avZ,MTime,MChrg,MInd0MaxHit,MInd1MaxHit,MInd0dT,MInd1dT,MInd0RecoY,MInd1RecoY,MRecZ,MPur;
  std::vector<int> MAdjClGen,TPart,MarleyPDGList,MarleyIDList,MarleyParentIDList;
  std::vector<float> MAdjClTime,MAdjClCharge,MAdjClNHit,MAdjClRecoY,MAdjClRecoZ,MAdjClR,MAdjClPur,MMarleyFrac,MGenFrac;
  std::vector<float> MAdjFlashTime,MAdjFlashPE,MAdjFlashNHit,MAdjFlashMaxPE,MAdjFlashRecoY,MAdjFlashRecoZ,MAdjFlashR,MAdjFlashPur;
  std::vector<float> MarleyEList, MarleyPList, MarleyXList, MarleyYList, MarleyZList;
  std::vector<std::map<int,simb::MCParticle>> Parts = {};
  
  // --- OpFlash Variables
  std::vector<float> OpFlashMarlPur,OpFlashPE,OpFlashMaxPE,OpFlashY,OpFlashZ,OpFlashT,OpFlashDeltaT,OpFlashNHit;
  
  // --- Histograms to fill about collection plane hits
  TH2F* hXTruth;
  TH2F* hYTruth;
  TH2F* hZTruth;
  TH1I* hAdjHits; 
  TH1F* hAdjHitsADCInt; 
  TH2F* hDriftTime;
  
  // --- Declare our services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
};
#endif

//......................................................
SolarNuAna::SolarNuAna(fhicl::ParameterSet const & p):EDAnalyzer(p){this->reconfigure(p);}

//......................................................
void SolarNuAna::reconfigure(fhicl::ParameterSet const & p){
  fLabels = p.get<std::vector<std::string>> ("ParticleLabelVector");
  
  fRawDigitLabel      = p.get<std::string>  ("RawDigitLabel");
  fHitLabel           = p.get<std::string>  ("HitLabel");
  fOpFlashLabel       = p.get<std::string>  ("OpFlashLabel");
  fOpHitLabel         = p.get<std::string>  ("OpHitLabel");
  fGEANTLabel         = p.get<std::string>  ("GEANT4Label");

  sPositionReco        = p.get<std::string> ("PositionRecoMode","DEFAULT");
  fMaxDetSizeY         = p.get<int>         ("DetectorSizeY");
  fClusterMatchMinNHit = p.get<int>         ("ClusterMatchMinNHit");
  fClusterMatchTime    = p.get<float>       ("ClusterMatchTime");
  fGoalInd0MatchTime   = p.get<float>       ("GoalInd0MatchTime");
  fGoalInd1MatchTime   = p.get<float>       ("GoalInd1MatchTime");
  fAdjClusterTime      = p.get<float>       ("AdjClusterTime");
  fAdjClusterRad       = p.get<float>       ("AdjClusterRad");
  fAdjOpFlashTime      = p.get<float>       ("AdjOpFlashTime");
  fAdjOpFlashRad       = p.get<float>       ("AdjOpFlashRad");
  fAdjOpFlashMaxPECut  = p.get<float>       ("AdjOpFlashMaxPECut");
  fTestNewClReco       = p.get<bool>        ("TestNewClReco",false);
  fDebug               = p.get<bool>        ("Debug",false);
} // Reconfigure

//......................................................
void SolarNuAna::beginJob(){   
  // --- Make our handle to the TFileService
  art::ServiceHandle<art::TFileService> tfs;
  fMCTruthTree = tfs->make<TTree>("MCTruthTree","MC Truth Tree");
  fSolarNuAnaTree = tfs->make<TTree>("SolarNuAnaTree","Solar Ana Tree");

  // MC Truth info.
  fMCTruthTree -> Branch("Run",               &Run,              "Run/I");
  fMCTruthTree -> Branch("SubRun",            &SubRun,           "SubRun/I");
  fMCTruthTree -> Branch("Event",             &Event,            "Event/I");
  fMCTruthTree -> Branch("Flag",              &Flag,             "Flag/I");          // Flag used to match truth with reco tree entries
  fMCTruthTree -> Branch("TruthPart",         &TPart);                               // Number particles per generator
  fMCTruthTree -> Branch("TNuQSqr",           &TNuQSqr,          "TruthNuQSqr/F");   // True neutrino transfer momentum [GeV]
  fMCTruthTree -> Branch("TNuE",              &TNuE,             "TruthNuE/F");      // True neutrino energy [GeV]
  fMCTruthTree -> Branch("TNuP",              &TNuP,             "TruthNuP/F");      // True neutrino momentum [GeV]
  fMCTruthTree -> Branch("TNuX",              &TNuX,             "TruthNuX/F");      // True neutrino X [cm] 
  fMCTruthTree -> Branch("TNuY",              &TNuY,             "TruthNuY/F");      // True neutrino Y [cm]
  fMCTruthTree -> Branch("TNuZ",              &TNuZ,             "TruthNuZ/F");      // True neutrino Z [cm]
  fMCTruthTree -> Branch("TMarleyPDG",        &MarleyPDGList);                       // PDG of marley marticles
  fMCTruthTree -> Branch("TMarleyE",          &MarleyEList);                         // Energy of marley particles [GeV]
  fMCTruthTree -> Branch("TMarleyP",          &MarleyPList);                         // Momentum of marley particles [GeV]
  fMCTruthTree -> Branch("TMarleyX",          &MarleyXList);                         // X of marley particles [cm] 
  fMCTruthTree -> Branch("TMarleyY",          &MarleyYList);                         // Y of marley particles [cm]
  fMCTruthTree -> Branch("TMarleyZ",          &MarleyZList);                         // Z of marley particles [cm]
  fMCTruthTree -> Branch("TMarleyID",         &MarleyIDList);                        // TrackID of marley particles
  fMCTruthTree -> Branch("TMarleyParentID",   &MarleyParentIDList);                  // ParentID of marley particles

  // Repeated Truth info.
  fSolarNuAnaTree -> Branch("Run",            &Run,              "Run/I");
  fSolarNuAnaTree -> Branch("SubRun",         &SubRun,           "SubRun/I");
  fSolarNuAnaTree -> Branch("Event",          &Event,            "Event/I");
  fSolarNuAnaTree -> Branch("Flag",           &Flag,             "Flag/I");          // Flag used to match truth with reco tree entries
  fSolarNuAnaTree -> Branch("TruthPart",      &TPart);                               // Number particles per generator
  fSolarNuAnaTree -> Branch("TNuQSqr",        &TNuQSqr,          "TruthNuQSqr/F");   // True neutrino transfer momentum [GeV]
  fSolarNuAnaTree -> Branch("TNuE",           &TNuE,             "TruthNuE/F");      // True neutrino energy
  fSolarNuAnaTree -> Branch("TNuP",           &TNuP,             "TruthNuP/F");      // True neutrino momentum
  fSolarNuAnaTree -> Branch("TNuX",           &TNuX,             "TruthNuX/F");      // True neutrino X
  fSolarNuAnaTree -> Branch("TNuY",           &TNuY,             "TruthNuY/F");      // True neutrino Y
  fSolarNuAnaTree -> Branch("TNuZ",           &TNuZ,             "TruthNuZ/F");      // True neutrino Z
  fSolarNuAnaTree -> Branch("TMarleyPDG",     &MarleyPDGList);                       // PDG of marley particles
  fSolarNuAnaTree -> Branch("TMarleyE",       &MarleyEList);                         // Energy of marley particles
  fSolarNuAnaTree -> Branch("TMarleyP",       &MarleyPList);                         // Momentum of marley particles
  fSolarNuAnaTree -> Branch("TMarleyX",       &MarleyXList);                         // X of marley particles
  fSolarNuAnaTree -> Branch("TMarleyY",       &MarleyYList);                         // Y of marley particles
  fSolarNuAnaTree -> Branch("TMarleyZ",       &MarleyZList);                         // Z of marley particles
  fSolarNuAnaTree -> Branch("TMarleyID",      &MarleyIDList);                        // TrackID of marley particles")
  fSolarNuAnaTree -> Branch("TMarleyParentID",&MarleyParentIDList);                  // ParentID of marley particles
  
  // Main Cluster info.
  fSolarNuAnaTree -> Branch("Generator",      &MGen,             "Generator/I");     // Main cluster generator idx  
  fSolarNuAnaTree -> Branch("Time",           &MTime,            "ColTime/F");       // Main cluster time [ticks]
  fSolarNuAnaTree -> Branch("Charge",         &MChrg,            "ColCharge/F");     // Main cluster charge [ADC*ticks]
  fSolarNuAnaTree -> Branch("NHits",          &MNHit,            "ColNHits/I");      // Main cluster #hits
  fSolarNuAnaTree -> Branch("Ind0TPC",        &MInd0TPC,         "Ind0TPC/I");       // Main cluster ind0 TPC
  fSolarNuAnaTree -> Branch("Ind1TPC",        &MInd1TPC,         "Ind1TPC/I");       // Main cluster ind1 TPC  
  fSolarNuAnaTree -> Branch("Ind0NHits",      &MInd0NHits,       "Ind0NHits/I");     // Main cluster ind0 Hits
  fSolarNuAnaTree -> Branch("Ind1NHits",      &MInd1NHits,       "Ind1NHits/I");     // Main cluster ind1 Hits
  fSolarNuAnaTree -> Branch("Ind0MaxHit",     &MInd0MaxHit,      "Ind0MaxHit/F");    // Main cluster ind0 MaxHit
  fSolarNuAnaTree -> Branch("Ind1MaxHit",     &MInd1MaxHit,      "Ind1MaxHit/F");    // Main cluster ind1 MaxHit
  fSolarNuAnaTree -> Branch("Ind0dT",         &MInd0dT,          "Ind0dT/F");        // Main cluster ind0 DT [Ticks]
  fSolarNuAnaTree -> Branch("Ind1dT",         &MInd1dT,          "Ind1dT/F");        // Main cluster ind1 DT [Ticks]
  fSolarNuAnaTree -> Branch("Ind0RecoY",      &MInd0RecoY,       "Ind0RecoY/F");     // Main cluster ind0 reco Y [cm]
  fSolarNuAnaTree -> Branch("Ind1RecoY",      &MInd1RecoY,       "Ind1RecoY/F");     // Main cluster ind1 reco Y [cm]
  fSolarNuAnaTree -> Branch("MainID",         &MMainID,          "MainID/I");        // Main cluster main track ID
  fSolarNuAnaTree -> Branch("RecoZ",          &MRecZ,            "RecoZ/F");         // Main cluster reco Z [cm]
  fSolarNuAnaTree -> Branch("Purity",         &MPur,             "Purity/F");        // Main cluster reco purity
  fSolarNuAnaTree -> Branch("GenFrac",        &MGenFrac);                            // Main cluster reco purity complete
  fSolarNuAnaTree -> Branch("MarleyFrac",     &MMarleyFrac);                         // Main cluster particle contribution (electron, gamma, neutron)
  // fSolarNuAnaTree -> Branch("Label",          &MGenLabel);                        // Main cluster generator label  
      
  // Adj. Cluster info.
  fSolarNuAnaTree -> Branch("AdjClGen",       &MAdjClGen);                           // Adj. clusters' generator idx
  fSolarNuAnaTree -> Branch("AdjClTime",      &MAdjClTime);                          // Adj. clusters' time [ticks]
  fSolarNuAnaTree -> Branch("AdjClCharge",    &MAdjClCharge);                        // Adj. clusters' charge [ADC*ticks]
  fSolarNuAnaTree -> Branch("AdjClNHit",      &MAdjClNHit);                          // Adj. clusters' #hits 
  fSolarNuAnaTree -> Branch("AdjClRecoY",     &MAdjClRecoY);                         // Adj. clusters' reco Y [cm]
  fSolarNuAnaTree -> Branch("AdjClRecoZ",     &MAdjClRecoZ);                         // Adj. clusters' reco Z [cm]
  fSolarNuAnaTree -> Branch("AdjClPur",       &MAdjClPur);                           // Adj. clusters' purity
  fSolarNuAnaTree -> Branch("AdjClR",         &MAdjClR);                             // Adj. clusters' distance [cm]

  // Adj. Flash info.
  fSolarNuAnaTree -> Branch("AdjOpFlashTime", &MAdjFlashTime);                       // Adj. flash' time [ticks]
  fSolarNuAnaTree -> Branch("AdjOpFlashPE",   &MAdjFlashPE);                         // Adj. flash' tot #PE [ADC*ticks]
  fSolarNuAnaTree -> Branch("AdjOpFlashNHit", &MAdjFlashNHit);                       // Adj. flash' #hits
  fSolarNuAnaTree -> Branch("AdjOpFlashMaxPE",&MAdjFlashMaxPE);                      // Adj. flash' max #PE [ADC*ticks]
  fSolarNuAnaTree -> Branch("AdjOpFlashRecoY",&MAdjFlashRecoY);                      // Adj. flash' reco Y [cm] 
  fSolarNuAnaTree -> Branch("AdjOpFlashRecoZ",&MAdjFlashRecoZ);                      // Adj. flash' reco Z [cm]
  fSolarNuAnaTree -> Branch("AdjOpFlashPur",  &MAdjFlashPur);                        // Adj. flash' purity
  fSolarNuAnaTree -> Branch("AdjOpFlashR",    &MAdjFlashR);                          // Adj. flash' reco distance [cm]

  // --- Our Histograms...
  hDriftTime      = tfs->make<TH2F>("hDriftTime", "hDriftTime"  , 100, -400., 400., 100, 0., 10000.);
  hXTruth         = tfs->make<TH2F>("hXTruth", "Missmatch in Y distance; Distance [cm]; True X position [cm]", 100, -600, 600, 100, -600, 600);
  hYTruth         = tfs->make<TH2F>("hYTruth", "Missmatch in Y distance; Distance [cm]; True Y position [cm]", 100, -600, 600, 100, -600, 600);
  hZTruth         = tfs->make<TH2F>("hZTruth", "Missmatch in Y distance; Distance [cm]; True Z position [cm]", 100, -600, 600, 100, 0, 1600);
  hAdjHits        = tfs->make<TH1I>("hAdjHits", "Number of adjacent collection plane hits; Number of adjacent collection plane hits; Number of events"  , 21, -0.5, 20.5 );
  hAdjHitsADCInt  = tfs->make<TH1F>("hAdjHitsADCInt", "Total summed ADC Integrals for clusters; Total summed ADC Integrals for clusters; Number of events"  , 1000, 0, 10000 );
} // BeginJob

//......................................................
void SolarNuAna::analyze(art::Event const & evt)
{ 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  //------------------------------------------------------------- Prepare everything for new event ----------------------------------------------------------------// 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  std::vector<std::set<int>>             trackids = {};
  std::map<int,simb::MCParticle>         ThisGeneratorParts;
  std::vector<recob::Hit>                ColHits0,ColHits1,ColHits2,ColHits3; 
  std::vector<std::vector<recob::Hit>>   ColHits = {ColHits0,ColHits1,ColHits2,ColHits3};
  std::vector<std::vector<recob::Hit>>   Clusters0, Clusters1, Clusters2, Clusters3;

  // --- We want to reset all of our previous run and TTree variables ---
  ResetVariables();
  ThisGeneratorParts.clear();
  
  Run = evt.run();SubRun = evt.subRun();Event = evt.event();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  
  Flag = rand() % 10000000000;
  std::cout << "\nTPC Frequency in [MHz]: " << clockData.TPCClock().Frequency() << std::endl;
  std::cout << "TPC Tick in [us]: " << clockData.TPCClock().TickPeriod() << std::endl;
  std::cout << "Event Flag: " << Flag << std::endl;
  std::cout << "\nSuccesfull reset of variables for evt " << Event << " run " << Run << std::endl << "########################################" << std::endl;
  
  std::cout << std::endl; //---------------------------------------------------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------- Create maps for ID tracking -----------------------------------------------------------------// 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  // --- Fill MC Truth IDs to tracking vectors. Get a list of all of my particles in one chunk. ---  
  const sim::ParticleList& PartList = pi_serv->ParticleList();
  std::cout << "There are a total of " << PartList.size() << " Particles in the event\n";
  
  // Loop over all signal+bkg handles and collect track IDs
  for ( size_t i = 0; i < fLabels.size(); i++){
    Parts.push_back(ThisGeneratorParts);                                    // For each label insert empty list
    
    art::Handle<std::vector<simb::MCTruth>> ThisHandle;
    evt.getByLabel(fLabels[i], ThisHandle);
    
    if(ThisHandle){
      auto ThisValidHanlde = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[i]); // Get generator handles
      art::FindManyP<simb::MCParticle> Assn(ThisValidHanlde,evt,fGEANTLabel);            // Assign labels to MCPArticles
      FillMyMaps( Parts[i], Assn, ThisValidHanlde);                                      // Fill empty list with previously assigned particles                                       
      if (Parts[i].size() < 1000){std::cout << "# of particles " << Parts[i].size() << "\tfrom " << fLabels[i] << std::endl;} // Print signal+bkg info to terminal
      else {std::cout << "# of particles " << Parts[i].size() << "\tfrom " << fLabels[i] << std::endl;}
      TPart.push_back(Parts[i].size()); // Insert #signal+bkg particles generated
      for ( std::map<int,simb::MCParticle>::iterator iter = Parts[i].begin(); iter != Parts[i].end(); iter++ ){
        std::set<int> ThisGeneratorIDs = {};
        trackids.push_back(ThisGeneratorIDs);
        trackids[i].insert( iter->first );// Contains a list of TrIDs
      }
    } 
    else{
      std::cout << "# of particles " << Parts[i].size() << "\tfrom " << fLabels[i] << " *not generated!" << std::endl;
      TPart.push_back(0);
      std::set<int> ThisGeneratorIDs = {};
      trackids.push_back(ThisGeneratorIDs);
    }
  }

  std::cout << std::endl;//----------------------------------------------------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------- Some MC Truth information -------------------------------------------------------------------// 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  std::set< int > signal_trackids;                                // Signal TrackIDs to be used in OpFlash matching
  std::vector<std::vector<int>> ClPartTrackIDs = {{},{},{},{}};   // Track IDs corresponding to each kind of MCTruth particle  {11,22,2112,else}
  art::Handle<std::vector<simb::MCTruth>> ThisHandle;
  evt.getByLabel(fLabels[0], ThisHandle);
  if (ThisHandle){
    auto MarlTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fLabels[0]);  // Get handle for MARLEY MCTruths
    // --- Loop over all neutrinos in the event ---
    for (auto const &MARLEYtruth : *MarlTrue){ 
      const simb::MCNeutrino &nue = MARLEYtruth.GetNeutrino();
      TNuQSqr = nue.QSqr();
      TNuE =    nue.Nu().E();
      TNuP =    nue.Nu().Pt();
      TNuX =    nue.Nu().Vx();
      TNuY =    nue.Nu().Vy();
      TNuZ =    nue.Nu().Vz();
      int N =   MARLEYtruth.NParticles(); 
      std::cout << "Number of Neutrino Daughters: " << N-2 << std::endl; 
      std::cout << "Neutrino energy: " << TNuE << " GeV; with transfer momentum: " << std::sqrt(TNuQSqr) << " GeV" << std::endl; 
      std::cout << "Position (" << TNuX << ", " << TNuY << ", " << TNuZ << ") cm" << std::endl;
      
      // Save information of each daughter particle of the marley process
      for ( int i = 0; i < N; i++) {
        const simb::MCParticle &MarleyParticle = MARLEYtruth.GetParticle(i);
        int MarleyParticlePDG = MarleyParticle.PdgCode(); MarleyPDGList.push_back(MarleyParticlePDG);
        int MarleyParticleID  = MarleyParticle.TrackId(); MarleyIDList.push_back(MarleyParticleID);
        int MarleyParentID    = MarleyParticle.Mother();  MarleyParentIDList.push_back(MarleyParentID);
        float MarleyParticleE = MarleyParticle.E();       MarleyEList.push_back(MarleyParticleE);
        float MarleyParticleP = MarleyParticle.P();       MarleyPList.push_back(MarleyParticleP);
        float MarleyParticleX = MarleyParticle.EndX();    MarleyXList.push_back(MarleyParticleX);
        float MarleyParticleY = MarleyParticle.EndY();    MarleyYList.push_back(MarleyParticleY);
        float MarleyParticleZ = MarleyParticle.EndZ();    MarleyZList.push_back(MarleyParticleZ);
      }
    }
    
    art::FindManyP<simb::MCParticle> MarlAssn(MarlTrue,evt,fGEANTLabel);
    std::cout << "\nGen.\t " << "PdgCode\t Energy\t\t TrackID" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    
    for ( size_t i = 0; i < MarlAssn.size(); i++) {
      auto parts = MarlAssn.at(i);
      for (auto part = parts.begin(); part != parts.end(); part++) {
        signal_trackids.emplace((*part)->TrackId());

        if ((*part)->PdgCode()<1000000){std::cout << fLabels[0] << "\t " << (*part)->PdgCode() << "\t\t " << (*part)->E() << "\t " << (*part)->TrackId() << std::endl;}
        else{std::cout << fLabels[0] << "\t " << (*part)->PdgCode() << "\t " << (*part)->E() << "\t " << (*part)->TrackId() << std::endl;}

        if ((*part)->PdgCode()==11){ // Electrons
          const TLorentzVector &v4_f = (*part)->EndPosition();
          auto x_f = v4_f.X();auto y_f = v4_f.Y();auto z_f = v4_f.Z();
          avX = x_f; avY = y_f; avZ = z_f;
          ClPartTrackIDs[0].push_back((*part)->TrackId());
          if (fDebug) std::cout << "\nMC Electron truth position x = " << avX << ", y = " << avY << ", z = " << avZ << std::endl;
          if (fDebug) std::cout << "Initial KE " << (*part)->E()-(*part)->Mass() << std::endl;
        }
        if ((*part)->PdgCode()==22){ClPartTrackIDs[1].push_back((*part)->TrackId());} // Gammas
        if ((*part)->PdgCode()==2112){ClPartTrackIDs[2].push_back((*part)->TrackId());} // Neutrons
        if ((*part)->PdgCode()!=11 && (*part)->PdgCode()!=22 && (*part)->PdgCode()!=2112){ClPartTrackIDs[3].push_back((*part)->TrackId());} // Others
        // else {ClPartTrackIDs[3].push_back((*part)->TrackId()); std::cout << "TrackID " << (*part)->TrackId() << " going to rest" << std::endl;}
      }
    }
    fMCTruthTree->Fill();
  }
  else{std::cout << "No MARLEY MCTruths found." << std::endl;}
  
  std::cout << std::endl;//----------------------------------------------------------------------------------------------------------------------------------------//
  //------------------------------------------------------------------- Optical Flash Analysis --------------------------------------------------------------------// 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  // Find OpHits and OpFlashes associated with the event
  art::Handle< std::vector< recob::OpHit >> OpHitHandle;
  art::Handle< std::vector< recob::OpFlash >> FlashHandle;
  std::vector<art::Ptr<recob::OpHit >> ophitlist;
  std::vector<art::Ptr<recob::OpFlash >> flashlist;
  if (evt.getByLabel(fOpHitLabel, OpHitHandle)){art::fill_ptr_vector(ophitlist, OpHitHandle);}
  if (evt.getByLabel(fOpFlashLabel, FlashHandle)){art::fill_ptr_vector(flashlist, FlashHandle);}
  
  // Grab assns with OpHits to get match to neutrino purity
  art::FindManyP< recob::OpHit > OpAssns(flashlist, evt, fOpFlashLabel);
  std::cout << "\nTotal number of flashes constructed: " << flashlist.size() << std::endl;

  // Loop over flashlist and assign OpHits to each flash
  for ( int i = 0; i < int(flashlist.size()); i++ ){
    recob::OpFlash TheFlash = *flashlist[i];
    std::vector< art::Ptr< recob::OpHit > > matchedHits = OpAssns.at(i);
    if (fDebug) std::cout << "Assigning OpHit to Flash" << std::endl;
    // Calculate the total PE of the flash and the time of the ophit with the highest PE 
    double totPE = 0; double MaxHitPE = 0;
    float OpHitT, OpHitPE;
    for (int j = 0; j < int(matchedHits.size()); j++){
      recob::OpHit ohit = *matchedHits[j];
      totPE += ohit.PE();
      OpHitPE = ohit.PE();
      if (OpHitPE > MaxHitPE){
        MaxHitPE = OpHitPE;
        OpHitT = ohit.PeakTimeAbs();
      }
    }

    if (fDebug) std::cout << "Evaluating Flash purity" << std::endl;
    double OpFlashPur = pbt->OpHitCollectionPurity(signal_trackids, matchedHits);
    if (fDebug) std::cout << "PE of this OpFlash " << totPE << " OpFlash time " << OpHitT << std::endl;

    // Calculate the flash purity, only for the Marley events
    if (MaxHitPE/totPE < fAdjOpFlashMaxPECut){
      OpFlashMarlPur.push_back(OpFlashPur);
      OpFlashPE.push_back(TheFlash.TotalPE());
      OpFlashMaxPE.push_back(MaxHitPE);
      OpFlashY.push_back(TheFlash.YCenter());
      OpFlashZ.push_back(TheFlash.ZCenter());
      OpFlashT.push_back(TheFlash.Time());
      OpFlashDeltaT.push_back(TheFlash.TimeWidth());
      OpFlashNHit.push_back(matchedHits.size());
    }
    if (fDebug && abs(OpHitT) < 30) std::cout << "OpFlash PE (max/tot) " << MaxHitPE << "/" << TheFlash.TotalPE() << " with purity " << OpFlashPur << " time " << TheFlash.Time() << std::endl;
  }

  std::cout << std::endl;//----------------------------------------------------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------- Hit collection and assignment ----------------------------------------------------------------// 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  // --- Lift out the reco hits:
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(fHitLabel);
  int NTotHits = reco_hits->size();

  for(int hit = 0; hit < NTotHits; ++hit){
    // --- Loop over the reconstructed hits to separate them among tpc planes according to view
    
    recob::Hit const& ThisHit = reco_hits->at(hit);
    if (fDebug) std::cout << "Hit " << hit << " has view " << ThisHit.View() << " and signal type " << ThisHit.SignalType() << std::endl;

    if      (ThisHit.SignalType() == 0 && ThisHit.View() == 0){ColHits0.push_back( ThisHit );} // SignalType = 0
    else if (ThisHit.SignalType() == 0 && ThisHit.View() == 1){ColHits1.push_back( ThisHit );} // SignalType = 0
    else if (ThisHit.SignalType() == 1)                       {ColHits2.push_back( ThisHit );} // SignalType = 1
    else    {ColHits3.push_back( ThisHit ); std::cout << "Hit was found with view out of scope" << std::endl;}
  } 

  std::cout << "# Hits in each view = " << ColHits0.size() << ", " << ColHits1.size() << ", " << ColHits2.size() << ", " << ColHits3.size() << std::endl;
  
  std::cout << std::endl;//----------------------------------------------------------------------------------------------------------------------------------------//
  //-------------------------------------------------------------- Cluster creation and analysis ------------------------------------------------------------------// 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  // --- Now calculate the clusters ...
  CalcAdjHits( ColHits0,Clusters0,hAdjHits,hAdjHitsADCInt,false);
  CalcAdjHits( ColHits1,Clusters1,hAdjHits,hAdjHitsADCInt,false);
  CalcAdjHits( ColHits2,Clusters2,hAdjHits,hAdjHitsADCInt,false);
  CalcAdjHits( ColHits3,Clusters3,hAdjHits,hAdjHitsADCInt,false);

  std::cout << "# Clusters from the hits = " << Clusters0.size() << ", " << Clusters1.size() << ", " << Clusters2.size() << ", " << Clusters3.size() << std::endl;
  
  std::vector< std::vector< std::vector<float>>> ClGenPur = {{},{},{}};
  std::vector< std::vector< std::vector<recob::Hit>>> AllPlaneClusters = {Clusters0,Clusters1,Clusters2};
  std::vector< std::vector< float>> ClTotChrg = {{},{},{}}, ClMaxChrg = {{},{},{}}, CldT      = {{},{},{}}, ClT        = {{},{},{}}, ClX   = {{},{},{}}, ClY    = {{},{},{}}, ClZ = {{},{},{}};
  std::vector< std::vector< float>> ClFracE   = {{},{},{}}, ClFracGa  = {{},{},{}}, ClFracNe  = {{},{},{}}, ClFracRest = {{},{},{}}, ClPur = {{},{},{}}, Cldzdy = {{},{},{}};
  std::vector< std::vector< int  >> ClMainID  = {{},{},{}}, ClTPC     = {{},{},{}}, ClNHits   = {{},{},{}}, ClGen      = {{},{},{}};

  //------------------------------------------------------------ First complete cluster analysis ------------------------------------------------------------------// 
  // --- Now loop over the planes and the clusters to calculate the cluster properties
  for (int idx = 0; idx < 3; idx++){ 
    int nhit, clustTPC;
    float FracE, FracGa, FracNe, FracRest, clustX, clustY, clustZ, clustT, ncharge, maxHit, dzdy;
    std::vector< std::vector<recob::Hit> > Clusters = AllPlaneClusters[idx];
    
    // --- Loop over the clusters
    for (int i = 0; i < int(Clusters.size()); i++){ 
      int MainTrID = 0;
      int Gen = 0; float Pur = 0;
      std::vector<float> thisdzdy = {};

      nhit = Clusters[i].size();
      ncharge = maxHit = clustT = FracE = FracGa = FracNe = FracRest = clustX = clustY = clustZ = clustTPC = dzdy = 0;
      std::vector<float> GenPur = {}; for (size_t genpur = 0; genpur < fLabels.size(); genpur++){GenPur.push_back(0);}
      
      for (recob::Hit hit : Clusters[i]){
        ncharge += hit.Integral();
        const geo::WireGeo* wire = geo->GeometryCore::WirePtr(hit.WireID()); // Wire directions should be the same for all hits of the same view (can be used to check)
        double hitCharge;
          
        geo::Point_t hXYZ = wire->GetCenter();
        geo::Point_t sXYZ = wire->GetStart();
        geo::Point_t eXYZ = wire->GetEnd();

        geo::Vector_t direction = eXYZ - sXYZ;
        auto dyds = direction.Y(), dzds = direction.Z();
        thisdzdy.push_back(dzds/dyds);
        
        // Choose the position reco method to use for the cluster reco
        if (sPositionReco == "DEFAULT") {if (fDebug) std::cout << "Using default position reco method" << std::endl;}
        else {if (fDebug) std::cout << "ERROR: Position reco method not recognised. Defaulting to standard." << std::endl;}

        int TPC = hit.WireID().TPC;
        clustTPC += hit.Integral() * TPC; 
        clustX += hit.Integral() * hXYZ.X(); clustY += hit.Integral() * hXYZ.Y();clustZ += hit.Integral() * hXYZ.Z();clustT += hit.Integral() * hit.PeakTime();

        if (hit.Integral()>maxHit) {maxHit = hit.Integral();} // Look for maxHit inside cluster
        
        MainTrID = 0; double TopEFrac = -DBL_MAX;
        std::vector< sim::TrackIDE > ThisHitIDE = bt_serv->HitToTrackIDEs(clockData, hit);
        
        for (size_t ideL=0; ideL < ThisHitIDE.size(); ++ideL){
          if (ThisHitIDE[ideL].energyFrac > TopEFrac){
            TopEFrac = ThisHitIDE[ideL].energyFrac;
            MainTrID = ThisHitIDE[ideL].trackID; 
            if (fDebug) std::cout << "This hit's IDE is: " << MainTrID << std::endl; 
          }
        }

        for (int frac=0; frac < int(ClPartTrackIDs.size()); ++frac){
          for (int trck=0; trck < int(ClPartTrackIDs[frac].size()); ++trck){
            if (abs(MainTrID) == ClPartTrackIDs[frac][trck]){
              if (frac == 0){FracE = FracE + hit.Integral();}
              if (frac == 1){FracGa = FracGa + hit.Integral();}
              if (frac == 2){FracNe = FracNe + hit.Integral();}
              if (frac == 3){FracRest = FracRest + hit.Integral();}
            }
          }    
        }
        
        long unsigned int ThisPType = WhichParType(abs(MainTrID));
        GenPur[int(ThisPType)] = GenPur[int(ThisPType)] + hit.Integral();
        if (fDebug) std::cout << "\nThis particle type " << ThisPType << "\nThis cluster's main track ID " << MainTrID;    
        if (ThisPType == 1){hitCharge = hit.Integral();Pur = Pur+hitCharge;}
      }

      float MaxGen = 0;
      for (size_t genpur = 0; genpur < GenPur.size(); genpur++){
        if (GenPur[genpur] > MaxGen){Gen = genpur+1;}
        GenPur[genpur] = GenPur[genpur]/ncharge;
      }

      for (size_t j = 0; j > thisdzdy.size(); j++) {if (thisdzdy[0] != thisdzdy[i]) std::cout << "MISSMATCH IN dzdy FOR CLUSTER " << idx << std::endl;}

      dzdy = thisdzdy[0]; thisdzdy.clear();      
      FracE /= ncharge; FracGa /= ncharge; FracNe /= ncharge; FracRest /= ncharge;
      clustTPC /= ncharge; clustX /= ncharge; clustY /= ncharge;clustZ /= ncharge;clustT /= ncharge;
      if (fDebug) std::cout << "\ndzdy " << dzdy << " for cluster " << " (" << clustY << ", " << clustZ << ") with track ID " << MainTrID <<  " in plane " << idx << std::endl; 
      
      ClTotChrg[idx].push_back(ncharge);
      ClMaxChrg[idx].push_back(maxHit);
      ClNHits[idx].push_back(nhit);
      ClT[idx].push_back(clustT);
      ClTPC[idx].push_back(int(clustTPC));
      ClX[idx].push_back(clustX);
      ClY[idx].push_back(clustY);
      ClZ[idx].push_back(clustZ);
      ClFracE[idx].push_back(FracE);
      ClFracGa[idx].push_back(FracGa);
      ClFracNe[idx].push_back(FracNe);
      ClFracRest[idx].push_back(FracRest);
      ClPur[idx].push_back(Pur/ncharge);
      ClGen[idx].push_back(Gen);
      Cldzdy[idx].push_back(dzdy);
      ClMainID[idx].push_back(MainTrID);
      ClGenPur[idx].push_back(GenPur);

      if (fDebug) std::cout << "\nCluster " << i << " in plane " << idx << " has " << nhit << " hits, " << ncharge << " charge, " << clustT << " time, " << clustX << " X, " << clustY << " Y, " << clustZ << " Z, " << FracE << " FracE, " << FracGa << " FracGa, " << FracNe << " FracNe, " << FracRest << " FracRest, " << Pur/ncharge << " Pur, " << Gen << " gen, " << dzdy << " dzdy, " << MainTrID << " MainTrID" << std::endl;
    }
  } // Finished first cluster processing
  std::cout << "\nLooking for matching clusters: " << std::endl;
  
  std::cout << std::endl;//-------------------------------------------------------------------- Cluster Matching -------------------------------------------------------------------------// 
  std::vector<std::vector<float>> MVecGenFrac = {};
  std::vector<int>   MVecNHit  = {}, MVecGen    = {}, MVecInd0NHits  = {}, MVecInd1NHits  = {}, MVecMainID = {}, MVecInd0TPC = {}, MVecInd1TPC = {};
  std::vector<float> MVecTime  = {}, MVecChrg   = {}, MVecInd0MaxHit = {}, MVecInd1MaxHit = {}, MVecInd0dT = {}, MVecInd1dT = {};
  std::vector<float> MVecInd0RecoY = {}, MVecInd1RecoY = {}, MVecRecY = {}, MVecRecZ = {};
  std::vector<float> MVecFracE = {}, MVecFracGa = {}, MVecFracNe = {}, MVecFracRest = {}, MVecPur = {};

  for (int ii = 0; ii < int(AllPlaneClusters[2].size()); ii++){    
    bool match = false;
    int ind0clustNHits = 0, ind1clustNHits = 0;
    double ind0clustY = -1e6, ind1clustY = -1e6, ind0clustMaxHit = 0, ind1clustMaxHit = 0;
    double ind0clustdT = fClusterMatchTime, ind1clustdT = fClusterMatchTime;
    // std::cout << "Test" << std::endl;
    if (!AllPlaneClusters[2][ii].empty()){
      if (fDebug) std::cout << "***Evaluating main cluster " << ii << " with charge " << ClTotChrg[2][ii] << "This cluster's purity is " << ClPur[2][ii]*100 << "%" << std::endl;
      if (fDebug) std::cout << "This cluster's position is (" << ClY[2][ii] << ", " << ClZ[2][ii] << ") and its time is " << ClT[2][ii] << std::endl;

      if (!AllPlaneClusters[0].empty()){
        for (int jj = 0; jj < int(AllPlaneClusters[0].size()); jj++){
          if (fDebug) std::cout << "Evaluating cluster " << jj << " from plane " << 0 << " with charge " << ClTotChrg[0][jj] << "This cluster's purity is " << ClPur[0][jj]*100 << "%" << std::endl;
          if (fDebug) std::cout << "This cluster's position is (" << ClY[0][jj] << ", " << ClZ[0][jj] << ") and its time is " << ClT[0][jj] << std::endl;
          if (abs(ClT[2][ii] - ClT[0][jj]) < fClusterMatchTime && abs(fGoalInd0MatchTime - abs(ClT[2][ii] - ClT[0][jj])) < abs(fGoalInd0MatchTime - ind0clustdT)){
            ind0clustY = ClY[0][jj] + (ClZ[2][ii] - ClZ[0][jj])/(Cldzdy[0][jj]);
            ind0clustdT = abs(ClT[2][ii] - ClT[0][jj]);
            ind0clustNHits = int(AllPlaneClusters[0][jj].size());
            ind0clustMaxHit = ClMaxChrg[0][jj];
            if (ind0clustY > -fMaxDetSizeY && ind0clustY < fMaxDetSizeY){match = true;}
            if (fDebug) std::cout << "¡¡¡ Matched cluster in plane 0 !!! --- Position x = " << ClX[0][jj] << ", y = " << ClY[0][jj] << ", z = " << ClZ[0][jj] << std::endl;
            if (fDebug) std::cout << "Reconstructed position y = " << ind0clustY << ", z = " << ClZ[2][ii] << std::endl;  
          }
        }
      }
      if (!AllPlaneClusters[1].empty()){
        for (int zz = 0; zz < int(AllPlaneClusters[1].size()); zz++){
          if (abs(ClT[2][ii] - ClT[1][zz]) < fClusterMatchTime && abs(fGoalInd1MatchTime - abs(ClT[2][ii] - ClT[1][zz])) < abs(fGoalInd1MatchTime - ind1clustdT)){
            ind1clustY = ClY[1][zz] + (ClZ[2][ii] - ClZ[1][zz])/(Cldzdy[1][zz]);
            ind1clustdT = abs(ClT[2][ii] - ClT[1][zz]);
            ind1clustNHits = int(AllPlaneClusters[1][zz].size());
            ind1clustMaxHit = ClMaxChrg[1][zz];
            if (ind1clustY > -fMaxDetSizeY && ind1clustY < fMaxDetSizeY){match = true;}
            if (fDebug) std::cout << "¡¡¡ Matched cluster in plane 1 !!! --- Position x = " << ClX[1][zz] << ", y = " << ClY[1][zz] << ", z = " << ClZ[1][zz] << std::endl;
            if (fDebug) std::cout << "Reconstructed position y = " << ind1clustY << ", z = " << ClZ[2][ii] << std::endl;
          }
        } // Loop over ind1 clusters
      }
    } // Loop over ind clusters
    else {std::cout << "Cluster " << ii << " in plane 2 has no hits" << std::endl;}
    
    //--------------------------------------------------------- Export Matched cluster vectors ------------------------------------------------------------------// 
    if (match == true){
      MVecTime.push_back(ClT[2][ii]);
      MVecChrg.push_back(ClTotChrg[2][ii]);
      MVecNHit.push_back(ClNHits[2][ii]);
      MVecInd0TPC.push_back(ClTPC[0][ii]);
      MVecInd1TPC.push_back(ClTPC[1][ii]);
      MVecInd0dT.push_back(ind0clustdT);
      MVecInd1dT.push_back(ind1clustdT);
      MVecInd0RecoY.push_back(ind0clustY);
      MVecInd1RecoY.push_back(ind1clustY);
      MVecInd0MaxHit.push_back(ind0clustMaxHit);
      MVecInd1MaxHit.push_back(ind1clustMaxHit);
      MVecInd0NHits.push_back(ind0clustNHits);
      MVecInd1NHits.push_back(ind1clustNHits);
      MVecRecZ.push_back(ClZ[2][ii]);
      MVecFracE.push_back(ClFracE[2][ii]);
      MVecFracGa.push_back(ClFracGa[2][ii]);
      MVecFracNe.push_back(ClFracNe[2][ii]);
      MVecFracRest.push_back(ClFracRest[2][ii]);
      MVecPur.push_back(ClPur[2][ii]);
      MVecGen.push_back(ClGen[2][ii]);
      MVecMainID.push_back(ClMainID[2][ii]);
      MVecGenFrac.push_back(ClGenPur[2][ii]);
      
      float buffer = 1;
      if ((ind0clustY > -buffer*fMaxDetSizeY && ind0clustY < buffer*fMaxDetSizeY) && (ind1clustY > -buffer*fMaxDetSizeY && ind1clustY < buffer*fMaxDetSizeY)){
        if (fDebug) std::cout << "BOTH IND RECO INSIDE OF DETECTOR" << std::endl;
        MVecRecY.push_back((ind0clustY+ind1clustY)/2);}
      else if (ind0clustY > -buffer*fMaxDetSizeY && ind0clustY < buffer*fMaxDetSizeY){
        if (fDebug) std::cout << "IND1 OUTSIDE OF DETECTOR" << std::endl;
        MVecRecY.push_back(ind0clustY);}
      else if (ind1clustY > -buffer*fMaxDetSizeY && ind1clustY < buffer*fMaxDetSizeY){
        if (fDebug) std::cout << "IND0 OUTSIDE OF DETECTOR" << std::endl;
        MVecRecY.push_back(ind1clustY);}
      else{
        if (fDebug) std::cout << "RECO OUTSIDE OF DETECTOR" << std::endl;
        MVecRecY.push_back((ind0clustY+ind1clustY)/2); 
        if (ClGen[2][ii] == 1){PrintInColor("Marley cluster recon structed outside of detector volume! RecoY = " + str((ind0clustY+ind1clustY)/2), GetColor("red"));}
      }

      // Print in color if the cluster is matched
      if (fDebug){
        std::string MatchedColor = "white";
        if (ClNHits[2][ii] > fClusterMatchMinNHit && (ind0clustNHits > fClusterMatchMinNHit || ind1clustNHits > fClusterMatchMinNHit)){MatchedColor = "blue";}
        else{MatchedColor = "white";}
        PrintInColor("¡¡¡ Matched clusters !!! ", GetColor(MatchedColor));
        PrintInColor(" - Cluster " + str(ClMainID[2][ii]) + " Gen " + str(ClGen[2][ii]) + " Purity " + str(ClPur[2][ii]) + " Hits " + str(ClNHits[2][ii]), GetColor(MatchedColor));
        PrintInColor(" - Hits(ind0, ind1, col) " + str(ind0clustNHits) + ", " + str(ind1clustNHits) + ", " + str(ClNHits[2][ii]), GetColor(MatchedColor));
        PrintInColor(" - Positions y(ind0, ind1) = " + str(ind0clustY) + ", " + str(ind1clustY) + ", z = " + str(ClZ[2][ii]) + "\n", GetColor(MatchedColor));
      }    
    } // if (match == true)
  } // Loop over collection plane clusters
  
  std::cout << std::endl;//-------------------------------------------------------------------- Cluster Tree Export -------------------------------------------------------------------------// 
  // Loop over matched clusters and export to tree if number of hits is above threshold
  for (int i = 0; i < int(MVecNHit.size()); i++){
    if(MVecNHit[i] > fClusterMatchMinNHit && (MVecInd0NHits[i] > fClusterMatchMinNHit || MVecInd1NHits[i] > fClusterMatchMinNHit)){
      MAdjClTime = {};MAdjClCharge = {};MAdjClNHit = {};MAdjClRecoY = {};MAdjClRecoZ = {};MAdjClR = {};MAdjClPur = {};MAdjClGen = {};
      MAdjFlashTime = {};MAdjFlashPE = {};MAdjFlashNHit = {};MAdjFlashMaxPE = {};MAdjFlashRecoY = {};MAdjFlashRecoZ = {};MAdjFlashR = {};MAdjFlashPur = {};
      
      std::string ResultColor = "white";
      PrintInColor("*** Matched preselection cluster: ",GetColor("blue"));
      if (abs(MVecRecY[i] - TNuY) < 10 && abs(MVecRecZ[i] - TNuZ) < 10) {ResultColor = "green";PrintInColor(" - Cluster is close to neutrino vertex!",GetColor(ResultColor));}
      else {ResultColor = "yellow";PrintInColor(" - Cluster is NOT close to neutrino vertex!",GetColor(ResultColor));}

      PrintInColor(" - Cluster " + str(MVecMainID[i]) + " Gen " + str(MVecGen[i]) + " Purity " + str(MVecPur[i]) + " Hits " + str(MVecNHit[i]), GetColor(ResultColor));
      PrintInColor(" - RecoY, Z (" + str(MVecRecY[i]) + ", " + str(MVecRecZ[i]) + ") Time " + str(MVecTime[i]) + "\n", GetColor(ResultColor));
      
      // Loop over collection plane clusters to find adjacent clusters with distance < fAdjClusterRad and time < fAdjClusterTime
      for (int j = 0; j < int(MVecNHit.size()); j++){
        if(j != i && sqrt(pow(MVecRecY[i]-MVecRecY[j],2)+pow(MVecRecZ[i]-MVecRecZ[j],2)) < fAdjClusterRad && abs(MVecTime[i]-MVecTime[j]) < fAdjClusterTime){
          MAdjClTime.push_back(MVecTime[j]);
          MAdjClCharge.push_back(MVecChrg[j]);
          MAdjClNHit.push_back(MVecNHit[j]);
          MAdjClRecoY.push_back(MVecRecY[j]);
          MAdjClRecoZ.push_back(MVecRecZ[j]);
          MAdjClR.push_back(sqrt(pow(MVecRecY[i]-MVecRecY[j],2)+pow(MVecRecZ[i]-MVecRecZ[j],2)));
          MAdjClPur.push_back(MVecPur[j]);
          MAdjClGen.push_back(MVecGen[j]);
        }
      }

      // Loop over optical flashes to find adjacent flashes with distance < fAdjOpFlashRad and time < fAdjOpFlashTime
      for (int j = 0; j < int(OpFlashPE.size()); j++){
        if(sqrt(pow(MVecRecY[i]-OpFlashY[j],2)+pow(MVecRecZ[i]-OpFlashZ[j],2)) < fAdjOpFlashRad  && MVecTime[i]/2 - OpFlashT[j] < fAdjOpFlashTime && MVecTime[i]/2 - OpFlashT[j] > 0){
          MAdjFlashTime.push_back(OpFlashT[j]);
          MAdjFlashPE.push_back(OpFlashPE[j]);
          MAdjFlashNHit.push_back(OpFlashNHit[j]);
          MAdjFlashMaxPE.push_back(OpFlashMaxPE[j]);
          MAdjFlashRecoY.push_back(OpFlashY[j]);
          MAdjFlashRecoZ.push_back(OpFlashZ[j]);
          MAdjFlashR.push_back(sqrt(pow(MVecRecY[i]-OpFlashY[j],2)+pow(MVecRecZ[i]-OpFlashZ[j],2)));
          MAdjFlashPur.push_back(OpFlashMarlPur[j]);
        }
      }

      // Fill the tree with the cluster information and the adjacent clusters and flashes
      MMarleyFrac =     {MVecFracE[i],MVecFracGa[i],MVecFracNe[i],MVecFracRest[i]};
      MGenFrac =        MVecGenFrac[i];
      MTime =           MVecTime[i];   
      MChrg =           MVecChrg[i];   
      MNHit =           MVecNHit[i];
      MInd0TPC =        MVecInd0TPC[i];
      MInd1TPC =        MVecInd1TPC[i];
      MInd0MaxHit =     MVecInd0MaxHit[i];   
      MInd1MaxHit =     MVecInd1MaxHit[i];
      MInd0NHits =      MVecInd0NHits[i];   
      MInd1NHits =      MVecInd1NHits[i];
      MInd0dT =         MVecInd0dT[i];   
      MInd1dT =         MVecInd1dT[i];   
      MInd0RecoY =      MVecInd0RecoY[i];   
      MInd1RecoY =      MVecInd1RecoY[i];   
      MRecZ =           MVecRecZ[i];   
      MPur =            MVecPur[i];
      MGen =            MVecGen[i];
      // std::cout << "MGenLabel = " << fLabels[MVecGen[i]] << std::endl;
      // MGenLabel =       fLabels[int(MGen)];
      MMainID =         MVecMainID[i];

      hDriftTime->      Fill(avX,MTime);
      fSolarNuAnaTree-> Fill();
      hXTruth->         Fill(MVecRecY[i]-TNuY,TNuX); 
      hYTruth->         Fill(MVecRecY[i]-TNuY,TNuY); 
      hZTruth->         Fill(MVecRecY[i]-TNuY,TNuZ); 
    }
  }
  std::cout << std::endl;
}

//########################################################################################################################################//
//_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_//
//########################################################################################################################################//

//......................................................
void SolarNuAna::ResetVariables()
// Reset variables for each event
{
  // Clear Marley MCTruth info.
  TNuE  = 0; TNuX = 0; TNuY = 0; TNuZ = 0;
  MarleyPDGList = {}; MarleyEList = {}; MarleyPList = {}; MarleyXList = {}; MarleyYList = {}; MarleyZList = {}, MarleyIDList = {}, MarleyParentIDList = {};
  TPart = {}; Parts = {};
  // ThisGeneratorParts.clear();
  OpFlashMarlPur.clear();OpFlashPE.clear();OpFlashMaxPE.clear();OpFlashY.clear();OpFlashZ.clear();OpFlashT.clear();OpFlashDeltaT.clear();OpFlashNHit.clear();

} // ResetVariables

void SolarNuAna::CalcAdjHits( std::vector< recob::Hit > MyVec,std::vector< std::vector<recob::Hit> >& Clusters,TH1I* MyHist, TH1F* ADCIntHist, bool HeavDebug ) 
/* 
Find adjacent hits in time and space:
- MyVec is the vector of hits to be clustered
- Clusters is the vector of clusters
- MyHist is the histogram to be filled with the number of hits in each cluster
- ADCIntHist is the histogram to be filled with the ADC integral of each cluster
- HeavDebug is a boolean to turn on/off debugging statements
*/
{
  const double TimeRange  = 10;
  const int    ChanRange  = 2;
  unsigned int FilledHits = 0;
  unsigned int NumOriHits = MyVec.size();

  while( NumOriHits != FilledHits ) 
  {
    if (HeavDebug) std::cerr << "\nStart of my while loop" << std::endl;
    std::vector< recob::Hit > AdjHitVec;
    AdjHitVec.push_back ( MyVec[0] );
    MyVec.erase( MyVec.begin()+0 );
    int LastSize = 0;
    int NewSize  = AdjHitVec.size();
    
    while ( LastSize != NewSize ) 
    {
      std::vector<int> AddNow;
      for (size_t aL=0; aL < AdjHitVec.size(); ++aL) 
      {
        for (size_t nL=0; nL < MyVec.size(); ++nL) 
        {
	        if (HeavDebug) 
          {
            std::cerr << "\t\tLooping though AdjVec " << aL << " and  MyVec " << nL
            << " AdjHitVec - " << AdjHitVec[aL].Channel() << " & " << AdjHitVec[aL].PeakTime()
            << " MVec - " << MyVec[nL].Channel() << " & " << MyVec[nL].PeakTime()
            << " Channel " << abs( (int)AdjHitVec[aL].Channel()  - (int)MyVec[nL].Channel()  )  << " bool " << (bool)(abs( (int)AdjHitVec[aL].Channel() - (int)MyVec[nL].Channel()  ) <= ChanRange)
            << " Time " << abs( AdjHitVec[aL].PeakTime() - MyVec[nL].PeakTime() ) << " bool " << (bool)(abs( (double)AdjHitVec[aL].PeakTime() - (double)MyVec[nL].PeakTime() ) <= TimeRange)
            << std::endl;
	        }
	        
          if ( abs( (int)AdjHitVec[aL].Channel()  - (int)MyVec[nL].Channel()  ) <= ChanRange &&
	        abs( (double)AdjHitVec[aL].PeakTime() - (double)MyVec[nL].PeakTime() ) <= TimeRange )
          {
	      
            if (HeavDebug) std::cerr << "\t\t\tFound a new thing!!!" << std::endl;
	          // --- Check that this element isn't already in AddNow.
	          bool AlreadyPres = false;
	          
            for (size_t zz=0; zz<AddNow.size(); ++zz) 
            {
	            if (AddNow[zz] == (int)nL) AlreadyPres = true;
	          }
	          
            if (!AlreadyPres)
	          AddNow.push_back( nL );
	        } // If this hit is within the window around one of my other hits.
	      } // Loop through my vector of colleciton plane hits.
      } // Loop through AdjHitVec

      // --- Now loop through AddNow and remove from Marley whilst adding to AdjHitVec
      std::sort(AddNow.begin(),AddNow.end());
      for (size_t aa=0; aa<AddNow.size(); ++aa) 
      { 
	      if (HeavDebug) 
        {
	        std::cerr << "\tRemoving element " << AddNow.size()-1-aa << " from MyVec ===> "
		      << MyVec[ AddNow[AddNow.size()-1-aa] ].Channel() << " & " << MyVec[ AddNow[AddNow.size()-1-aa] ].PeakTime()
		      << std::endl;
	      }

        AdjHitVec.push_back ( MyVec[ AddNow[AddNow.size()-1-aa] ] );
	      MyVec.erase( MyVec.begin() + AddNow[AddNow.size()-1-aa] ); // This line creates segmentation fault
	      // std::cout << "Erase works" << std::endl;
      }

      LastSize = NewSize;
      NewSize  = AdjHitVec.size();
      if (HeavDebug) 
      {
	      std::cerr << "\t---After that pass, AddNow was size " << AddNow.size() << " ==> LastSize is " << LastSize << ", and NewSize is " << NewSize
		    << "\nLets see what is in AdjHitVec...." << std::endl;
	      for (size_t aL=0; aL < AdjHitVec.size(); ++aL) 
        {
	        std::cout << "\tElement " << aL << " is ===> " << AdjHitVec[aL].Channel() << " & " << AdjHitVec[aL].PeakTime() << std::endl;
	      }
      }
    } // while ( LastSize != NewSize )

    int NumAdjColHits = AdjHitVec.size();
    float SummedADCInt = 0;
    for ( recob::Hit hit : AdjHitVec) SummedADCInt += hit.Integral();

    if (HeavDebug) std::cerr << "After that loop, I had " << NumAdjColHits << " adjacent collection plane hits." << std::endl;
    
    MyHist -> Fill( NumAdjColHits );
    ADCIntHist -> Fill( SummedADCInt );
    FilledHits += NumAdjColHits;
    
    if (AdjHitVec.size() > 0) Clusters.push_back(AdjHitVec);
  }

  if (HeavDebug)
  {
    std::vector<double> avgChannel;
    std::vector<double> avgTick;
    std::vector<double> summedADCInt;

    for (std::vector< recob::Hit > hits : Clusters)
    {
      double adcInt = 0;
      double channel = 0;
      double tick = 0;

      for (recob::Hit hit : hits)
      {
        tick += hit.Integral()*hit.PeakTime();
        channel += hit.Integral()*hit.Channel();
        adcInt += hit.Integral();
      }
      tick /= adcInt;
      channel /= adcInt;
      summedADCInt.push_back(adcInt);
      avgTick.push_back(tick);
      avgChannel.push_back(channel);
    }

    for (int i = 0; i < int(avgTick.size()-1); i++)
    {
      for (int j = i+1; j < int(avgTick.size()); j++)
      {
        std::cout << avgChannel[i] << " " << avgChannel[j] << "  " << std::abs(avgChannel[i]-avgChannel[j]) << std::endl;
        std::cout << avgTick[i] << " " << avgTick[j] << "  " << std::abs(avgTick[i]-avgTick[j]) << std::endl;
        std::cout << summedADCInt[i] << " " << summedADCInt[j] << std::endl;
      }
    }
  }
  return;
}

//...................................................... 
long unsigned int SolarNuAna::WhichParType( int TrID ) {
  for (long unsigned int i = 0; i < fLabels.size(); i++){
    if (InMyMap(TrID,Parts[i])) {return i+1;}
  }
  // If no match, then who knows???
  return 0;
}

//......................................................
// This function fills a map with the MCParticles from a given MCTruth
void SolarNuAna::FillMyMaps( std::map< int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand ){
  for ( size_t L1=0; L1 < Hand->size(); ++L1 ) {
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
      MyMap[ThisPar.TrackId()] = ThisPar;
      if (fDebug) std::cout << ThisPar.PdgCode() << " " << ThisPar.E() << std::endl;
    }
  }
  return;
}

//......................................................
// This function checks if a given TrackID is in a given map
bool SolarNuAna::InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap ){
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find( TrID );
  if ( ParIt != ParMap.end() ) 
  {return true;}
  else return false;
}

//......................................................
// This function creates a terminal color printout
void SolarNuAna::PrintInColor( std::string MyString, int Color ){
  std::cout << "\033[" << Color << "m" << MyString << "\033[0m" << std::endl;
  return;
}

// ......................................................
// This function returns an integer that corresponds to a given color name
int SolarNuAna::GetColor( std::string ColorName ){
  if (ColorName == "black") return 30;
  else if (ColorName == "red") return 31;
  else if (ColorName == "green") return 32;
  else if (ColorName == "yellow") return 33;
  else if (ColorName == "blue") return 34;
  else if (ColorName == "magenta") return 35;
  else if (ColorName == "cyan") return 36;
  else if (ColorName == "white") return 37;
  else {std::cout << "Color " << ColorName << " not recognized. Returning white." << std::endl; return 37;}
  return 0;
}

std::string SolarNuAna::str( int i ) {std::stringstream ss;ss << i;return ss.str();}
std::string SolarNuAna::str( double i ) {std::stringstream ss;ss << i;return ss.str();}
std::string SolarNuAna::str( float i ) {std::stringstream ss;ss << i;return ss.str();}

DEFINE_ART_MODULE(SolarNuAna)