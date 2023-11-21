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
#include <fcntl.h>

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
  void              PrintInColor  ( std::string MyString, int MyColor, std::string Type = "Info");
  int               GetColor      ( std::string MyString );
  std::string       str           ( int MyInt );
  std::string       str           ( float MyFloat );
  std::string       str           ( double MyDouble );
  std::string       str           ( std::vector<int> MyVec );
  std::string       str           ( std::vector<float> MyVec );
  std::string       str           ( std::vector<double> MyVec );
  int supress_stdout();
  void resume_stdout(int fd);
  // --- Our fcl parameter labels for the modules that made the data products
  std::string fRawDigitLabel,fHitLabel,fOpHitLabel,fOpDetWaveformLabel,fOpFlashLabel,fGEANTLabel; 

  // --- Input settings imported from the fcl
  std::string fGeometry;
  int fDetectorSizeX, fDetectorSizeY, fClusterInd0MatchTime, fClusterInd1MatchTime,fClusterPreselectionNHit;
  float fClusterMatchTime,fAdjClusterTime,fAdjClusterRad,fClusterMatchCharge,fAdjOpFlashRad,fAdjOpFlashTime,fAdjOpFlashMaxPECut,fAdjOpFlashMinPECut,fClusterMatchNHit;
  std::vector<std::string> fLabels;
  
  // --- Our TTrees, and its associated variables.
  TTree* fSolarNuAnaTree;
  TTree* fMCTruthTree;
  std::string MGenLabel;
  int Event,Flag,MNHit,MGen,MTPC,MInd0TPC,MInd1TPC,MInd0NHits,MInd1NHits,MMainID,MMainT,MMainPDG,MMainParentPDG;
  float TNuQSqr,TNuE,TNuP,TNuX,TNuY,TNuZ,avX,avY,avZ,MTime,MCharge,MMaxCharge,MInd0Charge,MInd1Charge,MInd0MaxCharge,MInd1MaxCharge,MInd0dT,MInd1dT,MInd0RecoY,MInd1RecoY,MRecZ,MPur,MMainE,MMainP,MMainParentE,MMainParentP,MMainParentT;
  std::vector<int> MAdjClGen,MAdjClMainID,TPart,MarleyPDGList,MarleyIDList,MarleyParentIDList,MAdjClMainPDG;
  std::vector<float> MAdjClTime,MAdjClCharge,MAdjClInd0Charge,MAdjClInd1Charge,MAdjClMaxCharge,MAdjClInd0MaxCharge,MAdjClInd1MaxCharge,MAdjClNHit,MAdjClInd0NHit,MAdjClInd1NHit,MAdjClRecoY,MAdjClRecoZ,MAdjClR,MAdjClPur,MAdjClMainE,MAdjClMainX,MAdjClMainY,MAdjClMainZ,MMarleyFrac,MGenFrac;
  std::vector<float> MAdjFlashTime,MAdjFlashPE,MAdjFlashNHit,MAdjFlashMaxPE,MAdjFlashRecoY,MAdjFlashRecoZ,MAdjFlashR,MAdjFlashPur;
  std::vector<float> MarleyEList,MarleyPList,MarleyXList,MarleyYList,MarleyZList;
  std::vector<double> MMainVertex,MEndVertex,MMainParentVertex;
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

  fGeometry                = p.get<std::string> ("Geometry");
  fDetectorSizeX           = p.get<int>         ("DetectorSizeX");
  fDetectorSizeY           = p.get<int>         ("DetectorSizeY");
  fClusterMatchNHit        = p.get<float>       ("ClusterMatchNHit");
  fClusterMatchCharge      = p.get<float>       ("ClusterMatchCharge");
  fClusterMatchTime        = p.get<float>       ("ClusterMatchTime");
  fClusterInd0MatchTime    = p.get<float>       ("ClusterInd0MatchTime");
  fClusterInd1MatchTime    = p.get<float>       ("ClusterInd1MatchTime");
  
  fClusterPreselectionNHit = p.get<int>         ("ClusterPreselectionNHit");
  
  fAdjClusterTime          = p.get<float>       ("AdjClusterTime");
  fAdjClusterRad           = p.get<float>       ("AdjClusterRad");
  
  fAdjOpFlashTime          = p.get<float>       ("AdjOpFlashTime");
  fAdjOpFlashRad           = p.get<float>       ("AdjOpFlashRad");
  fAdjOpFlashMaxPECut      = p.get<float>       ("AdjOpFlashMaxPECut");
  fAdjOpFlashMinPECut      = p.get<float>       ("AdjOpFlashMinPECut");
} // Reconfigure

//......................................................
void SolarNuAna::beginJob(){   
  // --- Make our handle to the TFileService
  art::ServiceHandle<art::TFileService> tfs;
  fMCTruthTree = tfs->make<TTree>("MCTruthTree","MC Truth Tree");
  fSolarNuAnaTree = tfs->make<TTree>("SolarNuAnaTree","Solar Ana Tree");

  // MC Truth info.
  // fMCTruthTree -> Branch("Run",               &Run,              "Run/I");
  // fMCTruthTree -> Branch("SubRun",            &SubRun,           "SubRun/I");
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
  // fSolarNuAnaTree -> Branch("Run",            &Run,              "Run/I");
  // fSolarNuAnaTree -> Branch("SubRun",         &SubRun,           "SubRun/I");
  fSolarNuAnaTree -> Branch("Event",           &Event,            "Event/I");
  fSolarNuAnaTree -> Branch("Flag",            &Flag,             "Flag/I");          // Flag used to match truth with reco tree entries
  fSolarNuAnaTree -> Branch("TruthPart",       &TPart);                               // Number particles per generator
  fSolarNuAnaTree -> Branch("TNuQSqr",         &TNuQSqr,          "TruthNuQSqr/F");   // True neutrino transfer momentum [GeV]
  fSolarNuAnaTree -> Branch("TNuE",            &TNuE,             "TruthNuE/F");      // True neutrino energy
  fSolarNuAnaTree -> Branch("TNuP",            &TNuP,             "TruthNuP/F");      // True neutrino momentum
  fSolarNuAnaTree -> Branch("TNuX",            &TNuX,             "TruthNuX/F");      // True neutrino X
  fSolarNuAnaTree -> Branch("TNuY",            &TNuY,             "TruthNuY/F");      // True neutrino Y
  fSolarNuAnaTree -> Branch("TNuZ",            &TNuZ,             "TruthNuZ/F");      // True neutrino Z
  fSolarNuAnaTree -> Branch("TMarleyPDG",      &MarleyPDGList);                       // PDG of marley particles
  fSolarNuAnaTree -> Branch("TMarleyE",        &MarleyEList);                         // Energy of marley particles
  fSolarNuAnaTree -> Branch("TMarleyP",        &MarleyPList);                         // Momentum of marley particles
  fSolarNuAnaTree -> Branch("TMarleyX",        &MarleyXList);                         // X of marley particles
  fSolarNuAnaTree -> Branch("TMarleyY",        &MarleyYList);                         // Y of marley particles
  fSolarNuAnaTree -> Branch("TMarleyZ",        &MarleyZList);                         // Z of marley particles
  fSolarNuAnaTree -> Branch("TMarleyID",       &MarleyIDList);                        // TrackID of marley particles")
  fSolarNuAnaTree -> Branch("TMarleyParentID", &MarleyParentIDList);                  // ParentID of marley particles
  
  // Main Cluster info.
  fSolarNuAnaTree -> Branch("Generator",        &MGen,             "Generator/I");     // Main cluster generator idx  
  fSolarNuAnaTree -> Branch("Purity",           &MPur,             "Purity/F");        // Main cluster reco purity
  fSolarNuAnaTree -> Branch("TPC",              &MTPC,             "ColTPC/I");        // Main cluster TPC
  fSolarNuAnaTree -> Branch("Time",             &MTime,            "ColTime/F");       // Main cluster time [ticks]
  fSolarNuAnaTree -> Branch("NHits",            &MNHit,            "ColNHits/I");      // Main cluster #hits
  fSolarNuAnaTree -> Branch("Charge",           &MCharge,          "ColCharge/F");     // Main cluster charge [ADC*ticks]
  fSolarNuAnaTree -> Branch("MaxCharge",        &MMaxCharge,       "ColCharge/F");     // Main cluster's max hit-charge [ADC*ticks]
  fSolarNuAnaTree -> Branch("RecoZ",            &MRecZ,            "RecoZ/F");         // Main cluster reco Z [cm]
  fSolarNuAnaTree -> Branch("Ind0TPC",          &MInd0TPC,         "Ind0TPC/I");       // Main cluster ind0 TPC
  fSolarNuAnaTree -> Branch("Ind1TPC",          &MInd1TPC,         "Ind1TPC/I");       // Main cluster ind1 TPC  
  fSolarNuAnaTree -> Branch("Ind0dT",           &MInd0dT,          "Ind0dT/F");        // Main cluster ind0 dT [Ticks]
  fSolarNuAnaTree -> Branch("Ind1dT",           &MInd1dT,          "Ind1dT/F");        // Main cluster ind1 dT [Ticks]
  fSolarNuAnaTree -> Branch("Ind0NHits",        &MInd0NHits,       "Ind0NHits/I");     // Main cluster ind0 Hits
  fSolarNuAnaTree -> Branch("Ind1NHits",        &MInd1NHits,       "Ind1NHits/I");     // Main cluster ind1 Hits
  fSolarNuAnaTree -> Branch("Ind0Charge",       &MInd0Charge,      "Ind0Charge/F");    // Main cluster ind0 MaxHit
  fSolarNuAnaTree -> Branch("Ind1Charge",       &MInd1Charge,      "Ind1Charge/F");    // Main cluster ind1 MaxHit
  fSolarNuAnaTree -> Branch("Ind0MaxCharge",    &MInd0MaxCharge,   "Ind0MaxCharge/F"); // Main cluster ind0 MaxHit
  fSolarNuAnaTree -> Branch("Ind1MaxCharge",    &MInd1MaxCharge,   "Ind1MaxCharge/F"); // Main cluster ind1 MaxHit
  fSolarNuAnaTree -> Branch("Ind0RecoY",        &MInd0RecoY,       "Ind0RecoY/F");     // Main cluster ind0 reco Y [cm]
  fSolarNuAnaTree -> Branch("Ind1RecoY",        &MInd1RecoY,       "Ind1RecoY/F");     // Main cluster ind1 reco Y [cm]
  fSolarNuAnaTree -> Branch("MainID",           &MMainID,          "MainID/I");        // Main cluster main track ID
  fSolarNuAnaTree -> Branch("MainT",            &MMainT,           "MainT/I");         // Main cluster main time [ticks]
  fSolarNuAnaTree -> Branch("MainE",            &MMainE,           "MainE/F");         // Main cluster main energy [GeV]
  fSolarNuAnaTree -> Branch("MainP",            &MMainP,           "MainP/F");         // Main cluster main momentum [GeV]
  fSolarNuAnaTree -> Branch("MainPDG",          &MMainPDG,         "MainPDG/I");       // Main cluster main pdg
  fSolarNuAnaTree -> Branch("MainParentPDG",    &MMainParentPDG,   "MainParentPDG/I"); // Main cluster main pdg
  fSolarNuAnaTree -> Branch("MainParentE",      &MMainParentE,     "MainParentE/F");   // Main cluster main parent energy [GeV]
  fSolarNuAnaTree -> Branch("MainParentP",      &MMainParentP,     "MainParentP/F");   // Main cluster main parent momentum [GeV]
  fSolarNuAnaTree -> Branch("MainParentT",      &MMainParentT,     "MainParentT/F");   // Main cluster main parent Time [ticks]
  fSolarNuAnaTree -> Branch("MainVertex",       &MMainVertex);                         // Main cluster main particle vertex [cm]
  fSolarNuAnaTree -> Branch("EndVertex",        &MEndVertex);                          // Main cluster end particle vertex [cm]
  fSolarNuAnaTree -> Branch("MainParentVertex", &MMainParentVertex);                   // Main cluster parent particle vertex [cm]
  fSolarNuAnaTree -> Branch("GenFrac",          &MGenFrac);                            // Main cluster reco purity complete
  fSolarNuAnaTree -> Branch("MarleyFrac",       &MMarleyFrac);                         // Main cluster particle contribution (electron, gamma, neutron)
  // fSolarNuAnaTree -> Branch("Label",          &MGenLabel);                          // Main cluster generator label  
      
  // Adj. Cluster info.
  fSolarNuAnaTree -> Branch("AdjClGen",           &MAdjClGen);                           // Adj. clusters' generator idx
  fSolarNuAnaTree -> Branch("AdjClNHit",          &MAdjClNHit);                          // Adj. clusters' #hits 
  fSolarNuAnaTree -> Branch("AdjClInd0NHit",      &MAdjClInd0NHit);                      // Adj. clusters' #hits 
  fSolarNuAnaTree -> Branch("AdjClInd1NHit",      &MAdjClInd1NHit);                      // Adj. clusters' #hits 
  fSolarNuAnaTree -> Branch("AdjClTime",          &MAdjClTime);                          // Adj. clusters' time [ticks]
  fSolarNuAnaTree -> Branch("AdjClCharge",        &MAdjClCharge);                        // Adj. clusters' charge [ADC*ticks]
  fSolarNuAnaTree -> Branch("AdjClInd0Charge",    &MAdjClInd0Charge);                    // Adj. clusters' charge [ADC*ticks]
  fSolarNuAnaTree -> Branch("AdjClInd1Charge",    &MAdjClInd1Charge);                    // Adj. clusters' charge [ADC*ticks]
  fSolarNuAnaTree -> Branch("AdjClMaxCharge",     &MAdjClMaxCharge);                     // Adj. clusters' charge [ADC*ticks]
  fSolarNuAnaTree -> Branch("AdjClInd0MaxCharge", &MAdjClInd0MaxCharge);                 // Adj. clusters' charge [ADC*ticks]
  fSolarNuAnaTree -> Branch("AdjClInd1MaxCharge", &MAdjClInd1MaxCharge);                 // Adj. clusters' charge [ADC*ticks]
  fSolarNuAnaTree -> Branch("AdjClRecoY",         &MAdjClRecoY);                         // Adj. clusters' reco Y [cm]
  fSolarNuAnaTree -> Branch("AdjClRecoZ",         &MAdjClRecoZ);                         // Adj. clusters' reco Z [cm]
  fSolarNuAnaTree -> Branch("AdjClR",             &MAdjClR);                             // Adj. clusters' distance [cm]
  fSolarNuAnaTree -> Branch("AdjClPur",           &MAdjClPur);                           // Adj. clusters' purity
  fSolarNuAnaTree -> Branch("AdjClMainID",        &MAdjClMainID);                        // Adj. clusters' main track ID
  fSolarNuAnaTree -> Branch("AdjClMainPDG",       &MAdjClMainPDG);                       // Adj. clusters' main PDG
  fSolarNuAnaTree -> Branch("AdjClMainE",         &MAdjClMainE);                         // Adj. clusters' main energy [GeV]
  fSolarNuAnaTree -> Branch("AdjClMainX",         &MAdjClMainX);                         // Adj. clusters' main X [cm]
  fSolarNuAnaTree -> Branch("AdjClMainY",         &MAdjClMainY);                         // Adj. clusters' main Y [cm]
  fSolarNuAnaTree -> Branch("AdjClMainZ",         &MAdjClMainZ);                         // Adj. clusters' main Z [cm]
  
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
  
  // Run = evt.run();SubRun = evt.subRun();
  Event = evt.event();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  
  Flag = rand() % 10000000000;
  std::string sHead = "";
  sHead = sHead + "\nTPC Frequency in [MHz]: " + str(clockData.TPCClock().Frequency());
  sHead = sHead + "\nTPC Tick in [us]: " + str(clockData.TPCClock().TickPeriod());
  sHead = sHead + "\nEvent Flag: " + str(Flag);
  sHead = sHead + "\nSuccesfull reset of variables for evt " + str(Event);
  sHead = sHead + "\n#########################################";
  PrintInColor(sHead,GetColor("magenta"));
  //---------------------------------------------------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------- Create maps for ID tracking -----------------------------------------------------------------// 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  // --- Fill MC Truth IDs to tracking vectors. Get a list of all of my particles in one chunk. ---  
  const sim::ParticleList& PartList = pi_serv->ParticleList();
  std::string sMcTruth = "";
  sMcTruth = sMcTruth + "\nThere are a total of " + str(int(PartList.size())) + " Particles in the event\n";
  
  // Loop over all signal+bkg handles and collect track IDs
  for ( size_t i = 0; i < fLabels.size(); i++){
    Parts.push_back(ThisGeneratorParts);                                    // For each label insert empty list
    
    art::Handle<std::vector<simb::MCTruth>> ThisHandle;
    evt.getByLabel(fLabels[i], ThisHandle);
    
    if(ThisHandle){
      auto ThisValidHanlde = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[i]); // Get generator handles
      art::FindManyP<simb::MCParticle> Assn(ThisValidHanlde,evt,fGEANTLabel);            // Assign labels to MCPArticles
      FillMyMaps( Parts[i], Assn, ThisValidHanlde);                                      // Fill empty list with previously assigned particles                                       
      if (Parts[i].size() < 1000){sMcTruth = sMcTruth + "\n# of particles " + str(int(Parts[i].size())) + "\tfrom " + fLabels[i];} // Print signal+bkg info to terminal
      else {sMcTruth = sMcTruth + "\n# of particles " + str(int(Parts[i].size())) + "\tfrom " + fLabels[i];}
      TPart.push_back(Parts[i].size()); // Insert #signal+bkg particles generated
      for ( std::map<int,simb::MCParticle>::iterator iter = Parts[i].begin(); iter != Parts[i].end(); iter++ ){
        std::set<int> ThisGeneratorIDs = {};
        trackids.push_back(ThisGeneratorIDs);
        trackids[i].insert( iter->first );// Contains a list of TrIDs
      }
    } 
    else{
      sMcTruth = sMcTruth + "\n# of particles " + str(int(Parts[i].size())) + "\tfrom " + fLabels[i] + " *not generated!";
      TPart.push_back(0);
      std::set<int> ThisGeneratorIDs = {};
      trackids.push_back(ThisGeneratorIDs);
    }
  }
  PrintInColor(sMcTruth,GetColor("yellow"));

  //----------------------------------------------------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------- Some MC Truth information -------------------------------------------------------------------// 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  std::set< int > signal_trackids;                                // Signal TrackIDs to be used in OpFlash matching
  std::vector<std::vector<int>> ClPartTrackIDs = {{},{},{},{}};   // Track IDs corresponding to each kind of MCTruth particle  {11,22,2112,else}
  art::Handle<std::vector<simb::MCTruth>> ThisHandle;
  std::string sNuTruth = "";
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
      
      sNuTruth = sNuTruth + "\nNumber of Neutrino Daughters: " + str(N-2);
      sNuTruth = sNuTruth + "\nNeutrino energy: " + str(TNuE) + " GeV"; 
      sNuTruth = sNuTruth + "\nMomentumTransfer: " + str(std::sqrt(TNuQSqr)) + " GeV";
      sNuTruth = sNuTruth + "\nPosition (" + str(TNuX) + ", " + str(TNuY) + ", " + str(TNuZ) + ") cm";
      
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
    sNuTruth = sNuTruth + "\nGen.\t PdgCode\t Energy\t\t TrackID \n------------------------------------------------";
    
    for ( size_t i = 0; i < MarlAssn.size(); i++) {
      auto parts = MarlAssn.at(i);
      for (auto part = parts.begin(); part != parts.end(); part++) {
        signal_trackids.emplace((*part)->TrackId());

        if ((*part)->PdgCode()<1000000){sNuTruth = sNuTruth + "\n" + fLabels[0] + "\t " + str((*part)->PdgCode()) + "\t\t " + str((*part)->E()) + "\t " + str((*part)->TrackId());}
        else{sNuTruth = sNuTruth + "\n" + fLabels[0] + "\t " + str((*part)->PdgCode()) + "\t " + str((*part)->E()) + "\t " + str((*part)->TrackId());}

        if ((*part)->PdgCode()==11){ // Electrons
          const TLorentzVector &v4_f = (*part)->EndPosition();
          auto x_f = v4_f.X();auto y_f = v4_f.Y();auto z_f = v4_f.Z();
          avX = x_f; avY = y_f; avZ = z_f;
          ClPartTrackIDs[0].push_back((*part)->TrackId());
          mf::LogDebug("SolarNuAna") << "\nMC Electron truth position x = " << avX << ", y = " << avY << ", z = " << avZ;
          mf::LogDebug("SolarNuAna") << "Initial KE " << (*part)->E()-(*part)->Mass();
        }
        if ((*part)->PdgCode()==22){ClPartTrackIDs[1].push_back((*part)->TrackId());} // Gammas
        if ((*part)->PdgCode()==2112){ClPartTrackIDs[2].push_back((*part)->TrackId());} // Neutrons
        if ((*part)->PdgCode()!=11 && (*part)->PdgCode()!=22 && (*part)->PdgCode()!=2112){ClPartTrackIDs[3].push_back((*part)->TrackId());} // Others
      }
    }
  }
  else{mf::LogWarning("SolarNuAna") << "No MARLEY MCTruths found.";}
  PrintInColor(sNuTruth,GetColor("blue"));
  fMCTruthTree->Fill();
  
  //----------------------------------------------------------------------------------------------------------------------------------------//
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

  // Loop over flashlist and assign OpHits to each flash
  for ( int i = 0; i < int(flashlist.size()); i++ ){

    recob::OpFlash TheFlash = *flashlist[i];
    // if (i%10 == 0) PrintInColor("Flash Time = " + str(TheFlash.Time()), GetColor("red"));
    std::vector< art::Ptr< recob::OpHit > > matchedHits = OpAssns.at(i);
    mf::LogDebug("SolarNuAna") << "Assigning OpHit to Flash";
    // Calculate the total PE of the flash and the time of the ophit with the highest PE 
    double totPE = 0; double MaxHitPE = 0;
    float OpHitT, OpHitPE;
    for (int j = 0; j < int(matchedHits.size()); j++){
      recob::OpHit ohit = *matchedHits[j];
      totPE += ohit.PE();
      OpHitPE = ohit.PE();
      if (OpHitPE > MaxHitPE){
        MaxHitPE = OpHitPE;
        OpHitT = ohit.PeakTime();
      }
    }

    mf::LogDebug("SolarNuAna") << "Evaluating Flash purity";
    // Calculate the flash purity
    int TerminalOutput = supress_stdout();
    double OpFlashPur = pbt->OpHitCollectionPurity(signal_trackids, matchedHits);
    resume_stdout(TerminalOutput);
    mf::LogDebug("SolarNuAna") << "PE of this OpFlash " << totPE << " OpFlash time " << OpHitT;

    // Calculate the flash purity, only for the Marley events
    if (MaxHitPE/totPE < fAdjOpFlashMaxPECut && totPE > fAdjOpFlashMinPECut){
      OpFlashMarlPur.push_back(OpFlashPur);
      OpFlashPE.push_back(TheFlash.TotalPE());
      OpFlashMaxPE.push_back(MaxHitPE);
      OpFlashY.push_back(TheFlash.YCenter());
      OpFlashZ.push_back(TheFlash.ZCenter());
      OpFlashT.push_back(TheFlash.Time());
      OpFlashDeltaT.push_back(TheFlash.TimeWidth());
      OpFlashNHit.push_back(matchedHits.size());
    }
    if (abs(OpHitT) < 30){
      mf::LogDebug("SolarNuAna") << "OpFlash PE (max/tot) " << MaxHitPE << "/" << TheFlash.TotalPE() << " with purity " << OpFlashPur << " time " << TheFlash.Time();
    }
  }

  //----------------------------------------------------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------- Hit collection and assignment ----------------------------------------------------------------// 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  // --- Lift out the reco hits:
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(fHitLabel);
  int NTotHits = reco_hits->size();

  for(int hit = 0; hit < NTotHits; ++hit){
    // --- Loop over the reconstructed hits to separate them among tpc planes according to view
    
    recob::Hit const& ThisHit = reco_hits->at(hit);
    if (ThisHit.PeakTime() < 0) PrintInColor("Negative Hit Time = " + str(ThisHit.PeakTime()), GetColor("red"));
    mf::LogDebug("SolarNuAna") << "Hit " << hit << " has view " << ThisHit.View() << " and signal type " << ThisHit.SignalType();

    if      (ThisHit.SignalType() == 0 && ThisHit.View() == 0){ColHits0.push_back( ThisHit );} // SignalType = 0
    else if (ThisHit.SignalType() == 0 && ThisHit.View() == 1){ColHits1.push_back( ThisHit );} // SignalType = 0
    else if (ThisHit.SignalType() == 1)                       {ColHits2.push_back( ThisHit );} // SignalType = 1
    else    {ColHits3.push_back( ThisHit ); mf::LogError("SolarNuAna") << "Hit was found with view out of scope";}
  } 

  //----------------------------------------------------------------------------------------------------------------------------------------//
  //-------------------------------------------------------------- Cluster creation and analysis ------------------------------------------------------------------// 
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
  // --- Now calculate the clusters ...
  CalcAdjHits(ColHits0,Clusters0,hAdjHits,hAdjHitsADCInt,false);
  CalcAdjHits(ColHits1,Clusters1,hAdjHits,hAdjHitsADCInt,false);
  CalcAdjHits(ColHits2,Clusters2,hAdjHits,hAdjHitsADCInt,false);
  CalcAdjHits(ColHits3,Clusters3,hAdjHits,hAdjHitsADCInt,false);

  std::vector< std::vector< std::vector<float>>> ClGenPur = {{},{},{}};
  std::vector< std::vector< std::vector<recob::Hit>>> AllPlaneClusters = {Clusters0,Clusters1,Clusters2};
  std::vector< std::vector< float>> ClCharge = {{},{},{}}, ClMaxCharge = {{},{},{}}, CldT = {{},{},{}}, ClT = {{},{},{}}, ClX = {{},{},{}}, ClY = {{},{},{}}, ClZ = {{},{},{}};
  std::vector< std::vector< float>> ClFracE = {{},{},{}}, ClFracGa  = {{},{},{}}, ClFracNe = {{},{},{}}, ClFracRest = {{},{},{}}, ClPur = {{},{},{}}, Cldzdy = {{},{},{}};
  std::vector< std::vector< int  >> ClMainID = {{},{},{}}, ClTPC = {{},{},{}}, ClNHits = {{},{},{}}, ClGen = {{},{},{}};

  std::string sRecoObjects = "";
  sRecoObjects = sRecoObjects + "\nTotal number of flashes constructed: " + str(int(flashlist.size()));
  sRecoObjects = sRecoObjects + "\n# Hits in each view = " + str(int(ColHits0.size())) + ", " + str(int(ColHits1.size())) + ", " + str(int(ColHits2.size())) + ", " + str(int(ColHits3.size()));
  sRecoObjects = sRecoObjects + "\n# Clusters from the hits = " + str(int(Clusters0.size())) + ", " + str(int(Clusters1.size())) + ", " + str(int(Clusters2.size())) + ", " + str(int(Clusters3.size()));
  PrintInColor(sRecoObjects,GetColor("cyan"));
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
        if (hit.PeakTime() < 0) PrintInColor("Negative Cluster Time = " + str(hit.PeakTime()), GetColor("red"));
        ncharge += hit.Integral();
        const geo::WireGeo* wire = geo->GeometryCore::WirePtr(hit.WireID()); // Wire directions should be the same for all hits of the same view (can be used to check)
        double hitCharge;
          
        geo::Point_t hXYZ = wire->GetCenter();
        geo::Point_t sXYZ = wire->GetStart();
        geo::Point_t eXYZ = wire->GetEnd();

        geo::Vector_t direction = eXYZ - sXYZ;
        auto dyds = direction.Y(), dzds = direction.Z();
        thisdzdy.push_back(dzds/dyds);

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
            mf::LogDebug("SolarNuAna") << "This hit's IDE is: " << MainTrID; 
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
        mf::LogDebug("SolarNuAna") << "\nThis particle type " << ThisPType << "\nThis cluster's main track ID " << MainTrID;    
        if (ThisPType == 1){hitCharge = hit.Integral();Pur = Pur+hitCharge;}
      }

      float MaxGen = 0;
      // PrintInColor("GenVector: "+ str(GenPur), GetColor("red"));
      for (size_t genpur = 0; genpur < GenPur.size(); genpur++){
        if (GenPur[genpur] > MaxGen){Gen = genpur; MaxGen = GenPur[genpur];}
        GenPur[genpur] = GenPur[genpur]/ncharge;
      }
      // PrintInColor("Gen: "+ str(Gen), GetColor("red"));

      for (size_t j = 0; j > thisdzdy.size(); j++) {if (thisdzdy[0] != thisdzdy[i]) mf::LogWarning("SolarNuAna") << "MISSMATCH IN dzdy FOR CLUSTER " << idx;}

      dzdy = thisdzdy[0]; thisdzdy.clear();      
      FracE /= ncharge; FracGa /= ncharge; FracNe /= ncharge; FracRest /= ncharge;
      clustTPC /= ncharge; clustX /= ncharge; clustY /= ncharge;clustZ /= ncharge;clustT /= ncharge;
      mf::LogDebug("SolarNuAna") << "\ndzdy " << dzdy << " for cluster " << " (" << clustY << ", " << clustZ << ") with track ID " << MainTrID <<  " in plane " << idx; 
      if (clustT < 0) PrintInColor("Negative Cluster Time = " + str(clustT), GetColor("red"));

      ClCharge[idx].push_back(ncharge);
      ClMaxCharge[idx].push_back(maxHit);
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

      mf::LogDebug("SolarNuAna") << "\nCluster " << i << " in plane " << idx << " has #hits" << nhit << " charge, " << ncharge << " time, " << clustT;
      mf::LogDebug("SolarNuAna") << " and position (" << clustY << ", " << clustZ << ") with main track ID " << MainTrID << " and purity " << Pur/ncharge;
    }
  } // Finished first cluster processing
  
  //-------------------------------------------------------------------- Cluster Matching -------------------------------------------------------------------------// 
  std::vector<std::vector<float>> MVecGenFrac = {};
  std::vector<int>   MVecNHit = {}, MVecGen = {}, MVecInd0NHits  = {}, MVecInd1NHits  = {}, MVecMainID = {}, MVecTPC = {}, MVecInd0TPC = {}, MVecInd1TPC = {};
  std::vector<float> MVecTime  = {}, MVecCharge = {}, MVecMaxCharge = {}, MVecInd0Charge = {}, MVecInd1Charge = {}, MVecInd0MaxCharge = {}, MVecInd1MaxCharge = {}, MVecInd0dT = {}, MVecInd1dT = {};
  std::vector<float> MVecInd0RecoY = {}, MVecInd1RecoY = {}, MVecRecY = {}, MVecRecZ = {};
  std::vector<float> MVecFracE = {}, MVecFracGa = {}, MVecFracNe = {}, MVecFracRest = {}, MVecPur = {};

  for (int ii = 0; ii < int(AllPlaneClusters[2].size()); ii++){
    bool match = false;
    int ind0clustNHits = 0, ind1clustNHits = 0;
    int ind0clustTPC = 0, ind1clustTPC = 0;
    double ind0clustY = -1e6, ind1clustY = -1e6, ind0clustMaxCharge = 0, ind1clustMaxCharge = 0, ind0clustCharge = 0, ind1clustCharge = 0;
    double ind0clustdT = fClusterMatchTime, ind1clustdT = fClusterMatchTime;
    if (!AllPlaneClusters[2][ii].empty()){
      if (!AllPlaneClusters[0].empty()){
        for (int jj = 0; jj < int(AllPlaneClusters[0].size()); jj++){
          if ( ClNHits[0][jj] < (1-fClusterMatchNHit)*ClNHits[2][ii] || ClNHits[0][jj] > (1+fClusterMatchNHit)*ClNHits[2][ii]) {continue;}
          if ( ClCharge[0][jj] < (1-fClusterMatchCharge)*ClCharge[2][ii] || ClCharge[0][jj] > (1+fClusterMatchCharge)*ClCharge[2][ii]) {continue;}
          if ( abs(ClT[2][ii] - ClT[0][jj]) < fClusterMatchTime && abs(fClusterInd0MatchTime - abs(ClT[2][ii] - ClT[0][jj])) < abs(fClusterInd0MatchTime - ind0clustdT)) {
            ind0clustY = ClY[0][jj] + (ClZ[2][ii] - ClZ[0][jj])/(Cldzdy[0][jj]);
            ind0clustdT = abs(ClT[2][ii] - ClT[0][jj]);
            ind0clustNHits = int(AllPlaneClusters[0][jj].size());
            ind0clustCharge = ClCharge[0][jj];
            ind0clustMaxCharge = ClMaxCharge[0][jj];
            ind0clustTPC = ClTPC[0][jj];
            if (ind0clustY > -fDetectorSizeY && ind0clustY < fDetectorSizeY){match = true;}
            mf::LogDebug("SolarNuAna") << "¡¡¡ Matched cluster in plane 0 !!! --- Position x = " << ClX[0][jj] << ", y = " << ClY[0][jj] << ", z = " << ClZ[0][jj];
            mf::LogDebug("SolarNuAna") << "Reconstructed position y = " << ind0clustY << ", z = " << ClZ[2][ii];  
          }
        }
      }
      if (!AllPlaneClusters[1].empty()){
        for (int zz = 0; zz < int(AllPlaneClusters[1].size()); zz++){
          if ( ClNHits[1][zz] < (1-fClusterMatchNHit)*ClNHits[2][ii] || ClNHits[1][zz] > (1+fClusterMatchNHit)*ClNHits[2][ii]){continue;}
          if ( ClCharge[1][zz] < (1-fClusterMatchCharge)*ClCharge[2][ii] || ClCharge[1][zz] > (1+fClusterMatchCharge)*ClCharge[2][ii]){continue;}
          if ( abs(ClT[2][ii] - ClT[1][zz]) < fClusterMatchTime && abs(fClusterInd1MatchTime - abs(ClT[2][ii] - ClT[1][zz])) < abs(fClusterInd1MatchTime - ind1clustdT)) {
            ind1clustY = ClY[1][zz] + (ClZ[2][ii] - ClZ[1][zz])/(Cldzdy[1][zz]);
            ind1clustdT = abs(ClT[2][ii] - ClT[1][zz]);
            ind1clustNHits = int(AllPlaneClusters[1][zz].size());
            ind1clustCharge = ClCharge[1][zz];
            ind1clustMaxCharge = ClMaxCharge[1][zz];
            ind1clustTPC = ClTPC[1][zz];
            if (ind1clustY > -fDetectorSizeY && ind1clustY < fDetectorSizeY){match = true;}
            mf::LogDebug("SolarNuAna") << "¡¡¡ Matched cluster in plane 1 !!! --- Position x = " << ClX[1][zz] << ", y = " << ClY[1][zz] << ", z = " << ClZ[1][zz];
            mf::LogDebug("SolarNuAna") << "Reconstructed position y = " << ind1clustY << ", z = " << ClZ[2][ii];
          }
        } // Loop over ind1 clusters
      }
    } // Loop over ind clusters
    else {mf::LogDebug("SolarNuAna") << "Cluster " << ii << " in plane 2 has no hits";}
    
    //--------------------------------------------------------- Export Matched cluster vectors ------------------------------------------------------------------// 
    if (match == true){
      // Cluster Charge
      MVecCharge.push_back(ClCharge[2][ii]);
      MVecMaxCharge.push_back(ClMaxCharge[2][ii]);
      MVecInd0Charge.push_back(ind0clustCharge);
      MVecInd1Charge.push_back(ind1clustCharge);
      MVecInd0MaxCharge.push_back(ind0clustMaxCharge);
      MVecInd1MaxCharge.push_back(ind1clustMaxCharge);
      // Cluster Hits
      MVecNHit.push_back(ClNHits[2][ii]);
      MVecInd0NHits.push_back(ind0clustNHits);
      MVecInd1NHits.push_back(ind1clustNHits);
      // Cluster TPC
      MVecTPC.push_back(ClTPC[2][ii]);
      MVecInd0TPC.push_back(ind0clustTPC);
      MVecInd1TPC.push_back(ind1clustTPC);
      // Cluster Time
      MVecTime.push_back(ClT[2][ii]);
      MVecInd0dT.push_back(ind0clustdT);
      MVecInd1dT.push_back(ind1clustdT);
      // Cluster RecoY
      MVecInd0RecoY.push_back(ind0clustY);
      MVecInd1RecoY.push_back(ind1clustY);
      // Cluster RecoZ	    
      MVecRecZ.push_back(ClZ[2][ii]);
      // Cluster Marley Fractions
      MVecFracE.push_back(ClFracE[2][ii]);
      MVecFracGa.push_back(ClFracGa[2][ii]);
      MVecFracNe.push_back(ClFracNe[2][ii]);
      MVecFracRest.push_back(ClFracRest[2][ii]);
      // Cluster Marley Purity
      MVecPur.push_back(ClPur[2][ii]);
      // Cluster Gen and GenFraction
      MVecMainID.push_back(ClMainID[2][ii]);
      MVecGen.push_back(ClGen[2][ii]);
      MVecGenFrac.push_back(ClGenPur[2][ii]);
      
      float buffer = 1;
      if ((ind0clustY > -buffer*fDetectorSizeY && ind0clustY < buffer*fDetectorSizeY) && (ind1clustY > -buffer*fDetectorSizeY && ind1clustY < buffer*fDetectorSizeY)){
        mf::LogDebug("SolarNuAna") << "BOTH IND RECO INSIDE OF DETECTOR";
        MVecRecY.push_back((ind0clustY+ind1clustY)/2);}
      else if (ind0clustY > -buffer*fDetectorSizeY && ind0clustY < buffer*fDetectorSizeY){
        mf::LogDebug("SolarNuAna") << "IND1 OUTSIDE OF DETECTOR";
        MVecRecY.push_back(ind0clustY);}
      else if (ind1clustY > -buffer*fDetectorSizeY && ind1clustY < buffer*fDetectorSizeY){
        mf::LogDebug("SolarNuAna") << "IND0 OUTSIDE OF DETECTOR";
        MVecRecY.push_back(ind1clustY);}
      else{
        mf::LogDebug("SolarNuAna") << "RECO OUTSIDE OF DETECTOR";
        MVecRecY.push_back((ind0clustY+ind1clustY)/2); 
        if (ClGen[2][ii] == 1){mf::LogWarning("SolarNuAna") << "Marley cluster reconstructed outside of detector volume! RecoY = " << str((ind0clustY+ind1clustY)/2);}
      }

      // Print in color if the cluster is matched
      mf::LogDebug("SolarNuAna") << "¡¡¡ Matched cluster !!! ";
      mf::LogDebug("SolarNuAna") << " - Cluster " << str(ClMainID[2][ii]) << " Gen " << str(ClGen[2][ii]) << " Purity " << str(ClPur[2][ii]) << " Hits " << str(ClNHits[2][ii]);
      mf::LogDebug("SolarNuAna") << " - Hits(ind0, ind1, col) " << str(ind0clustNHits) << ", " << str(ind1clustNHits) << ", " << str(ClNHits[2][ii]);
      mf::LogDebug("SolarNuAna") << " - Positions y(ind0, ind1) = " << str(ind0clustY) << ", " << str(ind1clustY) << ", z = " << str(ClZ[2][ii]) << "\n";
    } // if (match == true)
  } // Loop over collection plane clusters
  
  //-------------------------------------------------------------------- Cluster Tree Export -------------------------------------------------------------------------// 
  // Loop over matched clusters and export to tree if number of hits is above threshold
  for (int i = 0; i < int(MVecNHit.size()); i++){
    if(MVecNHit[i] > fClusterPreselectionNHit && (MVecInd0NHits[i] > fClusterPreselectionNHit || MVecInd1NHits[i] > fClusterPreselectionNHit)){
      MAdjClTime = {};MAdjClCharge = {};MAdjClInd0Charge = {};MAdjClInd1Charge = {};MAdjClMaxCharge = {};MAdjClInd0MaxCharge = {};MAdjClInd1MaxCharge = {};MAdjClNHit = {};;MAdjClInd0NHit = {};;MAdjClInd1NHit = {};
      MAdjClRecoY = {};MAdjClRecoZ = {};MAdjClR = {};MAdjClPur = {};MAdjClGen = {};MAdjClMainID = {};MAdjClMainPDG = {};
      MAdjClMainE = {}; MAdjClMainX = {};MAdjClMainY = {};MAdjClMainZ = {};
      MAdjFlashTime = {};MAdjFlashPE = {};MAdjFlashNHit = {};MAdjFlashMaxPE = {};MAdjFlashRecoY = {};MAdjFlashRecoZ = {};MAdjFlashR = {};MAdjFlashPur = {};
      
      std::string ResultColor = "white";
      if (abs(MVecRecY[i] - TNuY) < 10 && abs(MVecRecZ[i] - TNuZ) < 10) {ResultColor = "green";}
      else {ResultColor = "yellow";}

      PrintInColor("*** Matched preselection cluster: \n - MainCluster  " + str(MVecMainID[i]) + " Gen " + str(MVecGen[i]) + " Purity " + str(MVecPur[i]) + " Hits " + str(MVecNHit[i])
        + "\n - RecoY/RecoZ (" + str(MVecRecY[i]) + " / " + str(MVecRecZ[i]) + ") Time " + str(MVecTime[i]) + "\n", GetColor(ResultColor));

      // Loop over collection plane clusters to find adjacent clusters with distance < fAdjClusterRad and time < fAdjClusterTime
      for ( int j = 0; j < int(MVecNHit.size()); j++ ){
        if ( j == i ) {continue;}
        if ( fGeometry == "HD" && sqrt(pow(MVecRecY[i]-MVecRecY[j],2) + pow(MVecRecZ[i]-MVecRecZ[j],2) + pow((MVecTime[i]-MVecTime[j])*fDetectorSizeX/fAdjOpFlashTime,2)) > fAdjClusterRad ){continue;}
        if ( fGeometry == "VD" && sqrt(pow(MVecRecY[i]-MVecRecY[j],2) + pow(MVecRecZ[i]-MVecRecZ[j],2) + pow((MVecTime[i]-MVecTime[j])*fDetectorSizeX/(fAdjOpFlashTime/2),2)) > fAdjClusterRad ){continue;}
        MAdjClTime.push_back(MVecTime[j]);
        MAdjClCharge.push_back(MVecCharge[j]);
        MAdjClInd0Charge.push_back(MVecInd0Charge[j]);
        MAdjClInd1Charge.push_back(MVecInd1Charge[j]);
        MAdjClMaxCharge.push_back(MVecMaxCharge[j]);
        MAdjClInd0MaxCharge.push_back(MVecInd0MaxCharge[j]);
        MAdjClInd1MaxCharge.push_back(MVecInd1MaxCharge[j]);
        MAdjClNHit.push_back(MVecNHit[j]);
        MAdjClInd0NHit.push_back(MVecInd0NHits[j]);
        MAdjClInd1NHit.push_back(MVecInd1NHits[j]);
        MAdjClRecoY.push_back(MVecRecY[j]);
        MAdjClRecoZ.push_back(MVecRecZ[j]);
        MAdjClR.push_back(sqrt(pow(MVecRecY[i]-MVecRecY[j],2)+pow(MVecRecZ[i]-MVecRecZ[j],2)));
        MAdjClPur.push_back(MVecPur[j]);
        MAdjClGen.push_back(MVecGen[j]);
        MAdjClMainID.push_back(MVecMainID[j]);

        // If mother exists add the mother information
        const simb::MCParticle *MAdjClTruth;
        int TerminalOutput = supress_stdout();
        MAdjClTruth = pi_serv->TrackIdToParticle_P(MVecMainID[j]);
        resume_stdout(TerminalOutput);
        if (MAdjClTruth == 0) {
          MAdjClMainPDG.push_back(0);          
          MAdjClMainE.push_back(-1e6);
          MAdjClMainX.push_back(-1e6);
          MAdjClMainY.push_back(-1e6);
          MAdjClMainZ.push_back(-1e6);
        }
        else{
          MAdjClMainPDG.push_back(MAdjClTruth->PdgCode());          
          MAdjClMainE.push_back(MAdjClTruth->E());
          MAdjClMainX.push_back(MAdjClTruth->Vx());
          MAdjClMainY.push_back(MAdjClTruth->Vy());
          MAdjClMainZ.push_back(MAdjClTruth->Vz());
        }
      }

      for (int j = 0; j < int(OpFlashPE.size()); j++){
        if ( (MVecTime[i] - OpFlashT[j]) < 0) {continue;}
        if ( (MVecTime[i] - OpFlashT[j]) > fAdjOpFlashTime ) {continue;}
        if ( sqrt(pow(MVecRecY[i]-OpFlashY[j],2)+pow(MVecRecZ[i]-OpFlashZ[j],2)) > fAdjOpFlashRad) {continue;}  
        MAdjFlashTime.push_back(OpFlashT[j]);
        MAdjFlashPE.push_back(OpFlashPE[j]);
        MAdjFlashNHit.push_back(OpFlashNHit[j]);
        MAdjFlashMaxPE.push_back(OpFlashMaxPE[j]);
        MAdjFlashRecoY.push_back(OpFlashY[j]);
        MAdjFlashRecoZ.push_back(OpFlashZ[j]);
        MAdjFlashR.push_back(sqrt(pow(MVecRecY[i]-OpFlashY[j],2)+pow(MVecRecZ[i]-OpFlashZ[j],2)));
        MAdjFlashPur.push_back(OpFlashMarlPur[j]);
      }

      // Fill the tree with the cluster information and the adjacent clusters and flashes
      MMarleyFrac =    {MVecFracE[i],MVecFracGa[i],MVecFracNe[i],MVecFracRest[i]};
      MGenFrac =       MVecGenFrac[i];
      MTime =          MVecTime[i];   
      MCharge =        MVecCharge[i];   
      MMaxCharge =     MVecMaxCharge[i];   
      MNHit =          MVecNHit[i];
      // Cluster TPC
      MTPC =           MVecTPC[i];
      MInd0TPC =       MVecInd0TPC[i];
      MInd1TPC =       MVecInd1TPC[i];
      // Cluster MaxChargeHit
      MInd0Charge =    MVecInd0Charge[i];   
      MInd1Charge =    MVecInd1Charge[i];
      MInd0MaxCharge = MVecInd0MaxCharge[i];   
      MInd1MaxCharge = MVecInd1MaxCharge[i];
      MInd0NHits =     MVecInd0NHits[i];   
      MInd1NHits =     MVecInd1NHits[i];
      MInd0dT =        MVecInd0dT[i];   
      MInd1dT =        MVecInd1dT[i];   
      MInd0RecoY =     MVecInd0RecoY[i];   
      MInd1RecoY =     MVecInd1RecoY[i];   
      MMainID =        MVecMainID[i];
      MRecZ =          MVecRecZ[i];   
      MPur =           MVecPur[i];
      MGen =           MVecGen[i];

      // If mother exists add the mother information
      const simb::MCParticle *MClTruth;
      int TerminalOutput = supress_stdout();
      MClTruth = pi_serv->TrackIdToParticle_P(MVecMainID[i]);
      resume_stdout(TerminalOutput);
      if (MClTruth == 0){
        MMainVertex = {-1e6,-1e6,-1e6};
        MEndVertex = {-1e6,-1e6,-1e6};
        MMainPDG =       0;
        MMainE =      -1e6;
        MMainT =      -1e6;
        MMainP =      -1e6;
        
        MMainParentVertex = {-1e6,-1e6,-1e6};
        MMainParentPDG =       0;
        MMainParentE =      -1e6;
        MMainParentT =      -1e6;
        MMainParentP =      -1e6;
      }
      else{
        MMainVertex = {MClTruth->Vx(),MClTruth->Vy(),MClTruth->Vz()};
        MEndVertex =  {MClTruth->EndX(),MClTruth->EndY(),MClTruth->EndZ()};
        MMainPDG =    MClTruth->PdgCode();
        MMainE =      MClTruth->E();
        MMainT =      MClTruth->T();
        MMainP =      MClTruth->P();
        // If exists add the parent information
        const simb::MCParticle *MClParentTruth;
        int TerminalOutput = supress_stdout();
        MClParentTruth = pi_serv->TrackIdToParticle_P(MClTruth->Mother());
        resume_stdout(TerminalOutput);
        if (MClParentTruth == 0){
          MMainParentVertex = {-1e6,-1e6,-1e6};
          MMainParentPDG =       0;
          MMainParentE =      -1e6;
          MMainParentT =      -1e6;
          MMainParentP =      -1e6;
        }
        else{
          MMainParentVertex = {MClParentTruth->Vx(),MClParentTruth->Vy(),MClParentTruth->Vz()};
          MMainParentPDG =    MClParentTruth->PdgCode();
          MMainParentE =      MClParentTruth->E();
          MMainParentT =      MClParentTruth->T();
          MMainParentP =      MClParentTruth->P();
        }
      }

      hDriftTime->      Fill(avX,MTime);
      fSolarNuAnaTree-> Fill();
      hXTruth->         Fill(MVecRecY[i]-TNuY,TNuX); 
      hYTruth->         Fill(MVecRecY[i]-TNuY,TNuY); 
      hZTruth->         Fill(MVecRecZ[i]-TNuZ,TNuZ); 
      if (MVecTime[i]<0) mf::LogWarning("SolarNuAna") << "Negative Main Cluster Time = " << MVecTime[i];
    }
  }
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
      mf::LogDebug("SolarNuAna") << ThisPar.PdgCode() << " " << ThisPar.E();
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
void SolarNuAna::PrintInColor( std::string MyString, int Color, std::string Type ){
  if (Type == "Info"){mf::LogInfo("SolarNuAna") << "\033[" << Color << "m" << MyString << "\033[0m";}
  if (Type == "Degub"){mf::LogDebug("SolarNuAna") << "\033[" << Color << "m" << MyString << "\033[0m";}
  if (Type == "Error"){mf::LogError("SolarNuAna") << "\033[" << Color << "m" << MyString << "\033[0m";}
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
  else if (ColorName == "bright_black") return 90;
  else if (ColorName == "bright_red") return 91;
  else if (ColorName == "bright_green") return 92;
  else if (ColorName == "bright_yellow") return 93;
  else if (ColorName == "bright_blue") return 94;
  else if (ColorName == "bright_magenta") return 95;
  else if (ColorName == "bright_cyan") return 96;
  else if (ColorName == "bright_white") return 97;
  else {mf::LogError("SolarNuAna") << "Color " << ColorName << " not recognized. Returning white."; return 37;}
  return 0;
}

std::string SolarNuAna::str( int i ) {std::stringstream ss;ss << i;return ss.str();}
std::string SolarNuAna::str( double i ) {std::stringstream ss;ss << i;return ss.str();}
std::string SolarNuAna::str( float i ) {std::stringstream ss;ss << i;return ss.str();}
std::string SolarNuAna::str( std::vector<int> i ) {std::stringstream ss;for (int j = 0; j < int(i.size()); j++){ss << i[j] << " ";}return ss.str();}
std::string SolarNuAna::str( std::vector<double> i ) {std::stringstream ss;for (int j = 0; j < int(i.size()); j++){ss << i[j] << " ";}return ss.str();}
std::string SolarNuAna::str( std::vector<float> i ) {std::stringstream ss;for (int j = 0; j < int(i.size()); j++){ss << i[j] << " ";}return ss.str();}

int SolarNuAna::supress_stdout() {
  std::fflush(stdout);

  int ret = dup(1);
  int nullfd = open("/dev/null", O_WRONLY);
  // check nullfd for error omitted
  dup2(nullfd, 1);
  close(nullfd);

  return ret;
}

void SolarNuAna::resume_stdout(int fd) {
  std::fflush(stdout);
  dup2(fd, 1);
  close(fd);
}

DEFINE_ART_MODULE(SolarNuAna)