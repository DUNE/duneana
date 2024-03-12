////////////////////////////////////////////////////////////////////////////////////
// Class:       LowEAna                                                           //
// Module Type: analyzer                                                          //
// File:        LowEAna_module.cc                                                 //
//                                                                                //
// Written by Sergio Manthey Corchado with guidence of Daniel Pershey             //
// developed from Michael Baird's DAQSimAna_module                                //
////////////////////////////////////////////////////////////////////////////////////

// C++ includes
#ifndef LowEAna_h
#define LowEAna_h 

// C++ includes
// ROOT includes
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TVector3.h>
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

class LowEAna : public art::EDAnalyzer {

public:
  // --- Standard constructor and destructor for an ART module.
  explicit LowEAna(fhicl::ParameterSet const & p);
  LowEAna(LowEAna const &)               = delete;
  LowEAna(LowEAna &&)                    = delete;
  LowEAna & operator = (LowEAna const &) = delete;
  LowEAna & operator = (LowEAna &&)      = delete;
  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;

private:
  // --- Define a struct to hold the OpFlash info
  struct OpFlashes{
    std::vector<double> OpFlashPE;
    std::vector<double> OpFlashMaxPE;
    std::vector<double> OpFlashX;
    std::vector<double> OpFlashY;
    std::vector<double> OpFlashZ;
    std::vector<double> OpFlashT;
    std::vector<int> OpFlashNHit;
    std::vector<int> OpFlashGen;
    std::vector<double> OpFlashPur;
  };
  // --- Some of our own functions.
  void ResetEventVariables(detinfo::DetectorClocksData clockData);
  void ResetClusterVariables();
  
  void AssignOpFlashes(
    struct OpFlashes &EventOpFlashes,
    float &ClY,
    float &ClZ
  );

  void FillMCInteractionTree(
    std::map< int, simb::MCParticle> &MCParticleList,
    std::vector<std::string> ProcessList,
    bool HeavDebug );

  void FillClusterHitVectors(std::vector<recob::Hit> Cluster, 
    std::vector<int> &TPC,
    std::vector<int> &Channel,
    std::vector<double> &MotherX,
    std::vector<double> &MotherY,
    std::vector<double> &MotherZ,
    std::vector<double> &MotherE,
    std::vector<double> &MotherP,
    std::vector<int> &MotherPDG,
    std::vector<double> &AncestorX,
    std::vector<double> &AncestorY,
    std::vector<double> &AncestorZ,
    std::vector<double> &AncestorE,
    std::vector<double> &AncestorP,
    std::vector<int> &AncestorPDG, 
    std::vector<double> &Charge,
    std::vector<double> &Time, 
    std::vector<double> &Y,
    std::vector<double> &Z,
    std::vector<double> &Dir,
    detinfo::DetectorClocksData clockData,
    bool HeavDebug);

  void FillClusterVariables(
    std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
    std::vector<std::vector<int>>& ClNHits,
    std::vector<std::vector<float>>& ClT,
    std::vector<std::vector<float>>& ClCharge,
    bool HeavDebug);
  
  std::vector<std::vector<std::vector<recob::Hit>>> MatchClusters(
    std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
    std::vector<std::vector<int>>& ClNHits,
    std::vector<std::vector<float>>& ClT,
    std::vector<std::vector<float>>& ClCharge,
    bool HeavDebug);

  void FillMyMaps    ( std::map< int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand );
  void CalcAdjHits   ( std::vector< recob::Hit > MyVec,std::vector< std::vector<recob::Hit> >& Clusters,TH1I* MyHist, TH1F* MyADCIntHist, 
    detinfo::DetectorClocksData clockData,
    bool HeavDebug );

  std::string PrintInColor ( std::string InputString, std::string MyString, int MyColor );
  
  bool InMyMap(
    int TrID,
    std::map< int,
    simb::MCParticle> ParMap);

  long unsigned int WhichParType  ( int TrID );
  int  GetColor      ( std::string MyString );

  std::vector<double> ComputeInterpolationRecoY(
    int Event,
    std::vector<int> &IndTPC,
    std::vector<double> &Z, 
    std::vector<double> &T,
    std::vector<double> &IndZ,
    std::vector<double> &IndY,
    std::vector<double> &IndT,
    std::vector<double> &IndDir,
    bool HeavDebug);

  std::string str( int MyInt );
  std::string str( float MyFloat );
  std::string str( double MyDouble );
  std::string str( std::vector<int> MyVec );
  std::string str( std::vector<float> MyVec );
  std::string str( std::vector<double> MyVec );

  int supress_stdout();
  void resume_stdout(int fd);

  double Average(std::vector<double> &v);
  double Average(std::vector<float> &v);
  double Average(std::vector<int> &v);

  // --- Our fcl parameter labels for the modules that made the data products
  std::string fRawDigitLabel,fHitLabel,fOpHitLabel,fOpDetWaveformLabel,fOpFlashLabel,fGEANTLabel; 
  // --- Input settings imported from the fcl
  std::string fGeometry;
  std::vector<std::string> fLabels, fInteraction;
  int fClusterAlgoAdjChannel, fDetectorSizeX, fDetectorSizeY, fClusterInd0MatchTime, fClusterInd1MatchTime,fClusterPreselectionNHit;
  float fClusterMatchTime, fClusterMatchNHit, fClusterMatchCharge, fAdjClusterTime, fAdjClusterRad, fAdjOpFlashRad, fAdjOpFlashTime, fAdjOpFlashMaxPECut, fAdjOpFlashMinPECut;
  double fClusterAlgoTime;
  bool /*fTestNewClReco, */fDebug;
  // --- Our TTrees, and its associated variables.
  TTree* fMCTruthTree;
  TTree* fInteractionTree;
  TTree* fLowEAnaTree;
  // --- MC Truth Variables
  int /*Run,SubRun,*/Event,Flag;
  std::vector<int> TPart;
  std::vector<std::map<int,simb::MCParticle>> Parts = {};
  // --- MC Interaction Variables
  std::string Interaction;
  int PDG;
  float Energy;
  std::vector<int> DaughterPDG;
  std::vector<double> Momentum,StartVertex,EndVertex;
  std::vector<double> DaughterE,DaughterPx,DaughterPy,DaughterPz,DaughterStartVx,DaughterStartVy,DaughterStartVz,DaughterEndVx,DaughterEndVy,DaughterEndVz;
  // --- Cluster Hit Variables
  bool Main;
  int Idx,Match;
  std::vector<int> Gen,Ind0Gen,Ind1Gen,TPC,Ind0TPC,Ind1TPC,Channel,Ind0Channel,Ind1Channel,Ind0MaxHit,Ind1MaxHit;
  std::vector<double> Time,Ind0T,Ind1T,Charge,Ind0Charge,Ind1Charge,Y,Ind0Y,Ind1Y,Z,Ind0Z,Ind1Z; 
  // --- Backtracking Variables
  std::vector<int> MotherPDG,AncestorPDG;
  std::vector<int> Ind0MotherPDG,Ind0AncestorPDG,Ind1MotherPDG,Ind1AncestorPDG;
  std::vector<double> MotherE,AncestorE,MotherP,AncestorP,MotherX,AncestorX,MotherY,AncestorY,MotherZ,AncestorZ;
  std::vector<double> Ind0MotherE,Ind0AncestorE,Ind0MotherP,Ind0AncestorP,Ind0MotherX,Ind0AncestorX,Ind1MotherY,Ind0AncestorY,Ind0MotherZ,Ind0AncestorZ;
  std::vector<double> Ind1MotherE,Ind1AncestorE,Ind1MotherP,Ind1AncestorP,Ind1MotherX,Ind1AncestorX,Ind0MotherY,Ind1AncestorY,Ind1MotherZ,Ind1AncestorZ;
  // --- Adj. Flash Variables
  std::vector<int> AssOpFlashGen,AssOpFlashPDG,AssOpFlashNHit;
  std::vector<int> OpFlashGen,OpFlashPDG,OpFlashNHit;
  std::vector<double> AssOpFlashT,AssOpFlashPE,AssOpFlashMaxPE,AssOpFlashX,AssOpFlashY,AssOpFlashZ,AssOpFlashPur;
  std::vector<double> OpFlashT,OpFlashPE,OpFlashMaxPE,OpFlashX,OpFlashY,OpFlashZ,OpFlashPur;

  // --- Histograms to fill about collection plane hits
  TH1I* hAdjHits; 
  TH1F* hAdjHitsADCInt; 
  // --- Declare our services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
};

#endif

//......................................................
LowEAna::LowEAna(fhicl::ParameterSet const & p):EDAnalyzer(p){this->reconfigure(p);}

//......................................................
void LowEAna::reconfigure(fhicl::ParameterSet const & p){
  fLabels = p.get<std::vector<std::string>> ("ParticleLabelVector");
  fInteraction = p.get<std::vector<std::string>> ("InteractionLabelVector",{"nCapture"});
    
  fRawDigitLabel           = p.get<std::string> ("RawDigitLabel");
  fHitLabel                = p.get<std::string> ("HitLabel");
  fOpFlashLabel            = p.get<std::string> ("OpFlashLabel");
  fOpHitLabel              = p.get<std::string> ("OpHitLabel");
  fGEANTLabel              = p.get<std::string> ("GEANT4Label");

  fGeometry                = p.get<std::string> ("Geometry");
  fDetectorSizeX           = p.get<int>         ("DetectorSizeX");
  fDetectorSizeY           = p.get<int>         ("DetectorSizeY");
  fClusterMatchNHit        = p.get<float>       ("ClusterMatchNHit");
  fClusterMatchCharge      = p.get<float>       ("ClusterMatchCharge");
  fClusterMatchTime        = p.get<float>       ("ClusterMatchTime");
  fClusterInd0MatchTime    = p.get<float>       ("ClusterInd0MatchTime");
  fClusterInd1MatchTime    = p.get<float>       ("ClusterInd1MatchTime");
  
  fClusterPreselectionNHit = p.get<int>         ("ClusterPreselectionNHit");

  fClusterAlgoTime         = p.get<double>      ("ClusterAlgoTime");
  fClusterAlgoAdjChannel   = p.get<int>         ("ClusterAlgoAdjChannel");

  fAdjClusterTime          = p.get<float>       ("AdjClusterTime");
  fAdjClusterRad           = p.get<float>       ("AdjClusterRad");
  fAdjOpFlashTime          = p.get<float>       ("AdjOpFlashTime");
  fAdjOpFlashRad           = p.get<float>       ("AdjOpFlashRad");
  fAdjOpFlashMaxPECut      = p.get<float>       ("AdjOpFlashMaxPECut");
  fAdjOpFlashMinPECut      = p.get<float>       ("AdjOpFlashMinPECut");
  fDebug                   = p.get<bool>        ("Debug",true);
} // Reconfigure

//......................................................
void LowEAna::beginJob(){   
  // --- Make our handle to the TFileService
  art::ServiceHandle<art::TFileService> tfs;
  fMCTruthTree = tfs->make<TTree>("MCTruthTree","MC Truth Tree");
  fInteractionTree = tfs->make<TTree>("MCInteraction","MC Event Tree");
  fLowEAnaTree = tfs->make<TTree>("LowEAnaTree","Hit Ana Tree");
    
  // ----------------------------------------------- //
  // --------- Declare our TTree branches----------- //
  // ------------- MC Truth Variables -------------- //
  // ----------------------------------------------- //
    
  fMCTruthTree -> Branch("Event",                &Event,          "Event/I");  // Event number.
  fMCTruthTree -> Branch("Flag",                 &Flag,           "Flag/I");   // Flag used to match truth with reco tree entries.
  fMCTruthTree -> Branch("TruthPart",            &TPart);                      // Number particles per generator.

  fInteractionTree -> Branch("Event",            &Event,          "Event/I");  // Event number.
  fInteractionTree -> Branch("Flag",             &Flag,           "Flag/I");   // Flag used to match truth with reco tree entries.
  fInteractionTree -> Branch("PDG",              &PDG,            "PDG/I");    // Main interacting particle PDG.
  fInteractionTree -> Branch("Energy",           &Energy,         "Energy/F"); // Main interacting particle energy [GeV^2].
  fInteractionTree -> Branch("Interaction",      &Interaction);                // Type of interaction.
  fInteractionTree -> Branch("Momentum",         &Momentum);                   // Main interacting particle momentum [GeV^2].
  fInteractionTree -> Branch("StartVertex",      &StartVertex);                // Main interacting particle start vertex [cm].
  fInteractionTree -> Branch("EndVertex",        &EndVertex);                  // Main interacting particle end vertex [cm].
  fInteractionTree -> Branch("DaughterPDG",      &DaughterPDG);                // Main interacting particle daughter PDG.
  fInteractionTree -> Branch("DaughterE",        &DaughterE);                  // Main interacting particle daughter energy [GeV^2].
  fInteractionTree -> Branch("DaughterPx",       &DaughterPx);                 // Main interacting particle daughter momentum X [GeV^2].
  fInteractionTree -> Branch("DaughterPy",       &DaughterPy);                 // Main interacting particle daughter momentum Y [GeV^2].
  fInteractionTree -> Branch("DaughterPz",       &DaughterPz);                 // Main interacting particle daughter momentum Z [GeV^2].
  fInteractionTree -> Branch("DaughterStartVx",  &DaughterStartVx);            // Main interacting particle daughter start vertex X [cm].
  fInteractionTree -> Branch("DaughterStartVy",  &DaughterStartVy);            // Main interacting particle daughter start vertex Y [cm].
  fInteractionTree -> Branch("DaughterStartVz",  &DaughterStartVz);            // Main interacting particle daughter start vertex Z [cm].
  fInteractionTree -> Branch("DaughterEndVx",    &DaughterEndVx);              // Main interacting particle daughter end vertex X [cm].
  fInteractionTree -> Branch("DaughterEndVy",    &DaughterEndVy);              // Main interacting particle daughter end vertex Y [cm].
  fInteractionTree -> Branch("DaughterEndVz",    &DaughterEndVz);              // Main interacting particle daughter end vertex Z [cm].

  fLowEAnaTree -> Branch("Event",            &Event,          "Event/I");  // Event number.
  fLowEAnaTree -> Branch("Flag",             &Flag,           "Flag/I");   // Flag used to match truth with reco tree entries.
  fLowEAnaTree -> Branch("Index",            &Idx,            "Idx/I");    // Flag used to match truth with reco tree entries.
  fLowEAnaTree -> Branch("Main",             &Main,           "Main/O");   // Cluster bool (can be true or false main/adj cluster tagging).     
  fLowEAnaTree -> Branch("Match",            &Match,          "Match/I");  // Flag used to match truth with reco tree entries.
  fLowEAnaTree -> Branch("Gen",              &Gen);                        // Hit generator idx.
  fLowEAnaTree -> Branch("Ind0Gen",          &Ind0Gen);                    // Hit generator idx.  
  fLowEAnaTree -> Branch("Ind1Gen",          &Ind1Gen);                    // Hit generator idx.  
  fLowEAnaTree -> Branch("TPC",              &TPC);                        // Hit TPC.
  fLowEAnaTree -> Branch("Ind0TPC",          &Ind0TPC);                    // Hit ind0 TPC.
  fLowEAnaTree -> Branch("Ind1TPC",          &Ind1TPC);                    // Hit ind1 TPC.
  fLowEAnaTree -> Branch("Channel",          &Channel);                    // Hit Channel.
  fLowEAnaTree -> Branch("Ind0Channel",      &Ind0Channel);                // Hit ind0 Channel.
  fLowEAnaTree -> Branch("Ind1Channel",      &Ind1Channel);                // Hit ind1 Channel.  
  fLowEAnaTree -> Branch("Time",             &Time);                       // Hit time [ticks].
  fLowEAnaTree -> Branch("Ind0T",            &Ind0T);                      // Hit ind0 DT [Ticks].
  fLowEAnaTree -> Branch("Ind1T",            &Ind1T);                      // Hit ind1 DT [Ticks].
  fLowEAnaTree -> Branch("Charge",           &Charge);                     // Hit charge [ADC*ticks].
  fLowEAnaTree -> Branch("Ind0Charge",       &Ind0Charge);                 // Hit ind0 MaxHit.
  fLowEAnaTree -> Branch("Ind1Charge",       &Ind1Charge);                 // Hit ind1 MaxHit.
  fLowEAnaTree -> Branch("Y",                &Y);                          // Hit reco Y [cm]
  fLowEAnaTree -> Branch("Ind0Y",            &Ind0Y);                      // Hit ind0 reco Y [cm]
  fLowEAnaTree -> Branch("Ind1Y",            &Ind1Y);                      // Hit ind1 reco Y [cm]
  fLowEAnaTree -> Branch("Z",                &Z);                          // Hit reco Z [cm]
  fLowEAnaTree -> Branch("Ind0Z",            &Ind0Z);                      // Hit ind0 reco Z [cm]
  fLowEAnaTree -> Branch("Ind1Z",            &Ind1Z);                      // Hit ind1 reco Z [cm]
  
  // Per cluster backtracking info.
  fLowEAnaTree -> Branch("MotherPDG",        &MotherPDG);                  // Hit mother pdg
  fLowEAnaTree -> Branch("MotherE",          &MotherE);                    // Hit mother energy [GeV]
  fLowEAnaTree -> Branch("MotherP",          &MotherP);                    // Hit mother momentum [GeV]
  fLowEAnaTree -> Branch("MotherX",          &MotherX);                    // Hit mother vertex X [cm]
  fLowEAnaTree -> Branch("MotherY",          &MotherY);                    // Hit mother vertex Y [cm]
  fLowEAnaTree -> Branch("MotherZ",          &MotherZ);                    // Hit mother vertex Z [cm]
  fLowEAnaTree -> Branch("AncestorPDG",      &AncestorPDG);                // Hit ancestor pdg
  fLowEAnaTree -> Branch("AncestorE",        &AncestorE);                  // Hit ancestor energy [GeV]
  fLowEAnaTree -> Branch("AncestorP",        &AncestorP);                  // Hit ancestor momentum [GeV]
  fLowEAnaTree -> Branch("AncestorX",        &AncestorX);                  // Hit ancestor vertex X [cm]
  fLowEAnaTree -> Branch("AncestorY",        &AncestorY);                  // Hit ancestor vertex Y [cm]
  fLowEAnaTree -> Branch("AncestorZ",        &AncestorZ);                  // Hit ancestor vertex Z [cm]
    
  // Adj. Flash info.
  fLowEAnaTree -> Branch("OpFlashGen",       &AssOpFlashGen);              // Adj. OpFlash' generator idx.
  fLowEAnaTree -> Branch("OpFlashT",         &AssOpFlashT);                // Adj. OpFlash' time [ticks]
  fLowEAnaTree -> Branch("OpFlashNHit",      &AssOpFlashNHit);             // Adj. OpFlash' #hits
  fLowEAnaTree -> Branch("OpFlashPE",        &AssOpFlashPE);               // Adj. OpFlash' tot #PE [ADC*ticks]
  fLowEAnaTree -> Branch("OpFlashMaxPE",     &AssOpFlashMaxPE);            // Adj. OpFlash' max #PE [ADC*ticks]
  fLowEAnaTree -> Branch("OpFlashX",         &AssOpFlashX);                // Adj. OpFlash' reco X [cm]
  fLowEAnaTree -> Branch("OpFlashY",         &AssOpFlashY);                // Adj. OpFlash' reco Y [cm] 
  fLowEAnaTree -> Branch("OpFlashZ",         &AssOpFlashZ);                // Adj. OpFlash' reco Z [cm]
  fLowEAnaTree -> Branch("OpFlashPur",       &AssOpFlashPur);              // Adj. OpFlash' purity

  hAdjHits        = tfs->make<TH1I>("hAdjHits", "Number of hits; Number of hits; Number of events"  , 21, -0.5, 20.5 );
  hAdjHitsADCInt  = tfs->make<TH1F>("hAdjHitsADCInt","Total summed ADC Integrals; Total summed ADC Integrals; Number of events"  , 1000, 0, 10000 );
} // BeginJob

//......................................................
void LowEAna::analyze(art::Event const & evt){ 
  //--------------------------------------------------------------------------------------------//
  //------------------------------- Prepare everything for new event ---------------------------// 
  //--------------------------------------------------------------------------------------------//
  std::map<int,simb::MCParticle>                                               ThisGeneratorParts;
  std::vector<recob::Hit>                                                 Hits0,Hits1,Hits2,Hits3; 
  std::vector<std::vector<recob::Hit>>                           Hits = {Hits0,Hits1,Hits2,Hits3};
  std::vector<std::vector<recob::Hit>>                 Clusters0, Clusters1, Clusters2, Clusters3;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  // ---------- We want to reset all of our previous run and TTree variables -------------------// 
  // Run = evt.run();SubRun = evt.subRun();
  Event = evt.event();
  ResetEventVariables(clockData);
  
  //--------------------------------------------------------------------------------------------//
  //--------------------------------- Create maps for ID tracking ------------------------------//
  //--------------------------------------------------------------------------------------------//
  // -------- Fill MC Truth IDs to tracking vectors. Get a list of all of my particles ---------//
  mf::LogInfo lparticle("particles");
  std::string lparticlestr = "";
  const sim::ParticleList& PartList = pi_serv->ParticleList();
  lparticle << "\nTotal #particles: " << PartList.size();
  
  // Loop over all signal+bkg handles and collect track IDs
  for ( size_t i = 0; i < fLabels.size(); i++){
    Parts.push_back(ThisGeneratorParts); // For each label insert empty list
    
    art::Handle<std::vector<simb::MCTruth>> ThisHandle;
    evt.getByLabel(fLabels[i], ThisHandle);
    
    if(ThisHandle){
      lparticlestr = PrintInColor(lparticlestr,"\n"+fLabels[i]+" *is generated!",GetColor("green"));
      
      auto ThisValidHanlde = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[i]); 
      art::FindManyP<simb::MCParticle> Assn(ThisValidHanlde,evt,fGEANTLabel);            
      
      FillMyMaps( Parts[i], Assn, ThisValidHanlde);                                                                           
      FillMCInteractionTree(Parts[i], fInteraction, fDebug);
      TPart.push_back(Parts[i].size()); // Insert #signal+bkg particles generated
      lparticlestr = PrintInColor(lparticlestr,"\n-> #Particles: "+str(int(Parts[i].size())),GetColor("blue"));
    }
    else{
      TPart.push_back(0);
      lparticlestr = PrintInColor(lparticlestr,"\n"+fLabels[i]+": *not generated!",GetColor("yellow"));
    }
  }
  lparticle << lparticlestr;
  fMCTruthTree->Fill();
  std::cout << std::endl;//--------------------------------------------------------------------//
  //--------------------------------- Optical Flash Collection --------------------------------// 
  //-------------------------------------------------------------------------------------------//
  mf::LogInfo lflash("flashes");
  std::string lflashstr = "";
  // Find OpHits and OpFlashes associated with the event
  art::Handle< std::vector< recob::OpHit >> OpHitHandle;
  art::Handle< std::vector< recob::OpFlash >> FlashHandle;
  std::vector<art::Ptr<recob::OpHit >> OpHitList;
  std::vector<art::Ptr<recob::OpFlash >> OpFlashList;
  if (evt.getByLabel(fOpHitLabel, OpHitHandle)){art::fill_ptr_vector(OpHitList, OpHitHandle);}
  if (evt.getByLabel(fOpFlashLabel, FlashHandle)){art::fill_ptr_vector(OpFlashList, FlashHandle);}
  art::FindManyP< recob::OpHit > AssOpHits(OpFlashList, evt, fOpFlashLabel);
  
  // Grab assns with OpHits to get match to neutrino purity
  lflash << "\nTotal number of flashes constructed: " << OpFlashList.size();
  // Loop over flashlist and assign OpHits to each flash
  for ( int i = 0; i < int(OpFlashList.size()); i++ ){

    recob::OpFlash ThisFlash = *OpFlashList[i];
    // if (i%10 == 0) PrintInColor("Flash Time = " + str(TheFlash.Time()), GetColor("red"));
    std::vector< art::Ptr< recob::OpHit > > MOpHits = AssOpHits.at(i);
    
    // Calculate the total PE of the flash and the time of the ophit with the highest PE 
    double OpHitTotPE = 0; double MaxHitPE = 0;
    float OpHitT  = 0; 
    float OpHitPE = 0;
    for (int j = 0; j < int(MOpHits.size()); j++){
      recob::OpHit OpHit = *MOpHits[j];
      OpHitTotPE += OpHit.PE();
      OpHitPE = OpHit.PE();
      if (OpHitPE > MaxHitPE){
        MaxHitPE = OpHitPE;
        OpHitT = OpHit.PeakTime();
      }
    }
    
    lflash << "\nEvaluating Flash purity";
    lflash << "\nPE of this OpFlash " << OpHitTotPE; 
    lflash << "\nOpFlash time " << OpHitT;

    // Get trackID from Parts
    // Calculate the flash purity, only for the Marley events
    if (MaxHitPE/OpHitTotPE < fAdjOpFlashMaxPECut && OpHitTotPE > fAdjOpFlashMinPECut){
      // OpFlash Gen is ind of OpFlashPurVector with highest purity
      OpFlashPE.push_back(ThisFlash.TotalPE());
      OpFlashMaxPE.push_back(MaxHitPE);
      OpFlashX.push_back(ThisFlash.XCenter());
      OpFlashY.push_back(ThisFlash.YCenter());
      OpFlashZ.push_back(ThisFlash.ZCenter());
      OpFlashT.push_back(ThisFlash.Time());
      OpFlashNHit.push_back(MOpHits.size());
      std::vector<double> OpFlashPurVector = {};
      for (size_t i = 0; i < Parts.size(); i++){
        std::vector<int> GenTrackIDs = {};
        for ( std::map<int,simb::MCParticle>::iterator iter = Parts[i].begin(); iter != Parts[i].end(); iter++ ){
          if (iter->second.TrackId() != 0){
            GenTrackIDs.push_back(iter->second.TrackId());
          }
        }
        // Convert std::vector<int> to std::set<int> for signal_trackids
        std::set<int> SetGenTrackIDs(GenTrackIDs.begin(), GenTrackIDs.end());
        // Calculate the flash purity, for this generator's trackIDs
        int test = supress_stdout();
        double ThisOpFlashPur = pbt->OpHitCollectionPurity(SetGenTrackIDs, MOpHits);
        resume_stdout(test);
        OpFlashPurVector.push_back(ThisOpFlashPur);
      }
      double MaxOpFlashPur = 0; int Gen = 0; int MaxPur = 0;
      for (size_t pur = 0; pur < OpFlashPurVector.size(); pur++){
        if (OpFlashPurVector[pur] > MaxOpFlashPur){Gen = pur+1; MaxPur = OpFlashPurVector[pur];}
      }
      OpFlashGen.push_back(Gen);
      OpFlashPur.push_back(MaxPur);
    }
  }
  // Build a struct with all the flash vectors
  OpFlashes EventOpFlashes = {OpFlashPE,OpFlashMaxPE,OpFlashX,OpFlashY,OpFlashZ,OpFlashT,OpFlashNHit,OpFlashGen,OpFlashPur};
  lflash << lflashstr;
  std::cout << std::endl;//--------------------------------------------------------------------//
  //------------------------------- Hit collection and assignment -----------------------------// 
  //-------------------------------------------------------------------------------------------//
  mf::LogInfo lhit("hits");
  std::string lhitstr = "";
  // --- Lift out the reco hits:
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(fHitLabel);
  int NTotHits = reco_hits->size();

  for(int hit = 0; hit < NTotHits; ++hit){
    // --- Loop over the reconstructed hits to separate them among tpc planes according to view
    
    recob::Hit const& ThisHit = reco_hits->at(hit);
    if (ThisHit.PeakTime() < 0){
      lhitstr = PrintInColor(lhitstr,"Negative Hit Time = " + str(ThisHit.PeakTime()), GetColor("red"));
    }

    if      (ThisHit.SignalType() == 0 && ThisHit.View() == 0){Hits0.push_back( ThisHit );} // SignalType = 0
    else if (ThisHit.SignalType() == 0 && ThisHit.View() == 1){Hits1.push_back( ThisHit );} // SignalType = 0
    else if (ThisHit.SignalType() == 1)                       {Hits2.push_back( ThisHit );} // SignalType = 1
    else    {Hits3.push_back( ThisHit ); lhit << "Hit was found with view out of scope";}
  } 

  lhit << "\n# Hits per view: ";
  lhit << "\nInduction Plane 0:\t" << Hits0.size();
  lhit << "\nInduction Plane 1:\t" << Hits1.size();
  lhit << "\nCollection Plane: \t" << Hits2.size();
  lhit << lhitstr;
  std::cout << std::endl;//--------------------------------------------------------------------//
  //------------------------------------- Cluster Hit analysis --------------------------------// 
  //-------------------------------------------------------------------------------------------//
  mf::LogInfo lcluster("clusters");
  std::string lclusterstr = "";
  // --- Now calculate the clusters ...
  CalcAdjHits(Hits0,Clusters0,hAdjHits,hAdjHitsADCInt,clockData,false);
  CalcAdjHits(Hits1,Clusters1,hAdjHits,hAdjHitsADCInt,clockData,false);
  CalcAdjHits(Hits2,Clusters2,hAdjHits,hAdjHitsADCInt,clockData,false);
  CalcAdjHits(Hits3,Clusters3,hAdjHits,hAdjHitsADCInt,clockData,false);

  lcluster << "\n# Clusters per view: ";
  lcluster << "\nInduction Plane 0:\t" << Clusters0.size();
  lcluster << "\nInduction Plane 1:\t" << Clusters1.size();
  lcluster << "\nCollection Plane: \t" << Clusters2.size();
  // if (Clusters3.size() != 0) PrintInColor("Other: "+str(Clusters3.size()),GetColor("red"));
  
  std::vector< std::vector< std::vector<recob::Hit>>> AllClusters; 
  std::vector< std::vector< std::vector<recob::Hit>>> MatchedClusters; 
  AllClusters = {Clusters0,Clusters1,Clusters2};
  std::vector< std::vector<int>>  ClNHits = {{},{},{}};
  std::vector< std::vector<float>>  ClT = {{},{},{}};
  std::vector< std::vector<float>>  ClCharge = {{},{},{}};

  MatchedClusters = MatchClusters(AllClusters,ClNHits,ClT,ClCharge,fDebug);
  std::vector< std::vector<float>>  ClY = ClT;
  std::vector< std::vector<float>>  ClZ = ClT;
  // Reset all vector entries in ClY and ClZ to -1e6
  for (size_t i = 0; i < ClY.size(); i++){
    for (size_t j = 0; j < ClY[i].size(); j++){
      ClY[i][j] = -1e6;
      ClZ[i][j] = -1e6;
    }
  }
  lcluster << lclusterstr;
  std::cout << std::endl;//------------- Cluster Preselection ------------------------// 

  ResetClusterVariables();
  for (int ii = 0; ii < int(MatchedClusters[2].size()); ii++){
    mf::LogInfo lpreselection("preselection");
    std::string lpreselectionstr = "";
    // std::cout << "Evaluating cluster #" << ii << std::endl;    
    if (MatchedClusters[2][ii].empty()){continue;}
    if (ClNHits[2][ii] <= fClusterPreselectionNHit){continue;}
    // --- Remaining are MainCls
    Idx = ii;
    if (!MatchedClusters[0][ii].empty() && !MatchedClusters[1][ii].empty()){
      lpreselectionstr = PrintInColor(lpreselectionstr,"Found MainCl with match in Ind0 & Ind1 #"+str(ii),GetColor("green"));
      // --- Declare our vectors to fill
      std::vector<double> Dir, Ind0Dir, Ind1Dir;
      Dir = Ind0Dir = Ind1Dir = {};
      FillClusterHitVectors(MatchedClusters[0][ii], Ind0TPC, Ind0Channel, Ind0MotherX, Ind0MotherY, Ind0MotherZ, Ind0MotherE, Ind0MotherP, Ind0MotherPDG, Ind0AncestorX, Ind0AncestorY, Ind0AncestorZ, Ind0AncestorE, Ind0AncestorP, Ind0AncestorPDG, Ind0Charge, Ind0T, Ind0Y, Ind0Z, Ind0Dir, clockData, fDebug);
      FillClusterHitVectors(MatchedClusters[1][ii], Ind1TPC, Ind1Channel, Ind1MotherX, Ind1MotherY, Ind1MotherZ, Ind1MotherE, Ind1MotherP, Ind1MotherPDG, Ind1AncestorX, Ind1AncestorY, Ind1AncestorZ, Ind1AncestorE, Ind1AncestorP, Ind1AncestorPDG, Ind1Charge, Ind1T, Ind1Y, Ind1Z, Ind1Dir, clockData, fDebug);
      FillClusterHitVectors(MatchedClusters[2][ii], TPC, Channel,MotherX, MotherY, MotherZ, MotherE, MotherP, MotherPDG, AncestorX, AncestorY, AncestorZ, AncestorE, AncestorP, AncestorPDG, Charge, Time, Y, Z, Dir, clockData, fDebug);
      std::vector<double> RecoY0 = ComputeInterpolationRecoY(Event,Ind1TPC,Z,Time,Ind0Z,Ind0Y,Ind0T,Ind0Dir,fDebug); 
      std::vector<double> RecoY1 = ComputeInterpolationRecoY(Event,Ind0TPC,Z,Time,Ind1Z,Ind1Y,Ind1T,Ind1Dir,fDebug);
      
      for (size_t i = 0; i < Time.size(); i++){Y[i] = ((RecoY0[i]+RecoY1[i])/2);}

      ClY[2][ii] = Average(Y); ClY[1][ii] = Average(Ind1Y); ClY[0][ii] = Average(Ind0Y);
      ClZ[2][ii] = Average(Z); ClZ[1][ii] = Average(Ind1Z); ClZ[0][ii] = Average(Ind0Z);
      AssignOpFlashes(EventOpFlashes,ClY[2][ii],ClZ[2][ii]);
      Main = true;
      Match = 2;
      fLowEAnaTree->Fill();
      ResetClusterVariables();
    }
    else if (!MatchedClusters[0][ii].empty() && MatchedClusters[1][ii].empty()){
      // --- Declare our vectors to fill
      std::vector<double> Dir, Ind0Dir;
      Dir = Ind0Dir = {};
      lpreselectionstr = PrintInColor(lpreselectionstr,"Found MainCl with match in Ind0 #"+str(ii),GetColor("green"));
      FillClusterHitVectors(MatchedClusters[0][ii], Ind0TPC, Ind0Channel, Ind0MotherX, Ind0MotherY, Ind0MotherZ, Ind0MotherE, Ind0MotherP, Ind0MotherPDG, Ind0AncestorX, Ind0AncestorY, Ind0AncestorZ, Ind0AncestorE, Ind0AncestorP, Ind0AncestorPDG, Ind0Charge, Ind0T, Ind0Y, Ind0Z, Ind0Dir, clockData, fDebug);
      FillClusterHitVectors(MatchedClusters[2][ii], TPC, Channel, MotherX, MotherY, MotherZ, MotherE, MotherP, MotherPDG, AncestorX, AncestorY, AncestorZ, AncestorE, AncestorP, AncestorPDG, Charge, Time, Y, Z, Dir, clockData, fDebug);
      Y = ComputeInterpolationRecoY(Event,Ind0TPC,Z,Time,Ind0Z,Ind0Y,Ind0T,Ind0Dir,fDebug);
      ClY[2][ii] = Average(Y); ClY[0][ii] = Average(Ind0Y);
      ClZ[2][ii] = Average(Z); ClZ[0][ii] = Average(Ind0Z);
      AssignOpFlashes(EventOpFlashes,ClY[2][ii],ClZ[2][ii]);
      Main = true;
      Match = 0;
      fLowEAnaTree->Fill();
      ResetClusterVariables();
    }  
    else if (MatchedClusters[0][ii].empty() && !MatchedClusters[1][ii].empty()){
      // --- Declare our vectors to fill
      std::vector<double> Dir, Ind1Dir;
      Dir = Ind1Dir = {};
      lpreselectionstr = PrintInColor(lpreselectionstr,"Found MainCl with match in Ind1 #"+str(ii),GetColor("green"));
      FillClusterHitVectors(MatchedClusters[1][ii], Ind1TPC, Ind1Channel, Ind1MotherX, Ind1MotherY, Ind1MotherZ, Ind1MotherE, Ind1MotherP, Ind1MotherPDG, Ind1AncestorX, Ind1AncestorY, Ind1AncestorZ, Ind1AncestorE, Ind1AncestorP, Ind1AncestorPDG, Ind1Charge, Ind1T, Ind1Y, Ind1Z, Ind1Dir, clockData, fDebug);
      FillClusterHitVectors(MatchedClusters[2][ii], TPC, Channel, MotherX, MotherY, MotherZ, MotherE, MotherP, MotherPDG, AncestorX, AncestorY, AncestorZ, AncestorE, AncestorP, AncestorPDG, Charge, Time, Y, Z, Dir, clockData, fDebug);
      Y = ComputeInterpolationRecoY(Event,Ind1TPC,Z,Time,Ind1Z,Ind1Y,Ind1T,Ind1Dir,fDebug);
      ClY[2][ii] = Average(Y); ClY[1][ii] = Average(Ind1Y);
      ClZ[2][ii] = Average(Z); ClZ[1][ii] = Average(Ind1Z);
      AssignOpFlashes(EventOpFlashes,ClY[2][ii],ClZ[2][ii]);
      Main = true;
      Match = 1;
      fLowEAnaTree->Fill();
      ResetClusterVariables();
    }
    else{lpreselectionstr = PrintInColor(lpreselectionstr,"No match found! THIS SHOULD NOT HAPPEN, CLUSTERS HAVE ALREADY BEEN MATCHED!",GetColor("red"));}

    for (int jj = 0; jj < int(MatchedClusters[2].size()); jj++){
      if (jj == ii){continue;}
      if (MatchedClusters[2][jj].empty()){continue;}
      if (abs(ClT[2][ii] - ClT[2][jj]) > fAdjClusterTime){continue;}
      if (!MatchedClusters[0][jj].empty() && !MatchedClusters[1][jj].empty()){
        lpreselectionstr = PrintInColor(lpreselectionstr,"AdjCl cluster found!",GetColor("green"));
        // --- Declare our vectors to fill
        std::vector<double> Dir, Ind0Dir, Ind1Dir;
        Dir = Ind0Dir = Ind1Dir = {};
        FillClusterHitVectors(MatchedClusters[0][jj], Ind0TPC, Ind0Channel, Ind0MotherX, Ind0MotherY, Ind0MotherZ, Ind0MotherE, Ind0MotherP, Ind0MotherPDG, Ind0AncestorX, Ind0AncestorY, Ind0AncestorZ, Ind0AncestorE, Ind0AncestorP, Ind0AncestorPDG, Ind0Charge, Ind0T, Ind0Y, Ind0Z, Ind0Dir, clockData, fDebug);
        FillClusterHitVectors(MatchedClusters[1][jj], Ind1TPC, Ind1Channel, Ind1MotherX, Ind1MotherY, Ind1MotherZ, Ind1MotherE, Ind1MotherP, Ind1MotherPDG, Ind1AncestorX, Ind1AncestorY, Ind1AncestorZ, Ind1AncestorE, Ind1AncestorP, Ind1AncestorPDG, Ind1Charge, Ind1T, Ind1Y, Ind1Z, Ind1Dir, clockData, fDebug);
        FillClusterHitVectors(MatchedClusters[2][jj], TPC, Channel, MotherX, MotherY, MotherZ, MotherE, MotherP, MotherPDG, AncestorX, AncestorY, AncestorZ, AncestorE, AncestorP, AncestorPDG, Charge, Time, Y, Z, Dir, clockData, fDebug);
        std::vector<double> RecoY0 = ComputeInterpolationRecoY(Event,Ind1TPC,Z,Time,Ind0Z,Ind0Y,Ind0T,Ind0Dir,fDebug); 
        std::vector<double> RecoY1 = ComputeInterpolationRecoY(Event,Ind0TPC,Z,Time,Ind1Z,Ind1Y,Ind1T,Ind1Dir,fDebug);
        
        for (size_t i = 0; i < Time.size(); i++){Y[i] = ((RecoY0[i]+RecoY1[i])/2);}

        ClY[2][jj] = Average(Y); ClY[1][jj] = Average(Ind1Y); ClY[0][jj] = Average(Ind0Y);
        ClZ[2][jj] = Average(Z); ClZ[1][jj] = Average(Ind1Z); ClZ[0][jj] = Average(Ind0Z);
        if (sqrt(pow(ClY[2][jj]-ClY[2][ii],2)+pow(ClZ[2][jj]-ClZ[2][ii],2)) > fAdjClusterRad){continue;}
        AssignOpFlashes(EventOpFlashes,ClY[2][jj],ClZ[2][jj]);
        fLowEAnaTree->Fill();
        ResetClusterVariables();
      }
      else if (!MatchedClusters[0][jj].empty() && MatchedClusters[1][jj].empty()){
        lpreselectionstr = PrintInColor(lpreselectionstr,"AdjCl found!",GetColor("green"));
        // --- Declare our vectors to fill
        std::vector<double> Dir, Ind0Dir;
        Dir = Ind0Dir = {};
        FillClusterHitVectors(MatchedClusters[0][jj], Ind0TPC, Ind0Channel, Ind0MotherX, Ind0MotherY, Ind0MotherZ, Ind0MotherE, Ind0MotherP, Ind0MotherPDG, Ind0AncestorX, Ind0AncestorY, Ind0AncestorZ, Ind0AncestorE, Ind0AncestorP, Ind0AncestorPDG, Ind0Charge, Ind0T, Ind0Y, Ind0Z, Ind0Dir, clockData, fDebug);
        FillClusterHitVectors(MatchedClusters[2][jj], TPC, Channel, MotherX, MotherY, MotherZ, MotherE, MotherP, MotherPDG, AncestorX, AncestorY, AncestorZ, AncestorE, AncestorP, AncestorPDG, Charge, Time, Y, Z, Dir, clockData, fDebug);
        Y = ComputeInterpolationRecoY(Event,Ind0TPC,Z,Time,Ind0Z,Ind0Y,Ind0T,Ind0Dir,fDebug);
        ClY[2][jj] = Average(Y); ClY[0][jj] = Average(Ind0Y);
        ClZ[2][jj] = Average(Z); ClZ[0][jj] = Average(Ind0Z);
        if (sqrt(pow(ClY[2][jj]-ClY[2][ii],2)+pow(ClZ[2][jj]-ClZ[2][ii],2)) > fAdjClusterRad){continue;}
        AssignOpFlashes(EventOpFlashes,ClY[2][jj],ClZ[2][jj]);
        fLowEAnaTree->Fill();
        ResetClusterVariables();
      }
      else if (MatchedClusters[0][jj].empty() && !MatchedClusters[1][jj].empty()){
        // --- Declare our vectors to fill
        std::vector<double> Dir, Ind1Dir;
        Dir = Ind1Dir = {};
        lpreselectionstr = PrintInColor(lpreselectionstr,"AdjCl found!",GetColor("green"));
        FillClusterHitVectors(MatchedClusters[1][jj], Ind1TPC, Ind1Channel, Ind1MotherX, Ind1MotherY, Ind1MotherZ, Ind1MotherE, Ind1MotherP, Ind1MotherPDG, Ind1AncestorX, Ind1AncestorY, Ind1AncestorZ, Ind1AncestorE, Ind1AncestorP, Ind1AncestorPDG, Ind1Charge, Ind1T, Ind1Y, Ind1Z, Ind1Dir, clockData, fDebug);
        FillClusterHitVectors(MatchedClusters[2][jj], TPC, Channel, MotherX, MotherY, MotherZ, MotherE, MotherP, MotherPDG, AncestorX, AncestorY, AncestorZ, AncestorE, AncestorP, AncestorPDG, Charge, Time, Y, Z, Dir, clockData, fDebug);
        Y = ComputeInterpolationRecoY(Event,Ind1TPC,Z,Time,Ind1Z,Ind1Y,Ind1T,Ind1Dir,fDebug);
        ClY[2][jj] = Average(Y); ClY[1][jj] = Average(Ind1Y);
        ClZ[2][jj] = Average(Z); ClZ[1][jj] = Average(Ind1Z);
        if (sqrt(pow(ClY[2][jj]-ClY[2][ii],2)+pow(ClZ[2][jj]-ClZ[2][ii],2)) > fAdjClusterRad){continue;}
        AssignOpFlashes(EventOpFlashes,ClY[2][jj],ClZ[2][jj]);
        fLowEAnaTree->Fill();
        ResetClusterVariables();
      }
    } // End of loop over AdjCl
  } // Loop over MainCl
} // End of analyze

//########################################################################################################################################//
//_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_//
//########################################################################################################################################//
/* List all functions here:
- ResetEventVariables():        Reset variables for each event
- ResetClusterVariables():      Reset variables for each cluster
- FillClusterHitVectors():      Fill the cluster hit vectors 
- ComputeInterpolationRecoY():  Compute the reco Y position of a hit using interpolation
- FillClusterVariables():       Fill the cluster variables
- CalcAdjHits():                Calculate the adjacent hits
- FillMyMaps():                 Fill the maps particle maps for each generator
- FillMCInteractionTree():      Fill the MC interaction tree
- BeginJob():                   Begin job
*/
//########################################################################################################################################//
//_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_//
//########################################################################################################################################//

//......................................................
void LowEAna::ResetEventVariables( detinfo::DetectorClocksData clockData ) 
// Reset variables for each event
{
  Flag = rand() % 10000000000;
  mf::LogInfo lheader("header");
  lheader << "\n#####################################################################";
  lheader << "\n - TPC Frequency in [MHz]: " << clockData.TPCClock().Frequency()      ;
  lheader << "\n - TPC Tick in [us]: " << clockData.TPCClock().TickPeriod()           ;
  lheader << "\n - Event Flag: " << Flag                                              ;
  lheader << "\n - Succesfull reset of variables for Event " << Event << ": " << Flag ; 
  lheader << "\n#####################################################################";
  // Clear Marley MCTruth info.
  Parts = {}; TPart = {};
  Idx = 0;
  // Clear OpFlash Vectors
  OpFlashGen.clear();OpFlashPur.clear();OpFlashNHit.clear();OpFlashPE.clear();OpFlashMaxPE.clear();
  OpFlashT.clear();OpFlashX.clear();OpFlashY.clear();OpFlashZ.clear();
  return;
} // ResetVariables

//......................................................
void LowEAna::ResetClusterVariables()
{
  Main = false;
  Match = -1;
  Gen = Ind0Gen = Ind1Gen = TPC = Ind0TPC = Ind1TPC = {};
  Channel = Ind0Channel = Ind1Channel = {};

  Time = Ind0T = Ind1T = {};
  Charge = Ind0Charge = Ind1Charge = Y = Ind0Y = Ind1Y = Z = Ind0Z = Ind1Z = {};
  
  MotherPDG = AncestorPDG = {};
  Ind0MotherPDG = Ind0AncestorPDG = {};
  Ind1MotherPDG = Ind1AncestorPDG = {};

  MotherE = MotherP = AncestorE = AncestorP = {};
  Ind0MotherE = Ind0MotherP = Ind0AncestorE = Ind0AncestorP = {};
  Ind1MotherE = Ind1MotherP = Ind1AncestorE = Ind1AncestorP = {};
  
  MotherX = MotherY = MotherZ = AncestorX = AncestorY = AncestorZ = {};
  Ind0MotherX = Ind0MotherY = Ind0MotherZ = Ind0AncestorX = Ind0AncestorY = Ind0AncestorZ = {};
  Ind1MotherX = Ind1MotherY = Ind1MotherZ = Ind1AncestorX = Ind1AncestorY = Ind1AncestorZ = {};

  AssOpFlashGen.clear();AssOpFlashNHit.clear();
  AssOpFlashPur.clear();AssOpFlashPE.clear();AssOpFlashMaxPE.clear();
  AssOpFlashT.clear();AssOpFlashX.clear();AssOpFlashY.clear();AssOpFlashZ.clear();
  return;
}

void LowEAna::AssignOpFlashes(
  struct OpFlashes &EventOpFlashes,
  float &ClY,
  float &ClZ
)
/*
*/
{
  // Extract all vectors from the struct
  for(size_t zz = 0; zz < EventOpFlashes.OpFlashPE.size(); zz++){
    // Associate OpFlashes to the event according to Y and Z
    if(sqrt(std::pow(abs(EventOpFlashes.OpFlashY[zz] - ClY),2) + std::pow(abs(EventOpFlashes.OpFlashZ[zz] - ClZ),2)) > fAdjOpFlashRad){continue;}
    AssOpFlashGen.push_back(EventOpFlashes.OpFlashGen[zz]);
    AssOpFlashMaxPE.push_back(EventOpFlashes.OpFlashMaxPE[zz]);
    AssOpFlashNHit.push_back(EventOpFlashes.OpFlashNHit[zz]);
    AssOpFlashPur.push_back(EventOpFlashes.OpFlashPur[zz]);
    AssOpFlashPE.push_back(EventOpFlashes.OpFlashPE[zz]);
    AssOpFlashT.push_back(EventOpFlashes.OpFlashT[zz]);
    AssOpFlashX.push_back(EventOpFlashes.OpFlashX[zz]);
    AssOpFlashY.push_back(EventOpFlashes.OpFlashY[zz]);
    AssOpFlashZ.push_back(EventOpFlashes.OpFlashZ[zz]);
  } // Loop over OpFlash  
  return;
}

//......................................................
std::vector<double> LowEAna::ComputeInterpolationRecoY(
  int Event,
  std::vector<int> &HIndTPC,
  std::vector<double> &Z, 
  std::vector<double> &Time,
  std::vector<double> &IndZ,
  std::vector<double> &IndY,
  std::vector<double> &IndT,
  std::vector<double> &IndDir,
  bool HeavDebug)
/*
*/
{ 
  mf::LogInfo linterp("interpolation");
  std::string linterpstr = "";
  // Create the interpolator
  std::vector<double> RecoY = {};

  // TFile * f = new TFile(("Interpolation"+str(Event)+"_"+str(Time[0])+".root").c_str(),"RECREATE");
  TGraph * IndInterpZ = new TGraph();
  TGraph * IndInterpY = new TGraph();
  for (size_t i = 0; i < IndT.size(); i++){
    // Fill TGraphs
    IndInterpY->SetPoint(i, IndT[i], IndY[i]);
    IndInterpZ->SetPoint(i, IndT[i], IndZ[i]);
  }

  for (size_t i = 0; i < Z.size(); i++){
    float ThisHZ = Z[i];
    float ThisHT = Time[i];
    float ThisHIndDir = IndDir[0];
    float ThisHRefY = IndInterpY->Eval(ThisHT);
    float ThisHRefZ = IndInterpZ->Eval(ThisHT);
    float ThisRecoY = ThisHRefY + (ThisHZ - ThisHRefZ)/(ThisHIndDir);
    RecoY.push_back(ThisRecoY);
    linterpstr = PrintInColor(linterpstr,"\nReco Y = "+str(ThisRecoY)+" cm",GetColor("red"));
  }
  // IndInterpY->Write("InterpY");
  // IndInterpZ->Write("InterpZ");
  // f->Close();
  linterp << linterpstr;
  return RecoY;
}

//......................................................
void LowEAna::FillClusterHitVectors(std::vector<recob::Hit> Cluster, 
  std::vector<int> &TPC,
  std::vector<int> &Channel,
  std::vector<double> &MotherX,
  std::vector<double> &MotherY,
  std::vector<double> &MotherZ,
  std::vector<double> &MotherE,
  std::vector<double> &MotherP,
  std::vector<int> &MotherPDG,
  std::vector<double> &AncestorX,
  std::vector<double> &AncestorY,
  std::vector<double> &AncestorZ,
  std::vector<double> &AncestorE,
  std::vector<double> &AncestorP,
  std::vector<int> &AncestorPDG,
  std::vector<double> &Charge,
  std::vector<double> &Time,
  std::vector<double> &Y,
  std::vector<double> &Z,
  std::vector<double> &Dir,
  detinfo::DetectorClocksData ClockData,
  bool HeavDebug)
/*
*/
{ 
  // --- Declare our vectors to fill
  std::vector<int> TrackID;
  TrackID = {};
  for (recob::Hit ThisHit : Cluster){
    std::vector< sim::TrackIDE > ThisHitID = bt_serv->HitToTrackIDEs(ClockData, ThisHit);
    int MainTrID = 0;
    float TopEFrac = 0;
    
    for (size_t i=0; i < ThisHitID.size(); ++i){
      if (ThisHitID[i].energyFrac > TopEFrac){
        TopEFrac = ThisHitID[i].energyFrac;
        MainTrID = ThisHitID[i].trackID; 
      }
    }
    TrackID.push_back(MainTrID);
    
    // Call backtracker to get the particle from the trackID
    const simb::MCParticle *HTruth;
    int test = supress_stdout();
    HTruth = pi_serv->TrackIdToParticle_P(MainTrID);
    resume_stdout(test);
    if (HTruth != 0){
      MotherX.push_back(HTruth->Vx());
      MotherY.push_back(HTruth->Vy());
      MotherZ.push_back(HTruth->Vz());
      MotherE.push_back(HTruth->E());
      MotherP.push_back(HTruth->P());
      MotherPDG.push_back(HTruth->PdgCode());
      const simb::MCParticle *Ancestor;
      int test = supress_stdout();
      Ancestor = pi_serv->TrackIdToParticle_P(HTruth->Mother());
      resume_stdout(test);
      
      if (Ancestor != 0){
        AncestorX.push_back(Ancestor->Vx());
        AncestorY.push_back(Ancestor->Vy());
        AncestorZ.push_back(Ancestor->Vz());
        AncestorE.push_back(Ancestor->E());
        AncestorP.push_back(Ancestor->P());
        AncestorPDG.push_back(Ancestor->PdgCode());
      }
    }

    TPC.push_back(ThisHit.WireID().TPC);
    Time.push_back(ThisHit.PeakTime());
    Charge.push_back(ThisHit.Integral());
    const geo::WireGeo* ThisWire = geo->GeometryCore::WirePtr(ThisHit.WireID());
    geo::Point_t hXYZ = ThisWire->GetCenter();
    geo::Point_t sXYZ = ThisWire->GetStart();
    geo::Point_t eXYZ = ThisWire->GetEnd();
    Y.push_back(hXYZ.Y());
    Z.push_back(hXYZ.Z());
    geo::Vector_t Direction = eXYZ - sXYZ;
    float dyds = Direction.Y(); 
    float dzds = Direction.Z();
    Dir.push_back(dzds/dyds);
  }
}

void LowEAna::FillClusterVariables(std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
  std::vector<std::vector<int>> &ClNHits, 
  std::vector<std::vector<float>> &ClT, 
  std::vector<std::vector<float>> &ClCharge, 
  bool HeavDebug)
{
  for (size_t idx = 0; idx < Clusters.size(); idx++){ 
    std::vector<std::vector<recob::Hit>> TheseClusters = Clusters[idx];
    for (size_t i = 0; i < TheseClusters.size(); i++){
      float clustT = 0;
      float clustCharge = 0;
      std::vector<recob::Hit> ThisCluster = TheseClusters[i];
      for (recob::Hit hit : ThisCluster){
        clustCharge += hit.Integral();
        clustT += hit.PeakTime();
      }
      ClT[idx].push_back(clustT/clustCharge);
      ClCharge[idx].push_back(clustCharge);
      ClNHits[idx].push_back(ThisCluster.size());
    } // End of loop over clusters
  } // End of loop over planes
}

std::vector<std::vector<std::vector<recob::Hit>>> LowEAna::MatchClusters(
  std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
  std::vector<std::vector<int>> &ClNHits, 
  std::vector<std::vector<float>> &ClT, 
  std::vector<std::vector<float>> &ClCharge, 
  bool HeavDebug)
{
  mf::LogDebug lmatch("match");
  std::string lmatchstr = "";
  FillClusterVariables(Clusters,ClNHits,ClT,ClCharge,HeavDebug);
  std::vector<std::vector<std::vector<recob::Hit>>> MatchedClusters = {{},{},{}};
  std::vector<std::vector<int>> MatchedClNHits = {{},{},{}};
  std::vector<std::vector<float>> MatchedClT = {{},{},{}};
  std::vector<std::vector<float>> MatchedClCharge = {{},{},{}};
  
  // --- Declare our variables to fill
  int MatchInd0Idx = -1, MatchInd1Idx = -1;
  double Ind0ClustdT = fClusterMatchTime, Ind1ClustdT = fClusterMatchTime;
  bool MatchInd0 = false, MatchInd1 = false;

  for (size_t ii = 0; ii < Clusters[2].size(); ii++){
    if (Clusters[2][ii].empty()){continue;}
    // Reset variables for next match
    lmatch << "\nMatching cluster " << ii;
    MatchInd0 = false;
    MatchInd1 = false;
    Ind0ClustdT = fClusterMatchTime; 
    Ind1ClustdT = fClusterMatchTime;
    
    if (Clusters[0].empty()){continue;}
    for (int jj = 0; jj < int(Clusters[0].size()); jj++){
      if (ClNHits[0][jj] < (1-fClusterMatchNHit)*ClNHits[2][ii] || ClNHits[0][jj] > (1+fClusterMatchNHit)*ClNHits[2][ii]){continue;} // Cut on number of hits of Ind0 cluster
      if (ClCharge[0][jj] < (1-fClusterMatchCharge)*ClCharge[2][ii] || ClCharge[0][jj] > (1+fClusterMatchCharge)*ClCharge[2][ii]){continue;} // Cut on charge of Ind0 cluster
      if (abs(ClT[2][ii] - ClT[0][jj]) < fClusterMatchTime && abs(fClusterInd0MatchTime - abs(ClT[2][ii] - ClT[0][jj])) < abs(fClusterInd0MatchTime - Ind0ClustdT)){
        Ind0ClustdT = abs(ClT[2][ii] - ClT[0][jj]);
        MatchInd0 = true;
        MatchInd0Idx = jj;
        lmatch << "Matched Ind0 cluster " << jj;
      }
    }
    if (Clusters[1].empty()){continue;}
    for (int zz = 0; zz < int(Clusters[1].size()); zz++){
      if (ClNHits[1][zz] < (1-fClusterMatchNHit)*ClNHits[2][ii] || ClNHits[1][zz] > (1+fClusterMatchNHit)*ClNHits[2][ii]){continue;} // Cut on number of hits of Ind1 cluster
      if (ClCharge[1][zz] < (1-fClusterMatchCharge)*ClCharge[2][ii] || ClCharge[1][zz] > (1+fClusterMatchCharge)*ClCharge[2][ii]){continue;} // Cut on charge of Ind1 cluster
      if (abs(ClT[2][ii] - ClT[1][zz]) < fClusterMatchTime && abs(fClusterInd1MatchTime - abs(ClT[2][ii] - ClT[1][zz])) < abs(fClusterInd1MatchTime - Ind1ClustdT)){
        Ind1ClustdT = abs(ClT[2][ii] - ClT[1][zz]);
        MatchInd1 = true;
        MatchInd1Idx = zz;
        lmatch << "Matched Ind1 cluster " << zz;
      }
    } // Loop over ind1 clusters
    // Fill matched clusters according to the matching criteria
    if (MatchInd0 && MatchInd1){
      MatchedClusters[0].push_back(Clusters[0][MatchInd0Idx]);
      MatchedClNHits[0].push_back(ClNHits[0][MatchInd0Idx]);
      MatchedClT[0].push_back(ClT[0][MatchInd0Idx]);
      MatchedClCharge[0].push_back(ClCharge[0][MatchInd0Idx]);

      MatchedClusters[1].push_back(Clusters[1][MatchInd1Idx]);
      MatchedClNHits[1].push_back(ClNHits[1][MatchInd1Idx]);
      MatchedClT[1].push_back(ClT[1][MatchInd1Idx]);
      MatchedClCharge[1].push_back(ClCharge[1][MatchInd1Idx]);

      MatchedClusters[2].push_back(Clusters[2][ii]);
      MatchedClNHits[2].push_back(ClNHits[2][ii]);
      MatchedClT[2].push_back(ClT[2][ii]);
      MatchedClCharge[2].push_back(ClCharge[2][ii]);
    }
    else if (MatchInd0 && !MatchInd1){
      MatchedClusters[0].push_back(Clusters[0][MatchInd0Idx]);
      MatchedClNHits[0].push_back(ClNHits[0][MatchInd0Idx]);
      MatchedClT[0].push_back(ClT[0][MatchInd0Idx]);
      MatchedClCharge[0].push_back(ClCharge[0][MatchInd0Idx]);

      MatchedClusters[2].push_back(Clusters[2][ii]);
      MatchedClNHits[2].push_back(ClNHits[2][ii]);
      MatchedClT[2].push_back(ClT[2][ii]);
      MatchedClCharge[2].push_back(ClCharge[2][ii]);
      // Fill missing cluster with empty vector
      MatchedClusters[1].push_back({});
      MatchedClNHits[1].push_back(0);
      MatchedClT[1].push_back(0);
      MatchedClCharge[1].push_back(0);
    }
    else if (!MatchInd0 && MatchInd1){
      MatchedClusters[1].push_back(Clusters[1][MatchInd1Idx]);
      MatchedClNHits[1].push_back(ClNHits[1][MatchInd1Idx]);
      MatchedClT[1].push_back(ClT[1][MatchInd1Idx]);
      MatchedClCharge[1].push_back(ClCharge[1][MatchInd1Idx]);

      MatchedClusters[2].push_back(Clusters[2][ii]);
      MatchedClNHits[2].push_back(ClNHits[2][ii]);
      MatchedClT[2].push_back(ClT[2][ii]);
      MatchedClCharge[2].push_back(ClCharge[2][ii]);
      // Fill missing cluster with empty vector
      MatchedClusters[0].push_back({});
      MatchedClNHits[0].push_back(0);
      MatchedClT[0].push_back(0);
      MatchedClCharge[0].push_back(0);
    }  
  }
  ClNHits = MatchedClNHits;
  ClT = MatchedClT;
  ClCharge = MatchedClCharge;
  
  lmatch << "\n# Matched Clusters:\t" << MatchedClusters[0].size();
  lmatch << lmatchstr;
  
  return MatchedClusters;
}

//......................................................
void LowEAna::FillMCInteractionTree( std::map< int, simb::MCParticle> &MCParticleList, std::vector<std::string> ProcessList, bool HeavDebug )
/*
Fill MCInteraction Tree with information about the main interaction in the event:
- MCParticleList is the list of MCParticles with a given generator label in the event
- ProcessList is the list of processes to be considered as main interactions
- HeavDebug is a boolean to turn on/off debugging statements
*/
{ 
  mf::LogInfo lheader("header");
  std::string lheaderstr = "";
  // Make a copy of MCParticleList to be used for finding Daughter info
  std::map< int, simb::MCParticle> MCParticleListCopy = MCParticleList;
  bool FoundInteraction;

  if (ProcessList.empty()){
    lheaderstr = PrintInColor(lheaderstr,"\n-> No processes in the list!",GetColor("red"));
    lheader << lheaderstr;
    return;  
  }
  
  for (size_t j = 0; j < ProcessList.size(); j++){
    FoundInteraction = false;    
    for ( std::map<int,simb::MCParticle>::iterator mainiter = MCParticleList.begin(); mainiter != MCParticleList.end(); mainiter++ ){
      if ( mainiter->second.Process() != ProcessList[j] && mainiter->second.EndProcess() != ProcessList[j]){continue;}
      lheaderstr = lheaderstr+"\nFound a main interaction "+mainiter->second.EndProcess();
      FoundInteraction = true;
      simb::MCParticle MCParticle = mainiter->second;
      Interaction =  MCParticle.EndProcess();
      PDG =          MCParticle.PdgCode();
      Energy =       MCParticle.E();
      Momentum =    {MCParticle.Px(),MCParticle.Py(),MCParticle.Pz()};
      StartVertex = {MCParticle.Vx(),MCParticle.Vy(),MCParticle.Vz()};
      EndVertex =   {MCParticle.EndX(),MCParticle.EndY(),MCParticle.EndZ()};
      
      std::vector<int> DaughterList = {};
      for (int i = 0; i < MCParticle.NumberDaughters(); i++){
        DaughterList.push_back(MCParticle.Daughter(i));
      }
      // Print nice output with all the main interaction info
      lheaderstr = PrintInColor(lheaderstr,"\nMain interacting particle for process "+mainiter->second.Process()+": ",GetColor("magenta"));
      lheaderstr = PrintInColor(lheaderstr,"\nPDG ->\t"         + str(PDG),GetColor("cyan"));
      lheaderstr = PrintInColor(lheaderstr,"\nEnergy ->\t"      + str(Energy),GetColor("cyan"));
      lheaderstr = PrintInColor(lheaderstr,"\nMomentum ->\t"    + str(Momentum[0]) + " " + str(Momentum[1]) + " " + str(Momentum[2]),GetColor("cyan"));
      lheaderstr = PrintInColor(lheaderstr,"\nStartVertex ->\t" + str(StartVertex[0]) + " " + str(StartVertex[1]) + " " + str(StartVertex[2]),GetColor("cyan"));
      lheaderstr = PrintInColor(lheaderstr,"\nEndVertex ->\t"   + str(EndVertex[0]) + " " + str(EndVertex[1]) + " " + str(EndVertex[2]),GetColor("cyan"));

      for ( std::map<int,simb::MCParticle>::iterator daughteriter = MCParticleListCopy.begin(); daughteriter != MCParticleListCopy.end(); daughteriter++ ){
        for (size_t i = 0; i < DaughterList.size(); i++){
          if (daughteriter->first == MCParticle.Daughter(i)){
            DaughterPDG.push_back(daughteriter->second.PdgCode());
            DaughterE.push_back(daughteriter->second.E());
            DaughterPx.push_back(daughteriter->second.Px());
            DaughterPy.push_back(daughteriter->second.Py());
            DaughterPz.push_back(daughteriter->second.Pz());
            DaughterStartVx.push_back(daughteriter->second.Vx());
            DaughterStartVy.push_back(daughteriter->second.Vy());
            DaughterStartVz.push_back(daughteriter->second.Vz());
            DaughterEndVx.push_back(daughteriter->second.EndX());
            DaughterEndVy.push_back(daughteriter->second.EndY());
            DaughterEndVz.push_back(daughteriter->second.EndZ());
          } // If the particle is a daughter of the main interaction
        } // Loop over all daughters
      } // Loop over all particles in the map
      fInteractionTree -> Fill();
    } // Loop over all particles in the map
    if (!FoundInteraction) lheaderstr = PrintInColor(lheaderstr,"\n-> No main interaction found for process "+ProcessList[j]+"!",GetColor("yellow"));
    else lheaderstr = PrintInColor(lheaderstr,"\n-> Filled MCINteraction Tree for process "+ProcessList[j]+"!",GetColor("green"));
  } // Loop over all processes in the list
  lheader << lheaderstr;
  return;
} // FillMCInteractionTree

void LowEAna::CalcAdjHits( std::vector< recob::Hit > MyVec,std::vector< std::vector<recob::Hit> >& Clusters,TH1I* MyHist, TH1F* ADCIntHist, 
  detinfo::DetectorClocksData clockData,
  bool HeavDebug ) 
/* 
Find adjacent hits in time and space:
- MyVec is the vector of hits to be clustered
- Clusters is the vector of clusters
- MyHist is the histogram to be filled with the number of hits in each cluster
- ADCIntHist is the histogram to be filled with the ADC integral of each cluster
- HeavDebug is a boolean to turn on/off debugging statements
*/
{
  const double TimeRange  = fClusterAlgoTime;
  const int    ChanRange  = fClusterAlgoAdjChannel;
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
        std::vector< sim::TrackIDE > ThisHitID = bt_serv->HitToTrackIDEs(clockData, hit);
        int MainTrID = 0;
        float TopEFrac = 0;
        
        for (size_t i=0; i < ThisHitID.size(); ++i){
          if (ThisHitID[i].energyFrac > TopEFrac){
            TopEFrac = ThisHitID[i].energyFrac;
            MainTrID = ThisHitID[i].trackID; 
          }
        }
        const simb::MCParticle *HitMother;
        int test = supress_stdout();
        HitMother = pi_serv->TrackIdToParticle_P(MainTrID);
        resume_stdout(test);
        if (HitMother != 0){
          std::cout << "Hit PDG" << HitMother->PdgCode() << std::endl;
          std::cout << "Hit E" << HitMother->E() << std::endl;
        }
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
long unsigned int LowEAna::WhichParType( int TrID )
/*
For a given TrackID, return the index of the vector of MCParticles that it belongs to.
Essestially, this function tells you which generator a given particle can be associated with.
*/
{
  for (long unsigned int i = 0; i < fLabels.size(); i++){
    if (InMyMap(TrID,Parts[i])) {return i+1;}
  }
  // If no match, then who knows???
  return 0;
}

//......................................................
void LowEAna::FillMyMaps( std::map< int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand )
/*
This function fills a map with the MCParticles from a given MCTruth
*/ 
{
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
bool LowEAna::InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap ){
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find( TrID );
  if (ParIt != ParMap.end()) {return true;}
  else return false;
}

//......................................................
// This function creates a terminal color printout
std::string LowEAna::PrintInColor( std::string InputString, std::string MyString, int Color ){
  std::string OutputString = InputString + "\033[" + str(Color) + "m" + MyString + "\033[0m";
  return OutputString;
}

// ......................................................
// This function returns an integer that corresponds to a given color name
int LowEAna::GetColor( std::string ColorName ){
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

std::string LowEAna::str( int i ) {std::stringstream ss;ss << i;return ss.str();}
std::string LowEAna::str( double i ) {std::stringstream ss;ss << i;return ss.str();}
std::string LowEAna::str( float i ) {std::stringstream ss;ss << i;return ss.str();}
std::string LowEAna::str( std::vector<int> i ) {std::stringstream ss;for (int j = 0; j < int(i.size()); j++){ss << i[j] << " ";}return ss.str();}
std::string LowEAna::str( std::vector<double> i ) {std::stringstream ss;for (int j = 0; j < int(i.size()); j++){ss << i[j] << " ";}return ss.str();}
std::string LowEAna::str( std::vector<float> i ) {std::stringstream ss;for (int j = 0; j < int(i.size()); j++){ss << i[j] << " ";}return ss.str();}

int LowEAna::supress_stdout() {
  std::fflush(stdout);

  int ret = dup(1);
  int nullfd = open("/dev/null", O_WRONLY);
  // check nullfd for error omitted
  dup2(nullfd, 1);
  close(nullfd);

  return ret;
}

void LowEAna::resume_stdout(int fd) {
  std::fflush(stdout);
  dup2(fd, 1);
  close(fd);
}

//......................................................

// Function that returns the average of a vector
double LowEAna::Average(std::vector<double> &v)
{
  double sum = 0;
  for (size_t i = 0; i < v.size(); i++){sum += v[i];}
  return sum/v.size();
}
double LowEAna::Average(std::vector<int> &v)
{
  double sum = 0;
  for (size_t i = 0; i < v.size(); i++){sum += v[i];}
  return sum/v.size();
}
double LowEAna::Average(std::vector<float> &v)
{
  double sum = 0;
  for (size_t i = 0; i < v.size(); i++){sum += v[i];}
  return sum/v.size();
}

DEFINE_ART_MODULE(LowEAna)
