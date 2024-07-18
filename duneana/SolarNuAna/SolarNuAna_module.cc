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

// ROOT includes
#include <TApplication.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <fcntl.h>
#include "TRandom.h"
#include "TVector3.h"

// Framework includes (not all might be necessary)
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larcorealg/Geometry/GeometryCore.h"
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
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "AdjHitsUtils.h"
#include "AdjOpHitsUtils.h"

namespace solar
{
  class SolarNuAna : public art::EDAnalyzer
  {
  public:
    // --- Standard constructor and destructor for an ART module.
    explicit SolarNuAna(fhicl::ParameterSet const &p);
    SolarNuAna(SolarNuAna const &) = delete;
    SolarNuAna(SolarNuAna &&) = delete;
    SolarNuAna &operator=(SolarNuAna const &) = delete;
    SolarNuAna &operator=(SolarNuAna &&) = delete;
    void reconfigure(fhicl::ParameterSet const &p);
    void beginJob() override;
    void analyze(art::Event const &evt) override;

  private:
    // --- Some of our own functions.
    void ResetVariables();
    long unsigned int WhichParType(int TrID);
    bool InMyMap(int TrID, std::map<int, float> TrackIDMap);
    bool InMyMap(int TrID, std::map<int, simb::MCParticle> ParMap);
    void FillMyMaps(std::map<int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle<std::vector<simb::MCTruth>> Hand);
    void PrintInColor(std::string MyString, int MyColor, std::string Type = "Info");
    int GetColor(std::string MyString);
    std::string str(int MyInt);
    std::string str(float MyFloat, int MyPrecision = 2);
    std::string str(double MyDouble, int MyPrecision = 2);
    std::string str(std::vector<int> MyVec);
    std::string str(std::vector<float> MyVec, int MyPrecision = 2);
    std::string str(std::vector<double> MyVec, int MyPrecision = 2);
    int supress_stdout();
    void resume_stdout(int fd);

    // --- Our fcl parameter labels for the modules that made the data products
    std::string fRawDigitLabel, fHitLabel, fTrackLabel, fOpHitLabel, fOpFlashLabel, fGEANTLabel;

    // --- Input settings imported from the fcl
    std::string fGeometry;
    int fDetectorSizeX, fDetectorSizeY, fDetectorSizeZ,fDetectorDriftTime, fClusterAlgoAdjChannel, fClusterInd0MatchTime, fClusterInd1MatchTime, fClusterPreselectionNHit;
    float fClusterMatchTime, fAdjClusterRad, fMinClusterCharge, fClusterMatchCharge, fAdjOpFlashY, fAdjOpFlashZ, fAdjOpFlashTime, fAdjOpFlashMaxPERatioCut, fAdjOpFlashMinPECut, fClusterMatchNHit, fClusterAlgoTime;
    std::vector<std::string> fLabels;
    float fOpFlashAlgoTime, fOpFlashAlgoRad, fOpFlashAlgoPE, fOpFlashAlgoTriggerPE;
    bool fClusterPreselectionTrack, fGenerateAdjOpFlash, fSaveMarleyEDep, fSaveSignalOpHits, fOpFlashAlgoCentroid, fOpFlashAlgoDebug;

    // --- Our TTrees, and its associated variables.
    TTree *fConfigTree;
    TTree *fMCTruthTree;
    TTree *fSolarNuAnaTree;
    std::string TNuInteraction;
    int Event, Flag, MNHit, MGen, MTPC, MInd0TPC, MInd1TPC, MInd0NHits, MInd1NHits, MMainID, MMainPDG, MMainParentPDG, TrackNum, OpHitNum, OpFlashNum, MTrackNPoints;
    float TNuE, TNuX, TNuY, TNuZ, MTime, MCharge, MMaxCharge, MInd0Charge, MInd1Charge, MInd0MaxCharge, MInd1MaxCharge;
    float MInd0dT, MInd1dT, MInd0RecoY, MInd1RecoY, MRecY, MRecZ, MPur, MMainE, MMainP, MMainK, MMainT, MMainParentE, MMainParentP, MMainParentK, MMainParentT, MTrackChi2;
    std::vector<int> MAdjClGen, MAdjClMainID, TPart, MarleyPDGList, MarleyPDGDepList, MarleyIDList, MarleyMotherList, MarleyIDDepList, MAdjClMainPDG, HitNum, ClusterNum, MarleyElectronDepList;
    std::vector<float> MarleyEDepList, MarleyXDepList, MarleyYDepList, MarleyZDepList;
    std::vector<float> SOpHitPur, SOpHitPE, SOpHitX, SOpHitY, SOpHitZ, SOpHitT, SOpHitChannel, SOpHitFlashID;
    std::vector<float> MAdjClTime, MAdjClCharge, MAdjClInd0Charge, MAdjClInd1Charge, MAdjClMaxCharge, MAdjClInd0MaxCharge, MAdjClInd1MaxCharge;
    std::vector<float> MAdjClNHit, MAdjClInd0NHit, MAdjClInd1NHit, MAdjClRecoY, MAdjClRecoZ, MAdjClR, MAdjClPur, MAdjClMainE, MAdjClMainK;  
    std::vector<float> MAdjClMainX, MAdjClMainY, MAdjClMainZ, MAdjClEndX, MAdjClEndY, MAdjClEndZ, MMarleyFrac, MGenFrac;
    std::vector<float> MAdjFlashTime, MAdjFlashPE, MAdjFlashNHit, MAdjFlashMaxPE, MAdjFlashRecoX, MAdjFlashRecoY, MAdjFlashRecoZ, MAdjFlashR, MAdjFlashPur;
    std::vector<float> MarleyEList, MarleyPList, MarleyKList, MarleyTList, MarleyEndXList, MarleyEndYList, MarleyEndZList, MarleyMaxEDepList, MarleyMaxEDepXList, MarleyMaxEDepYList, MarleyMaxEDepZList;
    std::vector<double> MTrackStart, MTrackEnd;
    std::vector<double> MMainVertex, MEndVertex, MMainParentVertex;
    std::vector<std::map<int, simb::MCParticle>> Parts = {};
    bool MPrimary;

    // --- OpFlash Variables
    std::vector<float> OpFlashMarlPur, OpFlashPE, OpFlashMaxPE, OpFlashX, OpFlashY, OpFlashZ, OpFlashT, OpFlashDeltaT, OpFlashNHit;

    // --- Histograms to fill about collection plane hits
    float x_f;
    TH2F *hXTruth;
    TH2F *hYTruth;
    TH2F *hZTruth;
    TH1I *hAdjHits;
    TH1F *hAdjHitsADCInt;
    TH2F *hDriftTime;

    // --- Declare our services
    geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    std::unique_ptr<solar::AdjHitsUtils> adjhits;
    std::unique_ptr<solar::AdjOpHitsUtils> adjophits;

  };
  #endif

  //......................................................
  SolarNuAna::SolarNuAna(fhicl::ParameterSet const &p) 
    : EDAnalyzer(p),
    adjhits(new solar::AdjHitsUtils(p)),
    adjophits(new solar::AdjOpHitsUtils(p))
    {this->reconfigure(p);}

  //......................................................
  void SolarNuAna::reconfigure(fhicl::ParameterSet const &p)
  {
    fLabels                   = p.get<std::vector<std::string>>("ParticleLabelVector");
    fRawDigitLabel            = p.get<std::string>("RawDigitLabel");
    fHitLabel                 = p.get<std::string>("HitLabel");
    fOpFlashLabel             = p.get<std::string>("OpFlashLabel");
    fOpHitLabel               = p.get<std::string>("OpHitLabel");
    fTrackLabel               = p.get<std::string>("TrackLabel");
    fGEANTLabel               = p.get<std::string>("GEANT4Label");
    fGeometry                 = p.get<std::string>("Geometry");
    fDetectorSizeX            = p.get<int>("DetectorSizeX");
    fDetectorSizeY            = p.get<int>("DetectorSizeY");
    fDetectorSizeZ            = p.get<int>("DetectorSizeZ");
    fDetectorDriftTime        = p.get<int>("DetectorDriftTime");
    fClusterAlgoTime          = p.get<float>("ClusterAlgoTime");
    fClusterAlgoAdjChannel    = p.get<int>("ClusterAlgoAdjChannel");
    fClusterMatchNHit         = p.get<float>("ClusterMatchNHit");
    fClusterMatchCharge       = p.get<float>("ClusterMatchCharge");
    fClusterMatchTime         = p.get<float>("ClusterMatchTime");
    fClusterInd0MatchTime     = p.get<float>("ClusterInd0MatchTime");
    fClusterInd1MatchTime     = p.get<float>("ClusterInd1MatchTime");
    fClusterPreselectionNHit  = p.get<int>("ClusterPreselectionNHit");
    fClusterPreselectionTrack = p.get<bool>("ClusterPreselectionTrack");
    fAdjClusterRad            = p.get<float>("AdjClusterRad");
    fMinClusterCharge         = p.get<float>("MinClusterCharge");
    fGenerateAdjOpFlash       = p.get<bool>("GenerateAdjOpFlash");
    fOpFlashAlgoTime          = p.get<double>("OpFlashAlgoTime");
    fOpFlashAlgoRad           = p.get<double>("OpFlashAlgoRad");
    fOpFlashAlgoPE            = p.get<float>("OpFlashAlgoPE");
    fOpFlashAlgoTriggerPE     = p.get<float>("OpFlashAlgoTriggerPE");
    fOpFlashAlgoCentroid      = p.get<bool>("OpFlashAlgoCentroid");
    fOpFlashAlgoDebug         = p.get<bool>("OpFlashAlgoDebug");
    fAdjOpFlashTime           = p.get<float>("AdjOpFlashTime");
    fAdjOpFlashY              = p.get<float>("AdjOpFlashY");
    fAdjOpFlashZ              = p.get<float>("AdjOpFlashZ");
    fAdjOpFlashMaxPERatioCut  = p.get<float>("AdjOpFlashMaxPERatioCut");
    fAdjOpFlashMinPECut       = p.get<float>("AdjOpFlashMinPECut");
    fSaveMarleyEDep           = p.get<bool>("SaveMarleyEDep");
    fSaveSignalOpHits         = p.get<bool>("SaveSignalOpHits");
  } // Reconfigure

  //......................................................
  void SolarNuAna::beginJob()
  {
    // --- Make our handle to the TFileService
    art::ServiceHandle<art::TFileService> tfs;
    fConfigTree     = tfs->make<TTree>("ConfigTree", "Config Tree");
    fMCTruthTree    = tfs->make<TTree>("MCTruthTree", "MC Truth Tree");
    fSolarNuAnaTree = tfs->make<TTree>("SolarNuAnaTree", "Solar Ana Tree");

    // Larsoft Config info.
    fConfigTree->Branch("RawDigitLabel", &fRawDigitLabel);
    fConfigTree->Branch("HitLabel", &fHitLabel);
    fConfigTree->Branch("OpFlashLabel", &fOpFlashLabel);
    fConfigTree->Branch("OpHitLabel", &fOpHitLabel);
    fConfigTree->Branch("TrackLabel", &fTrackLabel);
    fConfigTree->Branch("GEANT4Label", &fGEANTLabel);
    fConfigTree->Branch("Geometry", &fGeometry);
    fConfigTree->Branch("DetectorSizeX", &fDetectorSizeX);
    fConfigTree->Branch("DetectorSizeY", &fDetectorSizeY);
    fConfigTree->Branch("DetectorSizeZ", &fDetectorSizeZ);
    fConfigTree->Branch("DetectorDriftTime", &fDetectorDriftTime);
    fConfigTree->Branch("ClusterAlgoTime", &fClusterAlgoTime);
    fConfigTree->Branch("ClusterAlgoAdjChannel", &fClusterAlgoAdjChannel);
    fConfigTree->Branch("ClusterMatchNHit", &fClusterMatchNHit);
    fConfigTree->Branch("ClusterMatchCharge", &fClusterMatchCharge);
    fConfigTree->Branch("ClusterMatchTime", &fClusterMatchTime);
    fConfigTree->Branch("ClusterInd0MatchTime", &fClusterInd0MatchTime);
    fConfigTree->Branch("ClusterInd1MatchTime", &fClusterInd1MatchTime);
    fConfigTree->Branch("ClusterPreselectionNHit", &fClusterPreselectionNHit);
    fConfigTree->Branch("ClusterPreselectionTrack", &fClusterPreselectionTrack);
    fConfigTree->Branch("AdjClusterRad", &fAdjClusterRad);
    fConfigTree->Branch("MinClusterCharge", &fMinClusterCharge);
    fConfigTree->Branch("GenerateAdjOpFlash", &fGenerateAdjOpFlash);
    fConfigTree->Branch("OpFlashAlgoTime", &fOpFlashAlgoTime);
    fConfigTree->Branch("OpFlashAlgoRad", &fOpFlashAlgoRad);
    fConfigTree->Branch("OpFlashAlgoPE", &fOpFlashAlgoPE);
    fConfigTree->Branch("OpFlashAlgoTriggerPE", &fOpFlashAlgoTriggerPE);
    fConfigTree->Branch("OpFlashAlgoCentroid", &fOpFlashAlgoCentroid);
    fConfigTree->Branch("OpFlashAlgoDebug", &fOpFlashAlgoDebug);
    fConfigTree->Branch("AdjOpFlashTime", &fAdjOpFlashTime);
    fConfigTree->Branch("AdjOpFlashY", &fAdjOpFlashY);
    fConfigTree->Branch("AdjOpFlashZ", &fAdjOpFlashZ);
    fConfigTree->Branch("AdjOpFlashMaxPERatioCut", &fAdjOpFlashMaxPERatioCut);
    fConfigTree->Branch("AdjOpFlashMinPECut", &fAdjOpFlashMinPECut);
    fConfigTree->Branch("SaveMarleyEDep", &fSaveMarleyEDep);
    fConfigTree->Branch("SaveSignalOpHits", &fSaveSignalOpHits);

    // MC Truth info.
    fMCTruthTree->Branch("Event", &Event, "Event/I");                // Event number
    fMCTruthTree->Branch("Flag", &Flag, "Flag/I");                   // Flag used to match truth with reco tree entries
    fMCTruthTree->Branch("TruthPart", &TPart);                       // Number particles per generator
    fMCTruthTree->Branch("Interaction", &TNuInteraction);            // True neutrino interaction process
    fMCTruthTree->Branch("TNuE", &TNuE, "TruthNuE/F");               // True neutrino energy [MeV]
    fMCTruthTree->Branch("TNuX", &TNuX, "TruthNuX/F");               // True neutrino X [cm]
    fMCTruthTree->Branch("TNuY", &TNuY, "TruthNuY/F");               // True neutrino Y [cm]
    fMCTruthTree->Branch("TNuZ", &TNuZ, "TruthNuZ/F");               // True neutrino Z [cm]
    fMCTruthTree->Branch("TMarleyPDG", &MarleyPDGList);              // PDG of marley marticles
    fMCTruthTree->Branch("TMarleyE", &MarleyEList);                  // Energy of marley particles [MeV]
    fMCTruthTree->Branch("TMarleyP", &MarleyPList);                  // Momentum of marley particles [MeV]
    fMCTruthTree->Branch("TMarleyK", &MarleyKList);                  // Kinetik Energy of marley particles [MeV]
    fMCTruthTree->Branch("TMarleyT", &MarleyTList);                  // Time of marley particles [ticks]
    fMCTruthTree->Branch("TMarleyEndX", &MarleyEndXList);            // X of marley particles [cm]
    fMCTruthTree->Branch("TMarleyEndY", &MarleyEndYList);            // Y of marley particles [cm]
    fMCTruthTree->Branch("TMarleyEndZ", &MarleyEndZList);            // Z of marley particles [cm]
    fMCTruthTree->Branch("TMarleyMaxEDep", &MarleyMaxEDepList);      // Energy of marley particles [MeV]
    fMCTruthTree->Branch("TMarleyX", &MarleyMaxEDepXList);           // X of marley particles [cm]
    fMCTruthTree->Branch("TMarleyY", &MarleyMaxEDepYList);           // Y of marley particles [cm]
    fMCTruthTree->Branch("TMarleyZ", &MarleyMaxEDepZList);           // Z of marley particles [cm]
    fMCTruthTree->Branch("TMarleyID", &MarleyIDList);                // TrackID of marley particles
    fMCTruthTree->Branch("TMarleyMother", &MarleyMotherList);        // TrackID of marley mother
    fMCTruthTree->Branch("TrackNum", &TrackNum, "TrackNum/I");       // Number of PMTracks
    fMCTruthTree->Branch("OpHitNum", &OpHitNum, "OpHitNum/I");       // Number of OpHits
    fMCTruthTree->Branch("OpFlashNum", &OpFlashNum, "OpFlashNum/I"); // Number of OpFlashes
    fMCTruthTree->Branch("HitNum", &HitNum);                         // Number of hits in each TPC plane
    fMCTruthTree->Branch("ClusterNum", &ClusterNum);                 // Number of clusters in each TPC plane
    if (fSaveMarleyEDep)
    { 
      fMCTruthTree->Branch("TMarleyPDGDepList", &MarleyPDGDepList);  // PDG for Energy deposited of marley particles
      fMCTruthTree->Branch("TMarleyEDepList", &MarleyEDepList);      // Energy deposited of marley particles [MeV]
      fMCTruthTree->Branch("TMarleyXDepList", &MarleyXDepList);      // X deposited of marley particles [cm]
      fMCTruthTree->Branch("TMarleyYDepList", &MarleyYDepList);      // Y deposited of marley particles [cm]
      fMCTruthTree->Branch("TMarleyZDepList", &MarleyZDepList);      // Z deposited of marley particles [cm]
      fMCTruthTree->Branch("TMarleyIDDepList", &MarleyIDDepList);    // ParentID of marley particles
      fMCTruthTree->Branch("TMarleyElectronDepList", &MarleyElectronDepList); // Number of electrons in the marley particles
    }
    if (fSaveSignalOpHits)
    {
      fMCTruthTree->Branch("OpHitPur", &SOpHitPur);         // OpHit Purity
      fMCTruthTree->Branch("OpHitPE", &SOpHitPE);           // OpHit PE
      fMCTruthTree->Branch("OpHitX", &SOpHitX);             // OpHit X
      fMCTruthTree->Branch("OpHitY", &SOpHitY);             // OpHit Y
      fMCTruthTree->Branch("OpHitZ", &SOpHitZ);             // OpHit Z
      fMCTruthTree->Branch("OpHitT", &SOpHitT);             // OpHit T
      fMCTruthTree->Branch("OpHitChannel", &SOpHitChannel); // OpHit Channel
      fMCTruthTree->Branch("OpHitFlashID", &SOpHitFlashID); // OpHit Area
    }

    // Repeated Truth info.
    fSolarNuAnaTree->Branch("Event", &Event, "Event/I");             // Event number
    fSolarNuAnaTree->Branch("Flag", &Flag, "Flag/I");                // Flag used to match truth with reco tree entries
    fSolarNuAnaTree->Branch("TruthPart", &TPart);                    // Number particles per generator
    fSolarNuAnaTree->Branch("Interaction", &TNuInteraction);         // True neutrino interaction process
    fSolarNuAnaTree->Branch("TNuE", &TNuE, "TruthNuE/F");            // True neutrino energy
    fSolarNuAnaTree->Branch("TNuX", &TNuX, "TruthNuX/F");            // True neutrino X
    fSolarNuAnaTree->Branch("TNuY", &TNuY, "TruthNuY/F");            // True neutrino Y
    fSolarNuAnaTree->Branch("TNuZ", &TNuZ, "TruthNuZ/F");            // True neutrino Z
    fSolarNuAnaTree->Branch("TMarleyPDG", &MarleyPDGList);           // PDG of marley particles
    fSolarNuAnaTree->Branch("TMarleyE", &MarleyEList);               // Energy of marley particles
    fSolarNuAnaTree->Branch("TMarleyK", &MarleyKList);               // Kinetik Energy of marley particles
    fSolarNuAnaTree->Branch("TMarleyEndX", &MarleyEndXList);         // X of marley particles
    fSolarNuAnaTree->Branch("TMarleyEndY", &MarleyEndYList);         // Y of marley particles
    fSolarNuAnaTree->Branch("TMarleyEndZ", &MarleyEndZList);         // Z of marley particles
    fSolarNuAnaTree->Branch("TMarleyMaxEDep", &MarleyMaxEDepList);   // Max Energy Deposition of marley particles
    fSolarNuAnaTree->Branch("TMarleyX", &MarleyMaxEDepXList);        // Max Energy Deposition X of marley particles
    fSolarNuAnaTree->Branch("TMarleyY", &MarleyMaxEDepYList);        // Max Energy Deposition Y of marley particles
    fSolarNuAnaTree->Branch("TMarleyZ", &MarleyMaxEDepZList);        // Max Energy Deposition Z of marley particles
    fSolarNuAnaTree->Branch("TMarleyID", &MarleyIDList);             // TrackID of marley particles")
    fSolarNuAnaTree->Branch("TMarleyMother", &MarleyMotherList);     // TrackID of marley particles")


    // Main Cluster info.
    fSolarNuAnaTree->Branch("Generator", &MGen, "Generator/I");                   // Main cluster generator idx
    fSolarNuAnaTree->Branch("Purity", &MPur, "Purity/F");                         // Main cluster reco purity
    fSolarNuAnaTree->Branch("TPC", &MTPC, "ColTPC/I");                            // Main cluster TPC
    fSolarNuAnaTree->Branch("Primary", &MPrimary);                                // Cluster hasn't any adjcl with AdjClCharge > MCharge (bool)
    fSolarNuAnaTree->Branch("Time", &MTime, "ColTime/F");                         // Main cluster time [ticks]
    fSolarNuAnaTree->Branch("NHits", &MNHit, "ColNHits/I");                       // Main cluster #hits
    fSolarNuAnaTree->Branch("Charge", &MCharge, "ColCharge/F");                   // Main cluster charge [ADC*ticks]
    fSolarNuAnaTree->Branch("MaxCharge", &MMaxCharge, "ColCharge/F");             // Main cluster's max TPCHit-charge [ADC*ticks]
    fSolarNuAnaTree->Branch("RecoZ", &MRecZ, "RecoZ/F");                          // Main cluster reco Z [cm]
    fSolarNuAnaTree->Branch("Ind0TPC", &MInd0TPC, "Ind0TPC/I");                   // Main cluster ind0 TPC
    fSolarNuAnaTree->Branch("Ind1TPC", &MInd1TPC, "Ind1TPC/I");                   // Main cluster ind1 TPC
    fSolarNuAnaTree->Branch("Ind0dT", &MInd0dT, "Ind0dT/F");                      // Main cluster ind0 dT [Ticks]
    fSolarNuAnaTree->Branch("Ind1dT", &MInd1dT, "Ind1dT/F");                      // Main cluster ind1 dT [Ticks]
    fSolarNuAnaTree->Branch("Ind0NHits", &MInd0NHits, "Ind0NHits/I");             // Main cluster ind0 Hits
    fSolarNuAnaTree->Branch("Ind1NHits", &MInd1NHits, "Ind1NHits/I");             // Main cluster ind1 Hits
    fSolarNuAnaTree->Branch("Ind0Charge", &MInd0Charge, "Ind0Charge/F");          // Main cluster ind0 MaxHit
    fSolarNuAnaTree->Branch("Ind1Charge", &MInd1Charge, "Ind1Charge/F");          // Main cluster ind1 MaxHit
    fSolarNuAnaTree->Branch("Ind0MaxCharge", &MInd0MaxCharge, "Ind0MaxCharge/F"); // Main cluster ind0 MaxHit
    fSolarNuAnaTree->Branch("Ind1MaxCharge", &MInd1MaxCharge, "Ind1MaxCharge/F"); // Main cluster ind1 MaxHit
    fSolarNuAnaTree->Branch("RecoY", &MRecY, "RecoY/F");                          // Main cluster ind0 reco Y [cm]
    fSolarNuAnaTree->Branch("Ind0RecoY", &MInd0RecoY, "Ind0RecoY/F");             // Main cluster ind0 reco Y [cm]
    fSolarNuAnaTree->Branch("Ind1RecoY", &MInd1RecoY, "Ind1RecoY/F");             // Main cluster ind1 reco Y [cm]
    fSolarNuAnaTree->Branch("MainID", &MMainID, "MainID/I");                      // Main cluster main track ID
    fSolarNuAnaTree->Branch("MainT", &MMainT, "MainT/I");                         // Main cluster main time [ticks]
    fSolarNuAnaTree->Branch("MainE", &MMainE, "MainE/F");                         // Main cluster main energy [MeV]
    fSolarNuAnaTree->Branch("MainP", &MMainP, "MainP/F");                         // Main cluster main momentum [MeV]
    fSolarNuAnaTree->Branch("MainK", &MMainK, "MainK/F");                         // Main cluster main kinetic energy [MeV]
    fSolarNuAnaTree->Branch("MainT", &MMainT, "MainT/F");                         // Main cluster main Time [ticks]
    fSolarNuAnaTree->Branch("MainPDG", &MMainPDG, "MainPDG/I");                   // Main cluster main pdg
    fSolarNuAnaTree->Branch("MainParentPDG", &MMainParentPDG, "MainParentPDG/I"); // Main cluster main pdg
    fSolarNuAnaTree->Branch("MainParentE", &MMainParentE, "MainParentE/F");       // Main cluster main parent energy [MeV]
    fSolarNuAnaTree->Branch("MainParentP", &MMainParentP, "MainParentP/F");       // Main cluster main parent momentum [MeV]
    fSolarNuAnaTree->Branch("MainParentK", &MMainParentK, "MainParentK/F");       // Main cluster main parent kinetic energy [MeV]
    fSolarNuAnaTree->Branch("MainParentT", &MMainParentT, "MainParentT/F");       // Main cluster main parent Time [ticks]
    fSolarNuAnaTree->Branch("MainVertex", &MMainVertex);                          // Main cluster main particle vertex [cm]
    fSolarNuAnaTree->Branch("EndVertex", &MEndVertex);                            // Main cluster end particle vertex [cm]
    fSolarNuAnaTree->Branch("MainParentVertex", &MMainParentVertex);              // Main cluster parent particle vertex [cm]
    fSolarNuAnaTree->Branch("GenFrac", &MGenFrac);                                // Main cluster reco purity complete
    fSolarNuAnaTree->Branch("MarleyFrac", &MMarleyFrac);                          // Main cluster particle contribution (electron, gamma, neutron)

    // Track info.
    fSolarNuAnaTree->Branch("MTrackNPoints", &MTrackNPoints, "TrackNPoints/I");   // Track #points
    fSolarNuAnaTree->Branch("MTrackStart", &MTrackStart);                         // Track start point
    fSolarNuAnaTree->Branch("MTrackEnd", &MTrackEnd);                             // Track end point
    fSolarNuAnaTree->Branch("MTrackChi2", &MTrackChi2);                           // Track chi2

    // Adj. Cluster info.
    fSolarNuAnaTree->Branch("AdjClGen", &MAdjClGen);                     // Adj. clusters' generator idx
    fSolarNuAnaTree->Branch("AdjClNHit", &MAdjClNHit);                   // Adj. clusters' #hits
    fSolarNuAnaTree->Branch("AdjClInd0NHit", &MAdjClInd0NHit);           // Adj. clusters' #hits
    fSolarNuAnaTree->Branch("AdjClInd1NHit", &MAdjClInd1NHit);           // Adj. clusters' #hits
    fSolarNuAnaTree->Branch("AdjClTime", &MAdjClTime);                   // Adj. clusters' time [ticks]
    fSolarNuAnaTree->Branch("AdjClCharge", &MAdjClCharge);               // Adj. clusters' charge [ADC*ticks]
    fSolarNuAnaTree->Branch("AdjClInd0Charge", &MAdjClInd0Charge);       // Adj. clusters' charge [ADC*ticks]
    fSolarNuAnaTree->Branch("AdjClInd1Charge", &MAdjClInd1Charge);       // Adj. clusters' charge [ADC*ticks]
    fSolarNuAnaTree->Branch("AdjClMaxCharge", &MAdjClMaxCharge);         // Adj. clusters' charge [ADC*ticks]
    fSolarNuAnaTree->Branch("AdjClInd0MaxCharge", &MAdjClInd0MaxCharge); // Adj. clusters' charge [ADC*ticks]
    fSolarNuAnaTree->Branch("AdjClInd1MaxCharge", &MAdjClInd1MaxCharge); // Adj. clusters' charge [ADC*ticks]
    fSolarNuAnaTree->Branch("AdjClRecoY", &MAdjClRecoY);                 // Adj. clusters' reco Y [cm]
    fSolarNuAnaTree->Branch("AdjClRecoZ", &MAdjClRecoZ);                 // Adj. clusters' reco Z [cm]
    fSolarNuAnaTree->Branch("AdjClR", &MAdjClR);                         // Adj. clusters' distance [cm]
    fSolarNuAnaTree->Branch("AdjClPur", &MAdjClPur);                     // Adj. clusters' purity
    fSolarNuAnaTree->Branch("AdjClMainID", &MAdjClMainID);               // Adj. clusters' main track ID
    fSolarNuAnaTree->Branch("AdjClMainPDG", &MAdjClMainPDG);             // Adj. clusters' main PDG
    fSolarNuAnaTree->Branch("AdjClMainE", &MAdjClMainE);                 // Adj. clusters' main energy [GeV]
    fSolarNuAnaTree->Branch("AdjClMainX", &MAdjClMainX);                 // Adj. clusters' main X [cm]
    fSolarNuAnaTree->Branch("AdjClMainY", &MAdjClMainY);                 // Adj. clusters' main Y [cm]
    fSolarNuAnaTree->Branch("AdjClMainZ", &MAdjClMainZ);                 // Adj. clusters' main Z [cm]
    fSolarNuAnaTree->Branch("AdjClEndX", &MAdjClEndX);                   // Adj. clusters' end X [cm]
    fSolarNuAnaTree->Branch("AdjClEndY", &MAdjClEndY);                   // Adj. clusters' end Y [cm]
    fSolarNuAnaTree->Branch("AdjClEndZ", &MAdjClEndZ);                   // Adj. clusters' end Z [cm]

    // Adj. Flash info.
    fSolarNuAnaTree->Branch("AdjOpFlashTime", &MAdjFlashTime);   // Adj. flash' time [ticks]
    fSolarNuAnaTree->Branch("AdjOpFlashPE", &MAdjFlashPE);       // Adj. flash' tot #PE [ADC*ticks]
    fSolarNuAnaTree->Branch("AdjOpFlashNHit", &MAdjFlashNHit);   // Adj. flash' #hits
    fSolarNuAnaTree->Branch("AdjOpFlashMaxPE", &MAdjFlashMaxPE); // Adj. flash' max #PE [ADC*ticks]
    fSolarNuAnaTree->Branch("AdjOpFlashRecoY", &MAdjFlashRecoX); // Adj. flash' reco X [cm]
    fSolarNuAnaTree->Branch("AdjOpFlashRecoY", &MAdjFlashRecoY); // Adj. flash' reco Y [cm]
    fSolarNuAnaTree->Branch("AdjOpFlashRecoZ", &MAdjFlashRecoZ); // Adj. flash' reco Z [cm]
    fSolarNuAnaTree->Branch("AdjOpFlashPur", &MAdjFlashPur);     // Adj. flash' purity
    fSolarNuAnaTree->Branch("AdjOpFlashR", &MAdjFlashR);         // Adj. flash' reco distance [cm]

    fConfigTree->AddFriend(fSolarNuAnaTree);
    fMCTruthTree->AddFriend(fSolarNuAnaTree);
    fConfigTree->Fill();

    // --- Our Histograms...
    hDriftTime      = tfs->make<TH2F>("hDriftTime",     "hDriftTime", 100, -400., 400., 100, 0., 10000.);
    hXTruth         = tfs->make<TH2F>("hXTruth",        "Missmatch in X distance; Distance [cm]; True X position [cm]", 100, -600, 600, 100, -600, 600);
    hYTruth         = tfs->make<TH2F>("hYTruth",        "Missmatch in Y distance; Distance [cm]; True Y position [cm]", 100, -600, 600, 100, -600, 600);
    hZTruth         = tfs->make<TH2F>("hZTruth",        "Missmatch in Z distance; Distance [cm]; True Z position [cm]", 100, -600, 600, 100, 0, 1600);
    hAdjHits        = tfs->make<TH1I>("hAdjHits",       "Number of adjacent collection plane hits; Number of adjacent collection plane hits; Number of events", 21, -0.5, 20.5);
    hAdjHitsADCInt  = tfs->make<TH1F>("hAdjHitsADCInt", "Total summed ADC Integrals for clusters; Total summed ADC Integrals for clusters; Number of events", 1000, 0, 10000);
  } // BeginJob

  //......................................................
  void SolarNuAna::analyze(art::Event const &evt)
  {
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------- Prepare everything for new event ----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    std::vector<std::set<int>> trackids = {};
    std::map<int, simb::MCParticle> ThisGeneratorParts;
    std::vector<recob::Hit> ColHits0, ColHits1, ColHits2, ColHits3;
    std::vector<std::vector<recob::Hit>> ColHits = {ColHits0, ColHits1, ColHits2, ColHits3};
    std::vector<std::vector<recob::Hit>> Clusters0, Clusters1, Clusters2, Clusters3;

    // --- We want to reset all of our previous run and TTree variables ---
    ResetVariables();
    ThisGeneratorParts.clear();
    Event = evt.event();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    Flag = rand() % 10000000000;
    std::string sHead = "";
    sHead = sHead + "\nTPC Frequency in [MHz]: " + str(clockData.TPCClock().Frequency());
    sHead = sHead + "\nTPC Tick in [us]: " + str(clockData.TPCClock().TickPeriod());
    sHead = sHead + "\nEvent Flag: " + str(Flag);
    sHead = sHead + "\nSuccesfull reset of variables for evt " + str(Event);
    sHead = sHead + "\n#########################################";
    PrintInColor(sHead, GetColor("magenta"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Create maps for ID tracking -----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // --- Fill MC Truth IDs to tracking vectors. Get a list of all of my particles in one chunk. ---
    const sim::ParticleList &PartList = pi_serv->ParticleList();
    std::string sMcTruth = "";
    sMcTruth = sMcTruth + "\nThere are a total of " + str(int(PartList.size())) + " Particles in the event\n";

    // Loop over all signal+bkg handles and collect track IDs
    for (size_t i = 0; i < fLabels.size(); i++)
    {
      Parts.push_back(ThisGeneratorParts); // For each label insert empty list

      art::Handle<std::vector<simb::MCTruth>> ThisHandle;
      evt.getByLabel(fLabels[i], ThisHandle);

      if (ThisHandle)
      {
        auto ThisValidHanlde = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[i]); // Get generator handles
        art::FindManyP<simb::MCParticle> Assn(ThisValidHanlde, evt, fGEANTLabel);          // Assign labels to MCPArticles
        FillMyMaps(Parts[i], Assn, ThisValidHanlde);                                       // Fill empty list with previously assigned particles
        if (Parts[i].size() < 1000)
        {
          sMcTruth = sMcTruth + "\n# of particles " + str(int(Parts[i].size())) + "\tfrom gen " + str(int(i)) + " " + fLabels[i];
        }
        else
        {
          sMcTruth = sMcTruth + "\n# of particles " + str(int(Parts[i].size())) + "\tfrom gen " + str(int(i)) + " " + fLabels[i];
        }
        TPart.push_back(Parts[i].size());
        for (std::map<int, simb::MCParticle>::iterator iter = Parts[i].begin(); iter != Parts[i].end(); iter++)
        {
          std::set<int> ThisGeneratorIDs = {};
          trackids.push_back(ThisGeneratorIDs);
          trackids[i].insert(iter->first);
        }
      }
      else
      {
        sMcTruth = sMcTruth + "\n# of particles " + str(int(Parts[i].size())) + "\tfrom gen " + str(int(i)) + " " + fLabels[i] + " *not generated!";
        TPart.push_back(0);
        std::set<int> ThisGeneratorIDs = {};
        trackids.push_back(ThisGeneratorIDs);
      }
    }
    PrintInColor(sMcTruth, GetColor("bright_red"));
    // PrintInColor("I'm here!", GetColor("red"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Some MC Truth information -------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    std::set<int> SignalTrackIDs;                                    // Signal TrackIDs to be used in OpFlash matching
    std::vector<std::vector<int>> ClPartTrackIDs = {{}, {}, {}, {}}; // Track IDs corresponding to each kind of MCTruth particle  {11,22,2112,else}
    art::Handle<std::vector<simb::MCTruth>> ThisHandle;
    std::string sNuTruth = "";
    evt.getByLabel(fLabels[0], ThisHandle);
    if (ThisHandle)
    {
      auto MarlTrue = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[0]); // Get handle for MARLEY MCTruths
      // --- Loop over all neutrinos in the event ---
      for (auto const &MARLEYtruth : *MarlTrue)
      {
        const simb::MCNeutrino &nue = MARLEYtruth.GetNeutrino();
        TNuInteraction = str(nue.InteractionType());
        TNuE = 1e3*nue.Nu().E();
        TNuX = nue.Nu().Vx();
        TNuY = nue.Nu().Vy();
        TNuZ = nue.Nu().Vz();
        int N = MARLEYtruth.NParticles();
        sNuTruth = sNuTruth + "\nNeutrino Interaction: " + TNuInteraction;
        sNuTruth = sNuTruth + "\nNumber of Producer Particles: " + str(N);
        sNuTruth = sNuTruth + "\nNeutrino energy: " + str(TNuE) + " MeV";
        sNuTruth = sNuTruth + "\nPosition (" + str(TNuX) + ", " + str(TNuY) + ", " + str(TNuZ) + ") cm";
      }
      art::FindManyP<simb::MCParticle> MarlAssn(MarlTrue, evt, fGEANTLabel);
      sNuTruth = sNuTruth + "\nGen.\tPdgCode\t\tEnergy\t\tEndPosition\t\tMother";
      sNuTruth = sNuTruth + "\n--------------------------------------------------------------------";

      for (size_t i = 0; i < MarlAssn.size(); i++)
      {
        auto parts = MarlAssn.at(i);
        for (auto part = parts.begin(); part != parts.end(); part++)
        {
          MarleyPDGList.push_back((*part)->PdgCode());
          MarleyEList.push_back(1e3*(*part)->E());
          MarleyPList.push_back(1e3*(*part)->P());
          MarleyKList.push_back(1e3*(*part)->E() - 1e3*(*part)->Mass());
          MarleyTList.push_back((*part)->T());
          MarleyEndXList.push_back((*part)->EndX());
          MarleyEndYList.push_back((*part)->EndY());
          MarleyEndZList.push_back((*part)->EndZ());
          MarleyIDList.push_back((*part)->TrackId());
          MarleyMotherList.push_back((*part)->Mother());
          std::map<int, float> MarleyMaxEDepMap, MarleyMaxEDepXMap, MarleyMaxEDepYMap, MarleyMaxEDepZMap;
          std::vector<const sim::IDE *> ides = bt_serv->TrackIdToSimIDEs_Ps((*part)->TrackId());
          for (auto const &ide : ides)
          {
            if (ide->numElectrons > 1e-6 && ide->energy > 1e-6 && abs(ide->x) > 1e-6 && abs(ide->y) > 1e-6 && abs(ide->z) > 1e-6){
              if (InMyMap((*part)->TrackId(), MarleyMaxEDepMap) == false)
              {
                MarleyMaxEDepMap[(*part)->TrackId()] = ide->energy;
                MarleyMaxEDepXMap[(*part)->TrackId()] = ide->x;
                MarleyMaxEDepYMap[(*part)->TrackId()] = ide->y;
                MarleyMaxEDepZMap[(*part)->TrackId()] = ide->z;
              }
              if (ide->energy > MarleyMaxEDepMap[(*part)->TrackId()])
              {
                MarleyMaxEDepMap[(*part)->TrackId()] = ide->energy;
                MarleyMaxEDepXMap[(*part)->TrackId()] = ide->x;
                MarleyMaxEDepYMap[(*part)->TrackId()] = ide->y;
                MarleyMaxEDepZMap[(*part)->TrackId()] = ide->z;
              }
              if (abs((*part)->PdgCode()) == 11 || abs((*part)->PdgCode()) == 22 || abs((*part)->PdgCode()) == 2112)
              {
                MarleyIDDepList.push_back((*part)->TrackId());
                MarleyEDepList.push_back(ide->energy);
                MarleyPDGDepList.push_back((*part)->PdgCode());
                MarleyXDepList.push_back(ide->x);
                MarleyYDepList.push_back(ide->y);
                MarleyZDepList.push_back(ide->z);
                MarleyElectronDepList.push_back(ide->numElectrons);
              }
            }
          }
          MarleyMaxEDepList.push_back(MarleyMaxEDepMap[(*part)->TrackId()]);
          MarleyMaxEDepXList.push_back(MarleyMaxEDepXMap[(*part)->TrackId()]);
          MarleyMaxEDepYList.push_back(MarleyMaxEDepYMap[(*part)->TrackId()]);
          MarleyMaxEDepZList.push_back(MarleyMaxEDepZMap[(*part)->TrackId()]);
          SignalTrackIDs.emplace((*part)->TrackId());

          if ((*part)->PdgCode() < 1000000)
          {
            sNuTruth = sNuTruth + "\n" + fLabels[0] + "\t" + str((*part)->PdgCode()) + "\t\t" + str(1e3*(*part)->E()) + "\t (" + str((*part)->EndX()) + ", " + str((*part)->EndY()) + ", " + str((*part)->EndZ()) + ")\t" + str((*part)->Mother());
          }
          else
          {
            sNuTruth = sNuTruth + "\n" + fLabels[0] + "\t" + str((*part)->PdgCode()) + "\t" + str(1e3*(*part)->E()) + " (" + str((*part)->EndX()) + ", " + str((*part)->EndY()) + ", " + str((*part)->EndZ()) + ")\t" + str((*part)->Mother());
          }

          if ((*part)->PdgCode() == 11) // Electrons
          {
            const TLorentzVector &v4_f = (*part)->EndPosition();
            x_f = v4_f.X();
            ClPartTrackIDs[0].push_back((*part)->TrackId());
            mf::LogDebug("SolarNuAna") << "\nMC Electron truth position x = " << v4_f.X() << ", y = " << v4_f.Y() << ", z = " << v4_f.Z();
            mf::LogDebug("SolarNuAna") << "Initial KE " << 1e3*(*part)->E() - 1e3*(*part)->Mass();
          } 
          if ((*part)->PdgCode() == 22) // Gammas
          {
            ClPartTrackIDs[1].push_back((*part)->TrackId());
          } 
          if ((*part)->PdgCode() == 2112) // Neutrons
          {
            ClPartTrackIDs[2].push_back((*part)->TrackId());
          } 
          if ((*part)->PdgCode() != 11 && (*part)->PdgCode() != 22 && (*part)->PdgCode() != 2112) // Others
          {
            ClPartTrackIDs[3].push_back((*part)->TrackId()); 
          } 
        }
      }
    }
    else
    {
      mf::LogWarning("SolarNuAna") << "No MARLEY MCTruths found.";
    }
    PrintInColor(sNuTruth, GetColor("yellow"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------- PMTrack Analysis -----------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    art::Handle<std::vector<recob::Track>> TrackHandle;
    std::vector<art::Ptr<recob::Track>> TrackList;
    if (evt.getByLabel(fTrackLabel, TrackHandle))
    {
      art::fill_ptr_vector(TrackList, TrackHandle);
    }
    TrackNum = int(TrackList.size());

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------- Optical Flash Analysis --------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // Find OpHits and OpFlashes associated with the event
    std::string sOpFlashTruth = "";
    art::Handle<std::vector<recob::OpHit>> OpHitHandle;
    std::vector<art::Ptr<recob::OpHit>> OpHitList;
    if (evt.getByLabel(fOpHitLabel, OpHitHandle))
    {
      art::fill_ptr_vector(OpHitList, OpHitHandle);
    }
    OpHitNum = int(OpHitList.size());
    if (fGenerateAdjOpFlash){
      std::vector<AdjOpHitsUtils::FlashInfo> FlashVec;
      std::vector<std::vector<art::Ptr<recob::OpHit>>> OpHitVec;
      adjophits->CalcAdjOpHitsFast(OpHitList, OpHitVec, fOpFlashAlgoDebug);
      adjophits->MakeFlashVector(FlashVec, OpHitVec, evt);
      OpFlashNum = int(FlashVec.size());
      for (int i = 0; i < int(FlashVec.size()); i++)
      {
        AdjOpHitsUtils::FlashInfo TheFlash = FlashVec[i];
        // struct FlashInfo{double NHit;double Time;double TimeWidth;double PE;double MaxPE;std::vector<double> PEperOpDet;double FastToTotal;double X;double Y;double Z;double YWidth;double ZWidth;};
        if (TheFlash.PE > fAdjOpFlashMinPECut)
        {
          double OpFlashPur = 0;
          OpFlashMaxPE.push_back(TheFlash.MaxPE);
          OpFlashPE.push_back(TheFlash.PE);
          OpFlashX.push_back(TheFlash.X);
          OpFlashY.push_back(TheFlash.Y);
          OpFlashZ.push_back(TheFlash.Z);
          OpFlashT.push_back(TheFlash.Time);
          OpFlashDeltaT.push_back(TheFlash.TimeWidth);
          OpFlashNHit.push_back(TheFlash.NHit);
          for (int j = 0; j < int(OpHitVec[i].size()); j++)
          {
            recob::OpHit OpHit = *OpHitVec[i][j];
            const std::vector< int > ThisOpHitTrackIds = pbt->OpHitToTrackIds(OpHit);
            float ThisOphitPurity = 0;
            for (auto const &ThisOpHitTrackId : ThisOpHitTrackIds)
            {
              if (SignalTrackIDs.find(ThisOpHitTrackId) != SignalTrackIDs.end())
              {
                ThisOphitPurity += 1;
              }
            }
            // Check if ThisOpHitTrackIds is empty
            if (ThisOpHitTrackIds.size() == 0){
              ThisOphitPurity = 0;
            }
            else{
              ThisOphitPurity /= int(ThisOpHitTrackIds.size());
            } 
            OpFlashPur += ThisOphitPurity;
            auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(OpHit.OpChannel()).GetCenter();
            SOpHitPur.push_back(ThisOphitPurity/int(ThisOpHitTrackIds.size()));
            SOpHitChannel.push_back(OpHit.OpChannel());
            SOpHitT.push_back(OpHit.PeakTime());
            SOpHitPE.push_back(OpHit.PE());
            SOpHitX.push_back(OpHitXYZ.X());
            SOpHitY.push_back(OpHitXYZ.Y());
            SOpHitZ.push_back(OpHitXYZ.Z());
            SOpHitFlashID.push_back(i);
          }
          // Check if OpHitVec[i] is empty
          if (OpHitVec[i].size() == 0){
            OpFlashPur = 0;
          }
          else{
            OpFlashPur /= int(OpHitVec[i].size());
          }
          OpFlashMarlPur.push_back(OpFlashPur); 
          if (abs(TheFlash.Time) < 20)
          {
            mf::LogDebug("SolarNuAna") << "Marley OpFlash PE (ratio/tot) " << TheFlash.MaxPE/TheFlash.PE << "/" << TheFlash.PE << " with purity " << OpFlashPur << " time " << TheFlash.Time; 
            sOpFlashTruth += "Marley OpFlash PE (fast/ratio/tot) " + str(TheFlash.FastToTotal) + "/" + str(TheFlash.MaxPE/TheFlash.PE) + "/" +  str(TheFlash.PE) + " with purity " + str(OpFlashPur) + " time " + str(TheFlash.Time) + " vertex (" + str(TheFlash.X) + ", " + str(TheFlash.Y) + ", " + str(TheFlash.Z) + ")\n";
          }
        }
      }
    }
    else{
      std::vector<art::Ptr<recob::OpFlash>> OpFlashList;
      art::Handle<std::vector<recob::OpFlash>> FlashHandle;
      if (evt.getByLabel(fOpFlashLabel, FlashHandle))
      {
        art::fill_ptr_vector(OpFlashList, FlashHandle);
      }
      OpFlashNum = int(OpFlashList.size());
      // Grab assns with OpHits to get match to neutrino purity
      art::FindManyP<recob::OpHit> OpAssns(OpFlashList, evt, fOpFlashLabel);

      // Loop over OpFlashList and assign OpHits to each flash
      for (int i = 0; i < int(OpFlashList.size()); i++)
      {
        recob::OpFlash TheFlash = *OpFlashList[i];
        std::vector<art::Ptr<recob::OpHit>> MatchedHits = OpAssns.at(i);
        mf::LogDebug("SolarNuAna") << "Assigning OpHit to Flash";
        double TotalFlashPE, MaxOpHitPE, FlashTime;
        TotalFlashPE = MaxOpHitPE = FlashTime = 0;
        int NMatchedHits = MatchedHits.size();
        for (int j = 0; j < NMatchedHits; j++)
        { // Loop over OpHits in the flash
          recob::OpHit OpHit = *MatchedHits[j];
          const std::vector< int > ThisOpHitTrackIds = pbt->OpHitToTrackIds(OpHit);
          float ThisOphitPurity = 0;
          for (auto const &ThisOpHitTrackId : ThisOpHitTrackIds)
          {
            if (SignalTrackIDs.find(ThisOpHitTrackId) != SignalTrackIDs.end())
            {
              ThisOphitPurity += 1;
            }
          } 
          auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(OpHit.OpChannel()).GetCenter();
          TotalFlashPE += OpHit.PE();
          FlashTime += OpHit.PeakTime()*OpHit.PE();
          if (OpHit.PE() > fAdjOpFlashMinPECut/NMatchedHits){
            SOpHitPur.push_back(ThisOphitPurity/int(ThisOpHitTrackIds.size()));
            if (OpHit.PE() > MaxOpHitPE){MaxOpHitPE = OpHit.PE();};
            SOpHitChannel.push_back(OpHit.OpChannel());
            SOpHitT.push_back(OpHit.PeakTime());
            SOpHitPE.push_back(OpHit.PE());
            SOpHitX.push_back(OpHitXYZ.X());
            SOpHitY.push_back(OpHitXYZ.Y());
            SOpHitZ.push_back(OpHitXYZ.Z());
            SOpHitFlashID.push_back(i);
          }
        } // End of OpHit loop
        FlashTime = FlashTime / TotalFlashPE;

        mf::LogDebug("SolarNuAna") << "Evaluating Flash purity";
        int TerminalOutput = supress_stdout();
        double OpFlashPur = pbt->OpHitCollectionPurity(SignalTrackIDs, MatchedHits);
        resume_stdout(TerminalOutput);
        mf::LogDebug("SolarNuAna") << "PE of this OpFlash " << TotalFlashPE << " OpFlash time " << FlashTime;

        // Calculate the flash purity, only for the Marley events
        if (MaxOpHitPE / TotalFlashPE < fAdjOpFlashMaxPERatioCut && TotalFlashPE > fAdjOpFlashMinPECut)
        {
          OpFlashMarlPur.push_back(OpFlashPur);
          OpFlashMaxPE.push_back(MaxOpHitPE);
          OpFlashPE.push_back(TheFlash.TotalPE());
          OpFlashX.push_back(TheFlash.XCenter());
          OpFlashY.push_back(TheFlash.YCenter());
          OpFlashZ.push_back(TheFlash.ZCenter());
          OpFlashT.push_back(TheFlash.Time());
          OpFlashDeltaT.push_back(TheFlash.TimeWidth());
          OpFlashNHit.push_back(MatchedHits.size());
        }
        if (abs(TheFlash.Time()) < 5)
        {
          mf::LogDebug("SolarNuAna") << "Marley OpFlash PE (max/tot) " << MaxOpHitPE << "/" << TheFlash.TotalPE() << " with purity " << OpFlashPur << " time " << TheFlash.Time(); 
          sOpFlashTruth += "Marley OpFlash PE (max/tot) " + str(MaxOpHitPE) + "/" +  str(TheFlash.TotalPE()) + " with purity " + str(OpFlashPur) + " time " + str(TheFlash.Time()) + " vertex (" + str(TheFlash.YCenter()) + ", " + str(TheFlash.ZCenter()) + ")\n";
        }
      }
    }
    PrintInColor(sOpFlashTruth, GetColor("blue"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------- Hit collection and assignment ----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // --- Lift out the reco hits:
    auto RecoHits = evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
    int NTotHits = RecoHits->size();

    for (int i = 0; i < NTotHits; ++i)
    {
      // --- Loop over the reconstructed hits to separate them among tpc planes according to view

      recob::Hit const &ThisHit = RecoHits->at(i);
      if (ThisHit.PeakTime() < 0)
        PrintInColor("Negative Hit Time = " + str(ThisHit.PeakTime()), GetColor("red"));
      mf::LogDebug("SolarNuAna") << "Hit " << i << " has view " << ThisHit.View() << " and signal type " << ThisHit.SignalType();

      if (ThisHit.SignalType() == 0 && ThisHit.View() == 0)
      {
        ColHits0.push_back(ThisHit);
      } // SignalType = 0
      else if (ThisHit.SignalType() == 0 && ThisHit.View() == 1)
      {
        ColHits1.push_back(ThisHit);
      } // SignalType = 0
      else if (ThisHit.SignalType() == 1)
      {
        ColHits2.push_back(ThisHit);
      } // SignalType = 1
      else
      {
        ColHits3.push_back(ThisHit);
        mf::LogError("SolarNuAna") << "Hit was found with view out of scope";
      }
    }

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //-------------------------------------------------------------- Cluster creation and analysis ------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // --- Now calculate the clusters ...
    adjhits->CalcAdjHits(ColHits0, Clusters0, hAdjHits, hAdjHitsADCInt, false);
    HitNum.push_back(ColHits0.size());
    ClusterNum.push_back(Clusters0.size());
    adjhits->CalcAdjHits(ColHits1, Clusters1, hAdjHits, hAdjHitsADCInt, false);
    HitNum.push_back(ColHits1.size());
    ClusterNum.push_back(Clusters1.size());
    adjhits->CalcAdjHits(ColHits2, Clusters2, hAdjHits, hAdjHitsADCInt, false);
    HitNum.push_back(ColHits2.size());
    ClusterNum.push_back(Clusters2.size());
    adjhits->CalcAdjHits(ColHits3, Clusters3, hAdjHits, hAdjHitsADCInt, false);
    HitNum.push_back(ColHits3.size());
    ClusterNum.push_back(Clusters3.size());
    fMCTruthTree->Fill();

    std::vector<std::vector<std::vector<float>>> ClGenPur = {{}, {}, {}};
    std::vector<std::vector<std::vector<recob::Hit>>> AllPlaneClusters = {Clusters0, Clusters1, Clusters2};
    std::vector<std::vector<int>> ClMainID = {{}, {}, {}}, ClTPC = {{}, {}, {}}, ClNHits = {{}, {}, {}}, ClGen = {{}, {}, {}};
    std::vector<std::vector<float>> ClCharge = {{}, {}, {}}, ClMaxCharge = {{}, {}, {}}, ClT = {{}, {}, {}}, ClX = {{}, {}, {}}, ClY = {{}, {}, {}}, ClZ = {{}, {}, {}};
    std::vector<std::vector<float>> ClFracE = {{}, {}, {}}, ClFracGa = {{}, {}, {}}, ClFracNe = {{}, {}, {}}, ClFracRest = {{}, {}, {}};
    std::vector<std::vector<float>> ClPur = {{}, {}, {}}, Cldzdy = {{}, {}, {}};

    std::string sRecoObjects = "";
    sRecoObjects += "\n# OpHits (" + fOpHitLabel + ") in full geometry: " + str(OpHitNum);
    sRecoObjects += "\n# OpFlashes (" + fOpFlashLabel + ") in full geometry: " + str(OpFlashNum);
    sRecoObjects += "\n# Hits (" + fHitLabel + ") in each view: " + str(int(ColHits0.size())) + ", " + str(int(ColHits1.size())) + ", " + str(int(ColHits2.size())) + ", " + str(int(ColHits3.size()));
    sRecoObjects += "\n# Cluster from the hits: " + str(int(Clusters0.size())) + ", " + str(int(Clusters1.size())) + ", " + str(int(Clusters2.size())) + ", " + str(int(Clusters3.size()));
    sRecoObjects += "\n# Tracks (" + fTrackLabel + ") in full geometry: " + str(TrackNum);
    PrintInColor(sRecoObjects, GetColor("cyan"));

    //------------------------------------------------------------ First complete cluster analysis ------------------------------------------------------------------//
    // --- Now loop over the planes and the clusters to calculate the cluster properties
    for (int idx = 0; idx < 3; idx++)
    {
      int nhit, clustTPC;
      float FracE, FracGa, FracNe, FracRest, clustX, clustY, clustZ, clustT, ncharge, maxHit, dzdy;
      std::vector<std::vector<recob::Hit>> Clusters = AllPlaneClusters[idx];

      // --- Loop over the clusters
      for (int i = 0; i < int(Clusters.size()); i++)
      {
        int MainTrID = 0;
        int Gen = 0;
        float Pur = 0;
        std::vector<float> thisdzdy = {};

        nhit = Clusters[i].size();
        ncharge = maxHit = clustT = FracE = FracGa = FracNe = FracRest = clustX = clustY = clustZ = clustTPC = dzdy = 0;
        std::vector<float> GenPur = {};
        for (size_t genpur = 0; genpur < fLabels.size(); genpur++)
        {
          GenPur.push_back(0);
        }

        for (recob::Hit TPCHit : Clusters[i])
        {
          if (TPCHit.PeakTime() < 0)
            PrintInColor("Negative Cluster Time = " + str(TPCHit.PeakTime()), GetColor("red"));
          ncharge += TPCHit.Integral();
          const geo::WireGeo *wire = wireReadout.WirePtr(TPCHit.WireID()); // Wire directions should be the same for all hits of the same view (can be used to check)
          double hitCharge;

          geo::Point_t hXYZ = wire->GetCenter();
          geo::Point_t sXYZ = wire->GetStart();
          geo::Point_t eXYZ = wire->GetEnd();
          geo::Vector_t direction = eXYZ - sXYZ;
          auto dyds = direction.Y(), dzds = direction.Z();
          thisdzdy.push_back(dzds / dyds);

          int TPC = TPCHit.WireID().TPC;
          clustTPC += TPCHit.Integral() * TPC;
          clustX += TPCHit.Integral() * hXYZ.X();
          clustY += TPCHit.Integral() * hXYZ.Y();
          clustZ += TPCHit.Integral() * hXYZ.Z();
          clustT += TPCHit.Integral() * TPCHit.PeakTime();

          if (TPCHit.Integral() > maxHit)
          {
            maxHit = TPCHit.Integral();
          } // Look for maxHit inside cluster

          MainTrID = 0;
          double TopEFrac = -DBL_MAX;
          std::vector<sim::TrackIDE> ThisHitIDE = bt_serv->HitToTrackIDEs(clockData, TPCHit);

          for (size_t ideL = 0; ideL < ThisHitIDE.size(); ++ideL)
          {
            if (ThisHitIDE[ideL].energyFrac > TopEFrac)
            {
              TopEFrac = ThisHitIDE[ideL].energyFrac;
              MainTrID = ThisHitIDE[ideL].trackID;
              mf::LogDebug("SolarNuAna") << "This TPCHit's IDE is: " << MainTrID;
            }
          }

          for (int frac = 0; frac < int(ClPartTrackIDs.size()); ++frac)
          {
            for (int trck = 0; trck < int(ClPartTrackIDs[frac].size()); ++trck)
            {
              if (abs(MainTrID) == ClPartTrackIDs[frac][trck])
              {
                if (frac == 0)
                {
                  FracE = FracE + TPCHit.Integral();
                }
                if (frac == 1)
                {
                  FracGa = FracGa + TPCHit.Integral();
                }
                if (frac == 2)
                {
                  FracNe = FracNe + TPCHit.Integral();
                }
                if (frac == 3)
                {
                  FracRest = FracRest + TPCHit.Integral();
                }
              }
            }
          }

          long unsigned int ThisPType = WhichParType(abs(MainTrID));
          GenPur[int(ThisPType)] = GenPur[int(ThisPType)] + TPCHit.Integral();
          mf::LogDebug("SolarNuAna") << "\nThis particle type " << ThisPType << "\nThis cluster's main track ID " << MainTrID;
          if (ThisPType == 1)
          {
            hitCharge = TPCHit.Integral();
            Pur = Pur + hitCharge;
          }
        }

        float MaxGen = 0;
        for (size_t genpur = 0; genpur < GenPur.size(); genpur++)
        {
          if (GenPur[genpur] > MaxGen)
          {
            Gen = genpur;
            MaxGen = GenPur[genpur];
          }
          GenPur[genpur] = GenPur[genpur] / ncharge;
        }

        for (size_t j = 0; j > thisdzdy.size(); j++)
        {
          if (thisdzdy[0] != thisdzdy[i])
            mf::LogWarning("SolarNuAna") << "MISSMATCH IN dzdy FOR CLUSTER " << idx;
        }

        dzdy = thisdzdy[0];
        thisdzdy.clear();
        FracE /= ncharge;
        FracGa /= ncharge;
        FracNe /= ncharge;
        FracRest /= ncharge;
        clustTPC /= ncharge;
        clustX /= ncharge;
        clustY /= ncharge;
        clustZ /= ncharge;
        clustT /= ncharge;
        mf::LogDebug("SolarNuAna") << "\ndzdy " << dzdy << " for cluster "
                                  << " (" << clustY << ", " << clustZ << ") with track ID " << MainTrID << " in plane " << idx;
        if (clustT < 0)
          PrintInColor("Negative Cluster Time = " + str(clustT), GetColor("red"));

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
        ClPur[idx].push_back(Pur / ncharge);
        ClGen[idx].push_back(Gen);
        Cldzdy[idx].push_back(dzdy);
        ClMainID[idx].push_back(MainTrID);
        ClGenPur[idx].push_back(GenPur);

        mf::LogDebug("SolarNuAna") << "\nCluster " << i << " in plane " << idx << " has #hits" << nhit << " charge, " << ncharge << " time, " << clustT;
        mf::LogDebug("SolarNuAna") << " and position (" << clustY << ", " << clustZ << ") with main track ID " << MainTrID << " and purity " << Pur / ncharge;
      }
    } // Finished first cluster processing

    //-------------------------------------------------------------------- Cluster Matching -------------------------------------------------------------------------//
    std::vector<std::vector<float>> MVecGenFrac = {};
    std::vector<int> MVecNHit = {}, MVecGen = {}, MVecInd0NHits = {}, MVecInd1NHits = {}, MVecMainID = {}, MVecTPC = {}, MVecInd0TPC = {}, MVecInd1TPC = {};
    std::vector<float> MVecTime = {}, MVecCharge = {}, MVecMaxCharge = {}, MVecInd0Charge = {}, MVecInd1Charge = {}, MVecInd0MaxCharge = {}, MVecInd1MaxCharge = {}, MVecInd0dT = {}, MVecInd1dT = {};
    std::vector<float> MVecInd0RecoY = {}, MVecInd1RecoY = {}, MVecRecY = {}, MVecRecZ = {};
    std::vector<float> MVecFracE = {}, MVecFracGa = {}, MVecFracNe = {}, MVecFracRest = {}, MVecPur = {};

    for (int ii = 0; ii < int(AllPlaneClusters[2].size()); ii++)
    {
      bool match = false;
      int ind0clustNHits = 0, ind1clustNHits = 0;
      int ind0clustTPC = 0, ind1clustTPC = 0;
      double ind0clustY = -1e6, ind1clustY = -1e6, ind0clustMaxCharge = 0, ind1clustMaxCharge = 0, ind0clustCharge = 0, ind1clustCharge = 0;
      double ind0clustdT = fClusterMatchTime, ind1clustdT = fClusterMatchTime;
      if (!AllPlaneClusters[2][ii].empty())
      {
        if (!AllPlaneClusters[0].empty())
        {
          for (int jj = 0; jj < int(AllPlaneClusters[0].size()); jj++)
          {
            if (ClNHits[0][jj] < (1 - fClusterMatchNHit) * ClNHits[2][ii] || ClNHits[0][jj] > (1 + fClusterMatchNHit) * ClNHits[2][ii])
            {
              continue;
            }
            if (ClCharge[0][jj] < (1 - fClusterMatchCharge) * ClCharge[2][ii] || ClCharge[0][jj] > (1 + fClusterMatchCharge) * ClCharge[2][ii])
            {
              continue;
            }
            if (abs(ClT[2][ii] - ClT[0][jj]) < fClusterMatchTime && abs(fClusterInd0MatchTime - abs(ClT[2][ii] - ClT[0][jj])) < abs(fClusterInd0MatchTime - ind0clustdT))
            {
              ind0clustY = ClY[0][jj] + (ClZ[2][ii] - ClZ[0][jj]) / (Cldzdy[0][jj]);
              ind0clustdT = abs(ClT[2][ii] - ClT[0][jj]);
              ind0clustNHits = int(AllPlaneClusters[0][jj].size());
              ind0clustCharge = ClCharge[0][jj];
              ind0clustMaxCharge = ClMaxCharge[0][jj];
              ind0clustTPC = ClTPC[0][jj];
              if (ind0clustY > -fDetectorSizeY && ind0clustY < fDetectorSizeY)
              {
                match = true;
              }
              mf::LogDebug("SolarNuAna") << "¡¡¡ Matched cluster in plane 0 !!! --- Position x = " << ClX[0][jj] << ", y = " << ClY[0][jj] << ", z = " << ClZ[0][jj];
              mf::LogDebug("SolarNuAna") << "Reconstructed position y = " << ind0clustY << ", z = " << ClZ[2][ii];
            }
          }
        }
        if (!AllPlaneClusters[1].empty())
        {
          for (int zz = 0; zz < int(AllPlaneClusters[1].size()); zz++)
          {
            if (ClNHits[1][zz] < (1 - fClusterMatchNHit) * ClNHits[2][ii] || ClNHits[1][zz] > (1 + fClusterMatchNHit) * ClNHits[2][ii])
            {
              continue;
            }
            if (ClCharge[1][zz] < (1 - fClusterMatchCharge) * ClCharge[2][ii] || ClCharge[1][zz] > (1 + fClusterMatchCharge) * ClCharge[2][ii])
            {
              continue;
            }
            if (abs(ClT[2][ii] - ClT[1][zz]) < fClusterMatchTime && abs(fClusterInd1MatchTime - abs(ClT[2][ii] - ClT[1][zz])) < abs(fClusterInd1MatchTime - ind1clustdT))
            {
              ind1clustY = ClY[1][zz] + (ClZ[2][ii] - ClZ[1][zz]) / (Cldzdy[1][zz]);
              ind1clustdT = abs(ClT[2][ii] - ClT[1][zz]);
              ind1clustNHits = int(AllPlaneClusters[1][zz].size());
              ind1clustCharge = ClCharge[1][zz];
              ind1clustMaxCharge = ClMaxCharge[1][zz];
              ind1clustTPC = ClTPC[1][zz];
              if (ind1clustY > -fDetectorSizeY && ind1clustY < fDetectorSizeY)
              {
                match = true;
              }
              mf::LogDebug("SolarNuAna") << "¡¡¡ Matched cluster in plane 1 !!! --- Position x = " << ClX[1][zz] << ", y = " << ClY[1][zz] << ", z = " << ClZ[1][zz];
              mf::LogDebug("SolarNuAna") << "Reconstructed position y = " << ind1clustY << ", z = " << ClZ[2][ii];
            }
          } // Loop over ind1 clusters
        }
      } // Loop over ind clusters
      else
      {
        mf::LogDebug("SolarNuAna") << "Cluster " << ii << " in plane 2 has no hits";
      }

      //--------------------------------------------------------- Export Matched cluster vectors ------------------------------------------------------------------//
      if (match == true)
      {
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
        if ((ind0clustY > -buffer * fDetectorSizeY && ind0clustY < buffer * fDetectorSizeY) && (ind1clustY > -buffer * fDetectorSizeY && ind1clustY < buffer * fDetectorSizeY))
        {
          mf::LogDebug("SolarNuAna") << "BOTH IND RECO INSIDE OF DETECTOR";
          MVecRecY.push_back((ind0clustY + ind1clustY) / 2);
        }
        else if (ind0clustY > -buffer * fDetectorSizeY && ind0clustY < buffer * fDetectorSizeY)
        {
          mf::LogDebug("SolarNuAna") << "IND1 OUTSIDE OF DETECTOR";
          MVecRecY.push_back(ind0clustY);
        }
        else if (ind1clustY > -buffer * fDetectorSizeY && ind1clustY < buffer * fDetectorSizeY)
        {
          mf::LogDebug("SolarNuAna") << "IND0 OUTSIDE OF DETECTOR";
          MVecRecY.push_back(ind1clustY);
        }
        else
        {
          mf::LogDebug("SolarNuAna") << "RECO OUTSIDE OF DETECTOR";
          MVecRecY.push_back((ind0clustY + ind1clustY) / 2);
          if (ClGen[2][ii] == 1)
          {
            mf::LogWarning("SolarNuAna") << "Marley cluster reconstructed outside of detector volume! RecoY = " << str((ind0clustY + ind1clustY) / 2);
          }
        }

        // Print in color if the cluster is matched
        mf::LogDebug("SolarNuAna") << "¡¡¡ Matched cluster !!! ";
        mf::LogDebug("SolarNuAna") << " - Cluster " << str(ClMainID[2][ii]) << " Gen " << str(ClGen[2][ii]) << " Purity " << str(ClPur[2][ii]) << " Hits " << str(ClNHits[2][ii]);
        mf::LogDebug("SolarNuAna") << " - Hits(ind0, ind1, col) " << str(ind0clustNHits) << ", " << str(ind1clustNHits) << ", " << str(ClNHits[2][ii]);
        mf::LogDebug("SolarNuAna") << " - Positions y(ind0, ind1) = " << str(ind0clustY) << ", " << str(ind1clustY) << ", z = " << str(ClZ[2][ii]) << "\n";
      } // if (match == true)
    }   // Loop over collection plane clusters

    //-------------------------------------------------------------------- Cluster Tree Export -------------------------------------------------------------------------//
    // Loop over matched clusters and export to tree if number of hits is above threshold
    for (int i = 0; i < int(MVecNHit.size()); i++)
    {
      std::string sClusterReco = "";
      std::string sResultColor = "white";
      bool TrackMatch = false;
      if (MVecCharge[i] < fMinClusterCharge)
      {
        continue;
      }
      if (MVecNHit[i] > fClusterPreselectionNHit && (MVecInd0NHits[i] > fClusterPreselectionNHit || MVecInd1NHits[i] > fClusterPreselectionNHit))
      {
        MPrimary = true;
        MAdjClTime = {};
        MAdjClCharge = {};
        MAdjClInd0Charge = {};
        MAdjClInd1Charge = {};
        MAdjClMaxCharge = {};
        MAdjClInd0MaxCharge = {};
        MAdjClInd1MaxCharge = {};
        MAdjClNHit = {};
        MAdjClInd0NHit = {};
        MAdjClInd1NHit = {};
        MAdjClRecoY = {};
        MAdjClRecoZ = {};
        MAdjClR = {};
        MAdjClPur = {};
        MAdjClGen = {};
        MAdjClMainID = {};
        MAdjClMainPDG = {};
        MAdjClMainE = {};
        MAdjClMainX = {};
        MAdjClMainY = {};
        MAdjClMainZ = {};
        MAdjClEndX = {};
        MAdjClEndY = {};
        MAdjClEndZ = {};
        MAdjFlashTime = {};
        MAdjFlashPE = {};
        MAdjFlashNHit = {};
        MAdjFlashMaxPE = {};
        MAdjFlashRecoX = {};
        MAdjFlashRecoY = {};
        MAdjFlashRecoZ = {};
        MAdjFlashR = {};
        MAdjFlashPur = {};
        MTrackStart = {-1e6, -1e6, -1e6};
        MTrackEnd = {-1e6, -1e6, -1e6};

        for (int j = 0; j < int(MVecNHit.size()); j++)
        {
          if (j == i)
          {
            continue;
          }
          if (fGeometry == "HD" && sqrt(pow(MVecRecY[i] - MVecRecY[j], 2) + pow(MVecRecZ[i] - MVecRecZ[j], 2) + pow((MVecTime[i] - MVecTime[j]) * fDetectorSizeX / fDetectorDriftTime, 2)) > fAdjClusterRad)
          {
            continue;
          }
          if (fGeometry == "VD" && sqrt(pow(MVecRecY[i] - MVecRecY[j], 2) + pow(MVecRecZ[i] - MVecRecZ[j], 2) + pow((MVecTime[i] - MVecTime[j]) * fDetectorSizeX / (fDetectorDriftTime/2), 2)) > fAdjClusterRad)
          {
            continue;
          }
          if (MVecCharge[j] < fMinClusterCharge)
          {
            continue;
          }
          if (MVecCharge[j] > MVecCharge[i])
          {
            MPrimary = false;
          }
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
          MAdjClR.push_back(sqrt(pow(MVecRecY[i] - MVecRecY[j], 2) + pow(MVecRecZ[i] - MVecRecZ[j], 2)));
          MAdjClPur.push_back(MVecPur[j]);
          MAdjClGen.push_back(MVecGen[j]);
          MAdjClMainID.push_back(MVecMainID[j]);

          // If mother exists add the mother information
          const simb::MCParticle *MAdjClTruth;
          int TerminalOutput = supress_stdout();
          MAdjClTruth = pi_serv->TrackIdToParticle_P(MVecMainID[j]);
          resume_stdout(TerminalOutput);
          if (MAdjClTruth == 0)
          {
            MAdjClMainPDG.push_back(0);
            MAdjClMainE.push_back(-1e6);
            MAdjClMainX.push_back(-1e6);
            MAdjClMainY.push_back(-1e6);
            MAdjClMainZ.push_back(-1e6);
            MAdjClEndX.push_back(-1e6);
            MAdjClEndY.push_back(-1e6);
            MAdjClEndZ.push_back(-1e6);
          }
          else
          {
            MAdjClMainPDG.push_back(MAdjClTruth->PdgCode());
            MAdjClMainE.push_back(MAdjClTruth->E());
            MAdjClMainX.push_back(MAdjClTruth->Vx());
            MAdjClMainY.push_back(MAdjClTruth->Vy());
            MAdjClMainZ.push_back(MAdjClTruth->Vz());
            MAdjClEndX.push_back(MAdjClTruth->EndX());
            MAdjClEndY.push_back(MAdjClTruth->EndY());
            MAdjClEndZ.push_back(MAdjClTruth->EndZ());
          }
        }

        if (sqrt(pow(MVecRecY[i] - TNuY,2) + pow(MVecRecZ[i] - TNuZ,2)) < 20)
        {
          sResultColor = "green";
        }
        else
        {
          sResultColor = "yellow";
        }
        sClusterReco += "*** Matched preselection cluster: \n - Primary  " + str(MPrimary) + " Gen " + str(MVecGen[i]) + " Purity " + str(MVecPur[i]) + " Hits " + str(MVecNHit[i]) + "\n - RecoY, RecoZ (" + str(MVecRecY[i]) + ", " + str(MVecRecZ[i]) + ") Time " + str(MVecTime[i]) + "\n";

        if (MPrimary){
          TVector3 ThisClVertex = {0, MVecRecY[i], MVecRecZ[i]};
          float MaxVertexDistance = 10; // if track is further away from ThisClVertex than
          for (int i = 0; i < TrackNum; i++) { // using index loop to get track idx
            recob::Track trk = *TrackList[i];
            TVector3 trk_start(0, trk.Start().Y(), trk.Start().Z());
            TVector3 trk_end(0, trk.End().Y(), trk.End().Z());
            // throw away bad tracks
            if ((trk_start - ThisClVertex).Mag() > MaxVertexDistance && (trk_end - ThisClVertex).Mag() > MaxVertexDistance)
              {continue;};
            MTrackNPoints = trk.NPoints();
            MTrackStart = {trk.Start().X(), trk.Start().Y(), trk.Start().Z()};
            MTrackEnd = {trk.End().X(), trk.End().Y(), trk.End().Z()};
            MTrackChi2 = trk.Chi2();
            sClusterReco += "*** Matched pmtrack: \n";
            sClusterReco += " - Track has start (" + str(trk.Start().X()) + ", " + str(trk.Start().Y()) + ", " + str(trk.Start().Z()) + ")\n";
            sClusterReco += " - Track has end   (" + str(trk.End().X()) + ", " + str(trk.End().Y()) + ", " + str(trk.End().Z()) + ")\n";
            TrackMatch = true;  
          }; // Loop over tracks
        };

        for (int j = 0; j < int(OpFlashPE.size()); j++)
        {
          if ((MVecTime[i] - OpFlashT[j]) < 0)
          {
            continue;
          }
          if ((MVecTime[i] - OpFlashT[j]) > fAdjOpFlashTime)
          {
            continue;
          }
          // if (sqrt(pow(MVecRecY[i] - OpFlashY[j], 2) + pow(MVecRecZ[i] - OpFlashZ[j], 2)) > fAdjOpFlashRad)
          // Instead of a circular cut, we apply an elliptical cut (dy/a)^2+(dz/b)^2<1
          if (pow(abs(MVecRecY[i] - OpFlashY[j])/fAdjOpFlashY, 2) + pow(abs(MVecRecZ[i] - OpFlashZ[j])/fAdjOpFlashZ, 2) > 1)
          {
            continue;
          }
          MAdjFlashTime.push_back(OpFlashT[j]);
          MAdjFlashPE.push_back(OpFlashPE[j]);
          MAdjFlashNHit.push_back(OpFlashNHit[j]);
          MAdjFlashMaxPE.push_back(OpFlashMaxPE[j]);
          MAdjFlashRecoX.push_back(OpFlashX[j]);
          MAdjFlashRecoY.push_back(OpFlashY[j]);
          MAdjFlashRecoZ.push_back(OpFlashZ[j]);
          MAdjFlashR.push_back(sqrt(pow(MVecRecY[i] - OpFlashY[j], 2) + pow(MVecRecZ[i] - OpFlashZ[j], 2)));
          MAdjFlashPur.push_back(OpFlashMarlPur[j]);
        }

        // Fill the tree with the cluster information and the adjacent clusters and flashes
        MPur = MVecPur[i];
        MGen = MVecGen[i];
        MMarleyFrac = {MVecFracE[i], MVecFracGa[i], MVecFracNe[i], MVecFracRest[i]};
        MGenFrac = MVecGenFrac[i];
        // Cluster TPC
        MTPC = MVecTPC[i];
        MInd0TPC = MVecInd0TPC[i];
        MInd1TPC = MVecInd1TPC[i];
        // Cluster Charge
        MCharge = MVecCharge[i];
        MMaxCharge = MVecMaxCharge[i];
        MInd0Charge = MVecInd0Charge[i];
        MInd1Charge = MVecInd1Charge[i];
        MInd0MaxCharge = MVecInd0MaxCharge[i];
        MInd1MaxCharge = MVecInd1MaxCharge[i];
        // Cluster Hits
        MNHit = MVecNHit[i];
        MInd0NHits = MVecInd0NHits[i];
        MInd1NHits = MVecInd1NHits[i];
        // Cluster Time
        MTime = MVecTime[i];
        MInd0dT = MVecInd0dT[i];
        MInd1dT = MVecInd1dT[i];
        // Cluster RecoY
        MRecY = MVecRecY[i];
        MInd0RecoY = MVecInd0RecoY[i];
        MInd1RecoY = MVecInd1RecoY[i];
        // Cluster RecoZ
        MRecZ = MVecRecZ[i];
        // Cluster MainID
        MMainID = MVecMainID[i];

        // If mother exists add the mother information
        const simb::MCParticle *MClTruth;
        int TerminalOutput = supress_stdout();
        MClTruth = pi_serv->TrackIdToParticle_P(MVecMainID[i]);
        resume_stdout(TerminalOutput);
        if (MClTruth == 0)
        {
          MMainVertex = {-1e6, -1e6, -1e6};
          MEndVertex = {-1e6, -1e6, -1e6};
          MMainPDG = 0;
          MMainE = -1e6;
          MMainP = -1e6;
          MMainK = -1e6;
          MMainT = -1e6;

          MMainParentVertex = {-1e6, -1e6, -1e6};
          MMainParentPDG = 0;
          MMainParentE = -1e6;
          MMainParentP = -1e6;
          MMainParentK = -1e6;
          MMainParentT = -1e6;
        }
        else
        {
          MMainVertex = {MClTruth->Vx(), MClTruth->Vy(), MClTruth->Vz()};
          MEndVertex = {MClTruth->EndX(), MClTruth->EndY(), MClTruth->EndZ()};
          MMainPDG = MClTruth->PdgCode();
          MMainE = 1e3*MClTruth->E();
          MMainP = 1e3*MClTruth->P();
          MMainK = MMainE - 1e3*MClTruth->Mass();
          MMainT = MClTruth->T();
          // If exists add the parent information
          const simb::MCParticle *MClParentTruth;
          int TerminalOutput = supress_stdout();
          MClParentTruth = pi_serv->TrackIdToParticle_P(MClTruth->Mother());
          resume_stdout(TerminalOutput);
          if (MClParentTruth == 0)
          {
            MMainParentVertex = {-1e6, -1e6, -1e6};
            MMainParentPDG = 0;
            MMainParentE = -1e6;
            MMainParentP = -1e6;
            MMainParentK = -1e6;
            MMainParentT = -1e6;
          }
          else
          {
            MMainParentVertex = {MClParentTruth->Vx(), MClParentTruth->Vy(), MClParentTruth->Vz()};
            MMainParentPDG = MClParentTruth->PdgCode();
            MMainParentE = 1e3*MClParentTruth->E();
            MMainParentP = 1e3*MClParentTruth->P();
            MMainParentK = MMainParentE - 1e3*MClParentTruth->Mass();
            MMainParentT = MClParentTruth->T();
          }
        }
        if (fClusterPreselectionTrack){if (TrackMatch == false){continue;}}
        fSolarNuAnaTree->Fill();
        hDriftTime->Fill(x_f, MTime);
        hXTruth->Fill(MVecRecY[i] - TNuY, TNuX);
        hYTruth->Fill(MVecRecY[i] - TNuY, TNuY);
        hZTruth->Fill(MVecRecZ[i] - TNuZ, TNuZ);
        if (MVecTime[i] < 0)
          mf::LogWarning("SolarNuAna") << "Negative Main Cluster Time = " << MVecTime[i];
      }
      // Check if the string sClusterReco is not empty and print it in color
      if (sClusterReco != "") {PrintInColor(sClusterReco, GetColor(sResultColor));}
    }
  }
  // ########################################################################################################################################//
  //_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_//
  // ########################################################################################################################################//

  //......................................................
  void SolarNuAna::ResetVariables()
  // Reset variables for each event
  {
    // Clear Marley MCTruth info.
    TNuE = 0;
    TNuX = 0;
    TNuY = 0;
    TNuZ = 0;
    TrackNum = 0;
    OpHitNum = 0;
    OpFlashNum = 0;
    MTrackNPoints = 0;
    MTrackChi2 = 0;
    MarleyElectronDepList = {};
    MarleyPDGList = {}; MarleyPDGDepList = {};
    MarleyIDList = {}, MarleyIDDepList = {};
    MarleyMotherList = {};
    MarleyEList = {}, MarleyPList = {}, MarleyKList = {}, MarleyTList = {}, MarleyEndXList = {}, MarleyEndYList = {}, MarleyEndZList = {};
    MarleyEDepList = {}, MarleyXDepList = {}, MarleyYDepList = {}, MarleyZDepList = {};
    MarleyMaxEDepList = {}, MarleyMaxEDepXList = {}, MarleyMaxEDepYList = {}, MarleyMaxEDepZList = {};
    SOpHitChannel = {}, SOpHitPur = {}, SOpHitPE = {}, SOpHitX = {}, SOpHitY = {}, SOpHitZ = {}, SOpHitT = {}, SOpHitFlashID = {};
    TPart = {}, Parts = {};
    HitNum = {};
    ClusterNum = {};
    // ThisGeneratorParts.clear();
    OpFlashMarlPur.clear();
    OpFlashPE.clear();
    OpFlashMaxPE.clear();
    OpFlashX.clear();
    OpFlashY.clear();
    OpFlashZ.clear();
    OpFlashT.clear();
    OpFlashDeltaT.clear();
    OpFlashNHit.clear();

  } // ResetVariables

  //......................................................
  long unsigned int SolarNuAna::WhichParType(int TrID)
  {
    for (long unsigned int i = 0; i < fLabels.size(); i++)
    {
      if (InMyMap(TrID, Parts[i]))
      {
        return i + 1;
      }
    }
    // If no match, then who knows???
    return 0;
  }

  //......................................................
  // This function fills a map with the MCParticles from a given MCTruth
  void SolarNuAna::FillMyMaps(std::map<int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle<std::vector<simb::MCTruth>> Hand)
  {
    for (size_t L1 = 0; L1 < Hand->size(); ++L1)
    {
      for (size_t L2 = 0; L2 < Assn.at(L1).size(); ++L2)
      {
        const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
        MyMap[ThisPar.TrackId()] = ThisPar;
        mf::LogDebug("SolarNuAna") << ThisPar.PdgCode() << " " << ThisPar.E();
      }
    }
    return;
  }

  //......................................................
  // This function checks if a given TrackID is in a given map
  bool SolarNuAna::InMyMap(int TrID, std::map<int, simb::MCParticle> ParMap)
  {
    std::map<int, simb::MCParticle>::iterator ParIt;
    ParIt = ParMap.find(TrID);
    if (ParIt != ParMap.end())
    {
      return true;
    }
    else
      return false;
  }

  bool SolarNuAna::InMyMap(int TrID, std::map<int, float> TrackIDMap)
  {
    std::map<int, float>::iterator ParIt;
    ParIt = TrackIDMap.find(TrID);
    if (ParIt != TrackIDMap.end())
    {
      return true;
    }
    else
      return false;
  }


  //......................................................
  // This function creates a terminal color printout
  void SolarNuAna::PrintInColor(std::string MyString, int Color, std::string Type)
  {
    if (Type == "Info")
    {
      mf::LogInfo("SolarNuAna") << "\033[" << Color << "m" << MyString << "\033[0m";
    }
    if (Type == "Degub")
    {
      mf::LogDebug("SolarNuAna") << "\033[" << Color << "m" << MyString << "\033[0m";
    }
    if (Type == "Error")
    {
      mf::LogError("SolarNuAna") << "\033[" << Color << "m" << MyString << "\033[0m";
    }
    return;
  }

  // ......................................................
  // This function returns an integer that corresponds to a given color name
  int SolarNuAna::GetColor(std::string ColorName)
  {
    if (ColorName == "black")
      return 30;
    else if (ColorName == "red")
      return 31;
    else if (ColorName == "green")
      return 32;
    else if (ColorName == "yellow")
      return 33;
    else if (ColorName == "blue")
      return 34;
    else if (ColorName == "magenta")
      return 35;
    else if (ColorName == "cyan")
      return 36;
    else if (ColorName == "white")
      return 37;
    else if (ColorName == "bright_black")
      return 90;
    else if (ColorName == "bright_red")
      return 91;
    else if (ColorName == "bright_green")
      return 92;
    else if (ColorName == "bright_yellow")
      return 93;
    else if (ColorName == "bright_blue")
      return 94;
    else if (ColorName == "bright_magenta")
      return 95;
    else if (ColorName == "bright_cyan")
      return 96;
    else if (ColorName == "bright_white")
      return 97;
    else
    {
      mf::LogError("SolarNuAna") << "Color " << ColorName << " not recognized. Returning white.";
      return 37;
    }
    return 0;
  }

  std::string SolarNuAna::str(int i)
  {
    std::stringstream ss;
    ss << i;
    return ss.str();
  }
  std::string SolarNuAna::str(double i, int prec)
  {
  // Use prec to define the precision of the output in terms of number of decimal places
    std::stringstream ss;
    ss << std::fixed << std::setprecision(prec) << i;
    return ss.str();
  }
  std::string SolarNuAna::str(float i, int prec)
  {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(prec) << i;
    return ss.str();
  }
  std::string SolarNuAna::str(std::vector<int> i)
  {
    std::stringstream ss;
    for (int j = 0; j < int(i.size()); j++)
    {
      ss << i[j] << " ";
    }
    return ss.str();
  }
  std::string SolarNuAna::str(std::vector<double> i, int prec)
  {
    std::stringstream ss;
    for (int j = 0; j < int(i.size()); j++)
    {
      ss << std::fixed << std::setprecision(prec) << i[j] << " ";
    }
    return ss.str();
  }
  std::string SolarNuAna::str(std::vector<float> i, int prec)
  {
    std::stringstream ss;
    for (int j = 0; j < int(i.size()); j++)
    {
      ss << std::fixed << std::setprecision(prec) << i[j] << " ";
    }
    return ss.str();
  }

  int SolarNuAna::supress_stdout()
  {
    std::fflush(stdout);

    int ret = dup(1);
    int nullfd = open("/dev/null", O_WRONLY);
    // check nullfd for error omitted
    dup2(nullfd, 1);
    close(nullfd);

    return ret;
  }

  void SolarNuAna::resume_stdout(int fd)
  {
    std::fflush(stdout);
    dup2(fd, 1);
    close(fd);
  }
} // namespace solar
DEFINE_ART_MODULE(solar::SolarNuAna)
