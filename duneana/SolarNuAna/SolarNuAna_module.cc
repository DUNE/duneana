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
#include "TVector3.h"
#include "TRandom.h"
#include <fcntl.h>

// Larsoft includes (not all might be necessary)
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

// Art includes and others
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
#include "fhiclcpp/ParameterSet.h"

#include "duneopdet/SolarNuUtils/SolarAuxUtils.h"
#include "duneopdet/SolarNuUtils/AdjOpHitsUtils.h"
#include "dunereco/LowEUtils/LowEUtils.h"

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
    void analyze(art::Event const &evt) override;
    void reconfigure(fhicl::ParameterSet const &p);
    void beginJob() override;

  private:
    // --- Some of our own functions.
    void ResetVariables();

    // --- Our fcl parameter labels for the modules that made the data products
    std::string fHitLabel, fTrackLabel, fOpHitLabel, fOpFlashLabel, fGEANTLabel;

    // --- Input settings imported from the fcl
    std::string fGeometry;
    int fDetectorSizeX, fDetectorSizeY, fDetectorSizeZ, fDetectorDriftTime, fClusterAlgoAdjChannel, fClusterInd0MatchTime, fClusterInd1MatchTime, fClusterPreselectionNHit, fAdjOpFlashMinNHitCut;
    float fClusterMatchTime, fAdjClusterRad, fMinClusterCharge, fClusterMatchCharge, fAdjOpFlashY, fAdjOpFlashZ, fAdjOpFlashTime, fAdjOpFlashMaxPERatioCut, fAdjOpFlashMinPECut, fClusterMatchNHit, fClusterAlgoTime;
    std::vector<std::string> fLabels;
    float fOpFlashAlgoMinTime, fOpFlashAlgoMaxTime, fOpFlashAlgoRad, fOpFlashAlgoPE, fOpFlashAlgoTriggerPE, fOpFlashAlgoHotVertexThld;
    bool fClusterPreselectionTrack, fClusterPreselectionPrimary, fGenerateAdjOpFlash, fSaveSignalEDep, fSaveSignalOpHits, fSaveOpFlashInfo, fSaveTrackInfo, fFlashMatchByResidual;
    // bool fOpFlashAlgoCentroid;

    // --- Our TTrees, and its associated variables.
    TTree *fConfigTree;
    TTree *fMCTruthTree;
    TTree *fSolarNuAnaTree;
    std::string TNuInteraction;
    std::vector<std::map<int, simb::MCParticle>> GeneratorParticles = {};
    int Event, Flag, MNHit, MGen, MTPC, MInd0TPC, MInd1TPC, MInd0NHits, MInd1NHits, MMainID, MMainPDG, MMainParentPDG, TrackNum, OpHitNum, OpFlashNum, MTrackNPoints;
    float TNuE, TNuX, TNuY, TNuZ, MTime, MCharge, MMaxCharge, MInd0Charge, MInd1Charge, MInd0MaxCharge, MInd1MaxCharge;
    float MInd0dT, MInd1dT, MInd0RecoY, MInd1RecoY, MRecX, MRecY, MRecZ, MPur, MGenPur, MMainE, MMainP, MMainK, MMainT, MMainParentE, MMainParentP, MMainParentK, MMainParentT, MTrackChi2;
    std::vector<int> MAdjClGen, MAdjClMainID, TPart, SignalPDGList, SignalPDGDepList, SignalIDList, SignalMotherList, SignalIDDepList, MAdjClMainPDG, HitNum, ClusterNum, SignalElectronDepList;
    std::vector<float> SignalEDepList, SignalXDepList, SignalYDepList, SignalZDepList;
    std::vector<float> SOpHitPur, SOpHitPE, SOpHitX, SOpHitY, SOpHitZ, SOpHitT, SOpHitChannel, SOpHitFlashID;
    std::vector<float> MAdjClTime, MAdjClCharge, MAdjClInd0Charge, MAdjClInd1Charge, MAdjClMaxCharge, MAdjClInd0MaxCharge, MAdjClInd1MaxCharge;
    std::vector<float> MAdjClNHit, MAdjClInd0NHit, MAdjClInd1NHit, MAdjClRecoY, MAdjClRecoZ, MAdjClR, MAdjClPur, MAdjClGenPur, MAdjClMainE, MAdjClMainK;
    std::vector<float> MAdjClMainX, MAdjClMainY, MAdjClMainZ, MAdjClEndX, MAdjClEndY, MAdjClEndZ, MSignalFrac, MGenFrac;
    std::vector<float> MAdjFlashTime, MAdjFlashResidual, MAdjFlashPE, MAdjFlashNHit, MAdjFlashMaxPE, MAdjFlashRecoX, MAdjFlashRecoY, MAdjFlashRecoZ, MAdjFlashR, MAdjFlashPur, MAdjFlashSTD, MAdjFlashFast;
    std::vector<float> SignalEList, SignalPList, SignalKList, SignalTList, SignalEndXList, SignalEndYList, SignalEndZList, SignalMaxEDepList, SignalMaxEDepXList, SignalMaxEDepYList, SignalMaxEDepZList;
    std::vector<double> MMainVertex, MEndVertex, MMainParentVertex;
    std::vector<double> MTrackStart, MTrackEnd;
    bool MPrimary;

    // --- OpFlash Variables
    std::vector<int> OpFlashID, OpFlashNHit;
    std::vector<float> OpFlashPur, OpFlashPE, OpFlashMaxPE, OpFlashX, OpFlashY, OpFlashZ, OpFlashT, OpFlashDeltaT, OpFlashSTD, OpFlashFast;

    // --- MatchedFlash Variables
    int MFlashNHit;
    float MFlashR, MFlashPE, MFlashMaxPE, MFlashPur, MFlashFast, MFlashTime, MFlashSTD, MFlashRecoX, MFlashRecoY, MFlashRecoZ, MFlashResidual;
    bool MFlashCorrect;

    // --- Histograms to fill about collection plane hits
    float MainElectronEndPointX;
    TH2F *hXTruth;
    TH2F *hYTruth;
    TH2F *hZTruth;
    TH1I *hAdjHits;
    TH1F *hAdjHitsADCInt;
    TH2F *hDriftTime;

    // --- Declare our services
    geo::WireReadoutGeom const &wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    std::unique_ptr<solar::SolarAuxUtils> solaraux;
    std::unique_ptr<solar::AdjOpHitsUtils> adjophits;
    std::unique_ptr<solar::LowEUtils> lowe;
  };
#endif

  //......................................................
  SolarNuAna::SolarNuAna(fhicl::ParameterSet const &p)
      : EDAnalyzer(p),
        solaraux(new solar::SolarAuxUtils(p)),
        adjophits(new solar::AdjOpHitsUtils(p)),
        lowe(new solar::LowEUtils(p))
  {
    this->reconfigure(p);
  }

  //......................................................
  void SolarNuAna::reconfigure(fhicl::ParameterSet const &p)
  {
    fLabels = p.get<std::vector<std::string>>("ParticleLabelVector");
    fHitLabel = p.get<std::string>("HitLabel");
    fOpFlashLabel = p.get<std::string>("OpFlashLabel");
    fOpHitLabel = p.get<std::string>("OpHitLabel");
    fTrackLabel = p.get<std::string>("TrackLabel");
    fGEANTLabel = p.get<std::string>("GEANT4Label");
    fGeometry = p.get<std::string>("Geometry");
    fDetectorSizeX = p.get<int>("DetectorSizeX");
    fDetectorSizeY = p.get<int>("DetectorSizeY");
    fDetectorSizeZ = p.get<int>("DetectorSizeZ");
    fDetectorDriftTime = p.get<int>("DetectorDriftTime");
    fClusterAlgoTime = p.get<float>("ClusterAlgoTime");
    fClusterAlgoAdjChannel = p.get<int>("ClusterAlgoAdjChannel");
    fClusterMatchNHit = p.get<float>("ClusterMatchNHit");
    fClusterMatchCharge = p.get<float>("ClusterMatchCharge");
    fClusterMatchTime = p.get<float>("ClusterMatchTime");
    fClusterInd0MatchTime = p.get<float>("ClusterInd0MatchTime");
    fClusterInd1MatchTime = p.get<float>("ClusterInd1MatchTime");
    fClusterPreselectionNHit = p.get<int>("ClusterPreselectionNHit");
    fClusterPreselectionTrack = p.get<bool>("ClusterPreselectionTrack");
    fClusterPreselectionPrimary = p.get<bool>("ClusterPreselectionPrimary");
    fAdjClusterRad = p.get<float>("AdjClusterRad");
    fMinClusterCharge = p.get<float>("MinClusterCharge");
    fGenerateAdjOpFlash = p.get<bool>("GenerateAdjOpFlash");
    fOpFlashAlgoMinTime = p.get<double>("OpFlashAlgoMinTime");
    fOpFlashAlgoMaxTime = p.get<double>("OpFlashAlgoMaxTime");
    fOpFlashAlgoRad = p.get<double>("OpFlashAlgoRad");
    fOpFlashAlgoPE = p.get<float>("OpFlashAlgoPE");
    fOpFlashAlgoTriggerPE = p.get<float>("OpFlashAlgoTriggerPE");
    fOpFlashAlgoHotVertexThld = p.get<float>("OpFlashAlgoHotVertexThld");
    // fOpFlashAlgoCentroid = p.get<bool>("OpFlashAlgoCentroid");
    fAdjOpFlashTime = p.get<float>("AdjOpFlashTime");
    fAdjOpFlashY = p.get<float>("AdjOpFlashY");
    fAdjOpFlashZ = p.get<float>("AdjOpFlashZ");
    fAdjOpFlashMaxPERatioCut = p.get<float>("AdjOpFlashMaxPERatioCut");
    fAdjOpFlashMinPECut = p.get<float>("AdjOpFlashMinPECut");
    fAdjOpFlashMinNHitCut = p.get<int>("AdjOpFlashMinNHitCut");
    fFlashMatchByResidual = p.get<bool>("FlashMatchByResidual");
    fSaveSignalEDep = p.get<bool>("SaveSignalEDep");
    fSaveSignalOpHits = p.get<bool>("SaveSignalOpHits");
    fSaveOpFlashInfo = p.get<bool>("SaveOpFlashInfo");
    fSaveTrackInfo = p.get<bool>("SaveTrackInfo");
  } // Reconfigure

  //......................................................
  void SolarNuAna::beginJob()
  {
    // --- Make our handle to the TFileService
    art::ServiceHandle<art::TFileService> tfs;
    fConfigTree = tfs->make<TTree>("ConfigTree", "Config Tree");
    fMCTruthTree = tfs->make<TTree>("MCTruthTree", "MC Truth Tree");
    fSolarNuAnaTree = tfs->make<TTree>("SolarNuAnaTree", "Solar Ana Tree");

    // Larsoft Config info.
    fConfigTree->Branch("GEANT4Label", &fGEANTLabel);
    fConfigTree->Branch("HitLabel", &fHitLabel);
    fConfigTree->Branch("TrackLabel", &fTrackLabel);
    fConfigTree->Branch("OpHitLabel", &fOpHitLabel);
    fConfigTree->Branch("OpFlashLabel", &fOpFlashLabel);
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
    fConfigTree->Branch("ClusterPreselectionPrimary", &fClusterPreselectionPrimary);
    fConfigTree->Branch("AdjClusterRad", &fAdjClusterRad);
    fConfigTree->Branch("MinClusterCharge", &fMinClusterCharge);
    fConfigTree->Branch("GenerateAdjOpFlash", &fGenerateAdjOpFlash);
    fConfigTree->Branch("OpFlashAlgoMinTime", &fOpFlashAlgoMinTime);
    fConfigTree->Branch("OpFlashAlgoMaxTime", &fOpFlashAlgoMaxTime);
    fConfigTree->Branch("OpFlashAlgoRad", &fOpFlashAlgoRad);
    fConfigTree->Branch("OpFlashAlgoPE", &fOpFlashAlgoPE);
    fConfigTree->Branch("OpFlashAlgoTriggerPE", &fOpFlashAlgoTriggerPE);
    fConfigTree->Branch("OpFlashAlgoHotVertexThld", &fOpFlashAlgoHotVertexThld);
    // fConfigTree->Branch("OpFlashAlgoCentroid", &fOpFlashAlgoCentroid);
    fConfigTree->Branch("AdjOpFlashTime", &fAdjOpFlashTime);
    fConfigTree->Branch("AdjOpFlashY", &fAdjOpFlashY);
    fConfigTree->Branch("AdjOpFlashZ", &fAdjOpFlashZ);
    fConfigTree->Branch("AdjOpFlashMaxPERatioCut", &fAdjOpFlashMaxPERatioCut);
    fConfigTree->Branch("AdjOpFlashMinPECut", &fAdjOpFlashMinPECut);
    fConfigTree->Branch("AdjOpFlashMinNHitCut", &fAdjOpFlashMinNHitCut);
    fConfigTree->Branch("FlashMatchByResidual", &fFlashMatchByResidual);
    fConfigTree->Branch("SaveSignalEDep", &fSaveSignalEDep);
    fConfigTree->Branch("SaveSignalOpHits", &fSaveSignalOpHits);
    fConfigTree->Branch("SaveOpFlashInfo", &fSaveOpFlashInfo);
    fConfigTree->Branch("SaveTrackInfo", &fSaveTrackInfo);

    // MC Truth info.
    fMCTruthTree->Branch("Event", &Event, "Event/I");                // Event number
    fMCTruthTree->Branch("Flag", &Flag, "Flag/I");                   // Flag used to match truth with reco tree entries
    fMCTruthTree->Branch("TruthPart", &TPart);                       // Number particles per generator
    fMCTruthTree->Branch("Interaction", &TNuInteraction);            // True neutrino interaction process
    fMCTruthTree->Branch("TNuE", &TNuE, "TruthNuE/F");               // True neutrino energy [MeV]
    fMCTruthTree->Branch("TNuX", &TNuX, "TruthNuX/F");               // True neutrino X [cm]
    fMCTruthTree->Branch("TNuY", &TNuY, "TruthNuY/F");               // True neutrino Y [cm]
    fMCTruthTree->Branch("TNuZ", &TNuZ, "TruthNuZ/F");               // True neutrino Z [cm]
    fMCTruthTree->Branch("TSignalPDG", &SignalPDGList);              // PDG of Signal marticles
    fMCTruthTree->Branch("TSignalE", &SignalEList);                  // Energy of Signal particles [MeV]
    fMCTruthTree->Branch("TSignalP", &SignalPList);                  // Momentum of Signal particles [MeV]
    fMCTruthTree->Branch("TSignalK", &SignalKList);                  // Kinetik Energy of Signal particles [MeV]
    fMCTruthTree->Branch("TSignalT", &SignalTList);                  // Time of Signal particles [ticks]
    fMCTruthTree->Branch("TSignalEndX", &SignalEndXList);            // X of Signal particles [cm]
    fMCTruthTree->Branch("TSignalEndY", &SignalEndYList);            // Y of Signal particles [cm]
    fMCTruthTree->Branch("TSignalEndZ", &SignalEndZList);            // Z of Signal particles [cm]
    fMCTruthTree->Branch("TSignalMaxEDep", &SignalMaxEDepList);      // Energy of Signal particles [MeV]
    fMCTruthTree->Branch("TSignalX", &SignalMaxEDepXList);           // X of Signal particles [cm]
    fMCTruthTree->Branch("TSignalY", &SignalMaxEDepYList);           // Y of Signal particles [cm]
    fMCTruthTree->Branch("TSignalZ", &SignalMaxEDepZList);           // Z of Signal particles [cm]
    fMCTruthTree->Branch("TSignalID", &SignalIDList);                // TrackID of Signal particles
    fMCTruthTree->Branch("TSignalMother", &SignalMotherList);        // TrackID of Signal mother
    fMCTruthTree->Branch("OpHitNum", &OpHitNum, "OpHitNum/I");       // Number of OpHits
    fMCTruthTree->Branch("OpFlashNum", &OpFlashNum, "OpFlashNum/I"); // Number of OpFlashes
    fMCTruthTree->Branch("HitNum", &HitNum);                         // Number of hits in each TPC plane
    fMCTruthTree->Branch("ClusterNum", &ClusterNum);                 // Number of clusters in each TPC plane
    fMCTruthTree->Branch("TrackNum", &TrackNum, "TrackNum/I");       // Number of PMTracks
    if (fSaveSignalEDep)
    {
      fMCTruthTree->Branch("TSignalPDGDepList", &SignalPDGDepList);           // PDG for Energy deposited of Signal particles
      fMCTruthTree->Branch("TSignalEDepList", &SignalEDepList);               // Energy deposited of Signal particles [MeV]
      fMCTruthTree->Branch("TSignalXDepList", &SignalXDepList);               // X deposited of Signal particles [cm]
      fMCTruthTree->Branch("TSignalYDepList", &SignalYDepList);               // Y deposited of Signal particles [cm]
      fMCTruthTree->Branch("TSignalZDepList", &SignalZDepList);               // Z deposited of Signal particles [cm]
      fMCTruthTree->Branch("TSignalIDDepList", &SignalIDDepList);             // ParentID of Signal particles
      fMCTruthTree->Branch("TSignalElectronDepList", &SignalElectronDepList); // Number of electrons in the Signal particles
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
    if (fSaveOpFlashInfo)
    {
      fMCTruthTree->Branch("OpFlashX", &OpFlashX);         // OpFlash X
      fMCTruthTree->Branch("OpFlashY", &OpFlashY);         // OpFlash Y
      fMCTruthTree->Branch("OpFlashZ", &OpFlashZ);         // OpFlash Z
      fMCTruthTree->Branch("OpFlashID", &OpFlashID);       // OpFlash ID
      fMCTruthTree->Branch("OpFlashPE", &OpFlashPE);       // OpFlash PE
      fMCTruthTree->Branch("OpFlashTime", &OpFlashT);      // OpFlash Time
      fMCTruthTree->Branch("OpFlashPur", &OpFlashPur);     // OpFlash Purity
      fMCTruthTree->Branch("OpFlashSTD", &OpFlashSTD);     // OpFlash STD
      fMCTruthTree->Branch("OpFlashNHit", &OpFlashNHit);   // OpFlash NHit
      fMCTruthTree->Branch("OpFlashMaxPE", &OpFlashMaxPE); // OpFlash Max PE
    }

    // Repeated Truth info.
    fSolarNuAnaTree->Branch("Event", &Event, "Event/I");           // Event number
    fSolarNuAnaTree->Branch("Flag", &Flag, "Flag/I");              // Flag used to match truth with reco tree entries
    fSolarNuAnaTree->Branch("TruthPart", &TPart);                  // Number particles per generator
    fSolarNuAnaTree->Branch("Interaction", &TNuInteraction);       // True neutrino interaction process
    fSolarNuAnaTree->Branch("TNuE", &TNuE, "TruthNuE/F");          // True neutrino energy
    fSolarNuAnaTree->Branch("TNuX", &TNuX, "TruthNuX/F");          // True neutrino X
    fSolarNuAnaTree->Branch("TNuY", &TNuY, "TruthNuY/F");          // True neutrino Y
    fSolarNuAnaTree->Branch("TNuZ", &TNuZ, "TruthNuZ/F");          // True neutrino Z
    fSolarNuAnaTree->Branch("TSignalPDG", &SignalPDGList);         // PDG of Signal particles
    fSolarNuAnaTree->Branch("TSignalE", &SignalEList);             // Energy of Signal particles
    fSolarNuAnaTree->Branch("TSignalK", &SignalKList);             // Kinetik Energy of Signal particles
    fSolarNuAnaTree->Branch("TSignalEndX", &SignalEndXList);       // X of Signal particles
    fSolarNuAnaTree->Branch("TSignalEndY", &SignalEndYList);       // Y of Signal particles
    fSolarNuAnaTree->Branch("TSignalEndZ", &SignalEndZList);       // Z of Signal particles
    fSolarNuAnaTree->Branch("TSignalMaxEDep", &SignalMaxEDepList); // Max Energy Deposition of Signal particles
    fSolarNuAnaTree->Branch("TSignalX", &SignalMaxEDepXList);      // Max Energy Deposition X of Signal particles
    fSolarNuAnaTree->Branch("TSignalY", &SignalMaxEDepYList);      // Max Energy Deposition Y of Signal particles
    fSolarNuAnaTree->Branch("TSignalZ", &SignalMaxEDepZList);      // Max Energy Deposition Z of Signal particles
    fSolarNuAnaTree->Branch("TSignalID", &SignalIDList);           // TrackID of Signal particles")
    fSolarNuAnaTree->Branch("TSignalMother", &SignalMotherList);   // TrackID of Signal particles")

    // Main Cluster info.
    fSolarNuAnaTree->Branch("Primary", &MPrimary);                                // Cluster hasn't any adjcl with AdjClCharge > MCharge (bool)
    fSolarNuAnaTree->Branch("Purity", &MPur, "Purity/F");                         // Main cluster reco signal purity
    fSolarNuAnaTree->Branch("Generator", &MGen, "Generator/I");                   // Main cluster generator idx
    fSolarNuAnaTree->Branch("GenPurity", &MGenPur, "GenPurity/F");                // Main cluster reco generator purity
    fSolarNuAnaTree->Branch("TPC", &MTPC, "ColTPC/I");                            // Main cluster TPC
    fSolarNuAnaTree->Branch("Time", &MTime, "ColTime/F");                         // Main cluster time [ticks]
    fSolarNuAnaTree->Branch("NHits", &MNHit, "ColNHits/I");                       // Main cluster #hits
    fSolarNuAnaTree->Branch("Charge", &MCharge, "ColCharge/F");                   // Main cluster charge [ADC*ticks]
    fSolarNuAnaTree->Branch("MaxCharge", &MMaxCharge, "ColCharge/F");             // Main cluster's max TPCHit-charge [ADC*ticks]
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
    fSolarNuAnaTree->Branch("Ind0RecoY", &MInd0RecoY, "Ind0RecoY/F");             // Main cluster ind0 reco Y [cm]
    fSolarNuAnaTree->Branch("Ind1RecoY", &MInd1RecoY, "Ind1RecoY/F");             // Main cluster ind1 reco Y [cm]
    fSolarNuAnaTree->Branch("RecoX", &MRecX, "RecoX/F");                          // Main cluster reco X [cm] (from matched Flash)
    fSolarNuAnaTree->Branch("RecoY", &MRecY, "RecoY/F");                          // Main cluster reco Y [cm]
    fSolarNuAnaTree->Branch("RecoZ", &MRecZ, "RecoZ/F");                          // Main cluster reco Z [cm]
    fSolarNuAnaTree->Branch("MainID", &MMainID, "MainID/I");                      // Main cluster main track ID
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
    fSolarNuAnaTree->Branch("SignalFrac", &MSignalFrac);                          // Main cluster particle contribution (electron, gamma, neutron)

    // Track info.
    if (fSaveTrackInfo){
      fSolarNuAnaTree->Branch("MTrackNPoints", &MTrackNPoints, "TrackNPoints/I"); // Track #points
      fSolarNuAnaTree->Branch("MTrackStart", &MTrackStart);                       // Track start point
      fSolarNuAnaTree->Branch("MTrackEnd", &MTrackEnd);                           // Track end point
      fSolarNuAnaTree->Branch("MTrackChi2", &MTrackChi2);                         // Track chi2
    }

    // Adj. Cluster info.
    fSolarNuAnaTree->Branch("AdjClGen", &MAdjClGen);                     // Adj. clusters' generator idx
    fSolarNuAnaTree->Branch("AdjClGenPur", &MAdjClGenPur);               // Adj. clusters' generator purity
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
    if (fSaveOpFlashInfo)
    {
      fSolarNuAnaTree->Branch("AdjOpFlashR", &MAdjFlashR);               // Adj. flash' reco distance [cm]
      fSolarNuAnaTree->Branch("AdjOpFlashPE", &MAdjFlashPE);             // Adj. flash' tot #PE [ADC*ticks]
      fSolarNuAnaTree->Branch("AdjOpFlashPur", &MAdjFlashPur);           // Adj. flash' purity
      fSolarNuAnaTree->Branch("AdjOpFlashSTD", &MAdjFlashSTD);           // Adj. flash' STD
      fSolarNuAnaTree->Branch("AdjOpFlashFast", &MAdjFlashFast);         // Adj. flash' Fast Component
      fSolarNuAnaTree->Branch("AdjOpFlashTime", &MAdjFlashTime);         // Adj. flash' time [ticks]
      fSolarNuAnaTree->Branch("AdjOpFlashNHit", &MAdjFlashNHit);         // Adj. flash' #hits
      fSolarNuAnaTree->Branch("AdjOpFlashMaxPE", &MAdjFlashMaxPE);       // Adj. flash' max #PE [ADC*ticks]
      fSolarNuAnaTree->Branch("AdjOpFlashRecoX", &MAdjFlashRecoX);       // Adj. flash' reco X [cm]
      fSolarNuAnaTree->Branch("AdjOpFlashRecoY", &MAdjFlashRecoY);       // Adj. flash' reco Y [cm]
      fSolarNuAnaTree->Branch("AdjOpFlashRecoZ", &MAdjFlashRecoZ);       // Adj. flash' reco Z [cm]
      fSolarNuAnaTree->Branch("AdjOpFlashResidual", &MAdjFlashResidual); // Adj. flash' residual wrt. cluster
    }

    // Matched Flash info.
    fSolarNuAnaTree->Branch("MatchedOpFlashR", &MFlashR, "MatchedOpFlashR/F");                      // Matched flash' reco distance [cm]
    fSolarNuAnaTree->Branch("MatchedOpFlashPE", &MFlashPE, "MatchedOpFlashPE/F");                   // Matched flash' tot #PE [ADC*ticks]
    fSolarNuAnaTree->Branch("MatchedOpFlashPur", &MFlashPur, "MatchedOpFlashPur/F");                // Matched flash' purity
    fSolarNuAnaTree->Branch("MatchedOpFlashSTD", &MFlashSTD, "MatchedOpFlashSTD/F");                // Matched flash' STD
    fSolarNuAnaTree->Branch("MatchedOpFlashFast", &MFlashFast, "MatchedOpFlashFast/F");             // Matched flash' Fast Component
    fSolarNuAnaTree->Branch("MatchedOpFlashNHit", &MFlashNHit, "MatchedOpFlashNHit/I");             // Matched flash' #hits
    fSolarNuAnaTree->Branch("MatchedOpFlashTime", &MFlashTime, "MatchedOpFlashTime/F");             // Matched flash' time [ticks]
    fSolarNuAnaTree->Branch("MatchedOpFlashMaxPE", &MFlashMaxPE, "MatchedOpFlashMaxPE/F");          // Matched flash' max #PE [ADC*ticks]
    fSolarNuAnaTree->Branch("MatchedOpFlashRecoX", &MFlashRecoX, "MatchedOpFlashRecoX/F");          // Matched flash' reco X [cm]
    fSolarNuAnaTree->Branch("MatchedOpFlashRecoY", &MFlashRecoY, "MatchedOpFlashRecoY/F");          // Matched flash' reco Y [cm]
    fSolarNuAnaTree->Branch("MatchedOpFlashRecoZ", &MFlashRecoZ, "MatchedOpFlashRecoZ/F");          // Matched flash' reco Z [cm]
    fSolarNuAnaTree->Branch("MatchedOpFlashResidual", &MFlashResidual, "MatchedOpFlashResidual/F"); // Matched flash' residual wrt. cluster
    fSolarNuAnaTree->Branch("MatchedOpFlashCorrectly", &MFlashCorrect);                             // Matched flash' correctnes (bool)

    fConfigTree->AddFriend(fSolarNuAnaTree);
    fMCTruthTree->AddFriend(fSolarNuAnaTree);
    fConfigTree->Fill();

    // --- Our Histograms...
    hDriftTime = tfs->make<TH2F>("hDriftTime", "hDriftTime", 100, -400., 400., 100, 0., 10000.);
    hXTruth = tfs->make<TH2F>("hXTruth", "Missmatch in X distance; Distance [cm]; True X position [cm]", 100, -600, 600, 100, -600, 600);
    hYTruth = tfs->make<TH2F>("hYTruth", "Missmatch in Y distance; Distance [cm]; True Y position [cm]", 100, -600, 600, 100, -600, 600);
    hZTruth = tfs->make<TH2F>("hZTruth", "Missmatch in Z distance; Distance [cm]; True Z position [cm]", 100, -600, 600, 100, 0, 1600);
    hAdjHits = tfs->make<TH1I>("hAdjHits", "Number of adjacent collection plane hits; Number of adjacent collection plane hits; Number of events", 21, -0.5, 20.5);
    hAdjHitsADCInt = tfs->make<TH1F>("hAdjHitsADCInt", "Total summed ADC Integrals for clusters; Total summed ADC Integrals for clusters; Number of events", 1000, 0, 10000);
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
    sHead = sHead + "\nTPC Frequency in [MHz]: " + SolarAuxUtils::str(clockData.TPCClock().Frequency());
    sHead = sHead + "\nTPC Tick in [us]: " + SolarAuxUtils::str(clockData.TPCClock().TickPeriod());
    sHead = sHead + "\nEvent Flag: " + SolarAuxUtils::str(Flag);
    sHead = sHead + "\nSuccesfull reset of variables for evt " + SolarAuxUtils::str(Event);
    sHead = sHead + "\n#########################################";
    solaraux->PrintInColor(sHead, SolarAuxUtils::GetColor("magenta"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Create maps for ID tracking -----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // --- Fill MC Truth IDs to tracking vectors. Get a list of all of my particles in one chunk. ---
    const sim::ParticleList &PartList = pi_serv->ParticleList();
    std::string sMcTruth = "";
    sMcTruth = sMcTruth + "\nThere are a total of " + SolarAuxUtils::str(int(PartList.size())) + " Particles in the event\n";

    // Loop over all signal+bkg handles and collect track IDs
    for (size_t i = 0; i < fLabels.size(); i++)
    {
      GeneratorParticles.push_back(ThisGeneratorParts); // For each label insert empty list

      art::Handle<std::vector<simb::MCTruth>> ThisHandle;
      evt.getByLabel(fLabels[i], ThisHandle);

      if (ThisHandle)
      {
        auto ThisValidHanlde = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[i]); // Get generator handles
        art::FindManyP<simb::MCParticle> Assn(ThisValidHanlde, evt, fGEANTLabel);          // Assign labels to MCPArticles
        solaraux->FillMyMaps(GeneratorParticles[i], Assn, ThisValidHanlde);                          // Fill empty list with previously assigned particles
        if (GeneratorParticles[i].size() < 1000)
        {
          sMcTruth = sMcTruth + "\n# of particles " + SolarAuxUtils::str(int(GeneratorParticles[i].size())) + "\tfrom gen " + SolarAuxUtils::str(int(i) + 1) + " " + fLabels[i];
        }
        else
        {
          sMcTruth = sMcTruth + "\n# of particles " + SolarAuxUtils::str(int(GeneratorParticles[i].size())) + "\tfrom gen " + SolarAuxUtils::str(int(i) + 1) + " " + fLabels[i];
        }
        TPart.push_back(GeneratorParticles[i].size());
        if (GeneratorParticles[i].size() > 0)
        {
          for (std::map<int, simb::MCParticle>::iterator iter = GeneratorParticles[i].begin(); iter != GeneratorParticles[i].end(); iter++)
          {
            std::set<int> ThisGeneratorIDs = {};
            trackids.push_back(ThisGeneratorIDs);
            trackids[i].insert(iter->first);
          }
        }
        else
        {
          std::set<int> ThisGeneratorIDs = {};
          trackids.push_back(ThisGeneratorIDs);
        }
      }
      else
      {
        sMcTruth = sMcTruth + "\n# of particles " + SolarAuxUtils::str(int(GeneratorParticles[i].size())) + "\tfrom gen " + SolarAuxUtils::str(int(i) + 1) + " " + fLabels[i] + " *not generated!";
        TPart.push_back(0);
        std::set<int> ThisGeneratorIDs = {};
        trackids.push_back(ThisGeneratorIDs);
      }
    }
    solaraux->PrintInColor(sMcTruth, SolarAuxUtils::GetColor("bright_red"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Some MC Truth information -------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    std::set<int> SignalTrackIDs;                                    // Signal TrackIDs to be used in OpFlash matching
    std::vector<std::vector<int>> ClPartTrackIDs = {{}, {}, {}, {}}; // Track IDs corresponding to each kind of MCTruth particle  {11,22,2112,else}
    art::Handle<std::vector<simb::MCTruth>> ThisHandle;
    std::string sSignalTruth = "";
    evt.getByLabel(fLabels[0], ThisHandle);
    if (ThisHandle)
    {
      auto Signal = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[0]); // Get handle for SIGNAL MCTruths
      // --- Loop over all neutrinos in the event ---
      for (auto const &SignalTruth : *Signal)
      {
        int NSignalParticles = SignalTruth.NParticles();
        sSignalTruth = sSignalTruth + "\nNumber of Signal Particles: " + SolarAuxUtils::str(NSignalParticles);
        if (fLabels[0] == "marley")
        {
          const simb::MCNeutrino &nue = SignalTruth.GetNeutrino();
          TNuInteraction = SolarAuxUtils::str(nue.InteractionType());
          TNuE = 1e3 * nue.Nu().E();
          TNuX = nue.Nu().Vx();
          TNuY = nue.Nu().Vy();
          TNuZ = nue.Nu().Vz();
          sSignalTruth = sSignalTruth + "\nNeutrino Interaction: " + TNuInteraction;
          sSignalTruth = sSignalTruth + "\nNeutrino Energy: " + SolarAuxUtils::str(TNuE) + " MeV";
          sSignalTruth = sSignalTruth + "\nPosition (" + SolarAuxUtils::str(TNuX) + ", " + SolarAuxUtils::str(TNuY) + ", " + SolarAuxUtils::str(TNuZ) + ") cm";
        }
        else
        {
          sSignalTruth = sSignalTruth + "\nMarley generator label not found in first position generator list! Defaulting to general truth information.";
          for (int i = 0; i < NSignalParticles; i++)
          {
            const simb::MCParticle &signalparticle = SignalTruth.GetParticle(i);
            float SignalParticleE = 1e3 * signalparticle.E();
            int SignalParticlePDG = signalparticle.PdgCode();
            float SignalParticleX = signalparticle.Vx();
            float SignalParticleY = signalparticle.Vy();
            float SignalParticleZ = signalparticle.Vz();
            sSignalTruth = sSignalTruth + "\nSignal PDG: " + SolarAuxUtils::str(SignalParticlePDG);
            sSignalTruth = sSignalTruth + "\nEnergy: " + SolarAuxUtils::str(SignalParticleE) + " MeV";
            sSignalTruth = sSignalTruth + "\nPosition (" + SolarAuxUtils::str(SignalParticleX) + ", " + SolarAuxUtils::str(SignalParticleY) + ", " + SolarAuxUtils::str(SignalParticleZ) + ") cm";
          }
        }
      }
      art::FindManyP<simb::MCParticle> SignalAssn(Signal, evt, fGEANTLabel);
      sSignalTruth = sSignalTruth + "\n\tGen.\tPdgCode\t\tEnergy\t\tEndPosition\t\tMother";
      sSignalTruth = sSignalTruth + "\n------------------------------------------------------------------------";

      for (size_t i = 0; i < SignalAssn.size(); i++)
      {
        auto SignalParticles = SignalAssn.at(i);
        for (auto SignalParticle = SignalParticles.begin(); SignalParticle != SignalParticles.end(); SignalParticle++)
        {
          SignalPDGList.push_back((*SignalParticle)->PdgCode());
          SignalEList.push_back(1e3 * (*SignalParticle)->E());
          SignalPList.push_back(1e3 * (*SignalParticle)->P());
          SignalKList.push_back(1e3 * (*SignalParticle)->E() - 1e3 * (*SignalParticle)->Mass());
          SignalTList.push_back((*SignalParticle)->T());
          SignalEndXList.push_back((*SignalParticle)->EndX());
          SignalEndYList.push_back((*SignalParticle)->EndY());
          SignalEndZList.push_back((*SignalParticle)->EndZ());
          SignalIDList.push_back((*SignalParticle)->TrackId());
          SignalMotherList.push_back((*SignalParticle)->Mother());
          std::map<int, float> SignalMaxEDepMap, SignalMaxEDepXMap, SignalMaxEDepYMap, SignalMaxEDepZMap;
          std::vector<const sim::IDE *> ides = bt_serv->TrackIdToSimIDEs_Ps((*SignalParticle)->TrackId());
          for (auto const &ide : ides)
          {
            if (ide->numElectrons > 1e-6 && ide->energy > 1e-6 && abs(ide->x) > 1e-6 && abs(ide->y) > 1e-6 && abs(ide->z) > 1e-6)
            {
              if (SolarAuxUtils::InMyMap((*SignalParticle)->TrackId(), SignalMaxEDepMap) == false)
              {
                SignalMaxEDepMap[(*SignalParticle)->TrackId()] = ide->energy;
                SignalMaxEDepXMap[(*SignalParticle)->TrackId()] = ide->x;
                SignalMaxEDepYMap[(*SignalParticle)->TrackId()] = ide->y;
                SignalMaxEDepZMap[(*SignalParticle)->TrackId()] = ide->z;
              }
              if (ide->energy > SignalMaxEDepMap[(*SignalParticle)->TrackId()])
              {
                SignalMaxEDepMap[(*SignalParticle)->TrackId()] = ide->energy;
                SignalMaxEDepXMap[(*SignalParticle)->TrackId()] = ide->x;
                SignalMaxEDepYMap[(*SignalParticle)->TrackId()] = ide->y;
                SignalMaxEDepZMap[(*SignalParticle)->TrackId()] = ide->z;
              }
              if (abs((*SignalParticle)->PdgCode()) == 11 || abs((*SignalParticle)->PdgCode()) == 22 || abs((*SignalParticle)->PdgCode()) == 2112)
              {
                SignalIDDepList.push_back((*SignalParticle)->TrackId());
                SignalEDepList.push_back(ide->energy);
                SignalPDGDepList.push_back((*SignalParticle)->PdgCode());
                SignalXDepList.push_back(ide->x);
                SignalYDepList.push_back(ide->y);
                SignalZDepList.push_back(ide->z);
                SignalElectronDepList.push_back(ide->numElectrons);
              }
            }
          }
          SignalMaxEDepList.push_back(SignalMaxEDepMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepXList.push_back(SignalMaxEDepXMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepYList.push_back(SignalMaxEDepYMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepZList.push_back(SignalMaxEDepZMap[(*SignalParticle)->TrackId()]);
          SignalTrackIDs.emplace((*SignalParticle)->TrackId());

          if ((*SignalParticle)->PdgCode() < 1000000)
          {
            sSignalTruth = sSignalTruth + "\n" + fLabels[0] + "\t" + SolarAuxUtils::str((*SignalParticle)->PdgCode()) + "\t\t" + SolarAuxUtils::str(1e3 * (*SignalParticle)->E()) + "\t (" + SolarAuxUtils::str((*SignalParticle)->EndX()) + ", " + SolarAuxUtils::str((*SignalParticle)->EndY()) + ", " + SolarAuxUtils::str((*SignalParticle)->EndZ()) + ")\t" + SolarAuxUtils::str((*SignalParticle)->Mother());
          }
          else
          {
            sSignalTruth = sSignalTruth + "\n" + fLabels[0] + "\t" + SolarAuxUtils::str((*SignalParticle)->PdgCode()) + "\t" + SolarAuxUtils::str(1e3 * (*SignalParticle)->E()) + " (" + SolarAuxUtils::str((*SignalParticle)->EndX()) + ", " + SolarAuxUtils::str((*SignalParticle)->EndY()) + ", " + SolarAuxUtils::str((*SignalParticle)->EndZ()) + ")\t" + SolarAuxUtils::str((*SignalParticle)->Mother());
          }

          if ((*SignalParticle)->PdgCode() == 11) // Electrons
          {
            const TLorentzVector &MainElectronEndPoint = (*SignalParticle)->EndPosition();
            MainElectronEndPointX = MainElectronEndPoint.X();
            ClPartTrackIDs[0].push_back((*SignalParticle)->TrackId());
            mf::LogDebug("SolarNuAna") << "\nMC Electron truth position x = " << MainElectronEndPoint.X() << ", y = " << MainElectronEndPoint.Y() << ", z = " << MainElectronEndPoint.Z();
            mf::LogDebug("SolarNuAna") << "Initial KE " << 1e3 * (*SignalParticle)->E() - 1e3 * (*SignalParticle)->Mass();
          }
          if ((*SignalParticle)->PdgCode() == 22) // Gammas
          {
            ClPartTrackIDs[1].push_back((*SignalParticle)->TrackId());
          }
          if ((*SignalParticle)->PdgCode() == 2112) // Neutrons
          {
            ClPartTrackIDs[2].push_back((*SignalParticle)->TrackId());
          }
          if ((*SignalParticle)->PdgCode() != 11 && (*SignalParticle)->PdgCode() != 22 && (*SignalParticle)->PdgCode() != 2112) // Others
          {
            ClPartTrackIDs[3].push_back((*SignalParticle)->TrackId());
          }
        }
      }
    }
    else
    {
      mf::LogWarning("SolarNuAna") << "No SIGNAL MCTruths found.";
    }
    sSignalTruth += "\nSignal Track IDs: ";
    for (auto const &SignalTrackID : SignalTrackIDs)
    {
      sSignalTruth += SolarAuxUtils::str(SignalTrackID) + "; ";
    }
    solaraux->PrintInColor(sSignalTruth, SolarAuxUtils::GetColor("yellow"));

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
    std::vector<art::Ptr<recob::OpHit>> OpHitList;
    art::Handle<std::vector<recob::OpHit>> OpHitHandle;
    std::vector<std::vector<art::Ptr<recob::OpHit>>> OpHitVec;
    std::vector<std::vector<int>> OpHitIdx;
    if (evt.getByLabel(fOpHitLabel, OpHitHandle))
    {
      art::fill_ptr_vector(OpHitList, OpHitHandle);
    }
    OpHitNum = int(OpHitList.size());
    if (fGenerateAdjOpFlash)
    {
      fOpFlashLabel = "solarflash";
      std::vector<AdjOpHitsUtils::FlashInfo> FlashVec;
      adjophits->CalcAdjOpHits(OpHitList, OpHitVec, OpHitIdx);
      adjophits->MakeFlashVector(FlashVec, OpHitVec, evt);
      OpFlashNum = int(FlashVec.size());
      for (int i = 0; i < int(FlashVec.size()); i++)
      {
        // if (FlashVec[i].MaxPE / FlashVec[i].PE < fAdjOpFlashMaxPERatioCut && FlashVec[i].PE < fAdjOpFlashMinPECut)
        // {
        //   continue;
        // }
        AdjOpHitsUtils::FlashInfo TheFlash = FlashVec[i];
        double ThisOpFlashPur = 0;
        OpFlashMaxPE.push_back(TheFlash.MaxPE);
        OpFlashSTD.push_back(TheFlash.STD);
        OpFlashFast.push_back(TheFlash.FastToTotal);
        OpFlashID.push_back(i);
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
          const std::vector<int> ThisOpHitTrackIds = pbt->OpHitToTrackIds(OpHit);
          float ThisOphitPurity = 0;
          for (auto const &ThisOpHitTrackId : ThisOpHitTrackIds)
          {
            if (SignalTrackIDs.find(ThisOpHitTrackId) != SignalTrackIDs.end())
            {
              ThisOphitPurity += 1;
            }
          }
          // Check if ThisOpHitTrackIds is empty
          if (ThisOpHitTrackIds.size() == 0)
          {
            ThisOphitPurity = 0;
          }
          else
          {
            ThisOphitPurity /= int(ThisOpHitTrackIds.size());
          }
          ThisOpFlashPur += ThisOphitPurity * OpHit.PE();
          auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(OpHit.OpChannel()).GetCenter();
          SOpHitPur.push_back(ThisOphitPurity);
          SOpHitChannel.push_back(OpHit.OpChannel());
          SOpHitT.push_back(OpHit.PeakTime());
          SOpHitPE.push_back(OpHit.PE());
          SOpHitX.push_back(OpHitXYZ.X());
          SOpHitY.push_back(OpHitXYZ.Y());
          SOpHitZ.push_back(OpHitXYZ.Z());
          SOpHitFlashID.push_back(i);
        }
        // Check if OpHitVec[i] is empty
        if (OpHitVec[i].size() == 0)
        {
          ThisOpFlashPur = 0;
        }
        else
        {
          ThisOpFlashPur /= TheFlash.PE;
        }
        OpFlashPur.push_back(ThisOpFlashPur);
        if (abs(TheFlash.Time) < 5)
        {
          mf::LogDebug("SolarNuAna") << "Signal OpFlash PE (fast/ratio/tot/STD) " << TheFlash.FastToTotal << "/" << TheFlash.MaxPE / TheFlash.PE << "/" << TheFlash.PE << "/" << TheFlash.STD << " with purity " << ThisOpFlashPur << " time " << TheFlash.Time;
          sOpFlashTruth += "OpFlash PE " + SolarAuxUtils::str(TheFlash.PE) + " with purity " + SolarAuxUtils::str(ThisOpFlashPur) + " time " + SolarAuxUtils::str(TheFlash.Time) + " vertex (" + SolarAuxUtils::str(TheFlash.X) + ", " + SolarAuxUtils::str(TheFlash.Y) + ", " + SolarAuxUtils::str(TheFlash.Z) + ")\n";
          sOpFlashTruth += "\t*** 1st Sanity check: Ratio " + SolarAuxUtils::str(TheFlash.MaxPE / TheFlash.PE) + " <= " + SolarAuxUtils::str(fAdjOpFlashMaxPERatioCut) + "\n";
          sOpFlashTruth += "\t*** 2nd Sanity check: #OpHits " + SolarAuxUtils::str(int(OpHitVec[i].size())) + " >= " + SolarAuxUtils::str(TheFlash.NHit) + "\n";
        }
      }
    }
    else
    {
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
        int NMatchedHits = MatchedHits.size();
        double FlashStdDev = 0.0, varY = 0.0, varZ = 0.0, TotalFlashPE = 0, MaxOpHitPE = 0, FlashTime = 0;

        for (int j = 0; j < NMatchedHits; j++)
        { // Loop over OpHits in the flash
          recob::OpHit OpHit = *MatchedHits[j];
          mf::LogDebug("SolarNuAna") << "Assigning OpHit to Flash";
          const std::vector<int> ThisOpHitTrackIds = pbt->OpHitToTrackIds(OpHit);
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
          varY += pow(TheFlash.YCenter() - OpHitXYZ.Y(), 2) * OpHit.PE();
          varZ += pow(TheFlash.ZCenter() - OpHitXYZ.Z(), 2) * OpHit.PE();
          FlashTime += OpHit.PeakTime() * OpHit.PE();
          // if (OpHit.PE() > fOpFlashAlgoPE)
          // {
          SOpHitPur.push_back(ThisOphitPurity / int(ThisOpHitTrackIds.size()));
          if (OpHit.PE() > MaxOpHitPE)
          {
            MaxOpHitPE = OpHit.PE();
          };
          SOpHitChannel.push_back(OpHit.OpChannel());
          SOpHitT.push_back(OpHit.PeakTime());
          SOpHitPE.push_back(OpHit.PE());
          SOpHitX.push_back(OpHitXYZ.X());
          SOpHitY.push_back(OpHitXYZ.Y());
          SOpHitZ.push_back(OpHitXYZ.Z());
          SOpHitFlashID.push_back(i);
          // }
        } // End of OpHit loop

        OpHitVec.push_back(MatchedHits);
        varY = varY / TotalFlashPE;
        varZ = varZ / TotalFlashPE;
        FlashTime = FlashTime / TotalFlashPE;
        FlashStdDev = sqrt(varY + varZ);

        mf::LogDebug("SolarNuAna") << "Evaluating Flash purity";
        int TerminalOutput = SolarAuxUtils::supress_stdout();
        double ThisOpFlashPur = pbt->OpHitCollectionPurity(SignalTrackIDs, MatchedHits);
        SolarAuxUtils::resume_stdout(TerminalOutput);
        mf::LogDebug("SolarNuAna") << "PE of this OpFlash " << TotalFlashPE << " OpFlash time " << FlashTime;

        // Calculate the flash purity, only for the Signal events
        // if (MaxOpHitPE / TotalFlashPE < fAdjOpFlashMaxPERatioCut && TotalFlashPE > fAdjOpFlashMinPECut)
        // {
        OpFlashID.push_back(i);
        OpFlashPur.push_back(ThisOpFlashPur);
        OpFlashMaxPE.push_back(MaxOpHitPE);
        OpFlashSTD.push_back(FlashStdDev);
        OpFlashT.push_back(TheFlash.Time());
        OpFlashX.push_back(TheFlash.XCenter());
        OpFlashY.push_back(TheFlash.YCenter());
        OpFlashZ.push_back(TheFlash.ZCenter());
        OpFlashPE.push_back(TheFlash.TotalPE());
        OpFlashFast.push_back(TheFlash.FastToTotal());
        OpFlashDeltaT.push_back(TheFlash.TimeWidth());
        OpFlashNHit.push_back(MatchedHits.size());
        // }
        if (abs(TheFlash.Time()) < 5)
        {
          mf::LogDebug("SolarNuAna") << "OpFlash PE " << TheFlash.TotalPE() << " with purity " << ThisOpFlashPur << " time " << TheFlash.Time();
          sOpFlashTruth += "OpFlash PE " + SolarAuxUtils::str(TheFlash.TotalPE()) + " with purity " + SolarAuxUtils::str(ThisOpFlashPur) + " time " + SolarAuxUtils::str(TheFlash.Time()) + " vertex (" + SolarAuxUtils::str(TheFlash.YCenter()) + ", " + SolarAuxUtils::str(TheFlash.ZCenter()) + ")\n";
          sOpFlashTruth += "\t*** 1st Sanity check: Ratio " + SolarAuxUtils::str(MaxOpHitPE / TotalFlashPE) + " <= " + SolarAuxUtils::str(fAdjOpFlashMaxPERatioCut) + "\n";
          sOpFlashTruth += "\t*** 2nd Sanity check: #OpHits " + SolarAuxUtils::str(int(NMatchedHits)) + " >= " + SolarAuxUtils::str(int(TheFlash.PEs().size())) + "\n";
        }
      }
    }
    solaraux->PrintInColor(sOpFlashTruth, SolarAuxUtils::GetColor("blue"));

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
        solaraux->PrintInColor("Negative Hit Time = " + SolarAuxUtils::str(ThisHit.PeakTime()), SolarAuxUtils::GetColor("red"));
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
    lowe->CalcAdjHits(ColHits0, Clusters0, hAdjHits, hAdjHitsADCInt, false);
    HitNum.push_back(ColHits0.size());
    ClusterNum.push_back(Clusters0.size());
    lowe->CalcAdjHits(ColHits1, Clusters1, hAdjHits, hAdjHitsADCInt, false);
    HitNum.push_back(ColHits1.size());
    ClusterNum.push_back(Clusters1.size());
    lowe->CalcAdjHits(ColHits2, Clusters2, hAdjHits, hAdjHitsADCInt, false);
    HitNum.push_back(ColHits2.size());
    ClusterNum.push_back(Clusters2.size());
    lowe->CalcAdjHits(ColHits3, Clusters3, hAdjHits, hAdjHitsADCInt, false);
    HitNum.push_back(ColHits3.size());
    ClusterNum.push_back(Clusters3.size());
    fMCTruthTree->Fill();

    std::vector<std::vector<std::vector<float>>> ClVecGenPur = {{}, {}, {}};
    std::vector<std::vector<std::vector<recob::Hit>>> AllPlaneClusters = {Clusters0, Clusters1, Clusters2};
    std::vector<std::vector<int>> ClMainID = {{}, {}, {}}, ClTPC = {{}, {}, {}}, ClNHits = {{}, {}, {}}, ClGen = {{}, {}, {}};
    std::vector<std::vector<float>> ClCharge = {{}, {}, {}}, ClMaxCharge = {{}, {}, {}}, ClT = {{}, {}, {}}, ClX = {{}, {}, {}}, ClY = {{}, {}, {}}, ClZ = {{}, {}, {}};
    std::vector<std::vector<float>> ClFracE = {{}, {}, {}}, ClFracGa = {{}, {}, {}}, ClFracNe = {{}, {}, {}}, ClFracRest = {{}, {}, {}};
    std::vector<std::vector<float>> ClPur = {{}, {}, {}}, Cldzdy = {{}, {}, {}}, ClGenPur = {{}, {}, {}};

    std::string sRecoObjects = "";
    sRecoObjects += "\n# OpHits (" + fOpHitLabel + ") in full geometry: " + SolarAuxUtils::str(OpHitNum);
    sRecoObjects += "\n# OpFlashes (" + fOpFlashLabel + ") in full geometry: " + SolarAuxUtils::str(OpFlashNum);
    sRecoObjects += "\n# Hits (" + fHitLabel + ") in each view: " + SolarAuxUtils::str(int(ColHits0.size())) + ", " + SolarAuxUtils::str(int(ColHits1.size())) + ", " + SolarAuxUtils::str(int(ColHits2.size())) + ", " + SolarAuxUtils::str(int(ColHits3.size()));
    sRecoObjects += "\n# Cluster from the hits: " + SolarAuxUtils::str(int(Clusters0.size())) + ", " + SolarAuxUtils::str(int(Clusters1.size())) + ", " + SolarAuxUtils::str(int(Clusters2.size())) + ", " + SolarAuxUtils::str(int(Clusters3.size()));
    sRecoObjects += "\n# Tracks (" + fTrackLabel + ") in full geometry: " + SolarAuxUtils::str(TrackNum);
    solaraux->PrintInColor(sRecoObjects, SolarAuxUtils::GetColor("cyan"));

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
        int MainTrID;
        int MainGenerator = 0;
        float Pur = 0;
        std::vector<float> thisdzdy = {};

        nhit = Clusters[i].size();
        ncharge = maxHit = clustT = FracE = FracGa = FracNe = FracRest = clustX = clustY = clustZ = clustTPC = dzdy = 0;
        // Define a vector of floats with size equal to the number of generators + 1
        std::vector<float> VecGenPur(fLabels.size() + 1, 0);

        for (recob::Hit TPCHit : Clusters[i])
        {
          if (TPCHit.PeakTime() < 0)
            solaraux->PrintInColor("Negative Cluster Time = " + SolarAuxUtils::str(TPCHit.PeakTime()), SolarAuxUtils::GetColor("red"));
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
          double TopEFrac = 0;
          std::vector<sim::TrackIDE> ThisHitIDE = bt_serv->HitToTrackIDEs(clockData, TPCHit);

          for (size_t ideL = 0; ideL < ThisHitIDE.size(); ++ideL)
          {
            if (ThisHitIDE[ideL].energyFrac > TopEFrac)
            {
              TopEFrac = ThisHitIDE[ideL].energyFrac;
              MainTrID = abs(ThisHitIDE[ideL].trackID);
              mf::LogDebug("SolarNuAna") << "This TPCHit's IDE is: " << MainTrID;
            }
          }

          for (int frac = 0; frac < int(ClPartTrackIDs.size()); ++frac)
          {
            for (int trck = 0; trck < int(ClPartTrackIDs[frac].size()); ++trck)
            {
              if (MainTrID == ClPartTrackIDs[frac][trck])
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

          long unsigned int GeneratorType = SolarAuxUtils::WhichGeneratorType(GeneratorParticles, MainTrID);
          VecGenPur[int(GeneratorType)] = VecGenPur[int(GeneratorType)] + TPCHit.Integral();
          mf::LogDebug("SolarNuAna") << "\nThis particle type " << GeneratorType << "\nThis cluster's main track ID " << MainTrID;
          if (SignalTrackIDs.find(MainTrID) != SignalTrackIDs.end())
          {
            hitCharge = TPCHit.Integral();
            Pur = Pur + hitCharge;
          }
        }
        // std::cout << SolarAuxUtils::str(VecGenPur) << "\n"
        //           << std::endl;
        float MainGenPurity = 0;
        for (size_t genpur = 0; genpur < VecGenPur.size(); genpur++)
        {
          VecGenPur[genpur] = VecGenPur[genpur] / ncharge;
          if (VecGenPur[genpur] > MainGenPurity)
          {
            MainGenerator = genpur;
            MainGenPurity = VecGenPur[genpur];
          }
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
        Pur /= ncharge;
        mf::LogDebug("SolarNuAna") << "\ndzdy " << dzdy << " for cluster "
                                   << " (" << clustY << ", " << clustZ << ") with track ID " << MainTrID << " in plane " << idx;
        if (clustT < 0)
          solaraux->PrintInColor("Negative Cluster Time = " + SolarAuxUtils::str(clustT), SolarAuxUtils::GetColor("red"));

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
        ClPur[idx].push_back(Pur);
        ClGen[idx].push_back(MainGenerator);
        ClGenPur[idx].push_back(MainGenPurity);
        Cldzdy[idx].push_back(dzdy);
        ClMainID[idx].push_back(MainTrID);
        ClVecGenPur[idx].push_back(VecGenPur);

        mf::LogDebug("SolarNuAna") << "\nCluster " << i << " in plane " << idx << " has #hits" << nhit << " charge, " << ncharge << " time, " << clustT;
        mf::LogDebug("SolarNuAna") << " and position (" << clustY << ", " << clustZ << ") with main track ID " << MainTrID << " and purity " << Pur;
      }
    } // Finished first cluster processing

    //-------------------------------------------------------------------- Cluster Matching -------------------------------------------------------------------------//
    std::vector<std::vector<float>> MVecGenFrac = {};
    std::vector<unsigned int> MVecGen = {};
    std::vector<int> MVecNHit = {}, MVecInd0NHits = {}, MVecInd1NHits = {}, MVecMainID = {}, MVecTPC = {}, MVecInd0TPC = {}, MVecInd1TPC = {};
    std::vector<float> MVecTime = {}, MVecCharge = {}, MVecMaxCharge = {}, MVecInd0Charge = {}, MVecInd1Charge = {}, MVecInd0MaxCharge = {}, MVecInd1MaxCharge = {}, MVecInd0dT = {}, MVecInd1dT = {};
    std::vector<float> MVecInd0RecoY = {}, MVecInd1RecoY = {}, MVecRecY = {}, MVecRecZ = {};
    std::vector<float> MVecFracE = {}, MVecFracGa = {}, MVecFracNe = {}, MVecFracRest = {}, MVecPur = {}, MVecGenPur = {};

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
            if (abs(ClNHits[0][jj] - ClNHits[2][ii]) / ClNHits[2][ii] > fClusterMatchNHit || abs(ClCharge[0][jj] - ClCharge[2][ii]) / ClCharge[2][ii] > fClusterMatchCharge)
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
            if (abs(ClNHits[1][zz] - ClNHits[2][ii]) / ClNHits[2][ii] > fClusterMatchNHit || abs(ClCharge[1][zz] - ClCharge[2][ii]) / ClCharge[2][ii] > fClusterMatchCharge)
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
        // Cluster Signal Fractions
        MVecFracE.push_back(ClFracE[2][ii]);
        MVecFracGa.push_back(ClFracGa[2][ii]);
        MVecFracNe.push_back(ClFracNe[2][ii]);
        MVecFracRest.push_back(ClFracRest[2][ii]);
        // Cluster Signal Purity
        MVecPur.push_back(ClPur[2][ii]);
        // Cluster Gen and GenFraction
        MVecMainID.push_back(ClMainID[2][ii]);
        MVecGen.push_back(ClGen[2][ii]);
        MVecGenPur.push_back(ClGenPur[2][ii]);
        MVecGenFrac.push_back(ClVecGenPur[2][ii]);

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
            mf::LogWarning("SolarNuAna") << "Signal cluster reconstructed outside of detector volume! RecoY = " << SolarAuxUtils::str((ind0clustY + ind1clustY) / 2);
          }
        }

        // Print in color if the cluster is matched
        mf::LogDebug("SolarNuAna") << "¡¡¡ Matched cluster !!! ";
        mf::LogDebug("SolarNuAna") << " - Cluster " << SolarAuxUtils::str(ClMainID[2][ii]) << " Gen " << SolarAuxUtils::str(ClGen[2][ii]) << " Purity " << SolarAuxUtils::str(ClGenPur[2][ii]) << " Hits " << SolarAuxUtils::str(ClNHits[2][ii]);
        mf::LogDebug("SolarNuAna") << " - Hits(ind0, ind1, col) " << SolarAuxUtils::str(ind0clustNHits) << ", " << SolarAuxUtils::str(ind1clustNHits) << ", " << SolarAuxUtils::str(ClNHits[2][ii]);
        mf::LogDebug("SolarNuAna") << " - Positions y(ind0, ind1) = " << SolarAuxUtils::str(ind0clustY) << ", " << SolarAuxUtils::str(ind1clustY) << ", z = " << SolarAuxUtils::str(ClZ[2][ii]) << "\n";
      } // if (match == true)
    } // Loop over collection plane clusters

    //-------------------------------------------------------------------- Cluster Tree Export -------------------------------------------------------------------------//
    // Loop over matched clusters and export to tree if number of hits is above threshold
    for (int i = 0; i < int(MVecNHit.size()); i++)
    {
      bool TrackMatch = false;
      std::string sFlashReco = "";
      std::string sClusterReco = "";
      std::string sResultColor = "white";
      float OpFlashResidual = 0;
      float MatchedOpFlashPE = -1e6;
      float MatchedOpFlashResidual = 1e6;
      float MatchedOpFlashX = -1e6;

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
        MAdjClGenPur = {};
        MAdjClMainID = {};
        MAdjClMainPDG = {};
        MAdjClMainE = {};
        MAdjClMainX = {};
        MAdjClMainY = {};
        MAdjClMainZ = {};
        MAdjClEndX = {};
        MAdjClEndY = {};
        MAdjClEndZ = {};
        MAdjFlashR = {};
        MAdjFlashPE = {};
        MAdjFlashPur = {};
        MAdjFlashSTD = {};
        MAdjFlashTime = {};
        MAdjFlashNHit = {};
        MAdjFlashFast = {};
        MAdjFlashMaxPE = {};
        MAdjFlashRecoX = {};
        MAdjFlashRecoY = {};
        MAdjFlashRecoZ = {};
        MAdjFlashResidual = {};
        MTrackStart = {-1e6, -1e6, -1e6};
        MTrackEnd = {-1e6, -1e6, -1e6};

        for (int j = 0; j < int(MVecNHit.size()); j++)
        {
          if (j == i)
          {
            continue;
          }

          double ClusterDistance = 0;
          solaraux->ComputeDistance3D(ClusterDistance, MVecTime[i], MVecRecY[i], MVecRecZ[i], MVecTime[j], MVecRecY[j], MVecRecZ[j]);
          if (ClusterDistance > fAdjClusterRad)
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
          MAdjClGenPur.push_back(MVecGenPur[j]);
          MAdjClMainID.push_back(MVecMainID[j]);

          // If mother exists add the mother information
          const simb::MCParticle *MAdjClTruth;
          int TerminalOutput = SolarAuxUtils::supress_stdout();
          MAdjClTruth = pi_serv->TrackIdToParticle_P(MVecMainID[j]);
          SolarAuxUtils::resume_stdout(TerminalOutput);
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

        sResultColor = "yellow";
        if (MVecPur[i] > 0.7)
        {
          sResultColor = "green";
        }

        if (MPrimary)
        {
          sClusterReco += "*** Matched preselection cluster: \n";
          sClusterReco += " - MainTrackID " + SolarAuxUtils::str(MVecMainID[i]) + "\n";
          sClusterReco += " - Purity " + SolarAuxUtils::str(MVecGenPur[i]) + " Hits " + SolarAuxUtils::str(MVecNHit[i]) + "\n";
          if (MVecGen[i] > 0 && int(MVecGen[i]) < (int(fLabels.size()) + 1))
          {
            sClusterReco += " - Gen " + SolarAuxUtils::str(int(MVecGen[i])) + " -> " + fLabels[MVecGen[i] - 1] + "\n";
          }
          else
          {
            sClusterReco += " - Gen ?? -> Unknown\n";
          }
          sClusterReco += " - RecoY, RecoZ (" + SolarAuxUtils::str(MVecRecY[i]) + ", " + SolarAuxUtils::str(MVecRecZ[i]) + ") Time " + SolarAuxUtils::str(MVecTime[i]) + "\n\n";
          TVector3 ThisClVertex = {0, MVecRecY[i], MVecRecZ[i]};
          float MaxVertexDistance = 10; // if track is further away from ThisClVertex than
          for (int i = 0; i < TrackNum; i++)
          { // using index loop to get track idx
            recob::Track trk = *TrackList[i];
            TVector3 trk_start(0, trk.Start().Y(), trk.Start().Z());
            TVector3 trk_end(0, trk.End().Y(), trk.End().Z());
            // throw away bad tracks
            if ((trk_start - ThisClVertex).Mag() > MaxVertexDistance && (trk_end - ThisClVertex).Mag() > MaxVertexDistance)
            {
              continue;
            }
            MTrackNPoints = trk.NPoints();
            MTrackStart = {trk.Start().X(), trk.Start().Y(), trk.Start().Z()};
            MTrackEnd = {trk.End().X(), trk.End().Y(), trk.End().Z()};
            MTrackChi2 = trk.Chi2();
            sClusterReco += "*** Matched pmtrack: \n";
            sClusterReco += " - Track has start (" + SolarAuxUtils::str(trk.Start().X()) + ", " + SolarAuxUtils::str(trk.Start().Y()) + ", " + SolarAuxUtils::str(trk.Start().Z()) + ")\n";
            sClusterReco += " - Track has end   (" + SolarAuxUtils::str(trk.End().X()) + ", " + SolarAuxUtils::str(trk.End().Y()) + ", " + SolarAuxUtils::str(trk.End().Z()) + ")\n";
            TrackMatch = true;
          }; // Loop over tracks
        };

        for (int j = 0; j < int(OpFlashPE.size()); j++)
        {
          // Skip flashes with time outside the cluster time window
          if ((MVecTime[i] - OpFlashT[j]) < 0 || (MVecTime[i] - OpFlashT[j]) > fAdjOpFlashTime)
          {
            continue;
          }
          // Make an eliptical cut on the flash position
          if (pow(MVecRecY[i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) + pow(MVecRecZ[i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1)
          {
            continue;
          }
          // Compute the time distance between the cluster and the flash. Use factor 2 to convert us to TPC tics
          double MAdjFlashX = 0;
          solaraux->ComputeDistanceX(MAdjFlashX, MVecTime[i], 2 * OpFlashT[j]);
          // For HD 1x2x6 (only 1 APA) the x coordinate is determined by the TPC number
          if (MVecTPC[i]%2 != 0)
          {
            MAdjFlashX = -MAdjFlashX;
          }
          float OpFlashR = sqrt(pow(MVecRecY[i] - OpFlashY[j], 2) + pow(MVecRecZ[i] - OpFlashZ[j], 2));
          MAdjFlashTime.push_back(OpFlashT[j]);
          MAdjFlashPE.push_back(OpFlashPE[j]);
          MAdjFlashNHit.push_back(OpFlashNHit[j]);
          MAdjFlashMaxPE.push_back(OpFlashMaxPE[j]);
          MAdjFlashSTD.push_back(OpFlashSTD[j]);
          MAdjFlashFast.push_back(OpFlashFast[j]);
          MAdjFlashRecoX.push_back(OpFlashX[j]);
          MAdjFlashRecoY.push_back(OpFlashY[j]);
          MAdjFlashRecoZ.push_back(OpFlashZ[j]);
          MAdjFlashR.push_back(OpFlashR);
          MAdjFlashPur.push_back(OpFlashPur[j]);
          // Compute the residual between the predicted cluster signal and the flash
          std::string sFlashMatching = "Testing flash " + SolarAuxUtils::str(j) + " with time " + SolarAuxUtils::str(OpFlashT[j]) + " and PE " + SolarAuxUtils::str(OpFlashPE[j]);
          solaraux->PrintInColor(sFlashMatching, SolarAuxUtils::GetColor(sResultColor), "Debug");
          adjophits->FlashMatchResidual(OpFlashResidual, OpHitVec[j], MAdjFlashX, double(MVecRecY[i]), double(MVecRecZ[i]));
          // Make a cut on the flash PE and MaxPE distributions
          if (OpFlashMaxPE[j] / OpFlashPE[j] > fAdjOpFlashMaxPERatioCut || OpFlashPE[j] < fAdjOpFlashMinPECut)
          {
            continue;
          }
          // Create a cut based on a sigmoid function of the opflash number of hits versus the associated drift time
          if (OpFlashNHit[j] < fAdjOpFlashMinNHitCut)
          {
            continue;
          }
          // If the residual is smaller than the minimum residual, update the minimum residual and the matched flash
          if ((fFlashMatchByResidual && OpFlashResidual < MatchedOpFlashResidual) || (!fFlashMatchByResidual && OpFlashPE[j] > MatchedOpFlashPE))
          {
            MFlashR = OpFlashR;
            MFlashPE = OpFlashPE[j];
            MFlashFast = OpFlashFast[j];
            MFlashNHit = OpFlashNHit[j];
            MFlashMaxPE = OpFlashMaxPE[j];
            MFlashPur = OpFlashPur[j];
            MFlashSTD = OpFlashSTD[j];
            MFlashTime = OpFlashT[j];
            MFlashRecoX = OpFlashX[j];
            MFlashRecoY = OpFlashY[j];
            MFlashRecoZ = OpFlashZ[j];
            MFlashResidual = OpFlashResidual;

            // Create an output string with the flash information
            sFlashReco = "*** Matched flash: \n - Purity " + SolarAuxUtils::str(OpFlashPur[j]) +
              " #Hits " + SolarAuxUtils::str(OpFlashNHit[j]) +
              " PE " + SolarAuxUtils::str(OpFlashPE[j]) +
              " MaxPE " + SolarAuxUtils::str(OpFlashMaxPE[j]) + "\n" +
              " Reco X,Y,Z (" + SolarAuxUtils::str(MAdjFlashX) + ", " + SolarAuxUtils::str(OpFlashY[j]) + ", " + SolarAuxUtils::str(OpFlashZ[j]) + ")" + "\n" +
              " Time " + SolarAuxUtils::str(OpFlashT[j]) +
              " Fast " + SolarAuxUtils::str(OpFlashFast[j]) +
              " Residual " + SolarAuxUtils::str(OpFlashResidual) + "\n";

            MatchedOpFlashX = MAdjFlashX;
            MatchedOpFlashResidual = OpFlashResidual;
            MatchedOpFlashPE = MFlashPE;
          }
          MAdjFlashResidual.push_back(OpFlashResidual);
        }
        sClusterReco += sFlashReco;

        // Fill the tree with the cluster information and the adjacent clusters and flashes
        MPur = MVecPur[i];
        MGen = MVecGen[i];
        MGenPur = MVecGenPur[i];
        MGenFrac = MVecGenFrac[i];
        MSignalFrac = {MVecFracE[i], MVecFracGa[i], MVecFracNe[i], MVecFracRest[i]};
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
        // Cluster RecoX
        MRecX = MatchedOpFlashX;
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
        int TerminalOutput = SolarAuxUtils::supress_stdout();
        MClTruth = pi_serv->TrackIdToParticle_P(MVecMainID[i]);
        SolarAuxUtils::resume_stdout(TerminalOutput);
        if (MClTruth == 0)
        {
          MMainVertex = {-1e6, -1e6, -1e6};
          MEndVertex = {-1e6, -1e6, -1e6};
          MMainPDG = 0;
          MMainE = -1e6;
          MMainT = -1e6;
          MMainP = -1e6;
          MMainK = -1e6;

          MMainParentVertex = {-1e6, -1e6, -1e6};
          MMainParentPDG = 0;
          MMainParentE = -1e6;
          MMainParentT = -1e6;
          MMainParentP = -1e6;
          MMainParentK = -1e6;
        }
        else
        {
          if (abs(MClTruth->Vx() - MatchedOpFlashX) < 10)
          {
            MFlashCorrect = true;
          };
          MMainVertex = {MClTruth->Vx(), MClTruth->Vy(), MClTruth->Vz()};
          MEndVertex = {MClTruth->EndX(), MClTruth->EndY(), MClTruth->EndZ()};
          MMainPDG = MClTruth->PdgCode();
          MMainE = 1e3 * MClTruth->E();
          MMainT = MClTruth->T();
          MMainP = 1e3 * MClTruth->P();
          MMainK = MMainE - 1e3 * MClTruth->Mass();
          // If exists add the parent information
          const simb::MCParticle *MClParentTruth;
          int TerminalOutput = SolarAuxUtils::supress_stdout();
          MClParentTruth = pi_serv->TrackIdToParticle_P(MClTruth->Mother());
          SolarAuxUtils::resume_stdout(TerminalOutput);
          if (MClParentTruth == 0)
          {
            MMainParentVertex = {-1e6, -1e6, -1e6};
            MMainParentPDG = 0;
            MMainParentE = -1e6;
            MMainParentT = -1e6;
            MMainParentP = -1e6;
            MMainParentK = -1e6;
          }
          else
          {
            MMainParentVertex = {MClParentTruth->Vx(), MClParentTruth->Vy(), MClParentTruth->Vz()};
            MMainParentPDG = MClParentTruth->PdgCode();
            MMainParentE = 1e3 * MClParentTruth->E();
            MMainParentT = MClParentTruth->T();
            MMainParentP = 1e3 * MClParentTruth->P();
            MMainParentK = MMainParentE - 1e3 * MClParentTruth->Mass();
          }
        }
        if ((fClusterPreselectionPrimary && !MPrimary) || (fClusterPreselectionTrack && !TrackMatch))
        {
          continue;
        }
        fSolarNuAnaTree->Fill();
        hDriftTime->Fill(MainElectronEndPointX, MTime);
        hXTruth->Fill(MVecRecY[i] - TNuY, TNuX);
        hYTruth->Fill(MVecRecY[i] - TNuY, TNuY);
        hZTruth->Fill(MVecRecZ[i] - TNuZ, TNuZ);
        if (MVecTime[i] < 0)
          mf::LogWarning("SolarNuAna") << "Negative Main Cluster Time = " << MVecTime[i];
      }
      // Check if the string sClusterReco is not empty and print it in color
      if (sClusterReco != "")
      {
        solaraux->PrintInColor(sClusterReco, SolarAuxUtils::GetColor(sResultColor));
      }
    }
  }

  // ########################################################################################################################################//
  // _FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_//
  // ########################################################################################################################################//

  //......................................................
  // Reset variables for each event
  void SolarNuAna::ResetVariables()
  {
    TNuE = 0;
    TNuX = 0;
    TNuY = 0;
    TNuZ = 0;
    TrackNum = 0;
    OpHitNum = 0;
    OpFlashNum = 0;
    MTrackNPoints = 0;
    MTrackChi2 = 0;
    MFlashPE = -1e6;
    MFlashFast = -1e6;
    MFlashNHit = -1e6;
    MFlashMaxPE = -1e6;
    MFlashPur = -1e6;
    MFlashSTD = -1e6;
    MFlashTime = -1e6;
    MFlashR = -1e6;
    MFlashRecoX = -1e6;
    MFlashRecoY = -1e6;
    MFlashRecoZ = -1e6;
    MFlashResidual = -1e6;
    MFlashCorrect = false;
    SignalElectronDepList = {};
    SignalPDGList = {};
    SignalPDGDepList = {};
    SignalIDList = {}, SignalIDDepList = {};
    SignalMotherList = {};
    SignalEList = {}, SignalPList = {}, SignalKList = {}, SignalTList = {}, SignalEndXList = {}, SignalEndYList = {}, SignalEndZList = {};
    SignalEDepList = {}, SignalXDepList = {}, SignalYDepList = {}, SignalZDepList = {};
    SignalMaxEDepList = {}, SignalMaxEDepXList = {}, SignalMaxEDepYList = {}, SignalMaxEDepZList = {};
    SOpHitChannel = {}, SOpHitPur = {}, SOpHitPE = {}, SOpHitX = {}, SOpHitY = {}, SOpHitZ = {}, SOpHitT = {}, SOpHitFlashID = {};
    TPart = {}, GeneratorParticles = {};
    HitNum = {};
    ClusterNum = {};
    OpFlashPur.clear();
    OpFlashID.clear();
    OpFlashPE.clear();
    OpFlashSTD.clear();
    OpFlashFast.clear();
    OpFlashMaxPE.clear();
    OpFlashX.clear();
    OpFlashY.clear();
    OpFlashZ.clear();
    OpFlashT.clear();
    OpFlashDeltaT.clear();
    OpFlashNHit.clear();
  }
} // namespace solar
DEFINE_ART_MODULE(solar::SolarNuAna)
