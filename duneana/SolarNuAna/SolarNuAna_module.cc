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

#include "duneopdet/LowEPDSUtils/AdjOpHitsUtils.h"
#include "dunecore/ProducerUtils/ProducerUtils.h"
#include "dunereco/LowEUtils/LowEUtils.h"

using namespace producer;

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
    std::vector<std::string> fLabels, fBackgroundLabels;
    std::string fSignalLabel, fClusterLabel, fSolarClusterLabel, fClusterChargeVariable, fOpHitTimeVariable;
    int fDetectorSizeY, fDetectorSizeZ, fClusterAlgoAdjChannel, fClusterInd0MatchTime, fClusterInd1MatchTime, fClusterPreselectionNHits, fAdjOpFlashMinNHitCut;
    float fClusterMatchTime, fAdjClusterRad, fMinClusterCharge, fClusterMatchCharge, fAdjOpFlashX, fAdjOpFlashY, fAdjOpFlashZ, fAdjOpFlashMaxPERatioCut, fAdjOpFlashMinPECut, fClusterMatchNHit, fClusterAlgoTime;
    float fOpFlashTimeOffset, fOpFlashAlgoMinTime, fOpFlashAlgoMaxTime, fOpFlashAlgoRad, fOpFlashAlgoPE, fOpFlashAlgoTriggerPE, fOpFlashAlgoHotVertexThld, fXACathodeX, fXAMembraneY, fXAStartCapZ, fXAFinalCapZ;
    bool fClusterPreselectionSignal, fClusterPreselectionPrimary, fClusterPreselectionTrack, fClusterPreselectionFlashMatch;
    bool fGenerateSolarCluster, fGenerateAdjCluster, fGenerateAdjOpFlash, fFlashMatchByResidual;
    bool fSaveSignalDaughters, fSaveSignalEDep, fSaveSignalOpHits, fSaveOpFlashInfo, fSaveTrackInfo;
    bool fAdjOpFlashMembraneProjection, fAdjOpFlashEndCapProjection; // If true, the TPC reco is projected to the membrane plane. If false, apply a 3D constraint dT, Y, Z.
    bool fOpFlashTime2us; // If true, the OpFlash time is in ticks, and we convert it to microseconds.

    // --- Our TTrees, and its associated variables.
    TTree *fConfigTree;
    TTree *fMCTruthTree;
    TTree *fSolarNuAnaTree;
    std::string sInteraction;
    std::vector<std::map<int, simb::MCParticle>> GeneratorParticles = {};
    int Event, Flag, MNHit, MGen, MTPC, MInd0TPC, MInd1TPC, MInd0NHits, MInd1NHits, MMainID, MMainPDG, MMainParentPDG, TrackNum, OpHitNum, OpFlashNum, MTrackNPoints, MAdjClNum, MSignalAdjClNum, SignalParticlePDG;
    float SignalParticleE, SignalParticleP, SignalParticleK, SignalParticleX, SignalParticleY, SignalParticleZ, SignalParticleTime, MTime, MCharge, MMaxCharge, MInd0Charge, MInd1Charge, MInd0MaxCharge, MInd1MaxCharge;
    float MInd0dTime, MInd1dTime, MInd0RecoY, MInd1RecoY, MRecX, MRecY, MRecZ, MPur, MInd0Pur, MInd1Pur, MGenPur, MMainE, MMainP, MMainK, MMainTime, MMainParentE, MMainParentP, MMainParentK, MMainParentTime, MTrackChi2;
    std::vector<int> MAdjClGen, MAdjClMainID, TPart, SignalPDGList, SignalPDGDepList, SignalIDList, SignalMotherList, SignalIDDepList, MAdjClMainPDG, HitNum, ClusterNum, SignalElectronDepList, SOpHitPlane;
    std::vector<float> SignalEDepList, SignalXDepList, SignalYDepList, SignalZDepList;
    std::vector<float> SOpHitPur, SOpHitPE, SOpHitX, SOpHitY, SOpHitZ, SOpHitTime, SOpHitChannel, SOpHitFlashID;
    std::vector<float> MAdjClTime, MAdjClCharge, MAdjClInd0Charge, MAdjClInd1Charge, MAdjClMaxCharge, MAdjClInd0MaxCharge, MAdjClInd1MaxCharge;
    std::vector<float> MAdjClNHits, MAdjClInd0NHits, MAdjClInd1NHits, MAdjClRecoY, MAdjClRecoZ, MAdjClR, MAdjClPur, MAdjClGenPur, MAdjClMainE, MAdjClMainP, MAdjClMainK;
    std::vector<float> MAdjClMainX, MAdjClMainY, MAdjClMainZ, MAdjClEndX, MAdjClEndY, MAdjClEndZ, MSignalFrac, MGenFrac;
    std::vector<int>  MAdjFlashPlane, MAdjFlashNHits;
    std::vector<float> MAdjFlashTime, MAdjFlashResidual, MAdjFlashPE, MAdjFlashMaxPE, MAdjFlashRecoX, MAdjFlashRecoY, MAdjFlashRecoZ, MAdjFlashR, MAdjFlashPur, MAdjFlashSTD, MAdjFlashFast;
    std::vector<float> SignalEList, SignalPList, SignalKList, SignalTimeList, SignalEndXList, SignalEndYList, SignalEndZList, SignalMaxEDepList, SignalMaxEDepXList, SignalMaxEDepYList, SignalMaxEDepZList;
    std::vector<double> MMainVertex, MEndVertex, MMainParentVertex;
    std::vector<double> MTrackStart, MTrackEnd;
    bool MPrimary;

    // --- OpFlash Variables
    std::vector<int> OpFlashID, OpFlashNHits, OpFlashPlane;
    std::vector<float> OpFlashPur, OpFlashPE, OpFlashMaxPE, OpFlashX, OpFlashY, OpFlashZ, OpFlashTime, OpFlashDeltaT, OpFlashSTD, OpFlashFast;

    // --- MatchedFlash Variables
    int MFlashNHits, MFlashPlane;
    float MFlashR, MFlashPE, MFlashMaxPE, MFlashPur, MFlashFast, MFlashTime, MFlashSTD, MFlashRecoX, MFlashRecoY, MFlashRecoZ, MFlashResidual;
    bool MFlashCorrect;

    // --- Maps to hold the geo::TPCID object for each TPCid
    std::map<unsigned int, geo::TPCID> TPCIDMap;    // Key is the TPC index, value is the TPCID object
    std::map<unsigned int, float> TPCIDdriftLength; // Key is the TPC index, value is the drift length in cm
    std::map<unsigned int, float> TPCIDdriftTime;   // Key is the TPC index, value is the drift time in us

    // --- Histograms to fill about collection plane hits
    float MainElectronEndPointX;
    TH2F *hDriftTime;
    TH2F *hXTruth;
    TH2F *hYTruth;
    TH2F *hZTruth;

    // --- Declare our services
    geo::WireReadoutGeom const &wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    std::unique_ptr<producer::ProducerUtils> producer;
    std::unique_ptr<solar::AdjOpHitsUtils> adjophits;
    std::unique_ptr<lowe::LowEUtils> lowe;
  };
#endif

  //......................................................
  SolarNuAna::SolarNuAna(fhicl::ParameterSet const &p)
      : EDAnalyzer(p),
        producer(new producer::ProducerUtils(p)),
        adjophits(new solar::AdjOpHitsUtils(p)),
        lowe(new lowe::LowEUtils(p))
  {
    this->reconfigure(p);
  }

  //......................................................
  void SolarNuAna::reconfigure(fhicl::ParameterSet const &p)
  {
    fSignalLabel = p.get<std::string>("SignalLabel", "marley");
    fBackgroundLabels = p.get<std::vector<std::string>>("BackgroundLabelVector", {""});
    fGEANTLabel = p.get<std::string>("GEANT4Label", "largeant");
    fHitLabel = p.get<std::string>("HitLabel", "hitfd");
    fClusterLabel = p.get<std::string>("ClusterLabel", "planecluster");
    fSolarClusterLabel = p.get<std::string>("SolarClusterLabel", "solarcluster");
    fClusterChargeVariable = p.get<std::string>("ClusterChargeVariable", "Integral");
    fTrackLabel = p.get<std::string>("TrackLabel", "pmtrack");
    fOpHitLabel = p.get<std::string>("OpHitLabel", "ophitspe");
    fOpHitTimeVariable = p.get<std::string>("OpHitTimeVariable", "PeakTime");
    fOpFlashLabel = p.get<std::string>("OpFlashLabel", "solarflash");
    fOpFlashTime2us = p.get<bool>("OpFlashTime2us", false);      // If true, the OpFlash time is in ticks, and we convert it to microseconds.
    fOpFlashTimeOffset = p.get<float>("OpFlashTimeOffset", 0.0); // Time offset to be applied to the OpFlash time (in us if fOpFlashTime2us is true, in ticks otherwise)
    fClusterAlgoTime = p.get<float>("ClusterAlgoTime", 12.5);    // Time window (in us) to look for hits to be clustered together
    fClusterAlgoAdjChannel = p.get<int>("ClusterAlgoAdjChannel");
    fGenerateSolarCluster = p.get<bool>("GenerateSolarCluster",true);
    fClusterMatchNHit = p.get<float>("ClusterMatchNHit", 2.0);
    fClusterMatchCharge = p.get<float>("ClusterMatchCharge", 0.6);
    fClusterMatchTime = p.get<float>("ClusterMatchTime", 10.0);
    fClusterInd0MatchTime = p.get<float>("ClusterInd0MatchTime", 0);
    fClusterInd1MatchTime = p.get<float>("ClusterInd1MatchTime", 0);
    fClusterPreselectionSignal = p.get<bool>("ClusterPreselectionSignal", true);
    fClusterPreselectionPrimary = p.get<bool>("ClusterPreselectionPrimary", true);
    fClusterPreselectionNHits = p.get<int>("ClusterPreselectionNHits", 0);
    fClusterPreselectionTrack = p.get<bool>("ClusterPreselectionTrack", false);
    fClusterPreselectionFlashMatch = p.get<bool>("ClusterPreselectionFlashMatch", false);
    fGenerateAdjCluster = p.get<bool>("GenerateAdjCluster", true);
    fAdjClusterRad = p.get<float>("AdjClusterRad");
    fMinClusterCharge = p.get<float>("MinClusterCharge");
    fGenerateAdjOpFlash = p.get<bool>("GenerateAdjOpFlash");
    fXACathodeX = p.get<float>("XACathodeX");
    fXAMembraneY = p.get<float>("XAMembraneY");
    fXAStartCapZ = p.get<float>("XAStartCapZ");
    fXAFinalCapZ = p.get<float>("XAFinalCapZ");
    fOpFlashAlgoMinTime = p.get<double>("OpFlashAlgoMinTime", 0.010); // 10 ns [0.6 tick]
    fOpFlashAlgoMaxTime = p.get<double>("OpFlashAlgoMaxTime", 0.016); // 16 ns [1 tick]
    fOpFlashAlgoRad = p.get<double>("OpFlashAlgoRad");
    fOpFlashAlgoPE = p.get<float>("OpFlashAlgoPE");
    fOpFlashAlgoTriggerPE = p.get<float>("OpFlashAlgoTriggerPE");
    fOpFlashAlgoHotVertexThld = p.get<float>("OpFlashAlgoHotVertexThld");
    fAdjOpFlashMembraneProjection = p.get<bool>("AdjOpFlashMembraneProjection");
    fAdjOpFlashEndCapProjection = p.get<bool>("AdjOpFlashEndCapProjection");
    fAdjOpFlashX = p.get<float>("AdjOpFlashX", 100.0);
    fAdjOpFlashY = p.get<float>("AdjOpFlashY", 100.0);
    fAdjOpFlashZ = p.get<float>("AdjOpFlashZ", 100.0);
    fAdjOpFlashMaxPERatioCut = p.get<float>("AdjOpFlashMaxPERatioCut");
    fAdjOpFlashMinPECut = p.get<float>("AdjOpFlashMinPECut");
    fAdjOpFlashMinNHitCut = p.get<int>("AdjOpFlashMinNHitCut");
    fFlashMatchByResidual = p.get<bool>("FlashMatchByResidual");
    fSaveSignalDaughters = p.get<bool>("SaveSignalDaughters");
    fSaveSignalEDep = p.get<bool>("SaveSignalEDep");
    fSaveSignalOpHits = p.get<bool>("SaveSignalOpHits");
    fSaveOpFlashInfo = p.get<bool>("SaveOpFlashInfo");
    fSaveTrackInfo = p.get<bool>("SaveTrackInfo");
    // Generate the list of labels to be used in the analysis
    fLabels.push_back(fSignalLabel);
    for (auto const &label : fBackgroundLabels)
    {
      fLabels.push_back(label);
    }
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
    fConfigTree->Branch("SignalLabel", &fSignalLabel);
    fConfigTree->Branch("GEANT4Label", &fGEANTLabel);
    fConfigTree->Branch("HitLabel", &fHitLabel);
    fConfigTree->Branch("ClusterLabel", &fClusterLabel);
    fConfigTree->Branch("SolarClusterLabel", &fSolarClusterLabel);
    fConfigTree->Branch("TrackLabel", &fTrackLabel);
    fConfigTree->Branch("OpHitLabel", &fOpHitLabel);
    fConfigTree->Branch("OpFlashLabel", &fOpFlashLabel);
    fConfigTree->Branch("OpFlashTime2us", &fOpFlashTime2us);
    fConfigTree->Branch("OpFlashTimeOffset", &fOpFlashTimeOffset);
    fConfigTree->Branch("ClusterAlgoTime", &fClusterAlgoTime);
    fConfigTree->Branch("ClusterAlgoAdjChannel", &fClusterAlgoAdjChannel);
    fConfigTree->Branch("GenerateSolarCluster", &fGenerateSolarCluster);
    fConfigTree->Branch("ClusterMatchNHit", &fClusterMatchNHit);
    fConfigTree->Branch("ClusterMatchCharge", &fClusterMatchCharge);
    fConfigTree->Branch("ClusterMatchTime", &fClusterMatchTime);
    fConfigTree->Branch("ClusterInd0MatchTime", &fClusterInd0MatchTime);
    fConfigTree->Branch("ClusterInd1MatchTime", &fClusterInd1MatchTime);
    fConfigTree->Branch("ClusterPreselectionSignal", &fClusterPreselectionSignal);
    fConfigTree->Branch("ClusterPreselectionPrimary", &fClusterPreselectionPrimary);
    fConfigTree->Branch("ClusterPreselectionNHits", &fClusterPreselectionNHits);
    fConfigTree->Branch("ClusterPreselectionTrack", &fClusterPreselectionTrack);
    fConfigTree->Branch("ClusterPreselectionFlashMatch", &fClusterPreselectionFlashMatch);
    fConfigTree->Branch("AdjClusterRad", &fAdjClusterRad);
    fConfigTree->Branch("MinClusterCharge", &fMinClusterCharge);
    fConfigTree->Branch("GenerateAdjOpFlash", &fGenerateAdjOpFlash);
    fConfigTree->Branch("XACathodeX", &fXACathodeX);
    fConfigTree->Branch("XAMembraneY", &fXAMembraneY);
    fConfigTree->Branch("XAStartCapZ", &fXAStartCapZ);
    fConfigTree->Branch("XAFinalCapZ", &fXAFinalCapZ);
    fConfigTree->Branch("OpFlashAlgoMinTime", &fOpFlashAlgoMinTime);
    fConfigTree->Branch("OpFlashAlgoMaxTime", &fOpFlashAlgoMaxTime);
    fConfigTree->Branch("OpFlashAlgoRad", &fOpFlashAlgoRad);
    fConfigTree->Branch("OpFlashAlgoPE", &fOpFlashAlgoPE);
    fConfigTree->Branch("OpFlashAlgoTriggerPE", &fOpFlashAlgoTriggerPE);
    fConfigTree->Branch("OpFlashAlgoHotVertexThld", &fOpFlashAlgoHotVertexThld);
    fConfigTree->Branch("AdjOpFlashMembraneProjection", &fAdjOpFlashMembraneProjection);
    fConfigTree->Branch("AdjOpFlashEndCapProjection", &fAdjOpFlashEndCapProjection);
    fConfigTree->Branch("AdjOpFlashX", &fAdjOpFlashX);
    fConfigTree->Branch("AdjOpFlashY", &fAdjOpFlashY);
    fConfigTree->Branch("AdjOpFlashZ", &fAdjOpFlashZ);
    fConfigTree->Branch("AdjOpFlashMaxPERatioCut", &fAdjOpFlashMaxPERatioCut);
    fConfigTree->Branch("AdjOpFlashMinPECut", &fAdjOpFlashMinPECut);
    fConfigTree->Branch("AdjOpFlashMinNHitCut", &fAdjOpFlashMinNHitCut);
    fConfigTree->Branch("FlashMatchByResidual", &fFlashMatchByResidual);
    fConfigTree->Branch("SaveSignalDaughters", &fSaveSignalDaughters);
    fConfigTree->Branch("SaveSignalEDep", &fSaveSignalEDep);
    fConfigTree->Branch("SaveSignalOpHits", &fSaveSignalOpHits);
    fConfigTree->Branch("SaveOpFlashInfo", &fSaveOpFlashInfo);
    fConfigTree->Branch("SaveTrackInfo", &fSaveTrackInfo);

    // MC Truth info.
    fMCTruthTree->Branch("Event", &Event, "Event/I");                                        // Event number
    fMCTruthTree->Branch("Flag", &Flag, "Flag/I");                                           // Flag used to match truth with reco tree entries
    fMCTruthTree->Branch("TruthPart", &TPart);                                               // Number particles per generator
    fMCTruthTree->Branch("Interaction", &sInteraction);                                      // True signal interaction process
    fMCTruthTree->Branch("SignalParticleE", &SignalParticleE, "SignalParticleE/F");          // True signal energy [MeV]
    fMCTruthTree->Branch("SignalParticleP", &SignalParticleP, "SignalParticleP/F");          // True signal momentum [MeV]
    fMCTruthTree->Branch("SignalParticleK", &SignalParticleK, "SignalParticleK/F");          // True signal K.E. [MeV]
    fMCTruthTree->Branch("SignalParticleX", &SignalParticleX, "SignalParticleX/F");          // True signal X [cm]
    fMCTruthTree->Branch("SignalParticleY", &SignalParticleY, "SignalParticleY/F");          // True signal Y [cm]
    fMCTruthTree->Branch("SignalParticleZ", &SignalParticleZ, "SignalParticleZ/F");          // True signal Z [cm]
    fMCTruthTree->Branch("SignalParticlePDG", &SignalParticlePDG, "SignalParticlePDG/I");    // True signal PDG
    fMCTruthTree->Branch("SignalParticleTime", &SignalParticleTime, "SignalParticleTime/F"); // True signal time [tick]
    fMCTruthTree->Branch("OpHitNum", &OpHitNum, "OpHitNum/I");                               // Number of OpHits
    fMCTruthTree->Branch("OpFlashNum", &OpFlashNum, "OpFlashNum/I");                         // Number of OpFlashes
    fMCTruthTree->Branch("HitNum", &HitNum, "HitNum/I");                                     // Number of hits in each TPC plane
    fMCTruthTree->Branch("ClusterNum", &ClusterNum, "ClusterNum/I");                         // Number of clusters in each TPC plane
    fMCTruthTree->Branch("TrackNum", &TrackNum, "TrackNum/I");                               // Number of PMTracks
    
    if (fSaveSignalDaughters)
    { // Save Signal Daughters.
      fMCTruthTree->Branch("TSignalPDG", &SignalPDGList);         // PDG of Signal marticles
      fMCTruthTree->Branch("TSignalE", &SignalEList);             // Energy of Signal particles [MeV]
      fMCTruthTree->Branch("TSignalP", &SignalPList);             // Energy of Signal momentum [MeV]
      fMCTruthTree->Branch("TSignalK", &SignalKList);             // Kinetik Energy of Signal particles [MeV]
      fMCTruthTree->Branch("TSignalT", &SignalTimeList);          // Time of Signal particles [ticks]
      fMCTruthTree->Branch("TSignalEndX", &SignalEndXList);       // X of Signal particles [cm]
      fMCTruthTree->Branch("TSignalEndY", &SignalEndYList);       // Y of Signal particles [cm]
      fMCTruthTree->Branch("TSignalEndZ", &SignalEndZList);       // Z of Signal particles [cm]
      fMCTruthTree->Branch("TSignalMaxEDep", &SignalMaxEDepList); // Energy of Signal particles [MeV]
      fMCTruthTree->Branch("TSignalX", &SignalMaxEDepXList);      // X of Signal particles [cm]
      fMCTruthTree->Branch("TSignalY", &SignalMaxEDepYList);      // Y of Signal particles [cm]
      fMCTruthTree->Branch("TSignalZ", &SignalMaxEDepZList);      // Z of Signal particles [cm]
      fMCTruthTree->Branch("TSignalID", &SignalIDList);           // TrackID of Signal particles
      fMCTruthTree->Branch("TSignalMother", &SignalMotherList);   // TrackID of Signal mother
    }
    if (fSaveSignalEDep)
    { // Save Energy deposited by Signal particles. (Can be very heavy for any production)
      fMCTruthTree->Branch("TSignalPDGDepList", &SignalPDGDepList);           // PDG for Energy deposited of Signal particles
      fMCTruthTree->Branch("TSignalEDepList", &SignalEDepList);               // Energy deposited of Signal particles [MeV]
      fMCTruthTree->Branch("TSignalXDepList", &SignalXDepList);               // X deposited of Signal particles [cm]
      fMCTruthTree->Branch("TSignalYDepList", &SignalYDepList);               // Y deposited of Signal particles [cm]
      fMCTruthTree->Branch("TSignalZDepList", &SignalZDepList);               // Z deposited of Signal particles [cm]
      fMCTruthTree->Branch("TSignalIDDepList", &SignalIDDepList);             // ParentID of Signal particles
      fMCTruthTree->Branch("TSignalElectronDepList", &SignalElectronDepList); // Number of electrons in the Signal particles
    }
    if (fSaveSignalOpHits)
    { // Save OpHits. (Can be very heavy for background productions)
      fMCTruthTree->Branch("OpHitPur", &SOpHitPur);         // OpHit Purity
      fMCTruthTree->Branch("OpHitPlane", &SOpHitPlane);     // OpHit Plane
      fMCTruthTree->Branch("OpHitPE", &SOpHitPE);           // OpHit PE
      fMCTruthTree->Branch("OpHitX", &SOpHitX);             // OpHit X
      fMCTruthTree->Branch("OpHitY", &SOpHitY);             // OpHit Y
      fMCTruthTree->Branch("OpHitZ", &SOpHitZ);             // OpHit Z
      fMCTruthTree->Branch("OpHitTime", &SOpHitTime);       // OpHit Time
      fMCTruthTree->Branch("OpHitChannel", &SOpHitChannel); // OpHit Channel
      fMCTruthTree->Branch("OpHitFlashID", &SOpHitFlashID); // OpHit Area
    }
    if (fSaveOpFlashInfo)
    {
      fMCTruthTree->Branch("OpFlashPur", &OpFlashPur);     // OpFlash Purity
      fMCTruthTree->Branch("OpFlashID", &OpFlashID);       // OpFlash ID
      fMCTruthTree->Branch("OpFlashPE", &OpFlashPE);       // OpFlash PE
      fMCTruthTree->Branch("OpFlashX", &OpFlashX);         // OpFlash X
      fMCTruthTree->Branch("OpFlashY", &OpFlashY);         // OpFlash Y
      fMCTruthTree->Branch("OpFlashZ", &OpFlashZ);         // OpFlash Z
      fMCTruthTree->Branch("OpFlashTime", &OpFlashTime);   // OpFlash Time
      fMCTruthTree->Branch("OpFlashSTD", &OpFlashSTD);     // OpFlash STD
      fMCTruthTree->Branch("OpFlashNHits", &OpFlashNHits); // OpFlash NHit
      fMCTruthTree->Branch("OpFlashPlane", &OpFlashPlane); // OpFlash Plane
      fMCTruthTree->Branch("OpFlashMaxPE", &OpFlashMaxPE); // OpFlash Max PE
    }

    // Repeated Truth info.
    fSolarNuAnaTree->Branch("Event", &Event, "Event/I");                                        // Event number
    fSolarNuAnaTree->Branch("Flag", &Flag, "Flag/I");                                           // Flag used to match truth with reco tree entries
    fSolarNuAnaTree->Branch("TruthPart", &TPart);                                               // Number particles per generator
    fSolarNuAnaTree->Branch("Interaction", &sInteraction);                                      // True signal interaction process
    fSolarNuAnaTree->Branch("SignalParticleE", &SignalParticleE, "SignalParticleE/F");          // True signal energy [MeV]
    fSolarNuAnaTree->Branch("SignalParticleP", &SignalParticleP, "SignalParticleP/F");          // True signal momentum [MeV]
    fSolarNuAnaTree->Branch("SignalParticleK", &SignalParticleK, "SignalParticleK/F");          // True signal K.E. [MeV]
    fSolarNuAnaTree->Branch("SignalParticleX", &SignalParticleX, "SignalParticleX/F");          // True signal X [cm]
    fSolarNuAnaTree->Branch("SignalParticleY", &SignalParticleY, "SignalParticleY/F");          // True signal Y [cm]
    fSolarNuAnaTree->Branch("SignalParticleZ", &SignalParticleZ, "SignalParticleZ/F");          // True signal Z [cm]
    fSolarNuAnaTree->Branch("SignalParticlePDG", &SignalParticlePDG, "SignalParticlePDG/I");    // True signal PDG
    fSolarNuAnaTree->Branch("SignalParticleTime", &SignalParticleTime, "SignalParticleTime/F"); // True signal Time [tick]
    
    if (fSaveSignalDaughters)
    { // Save Signal Daughters. (Only makes sense for marley)
      fSolarNuAnaTree->Branch("TSignalPDG", &SignalPDGList);         // PDG of Signal particles
      fSolarNuAnaTree->Branch("TSignalE", &SignalEList);             // Energy of Signal particles
      fSolarNuAnaTree->Branch("TSignalP", &SignalPList);             // Momentum of Signal particles
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
    }

    // Main Cluster info.
    fSolarNuAnaTree->Branch("Primary", &MPrimary);                                   // Cluster hasn't any adjcl with AdjClCharge > MCharge (bool)
    fSolarNuAnaTree->Branch("Ind0Purity", &MInd0Pur, "Ind0Purity/F");                // Main cluster ind0 reco signal purity
    fSolarNuAnaTree->Branch("Ind1Purity", &MInd1Pur, "Ind1Purity/F");                // Main cluster ind1 reco signal purity
    fSolarNuAnaTree->Branch("Purity", &MPur, "Purity/F");                            // Main cluster reco signal purity
    fSolarNuAnaTree->Branch("Generator", &MGen, "Generator/I");                      // Main cluster generator idx
    fSolarNuAnaTree->Branch("GenPurity", &MGenPur, "GenPurity/F");                   // Main cluster reco generator purity
    fSolarNuAnaTree->Branch("TPC", &MTPC, "ColTPC/I");                               // Main cluster TPC
    fSolarNuAnaTree->Branch("Time", &MTime, "ColTime/F");                            // Main cluster time [ticks]
    fSolarNuAnaTree->Branch("NHits", &MNHit, "ColNHits/I");                          // Main cluster #hits
    fSolarNuAnaTree->Branch("Charge", &MCharge, "ColCharge/F");                      // Main cluster charge [ADC*ticks]
    fSolarNuAnaTree->Branch("MaxCharge", &MMaxCharge, "ColCharge/F");                // Main cluster's max TPCHit-charge [ADC*ticks]
    fSolarNuAnaTree->Branch("Ind0TPC", &MInd0TPC, "Ind0TPC/I");                      // Main cluster ind0 TPC
    fSolarNuAnaTree->Branch("Ind1TPC", &MInd1TPC, "Ind1TPC/I");                      // Main cluster ind1 TPC
    fSolarNuAnaTree->Branch("Ind0dTime", &MInd0dTime, "Ind0dTime/F");                // Main cluster ind0 dT [Ticks]
    fSolarNuAnaTree->Branch("Ind1dTime", &MInd1dTime, "Ind1dTime/F");                // Main cluster ind1 dT [Ticks]
    fSolarNuAnaTree->Branch("Ind0NHits", &MInd0NHits, "Ind0NHits/I");                // Main cluster ind0 Hits
    fSolarNuAnaTree->Branch("Ind1NHits", &MInd1NHits, "Ind1NHits/I");                // Main cluster ind1 Hits
    fSolarNuAnaTree->Branch("Ind0Charge", &MInd0Charge, "Ind0Charge/F");             // Main cluster ind0 MaxHit
    fSolarNuAnaTree->Branch("Ind1Charge", &MInd1Charge, "Ind1Charge/F");             // Main cluster ind1 MaxHit
    fSolarNuAnaTree->Branch("Ind0MaxCharge", &MInd0MaxCharge, "Ind0MaxCharge/F");    // Main cluster ind0 MaxHit
    fSolarNuAnaTree->Branch("Ind1MaxCharge", &MInd1MaxCharge, "Ind1MaxCharge/F");    // Main cluster ind1 MaxHit
    fSolarNuAnaTree->Branch("Ind0RecoY", &MInd0RecoY, "Ind0RecoY/F");                // Main cluster ind0 reco Y [cm]
    fSolarNuAnaTree->Branch("Ind1RecoY", &MInd1RecoY, "Ind1RecoY/F");                // Main cluster ind1 reco Y [cm]
    fSolarNuAnaTree->Branch("RecoX", &MRecX, "RecoX/F");                             // Main cluster reco X [cm] (from matched Flash)
    fSolarNuAnaTree->Branch("RecoY", &MRecY, "RecoY/F");                             // Main cluster reco Y [cm]
    fSolarNuAnaTree->Branch("RecoZ", &MRecZ, "RecoZ/F");                             // Main cluster reco Z [cm]
    fSolarNuAnaTree->Branch("MainID", &MMainID, "MainID/I");                         // Main cluster main track ID
    fSolarNuAnaTree->Branch("MainE", &MMainE, "MainE/F");                            // Main cluster main energy [MeV]
    fSolarNuAnaTree->Branch("MainP", &MMainP, "MainP/F");                            // Main cluster main momentum [MeV]
    fSolarNuAnaTree->Branch("MainK", &MMainK, "MainK/F");                            // Main cluster main kinetic energy [MeV]
    fSolarNuAnaTree->Branch("MainTime", &MMainTime, "MainTime/F");                   // Main cluster main Time [ticks]
    fSolarNuAnaTree->Branch("MainPDG", &MMainPDG, "MainPDG/I");                      // Main cluster main pdg
    fSolarNuAnaTree->Branch("MainParentPDG", &MMainParentPDG, "MainParentPDG/I");    // Main cluster main pdg
    fSolarNuAnaTree->Branch("MainParentE", &MMainParentE, "MainParentE/F");          // Main cluster main parent energy [MeV]
    fSolarNuAnaTree->Branch("MainParentP", &MMainParentP, "MainParentP/F");          // Main cluster main parent momentum [MeV]
    fSolarNuAnaTree->Branch("MainParentK", &MMainParentK, "MainParentK/F");          // Main cluster main parent kinetic energy [MeV]
    fSolarNuAnaTree->Branch("MainParentTime", &MMainParentTime, "MainParentTime/F"); // Main cluster main parent Time [ticks]
    fSolarNuAnaTree->Branch("MainVertex", &MMainVertex);                             // Main cluster main particle vertex [cm]
    fSolarNuAnaTree->Branch("EndVertex", &MEndVertex);                               // Main cluster end particle vertex [cm]
    fSolarNuAnaTree->Branch("MainParentVertex", &MMainParentVertex);                 // Main cluster parent particle vertex [cm]
    fSolarNuAnaTree->Branch("GenFrac", &MGenFrac);                                   // Main cluster reco purity complete
    fSolarNuAnaTree->Branch("SignalFrac", &MSignalFrac);                             // Main cluster particle contribution (electron, gamma, neutron)

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
    fSolarNuAnaTree->Branch("AdjClNHits", &MAdjClNHits);                 // Adj. clusters' #hits
    fSolarNuAnaTree->Branch("AdjClInd0NHits", &MAdjClInd0NHits);         // Adj. clusters' #hits
    fSolarNuAnaTree->Branch("AdjClInd1NHits", &MAdjClInd1NHits);         // Adj. clusters' #hits
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
    fSolarNuAnaTree->Branch("AdjClMainE", &MAdjClMainE);                 // Adj. clusters' main energy [MeV]
    fSolarNuAnaTree->Branch("AdjClMainP", &MAdjClMainP);                 // Adj. clusters' main momentum [MeV]
    fSolarNuAnaTree->Branch("AdjClMainK", &MAdjClMainK);                 // Adj. clusters' main K.E. [MeV]
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
      fSolarNuAnaTree->Branch("AdjOpFlashNHits", &MAdjFlashNHits);       // Adj. flash' #hits
      fSolarNuAnaTree->Branch("AdjOpFlashPlane", &MAdjFlashPlane);       // Adj. flash' Plane
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
    fSolarNuAnaTree->Branch("MatchedOpFlashTime", &MFlashTime, "MatchedOpFlashTime/F");             // Matched flash' time [ticks]
    fSolarNuAnaTree->Branch("MatchedOpFlashNHits", &MFlashNHits, "MatchedOpFlashNHits/I");          // Matched flash' #hits
    fSolarNuAnaTree->Branch("MatchedOpFlashPlane", &MFlashPlane, "MatchedOpFlashPlane/I");          // Matched flash' Plane
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
  } // BeginJob

  //......................................................
  void SolarNuAna::analyze(art::Event const &evt)
  {
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------- Prepare everything for new event ----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    std::string fGeometry = "";
    std::string geoName = geom->DetectorName();
    std::vector<std::set<int>> trackids = {};
    std::map<int, simb::MCParticle> ThisGeneratorParts;
    std::vector<recob::Hit> Ind0Hits, Ind1Hits, ColHits, GhostHits;
    std::vector<std::vector<recob::Hit>> Clusters0, Clusters1, Clusters2, Clusters3;
    // --- We want to reset all of our previous run and TTree variables ---
    ResetVariables();
    ThisGeneratorParts.clear();
    Event = evt.event();
    Flag = rand() % 10000000000;
    // Initialize the services we need
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    // Get the drift properties of the detector
    geo::CryostatID c(0);
    const geo::CryostatGeo& cryostat = geom->Cryostat(c);
    const geo::TPCGeo &tpcg = cryostat.TPC(0);
    float driftLength = 0.;
    float driftTime = 0.;
    // Get size of the fiducial LAr volume
    const double fidVolX = cryostat.BoundingBox().MaxX() - cryostat.BoundingBox().MinX();
    const double fidVolY = cryostat.BoundingBox().MaxY() - cryostat.BoundingBox().MinY();
    const double fidVolZ = cryostat.BoundingBox().MaxZ() - cryostat.BoundingBox().MinZ();
    // Determine geometry type and calculate drift length and time
    if (geoName.find("dune10kt") != std::string::npos) {
        driftLength = tpcg.DriftDistance();
        driftTime = driftLength / detProp.DriftVelocity(); // in us
        fGeometry = "HD";
    }
    else if(geoName.find("dunevd10kt") != std::string::npos) {
        driftLength = tpcg.DriftDistance();
        driftTime = driftLength / detProp.DriftVelocity();
        fGeometry = "VD";
    }
    else {
        throw cet::exception("SolarNuAna") << "Geometry " << geoName << " not supported. Only dune10kt and dunevd10kt are supported.\n";
    }
    // Loop over all TPCs in the cryostat and fill the map
    std::string sTPCMap = "";
    unsigned int maxTPC = 0;
    for (auto const& tpcid : geom->Iterate<geo::TPCID>()) {
      if (tpcid.isValid) { // Fill the TPC map with the TPC ID
        TPCIDMap[tpcid.TPC] = tpcid;
        TPCIDdriftLength[tpcid.TPC] = cryostat.TPC(tpcid).DriftDistance();
        TPCIDdriftTime[tpcid.TPC] = driftLength / detProp.DriftVelocity();
        if (tpcid.TPC > maxTPC) { // Keep track of the maximum TPC ID
          maxTPC = tpcid.TPC; 
        }
        sTPCMap += "Found TPC ID: " + std::to_string(tpcid.TPC) + " in Cryostat: " + std::to_string(c.Cryostat) + 
                   " with Drift Length: " + ProducerUtils::str(TPCIDdriftLength[tpcid.TPC]) + 
                   " cm and Drift Time: " + ProducerUtils::str(TPCIDdriftTime[tpcid.TPC]) + " us\n";
      }
    }
    // Add extra TPCID entry -1 for all clusters that are not associated with a TPC
    TPCIDMap[-1] = geo::TPCID(); // Invalid TPC ID
    TPCIDdriftLength[-1] = TPCIDdriftLength.begin()->second; // Use the first valid TPC's drift length
    TPCIDdriftTime[-1] = TPCIDdriftTime.begin()->second;     // Use the first valid TPC's drift time
    producer->PrintInColor(sTPCMap, ProducerUtils::GetColor("yellow"), "Debug");

    std::string sHead = "";
    sHead = sHead + "\n#########################################";
    sHead = sHead + "\nEvent: " + ProducerUtils::str(Event) + " Flag: " + ProducerUtils::str(Flag);
    sHead = sHead + "\nGeometry: " + geoName + " (" + fGeometry + ")";
    sHead = sHead + "\nSignal Label: " + fLabels[0];
    sHead = sHead + "\nPDS Frequency in [MHz]: " + ProducerUtils::str(clockData.OpticalClock().Frequency());
    sHead = sHead + "\nPDS Tick in [us]: " + ProducerUtils::str(clockData.OpticalClock().TickPeriod(), 3);
    sHead = sHead + "\nTPC Frequency in [MHz]: " + ProducerUtils::str(clockData.TPCClock().Frequency());
    sHead = sHead + "\nTPC Tick in [us]: " + ProducerUtils::str(clockData.TPCClock().TickPeriod());
    sHead = sHead + "\nTPC Map: " + ProducerUtils::str(TPCIDMap.size()) + " TPCs found";
    sHead = sHead + "\nTPC DriftLength in [cm]: " + ProducerUtils::str(driftLength);
    sHead = sHead + "\nTPC DriftTime in [us]: " + ProducerUtils::str(driftTime);
    sHead = sHead + "\nFiducial Volume in [cm]: " + ProducerUtils::str(fidVolX) + ", " + ProducerUtils::str(fidVolY) + ", " + ProducerUtils::str(fidVolZ);
    sHead = sHead + "\n#########################################";
    producer->PrintInColor(sHead, ProducerUtils::GetColor("magenta"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Create maps for ID tracking -----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // --- Fill MC Truth IDs to tracking vectors. Get a list of all of my particles in one chunk. ---
    std::string sMcTruth = "";
    const sim::ParticleList &PartList = pi_serv->ParticleList();
    sMcTruth = sMcTruth + "\nThere are a total of " + ProducerUtils::str(int(PartList.size())) + " Particles in the event\n";
    // Loop over all signal+bkg handles and collect track IDs
    for (size_t i = 0; i < fLabels.size(); i++)
    {
      GeneratorParticles.push_back(ThisGeneratorParts); // For each label insert empty list

      art::Handle<std::vector<simb::MCTruth>> ThisHandle;
      evt.getByLabel(fLabels[i], ThisHandle);

      if (ThisHandle) {
        auto ThisValidHanlde = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[i]); // Get generator handles
        art::FindManyP<simb::MCParticle> Assn(ThisValidHanlde, evt, fGEANTLabel);          // Assign labels to MCPArticles
        producer->FillMyMaps(GeneratorParticles[i], Assn, ThisValidHanlde);                // Fill empty list with previously assigned particles

        sMcTruth = sMcTruth + "\n# of particles " + ProducerUtils::str(int(GeneratorParticles[i].size())) + "\tfrom gen " + ProducerUtils::str(int(i) + 1) + " " + fLabels[i];
        TPart.push_back(GeneratorParticles[i].size());
        if (GeneratorParticles[i].size() > 0) {
          for (std::map<int, simb::MCParticle>::iterator iter = GeneratorParticles[i].begin(); iter != GeneratorParticles[i].end(); iter++)
          {
            std::set<int> ThisGeneratorIDs = {};
            trackids.push_back(ThisGeneratorIDs);
            trackids[i].insert(iter->first);
          }
        }
        else {
          std::set<int> ThisGeneratorIDs = {};
          trackids.push_back(ThisGeneratorIDs);
        }
      }
      else {
        sMcTruth = sMcTruth + "\n# of particles " + ProducerUtils::str(int(GeneratorParticles[i].size())) + "\tfrom gen " + ProducerUtils::str(int(i) + 1) + " " + fLabels[i] + " *not generated!";
        TPart.push_back(0);
        std::set<int> ThisGeneratorIDs = {};
        trackids.push_back(ThisGeneratorIDs);
      }
    }

    producer->PrintInColor(sMcTruth, ProducerUtils::GetColor("bright_red"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Some MC Truth information -------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    std::set<int> SignalPrimaryTrackIDs, SignalTrackIDs;                // Signal TrackIDs to be used in OpFlash matching 
    std::map<int, int> SignalPDGMap, SignalMotherMap, SignalPrimaryMap; // Count of different PDGs in the event
    std::map<int, float> SignalEMap, SignalPMap, SignalKMap, SignalStartXMap, SignalStartYMap, SignalStartZMap, SignalFinalXMap, SignalFinalYMap, SignalFinalZMap, SignalMaxEDepMap, SignalMaxEDepXMap, SignalMaxEDepYMap, SignalMaxEDepZMap, SignalTimeMap;
    std::vector<std::vector<int>> ClPartTrackIDs = {{}, {}, {}, {}};    // Track IDs corresponding to each kind of MCTruth particle  {11,22,2112,else}
    art::Handle<std::vector<simb::MCTruth>> ThisHandle;
    std::string sSignalParticle = "";
    std::string sSignalTruth = "";

    evt.getByLabel(fLabels[0], ThisHandle);
    if (ThisHandle) {
      auto Signal = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[0]); // Get handle for SIGNAL MCTruths
      // --- Loop over all neutrinos in the event ---
      for (auto const &SignalTruth : *Signal)
      {
        int NSignalParticles = SignalTruth.NParticles();
        sSignalTruth = sSignalTruth + "\nNumber of Signal Particles: " + ProducerUtils::str(NSignalParticles);
        
        if (fLabels[0] == "marley") {
          sSignalParticle = "Neutrino";
          const simb::MCNeutrino &nue = SignalTruth.GetNeutrino();
          sInteraction = ProducerUtils::str(nue.InteractionType());
          SignalParticleE = 1e3 * nue.Nu().E();
          SignalParticleP = 1e3 * nue.Nu().P();
          SignalParticleK = 1e3 * nue.Nu().E() - 1e3 * nue.Nu().Mass();
          SignalParticleX = nue.Nu().Vx();
          SignalParticleY = nue.Nu().Vy();
          SignalParticleZ = nue.Nu().Vz();
          SignalParticlePDG = nue.Nu().PdgCode();
          SignalParticleTime = nue.Nu().T();
          sSignalTruth = sSignalTruth + "\nNeutrino Interaction: " + sInteraction;
          sSignalTruth = sSignalTruth + "\nNeutrino Energy: " + ProducerUtils::str(SignalParticleE) + " MeV";
          sSignalTruth = sSignalTruth + "\nPosition (" + ProducerUtils::str(SignalParticleX) + ", " + ProducerUtils::str(SignalParticleY) + ", " + ProducerUtils::str(SignalParticleZ) + ") cm";
        }

        else if (fLabels[0] == "generator") {
          sInteraction = "Background";
          sSignalTruth = sSignalTruth + "\nFound generator label: " + fLabels[0] + ". Using single generator config.\n";
          
          if (NSignalParticles > 1) {
            sSignalTruth = sSignalTruth + "\n[WARNING] Multiple particles found in the Signal MCTruth.\n";
          }

          for (int i = 0; i < NSignalParticles; i++)
          {
            const simb::MCParticle &SignalParticle = SignalTruth.GetParticle(i);
            SignalParticleE = 1e3 * SignalParticle.E();
            SignalParticleP = 1e3 * SignalParticle.P();
            SignalParticleK = 1e3 * SignalParticle.E() - 1e3 * SignalParticle.Mass();
            SignalParticleX = SignalParticle.Vx();
            SignalParticleY = SignalParticle.Vy();
            SignalParticleZ = SignalParticle.Vz();
            SignalParticlePDG = SignalParticle.PdgCode();
            SignalParticleTime = SignalParticle.T();
            
            if (abs(SignalParticle.PdgCode()) == 12) {
              sSignalParticle = "Neutrino";
            }
            else if (abs(SignalParticle.PdgCode()) == 11) {
              sSignalParticle = "Electron";
            }
            else if (abs(SignalParticle.PdgCode()) == 22) {
              sSignalParticle = "Photon";
            }
            else if (abs(SignalParticle.PdgCode()) == 2112) {
              sSignalParticle = "Neutron";
            }
            else {
              sSignalParticle = "Other";
            }

            sSignalTruth = sSignalTruth + sSignalParticle + " Energy: " + ProducerUtils::str(SignalParticleK) + " MeV\n"; 
            sSignalTruth += "\t- Strat Position: (" + ProducerUtils::str(SignalParticleX) + ", " + ProducerUtils::str(SignalParticleY) + ", " + ProducerUtils::str(SignalParticleZ) + ") cm\n";
            sSignalTruth += "\t- Final Position: (" + ProducerUtils::str(SignalParticle.EndX()) + ", " + ProducerUtils::str(SignalParticle.EndY()) + ", " + ProducerUtils::str(SignalParticle.EndZ()) + ") cm\n";
          }
        }
        else {
          sSignalTruth = sSignalTruth + "\n[WARNING] Unknown generator label: " + fLabels[0] + ". No truth information saved.\n";
          sInteraction = "Unknown";
        }
      }
      art::FindManyP<simb::MCParticle> SignalAssn(Signal, evt, fGEANTLabel);
      sSignalTruth = sSignalTruth + "\nGen.\tPdgCode\t\tEnergy\t\tEndPosition\t\tTrackID\tMotherID";
      sSignalTruth = sSignalTruth + "\n--------------------------------------------------------------------------------";

      for (size_t i = 0; i < SignalAssn.size(); i++)
      {
        auto SignalParticles = SignalAssn.at(i);
        for (auto SignalParticle = SignalParticles.begin(); SignalParticle != SignalParticles.end(); SignalParticle++)
        {
          SignalPDGList.push_back((*SignalParticle)->PdgCode());
          SignalEList.push_back(1e3 * (*SignalParticle)->E());
          SignalPList.push_back(1e3 * (*SignalParticle)->P());
          SignalKList.push_back(1e3 * (*SignalParticle)->E() - 1e3 * (*SignalParticle)->Mass());
          SignalTimeList.push_back((*SignalParticle)->T());
          SignalEndXList.push_back((*SignalParticle)->EndX());
          SignalEndYList.push_back((*SignalParticle)->EndY());
          SignalEndZList.push_back((*SignalParticle)->EndZ());
          SignalIDList.push_back((*SignalParticle)->TrackId());
          SignalMotherList.push_back((*SignalParticle)->Mother());
          SignalMotherMap[(*SignalParticle)->TrackId()] = (*SignalParticle)->Mother();
          
          if ((*SignalParticle)->Mother() == 0) {
            SignalPrimaryTrackIDs.insert((*SignalParticle)->TrackId());
            SignalPDGMap[(*SignalParticle)->TrackId()] = (*SignalParticle)->PdgCode();
            SignalStartXMap[(*SignalParticle)->TrackId()] = (*SignalParticle)->Vx();
            SignalStartYMap[(*SignalParticle)->TrackId()] = (*SignalParticle)->Vy();
            SignalStartZMap[(*SignalParticle)->TrackId()] = (*SignalParticle)->Vz();
            SignalFinalXMap[(*SignalParticle)->TrackId()] = (*SignalParticle)->EndX();
            SignalFinalYMap[(*SignalParticle)->TrackId()] = (*SignalParticle)->EndY();
            SignalFinalZMap[(*SignalParticle)->TrackId()] = (*SignalParticle)->EndZ();
            SignalEMap[(*SignalParticle)->TrackId()] = 1e3 * (*SignalParticle)->E();
            SignalPMap[(*SignalParticle)->TrackId()] = 1e3 * (*SignalParticle)->P();
            SignalKMap[(*SignalParticle)->TrackId()] = 1e3 * (*SignalParticle)->E() - 1e3 * (*SignalParticle)->Mass();
            SignalTimeMap[(*SignalParticle)->TrackId()] = (*SignalParticle)->T();
          }

          std::vector<const sim::IDE *> ides = bt_serv->TrackIdToSimIDEs_Ps((*SignalParticle)->TrackId());
          for (auto const &ide : ides)
          {
            if (ide->numElectrons < 1 || ide->energy < 1e-6 || abs(ide->x) > TPCIDdriftLength[0] || abs(ide->y) > fidVolY/2 || abs(ide->z) > fidVolZ) {
              continue;
            } 

            if (ProducerUtils::InMyMap((*SignalParticle)->TrackId(), SignalMaxEDepMap) == false) {
              SignalMaxEDepMap[(*SignalParticle)->TrackId()] = ide->energy;
              SignalMaxEDepXMap[(*SignalParticle)->TrackId()] = ide->x;
              SignalMaxEDepYMap[(*SignalParticle)->TrackId()] = ide->y;
              SignalMaxEDepZMap[(*SignalParticle)->TrackId()] = ide->z;
            }
            
            if (ide->energy > SignalMaxEDepMap[(*SignalParticle)->TrackId()]) {
              SignalMaxEDepMap[(*SignalParticle)->TrackId()] = ide->energy;
              SignalMaxEDepXMap[(*SignalParticle)->TrackId()] = ide->x;
              SignalMaxEDepYMap[(*SignalParticle)->TrackId()] = ide->y;
              SignalMaxEDepZMap[(*SignalParticle)->TrackId()] = ide->z;
            }

            if (abs((*SignalParticle)->PdgCode()) == 11 || abs((*SignalParticle)->PdgCode()) == 22 || abs((*SignalParticle)->PdgCode()) == 2112) {
              SignalIDDepList.push_back((*SignalParticle)->TrackId());
              SignalEDepList.push_back(ide->energy);
              SignalPDGDepList.push_back((*SignalParticle)->PdgCode());
              SignalXDepList.push_back(ide->x);
              SignalYDepList.push_back(ide->y);
              SignalZDepList.push_back(ide->z);
              SignalElectronDepList.push_back(ide->numElectrons);
            }
          }
          // Keep track of all signal track IDs for later
          SignalTrackIDs.emplace((*SignalParticle)->TrackId());
          SignalMaxEDepList.push_back(SignalMaxEDepMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepXList.push_back(SignalMaxEDepXMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepYList.push_back(SignalMaxEDepYMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepZList.push_back(SignalMaxEDepZMap[(*SignalParticle)->TrackId()]);
          std::string sLabel = "";
          // Special case for marley and generator labels
          if (fLabels[0] == "marley") {
            sLabel = "marley";
          }
          else if (fLabels[0] == "generator") {
            sLabel = "single";
          }
          else {
            sLabel = fLabels[0];
          }
          // Print differently some information for better alignment
          if ((*SignalParticle)->PdgCode() < 1000000) {
            sSignalTruth = sSignalTruth + "\n" + sLabel + "\t" + ProducerUtils::str((*SignalParticle)->PdgCode()) + "\t\t" + ProducerUtils::str(1e3 * (*SignalParticle)->E() - 1e3 * (*SignalParticle)->Mass()) + "\t(" + ProducerUtils::str((*SignalParticle)->Vx()) + ", " + ProducerUtils::str((*SignalParticle)->Vy()) + ", " + ProducerUtils::str((*SignalParticle)->Vz()) + ")\t" + ProducerUtils::str((*SignalParticle)->TrackId()) + "\t" + ProducerUtils::str((*SignalParticle)->Mother());
          }
          else {
            sSignalTruth = sSignalTruth + "\n" + sLabel + "\t" + ProducerUtils::str((*SignalParticle)->PdgCode()) + "\t" + ProducerUtils::str(1e3 * (*SignalParticle)->E() - 1e3 * (*SignalParticle)->Mass()) + "\t(" + ProducerUtils::str((*SignalParticle)->Vx()) + ", " + ProducerUtils::str((*SignalParticle)->Vy()) + ", " + ProducerUtils::str((*SignalParticle)->Vz()) + ")\t" + ProducerUtils::str((*SignalParticle)->TrackId()) + "\t" + ProducerUtils::str((*SignalParticle)->Mother());
          }
          // Save track IDs of certain particles for cluster association
          if ((*SignalParticle)->PdgCode() == 11) {
            const TLorentzVector &MainElectronEndPoint = (*SignalParticle)->EndPosition();
            MainElectronEndPointX = MainElectronEndPoint.X();
            ClPartTrackIDs[0].push_back((*SignalParticle)->TrackId());
          }
          else if ((*SignalParticle)->PdgCode() == 22) {
            ClPartTrackIDs[1].push_back((*SignalParticle)->TrackId());
          }
          else if ((*SignalParticle)->PdgCode() == 2112) {
            ClPartTrackIDs[2].push_back((*SignalParticle)->TrackId());
          }
          else {
            ClPartTrackIDs[3].push_back((*SignalParticle)->TrackId());
          }
        }
      }
    }
    
    else {
      mf::LogWarning("SolarNuAna") << "No SIGNAL MCTruths found.";
      sSignalTruth = sSignalTruth + "\nNo SIGNAL MCTruths found.\n";
    }

    sSignalTruth += "\nSignal Track IDs: ";
    // Backtrace each SignalParticle until find a trackID that is in the SignalPrimaryTrackIDs set (i.e. primary generated particles)
    for (auto const &SignalTrackID : SignalTrackIDs)
    {
      int MotherID = SignalMotherMap[SignalTrackID];
      sSignalTruth += ProducerUtils::str(SignalTrackID) + "; ";
      while (MotherID != 0 && SignalPrimaryTrackIDs.find(MotherID) == SignalPrimaryTrackIDs.end()) { MotherID = SignalMotherMap[MotherID]; }

      if (MotherID != 0) {
        SignalPrimaryMap[SignalTrackID] = MotherID;
      }
      else {
        SignalPrimaryMap[SignalTrackID] = SignalTrackID; // If no mother found, assign itself as primary
      }
    }
    producer->PrintInColor(sSignalTruth, ProducerUtils::GetColor("yellow"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------- PMTrack Analysis -----------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    art::Handle<std::vector<recob::Track>> TrackHandle;
    std::vector<art::Ptr<recob::Track>> TrackList;
    if (evt.getByLabel(fTrackLabel, TrackHandle)) {
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
    if (evt.getByLabel(fOpHitLabel, OpHitHandle)) {
      art::fill_ptr_vector(OpHitList, OpHitHandle);
    }
    // Grab assns with OpHits to get match to neutrino purity
    OpHitNum = int(OpHitList.size());
    if (fGenerateAdjOpFlash) {
      fOpFlashLabel = "solarflash";
      std::vector<AdjOpHitsUtils::FlashInfo> FlashVec;
      adjophits->CalcAdjOpHits(OpHitList, OpHitVec, OpHitIdx, evt);
      adjophits->MakeFlashVector(FlashVec, OpHitVec, evt);
      OpFlashNum = int(FlashVec.size());
      
      for (int i = 0; i < int(FlashVec.size()); i++)
      {
        AdjOpHitsUtils::FlashInfo TheFlash = FlashVec[i];
        double ThisOpFlashPur = 0;
        OpFlashPlane.push_back(TheFlash.Plane);
        OpFlashNHits.push_back(TheFlash.NHit);
        OpFlashTime.push_back(TheFlash.Time - fOpFlashTimeOffset); // Convert to microseconds happens in AdjOpHits
        OpFlashDeltaT.push_back(TheFlash.TimeWidth); // Convert to microseconds
        OpFlashPE.push_back(TheFlash.PE);
        OpFlashMaxPE.push_back(TheFlash.MaxPE);
        OpFlashFast.push_back(TheFlash.FastToTotal);
        OpFlashID.push_back(i);
        OpFlashX.push_back(TheFlash.X);
        OpFlashY.push_back(TheFlash.Y);
        OpFlashZ.push_back(TheFlash.Z);
        OpFlashSTD.push_back(TheFlash.STD);
        for (int j = 0; j < int(OpHitVec[i].size()); j++)
        {
          recob::OpHit OpHit = *OpHitVec[i][j];
          const std::vector<int> ThisOpHitTrackIds = pbt->OpHitToTrackIds(OpHit);
          float ThisOphitPurity = 0;
          for (auto const &ThisOpHitTrackId : ThisOpHitTrackIds)
          {
            if (SignalTrackIDs.find(ThisOpHitTrackId) != SignalTrackIDs.end())
              ThisOphitPurity += 1;
          }
          // Check if ThisOpHitTrackIds is empty
          if (ThisOpHitTrackIds.size() == 0)
            ThisOphitPurity = 0;
          else
            ThisOphitPurity /= int(ThisOpHitTrackIds.size());

          ThisOpFlashPur += ThisOphitPurity * OpHit.PE();
          auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(OpHit.OpChannel()).GetCenter();
          SOpHitPur.push_back(ThisOphitPurity);
          SOpHitChannel.push_back(OpHit.OpChannel());

          if (fOpHitTimeVariable == "StartTime")
            SOpHitTime.push_back(OpHit.StartTime() * clockData.OpticalClock().TickPeriod() - fOpFlashTimeOffset); // Convert to microseconds
          else // Default to PeakTime
            SOpHitTime.push_back(OpHit.PeakTime() * clockData.OpticalClock().TickPeriod() - fOpFlashTimeOffset); // Convert to microseconds

          SOpHitPE.push_back(OpHit.PE());
          SOpHitX.push_back(OpHitXYZ.X());
          SOpHitY.push_back(OpHitXYZ.Y());
          SOpHitZ.push_back(OpHitXYZ.Z());
          SOpHitFlashID.push_back(i);
          SOpHitPlane.push_back(TheFlash.Plane);
        }
        // Check if OpHitVec[i] is empty
        if (OpHitVec[i].size() == 0)
          ThisOpFlashPur = 0;
        else
          ThisOpFlashPur /= TheFlash.PE;

        OpFlashPur.push_back(ThisOpFlashPur);
        if (ThisOpFlashPur > 0) {
          sOpFlashTruth += "OpFlash PE " + ProducerUtils::str(TheFlash.PE) + " with purity " + ProducerUtils::str(ThisOpFlashPur) + " time " + ProducerUtils::str(TheFlash.Time) + " plane " + ProducerUtils::str(TheFlash.Plane) + "\n";
          sOpFlashTruth += " - Vertex (" + ProducerUtils::str(TheFlash.X) + ", " + ProducerUtils::str(TheFlash.Y) + ", " + ProducerUtils::str(TheFlash.Z) + ")\n";
          sOpFlashTruth += "\t*** 1st Sanity check: Ratio " + ProducerUtils::str(TheFlash.MaxPE / TheFlash.PE) + " <= " + ProducerUtils::str(fAdjOpFlashMaxPERatioCut) + "\n";
          sOpFlashTruth += "\t*** 2nd Sanity check: #OpHits " + ProducerUtils::str(int(OpHitVec[i].size())) + " >= " + ProducerUtils::str(TheFlash.NHit) + "\n";
        }
      }
    }
    else {
      std::vector<art::Ptr<recob::OpFlash>> OpFlashList;
      art::Handle<std::vector<recob::OpFlash>> FlashHandle;
      
      if (evt.getByLabel(fOpFlashLabel, FlashHandle)){
        art::fill_ptr_vector(OpFlashList, FlashHandle);
      }
      // Grab assns with OpHits to get match to neutrino purity
      OpFlashNum = int(OpFlashList.size());
      art::FindManyP<recob::OpHit> OpAssns(OpFlashList, evt, fOpFlashLabel);
      // Loop over OpFlashList and assign OpHits to each flash
      for (int i = 0; i < int(OpFlashList.size()); i++)
      {
        recob::OpFlash TheFlash = *OpFlashList[i];
        std::vector<art::Ptr<recob::OpHit>> MatchedHits = OpAssns.at(i);
        int NMatchedHits = MatchedHits.size();
        double FlashStdDev = 0.0, TotalFlashPE = 0, MaxOpHitPE = 0;
        std::vector<float> varXY, varYZ, varXZ;
        varXY = varYZ = varXZ = {};

        for (int j = 0; j < NMatchedHits; j++)
        { // Loop over OpHits in the flash
          recob::OpHit OpHit = *MatchedHits[j];
          art::Ptr<recob::OpHit> OpHitPtr = MatchedHits[j];
          mf::LogDebug("SolarNuAna") << "Assigning OpHit to Flash";
          const std::vector<int> ThisOpHitTrackIds = pbt->OpHitToTrackIds(OpHit);
          float ThisOphitPurity = 0;
          
          for (auto const &ThisOpHitTrackId : ThisOpHitTrackIds)
          {
            if (SignalTrackIDs.find(ThisOpHitTrackId) != SignalTrackIDs.end()) {
              ThisOphitPurity += 1;
            }
          }

          auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(OpHit.OpChannel()).GetCenter();
          TotalFlashPE += OpHit.PE();
          
          varXY.push_back(sqrt(pow(TheFlash.XCenter() - OpHitXYZ.X(), 2) +  pow(TheFlash.YCenter() - OpHitXYZ.Y(), 2)) * OpHit.PE());
          varYZ.push_back(sqrt(pow(TheFlash.YCenter() - OpHitXYZ.Y(), 2) +  pow(TheFlash.ZCenter() - OpHitXYZ.Z(), 2)) * OpHit.PE());
          varXZ.push_back(sqrt(pow(TheFlash.XCenter() - OpHitXYZ.X(), 2) +  pow(TheFlash.ZCenter() - OpHitXYZ.Z(), 2)) * OpHit.PE());
          SOpHitPur.push_back(ThisOphitPurity / int(ThisOpHitTrackIds.size()));
          
          if (OpHit.PE() > MaxOpHitPE) {
            MaxOpHitPE = OpHit.PE();
          }

          SOpHitFlashID.push_back(i);
          SOpHitPE.push_back(OpHit.PE());
          SOpHitX.push_back(OpHitXYZ.X());
          SOpHitY.push_back(OpHitXYZ.Y());
          SOpHitZ.push_back(OpHitXYZ.Z());
          
          if (fOpHitTimeVariable == "StartTime")
            SOpHitTime.push_back(OpHit.StartTime() * clockData.OpticalClock().TickPeriod() - fOpFlashTimeOffset); // Convert to microseconds
          else // Default to PeakTime
            SOpHitTime.push_back(OpHit.PeakTime() * clockData.OpticalClock().TickPeriod() - fOpFlashTimeOffset); // Convert to microseconds
            
          SOpHitChannel.push_back(OpHit.OpChannel());
          SOpHitPlane.push_back(adjophits->GetOpHitPlane(OpHitPtr, 0.01)); // Get plane assignment for the OpHit
        } // End of OpHit loop

        OpHitVec.push_back(MatchedHits);
        FlashStdDev = adjophits->GetOpFlashPlaneSTD(TheFlash.Frame(), varXY, varYZ, varXZ);
        int TerminalOutput = ProducerUtils::supress_stdout();
        double ThisOpFlashPur = pbt->OpHitCollectionPurity(SignalTrackIDs, MatchedHits);
        ProducerUtils::resume_stdout(TerminalOutput);

        // Calculate the flash purity, only for the Signal events
        OpFlashID.push_back(i);
        OpFlashPlane.push_back(TheFlash.Frame());
        OpFlashPur.push_back(ThisOpFlashPur);
        OpFlashMaxPE.push_back(MaxOpHitPE);
        OpFlashSTD.push_back(FlashStdDev);
        OpFlashX.push_back(TheFlash.XCenter());
        OpFlashY.push_back(TheFlash.YCenter());
        OpFlashZ.push_back(TheFlash.ZCenter());
        OpFlashPE.push_back(TheFlash.TotalPE());
        OpFlashNHits.push_back(MatchedHits.size());
        OpFlashFast.push_back(TheFlash.FastToTotal());
        
        if (fOpFlashTime2us) {
          OpFlashTime.push_back(TheFlash.Time() * clockData.OpticalClock().TickPeriod() - fOpFlashTimeOffset); // Expected flash to provide time in ticks, convert to microseconds
          OpFlashDeltaT.push_back(TheFlash.TimeWidth() * clockData.OpticalClock().TickPeriod()); // Expected flash to provide time width in ticks, convert to microseconds
        }
        else {
          OpFlashTime.push_back(TheFlash.Time()); // Expected flash to provide time in microseconds
          OpFlashDeltaT.push_back(TheFlash.TimeWidth()); // Expected flash to provide time width in microseconds
        }

        if (ThisOpFlashPur > 0) {
          mf::LogDebug("SolarNuAna") << "OpFlash PE " << TheFlash.TotalPE() << " with purity " << ThisOpFlashPur << " time " << TheFlash.Time();
          sOpFlashTruth += "OpFlash PE " + ProducerUtils::str(TheFlash.TotalPE()) + " with purity " + ProducerUtils::str(ThisOpFlashPur) + " time " + ProducerUtils::str(TheFlash.Time()) + " plane " + ProducerUtils::str(int(TheFlash.Frame())) + "\n";
          sOpFlashTruth += " - Vertex (" + ProducerUtils::str(TheFlash.XCenter()) + ", " + ProducerUtils::str(TheFlash.YCenter()) + ", " + ProducerUtils::str(TheFlash.ZCenter()) + ")\n";
          sOpFlashTruth += "\t*** 1st Sanity check: Ratio " + ProducerUtils::str(MaxOpHitPE / TotalFlashPE) + " <= " + ProducerUtils::str(fAdjOpFlashMaxPERatioCut) + "\n";
          sOpFlashTruth += "\t*** 2nd Sanity check: #OpHits " + ProducerUtils::str(int(NMatchedHits)) + " >= " + ProducerUtils::str(int(TheFlash.PEs().size())) + "\n";
        }
      }
    }
    sOpFlashTruth = sOpFlashTruth + "\n# of OpFlashes (" + fOpFlashLabel + ") in full geometry: " + ProducerUtils::str(OpFlashNum) + "\n";
    producer->PrintInColor(sOpFlashTruth, ProducerUtils::GetColor("blue"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------- Hit collection and assignment ----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // --- Lift out the reco hits:
    auto RecoHits = evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
    std::vector<art::Ptr<recob::Hit>> RecoHitsPtr;
    int NTotHits = RecoHits->size();

    for (int i = 0; i < NTotHits; ++i)
    {
      // --- Loop over the reconstructed hits to separate them among tpc planes according to view and signal type
      recob::Hit const &ThisHit = RecoHits->at(i);
      // Add to RecoHitsPtr
      RecoHitsPtr.push_back(art::Ptr<recob::Hit>(RecoHits, i));
      if (ThisHit.View() == 0)
      {
        Ind0Hits.push_back(ThisHit);
      } // SignalType = 0
      else if (ThisHit.View() == 1)
      {
        Ind1Hits.push_back(ThisHit);
      } // SignalType = 0
      else if (ThisHit.View() == 2)
      {
        ColHits.push_back(ThisHit);
      } // SignalType = 1
      else
      {
        GhostHits.push_back(ThisHit);
        mf::LogError("SolarNuAna") << "Hit was found with view out of scope";
      }
    }

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //-------------------------------------------------------------- Cluster creation and analysis ------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    std::string sRecoObjects = "";
    std::vector<std::vector<art::Ptr<recob::Hit>>> ClustersPtr;
    std::vector<lowe::LowEUtils::RawPerPlaneCluster> PerPlaneClusters;
    std::vector<std::vector<std::vector<recob::Hit>>> AllPlaneClusters;
    std::vector<std::vector<int>> ClustersIdx = {{}, {}, {}};
    std::vector<std::vector<int>> RecoHitIdx;
    // Map to associate the ClusterIdx with the position in the ClVectors
    std::map<int, std::vector<int>> ClIdxMap;

    if (fGenerateSolarCluster == false) {
      // Get clusters from event recob::Cluster with label "planecluster"
      // ...
    }
    else {
      lowe->CalcAdjHits(RecoHitsPtr, ClustersPtr, RecoHitIdx, evt);
      for (int i = 0; i < int(ClustersPtr.size()); i++)
      {
        std::vector<recob::Hit> ThisHitVector = {}; // Convert pointer to vector
        for (int j = 0; j < int(ClustersPtr[i].size()); j++)
        {
          ThisHitVector.push_back(*ClustersPtr[i][j]);
        }
        int ThisIdx = RecoHitIdx[i][0];
        if (RecoHitsPtr[ThisIdx]->View() == 0)
        {
          Clusters0.push_back(ThisHitVector);
          ClustersIdx[0].push_back(i);
        }
        else if (RecoHitsPtr[ThisIdx]->View() == 1)
        {
          Clusters1.push_back(ThisHitVector);
          ClustersIdx[1].push_back(i);
        }
        else if (RecoHitsPtr[ThisIdx]->View() == 2)
        {
          Clusters2.push_back(ThisHitVector);
          ClustersIdx[2].push_back(i);
        }
        else if (RecoHitsPtr[ThisIdx]->View() == 3)
        {
          Clusters3.push_back(ThisHitVector);
        }
      }
      lowe->MakeClusterVector(PerPlaneClusters, ClustersPtr, evt);
      AllPlaneClusters = {Clusters0, Clusters1, Clusters2};
    }

    std::vector<std::vector<std::vector<float>>> ClVecGenPur = {{}, {}, {}};
    std::vector<std::vector<int>> ClMainID = {{}, {}, {}}, ClTPC = {{}, {}, {}}, ClNHits = {{}, {}, {}}, ClGen = {{}, {}, {}};
    std::vector<std::vector<float>> ClCharge = {{}, {}, {}}, ClMaxCharge = {{}, {}, {}}, ClT = {{}, {}, {}}, ClX = {{}, {}, {}}, ClY = {{}, {}, {}}, ClZ = {{}, {}, {}};
    std::vector<std::vector<float>> ClFracE = {{}, {}, {}}, ClFracGa = {{}, {}, {}}, ClFracNe = {{}, {}, {}}, ClFracRest = {{}, {}, {}};
    std::vector<std::vector<float>> ClPur = {{}, {}, {}}, Cldzdy = {{}, {}, {}}, ClGenPur = {{}, {}, {}};

    sRecoObjects += "\n# OpHits (" + fOpHitLabel + ") in full geometry: " + ProducerUtils::str(OpHitNum);
    sRecoObjects += "\n# OpFlashes (" + fOpFlashLabel + ") in full geometry: " + ProducerUtils::str(OpFlashNum);
    sRecoObjects += "\n# Hits (" + fHitLabel + ") in each view: " + ProducerUtils::str(int(Ind0Hits.size())) + ", " + ProducerUtils::str(int(Ind1Hits.size())) + ", " + ProducerUtils::str(int(ColHits.size())) + ", " + ProducerUtils::str(int(GhostHits.size()));
    sRecoObjects += "\n# Cluster from the hits: " + ProducerUtils::str(int(Clusters0.size())) + ", " + ProducerUtils::str(int(Clusters1.size())) + ", " + ProducerUtils::str(int(Clusters2.size())) + ", " + ProducerUtils::str(int(Clusters3.size()));
    sRecoObjects += "\n# Tracks (" + fTrackLabel + ") in full geometry: " + ProducerUtils::str(TrackNum);
    producer->PrintInColor(sRecoObjects, ProducerUtils::GetColor("cyan"));

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
        ncharge = maxHit = clustT = FracE = FracGa = FracNe = FracRest = clustX = clustY = clustZ = dzdy = 0;
        clustTPC = -1; // Initialize clustTPC to -1
        // Define a vector of floats with size equal to the number of generators + 1
        std::vector<float> VecGenPur(fLabels.size() + 1, 0);

        for (recob::Hit TPCHit : Clusters[i])
        {
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
          clustX += TPCHit.Integral() * hXYZ.X();
          clustY += TPCHit.Integral() * hXYZ.Y();
          clustZ += TPCHit.Integral() * hXYZ.Z();
          clustT += TPCHit.Integral() * TPCHit.PeakTime() * clockData.TPCClock().TickPeriod(); // Convert to microseconds

          if (TPCHit.Integral() > maxHit) { // If clusterTPC not in TPCIDMap, set it to -1
            if (TPCIDMap.find(TPC) == TPCIDMap.end()) {
              clustTPC = -1;
            }
            else {
              clustTPC = TPC;
            }
            // Look for maxHit inside cluster
            maxHit = TPCHit.Integral();            
          }
          else {
            clustTPC = TPC;
          }

          MainTrID = 0;
          double TopEFrac = 0;
          std::vector<sim::TrackIDE> ThisHitIDE = bt_serv->HitToTrackIDEs(clockData, TPCHit);

          for (size_t ideL = 0; ideL < ThisHitIDE.size(); ++ideL)
          {
            if (ThisHitIDE[ideL].energyFrac > TopEFrac) {
              TopEFrac = ThisHitIDE[ideL].energyFrac;
              MainTrID = abs(ThisHitIDE[ideL].trackID);
            }
          }

          for (int frac = 0; frac < int(ClPartTrackIDs.size()); ++frac)
          {
            for (int trck = 0; trck < int(ClPartTrackIDs[frac].size()); ++trck)
            {
              if (MainTrID == ClPartTrackIDs[frac][trck]) {
                if (frac == 0) {
                  FracE = FracE + TPCHit.Integral();
                }
                else if (frac == 1) {
                  FracGa = FracGa + TPCHit.Integral();
                }
                else if (frac == 2) {
                  FracNe = FracNe + TPCHit.Integral();
                }
                else {
                  FracRest = FracRest + TPCHit.Integral();
                }
              }
            }
          }

          long unsigned int GeneratorType = ProducerUtils::WhichGeneratorType(GeneratorParticles, MainTrID);
          VecGenPur[int(GeneratorType)] = VecGenPur[int(GeneratorType)] + TPCHit.Integral();
          if (SignalTrackIDs.find(MainTrID) != SignalTrackIDs.end()) {
            hitCharge = TPCHit.Integral();
            Pur = Pur + hitCharge;
          }
        }

        float MainGenPurity = 0;
        for (size_t genpur = 0; genpur < VecGenPur.size(); genpur++)
        {
          VecGenPur[genpur] = VecGenPur[genpur] / ncharge;
          if (VecGenPur[genpur] > MainGenPurity) {
            MainGenerator = genpur;
            MainGenPurity = VecGenPur[genpur];
          }
        }

        dzdy = thisdzdy[0];
        thisdzdy.clear();
        FracE /= ncharge;
        FracGa /= ncharge;
        FracNe /= ncharge;
        FracRest /= ncharge;
        clustX /= ncharge;
        clustY /= ncharge;
        clustZ /= ncharge;
        clustT /= ncharge;
        Pur /= ncharge;
        
        ClIdxMap[ClustersIdx[idx][i]] = {idx, i}; // Map the cluster index to the plane and cluster number
        ClNHits[idx].push_back(nhit);
        ClCharge[idx].push_back(ncharge);
        ClMaxCharge[idx].push_back(maxHit);
        ClT[idx].push_back(clustT);
        ClTPC[idx].push_back(clustTPC);
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
      }
    } // Finished first cluster processing

    //-------------------------------------------------------------------- Cluster Matching -------------------------------------------------------------------------//
    std::vector<unsigned int> MVecGen = {};
    std::vector<std::vector<float>> MVecGenFrac = {};
    std::vector<float> MVecFracE = {}, MVecFracGa = {}, MVecFracNe = {}, MVecFracRest = {}, MVecGenPur = {};
    std::vector<std::vector<int>>  MatchedClustersIdx = {{}, {}, {}};
    std::vector<std::vector<int>> MVecMainID = {{}, {}, {}}, MVecNHits = {{}, {}, {}}, MVecTPC = {{}, {}, {}}, MVecChannel = {{}, {}, {}};
    std::vector<std::vector<float>> MVecPur = {{}, {}, {}}, MVecMaxCharge = {{}, {}, {}}, MVecCharge = {{}, {}, {}}, MVecTime = {{}, {}, {}}, MVecRecoX = {{}, {}, {}}, MVecRecoY = {{}, {}, {}}, MVecRecoZ = {{}, {}, {}};
    std::vector<std::vector<float>> MVecDirDir = {{}, {}, {}}, MatchedClCompleteness = {{}, {}, {}}, MVecdT = {{}, {}, {}};
    std::vector<solar::LowECluster> SolarClusters;
    std::vector<art::Ptr<solar::LowECluster>> SolarClustersPtr;
    
    std::vector<std::vector<std::vector<recob::Hit>>> MatchedClusters = {{}, {}, {}};
  
    // If present, grab the SolarClusters from the event
    if (fGenerateSolarCluster == false) {
      art::Handle<std::vector<solar::LowECluster>> SolarClusterHandle;
      evt.getByLabel(fSolarClusterLabel, SolarClusterHandle);
      if (SolarClusterHandle.isValid()) {
        for (size_t i = 0; i < SolarClusterHandle->size(); i++) {
          SolarClustersPtr.push_back(art::Ptr<solar::LowECluster>(SolarClusterHandle, i));
        }
      }
      // Requires further implementation to match the previous "planecluster"
      // ...
    }   
    else {
      std::string SolarClusterInfo = "SolarClusterInfo: ";
      SolarClusterInfo = SolarClusterInfo + "(" + ProducerUtils::str(Clusters0.size()) + "," + ProducerUtils::str(Clusters1.size()) + "," + ProducerUtils::str(Clusters2.size()) + ")";
      lowe->MatchClusters(SignalTrackIDs, MatchedClustersIdx, MatchedClusters, ClustersIdx, AllPlaneClusters, MVecMainID, MVecNHits, MVecChannel, MVecTime, MVecRecoY, MVecRecoZ, MVecDirDir, MVecCharge, MVecPur, MatchedClCompleteness, clockData, true);

      SolarClusterInfo = SolarClusterInfo + "\nFound " + ProducerUtils::str(int(MatchedClustersIdx[2].size())) + " MatchedClusters (from col. plane loop)!";
      for (int ThisClIdx = 0; ThisClIdx < int(MatchedClustersIdx[2].size()); ThisClIdx++)
      {
        // MVecTime[2][ThisClIdx] *= clockData.TPCClock().TickPeriod(); // Convert to microseconds
        for (int plane = 0; plane < 2; plane++)
        {
          int RefClIdx = ClIdxMap[MatchedClustersIdx[plane][ThisClIdx]][1]; // Get the cluster index in the plane
          MVecMainID[plane].push_back(ClMainID[plane][RefClIdx]); 
          if (MVecTime[plane][ThisClIdx] > -1e6) {
            // MVecTime[plane][ThisClIdx] *= clockData.TPCClock().TickPeriod(); // Convert to microseconds
            MVecdT[plane].push_back(abs(MVecTime[2][ThisClIdx] - MVecTime[plane][ThisClIdx]));
            MVecTPC[plane].push_back(ClTPC[plane][RefClIdx]);
            MVecMaxCharge[plane].push_back(ClMaxCharge[plane][RefClIdx]);
            SolarClusterInfo = SolarClusterInfo + "\nMatched Cluster in plane " + ProducerUtils::str(plane) + " with time " + ProducerUtils::str(MVecTime[plane][ThisClIdx]) + " and charge " + ProducerUtils::str(MVecCharge[plane][ThisClIdx]) + " with TPC " + ProducerUtils::str(MVecTPC[plane][ThisClIdx]);
          }
          
          else {
            MVecdT[plane].push_back(-1e6);
            MVecTPC[plane].push_back(-1);
            MVecMaxCharge[plane].push_back(-1e6);
            SolarClusterInfo = SolarClusterInfo + "\nMatched Cluster in plane " + ProducerUtils::str(plane) + " with time -1e6 and charge -1e6 with TPC -1";
          }
        }
        int RefClIdx = ClIdxMap[MatchedClustersIdx[2][ThisClIdx]][1]; // Get the plane index of the matched cluster
        MVecMainID[2].push_back(ClMainID[2][RefClIdx]);
        MVecRecoX[2].push_back(ClT[2][RefClIdx] *driftLength/driftTime); // Convert to microseconds and then to cm
        MVecTPC[2].push_back(ClTPC[2][RefClIdx]);
        MVecMaxCharge[2].push_back(ClMaxCharge[2][RefClIdx]);
        MVecGenPur.push_back(ClGenPur[2][RefClIdx]);
        MVecGen.push_back(ClGen[2][RefClIdx]);
        MVecFracE.push_back(ClFracE[2][RefClIdx]);
        MVecFracGa.push_back(ClFracGa[2][RefClIdx]);
        MVecFracNe.push_back(ClFracNe[2][RefClIdx]);
        MVecFracRest.push_back(ClFracRest[2][RefClIdx]);  
        MVecGenFrac.push_back(ClVecGenPur[2][RefClIdx]);
        SolarClusterInfo = SolarClusterInfo + "\nMatched Cluster in plane 2 with time " + ProducerUtils::str(MVecTime[2][ThisClIdx]) + " and charge " + ProducerUtils::str(MVecCharge[2][ThisClIdx]) + " with TPC " + ProducerUtils::str(MVecTPC[2][ThisClIdx]);
      }
      producer->PrintInColor(SolarClusterInfo, ProducerUtils::GetColor("yellow"), "Debug");

      for (int i = 0; i < int(MVecNHits[2].size()); i++)
      {
        if (fClusterPreselectionSignal && MVecPur[2][i] == 0)
        {
          continue;
        }
        std::vector<float> clustPos = {MVecRecoX[2][i], MVecRecoY[2][i], MVecRecoZ[2][i]};
        int clustMainID = MVecMainID[2][i];
        int clustNHits = MVecNHits[2][i];
        int clustChannel = MVecChannel[2][i];
        float clustCharge = MVecCharge[2][i];
        float clustTime = MVecTime[2][i];
        float clustPurity = MVecPur[2][i];
        float clustCompleteness = MatchedClCompleteness[2][i];
        std::vector<recob::Cluster> clustVector = {}; // Vector of recob::Cluster
        // Add clusters according to the indices in MatchedClustersIdx
        for (int plane = 0; plane < 3; plane++)
        {
          int clustIdx = ClIdxMap[MatchedClustersIdx[plane][i]][1];
          if (clustIdx >= 0 && clustIdx < int(AllPlaneClusters[plane].size()))
          {
            // Cluster(float start_wire,float sigma_start_wire,float start_tick,float sigma_start_tick,float start_charge,float start_angle,float start_opening,float end_wire,float sigma_end_wire,float end_tick,float sigma_end_tick,float end_charge,float end_angle,float end_opening,float integral,float integral_stddev,float summedADC,float summedADC_stddev,unsigned int n_hits,float multiple_hit_density,float width,ID_t ID,geo::View_t view,geo::PlaneID const& plane,SentryArgument_t sentry = Sentry);
            recob::Cluster thisCluster(
              ClY[plane][clustIdx], 0, ClT[plane][clustIdx], 0, 
              ClCharge[plane][clustIdx], 0, 0, 
              ClY[plane][clustIdx], 0, ClT[plane][clustIdx], 0, 
              ClCharge[plane][clustIdx], 0, 0, 
              ClCharge[plane][clustIdx], 0, 
              ClCharge[plane][clustIdx], 0,
              int(AllPlaneClusters[plane][clustIdx].size()), 
              0,
              0,
              clustIdx,
              geo::View_t(plane), 
              geo::PlaneID(0, 0, plane), 
              {} // Assuming default for SentryArgument_t
            );
            clustVector.push_back(thisCluster);
          }
        }

        solar::LowECluster ThisSolarCluster(clustPos, clustMainID, clustNHits, clustChannel, clustCharge, clustTime, clustPurity, clustCompleteness, clustVector);
        SolarClusters.push_back(ThisSolarCluster);
      }    
    }

    //-------------------------------------------------------------------- Cluster Tree Export -------------------------------------------------------------------------//
    // Need to implement the primary cluster finding based on external algorithm that uses the solar::LowECluster
    // std::vector<bool> EventCandidateFound = {};
    // std::vector<std::vector<art::Ptr<solar::LowECluster>>> EventCandidateVector;
    // std::vector<std::vector<int>> EventCandidateIdx;
    // lowe->FindPrimaryClusters(SolarClustersPtr, EventCandidateFound, EventCandidateVector, EventCandidateIdx, clockData, evt);

    // For now, loop over matched clusters and export to tree if all conditions are satisfied
    std::string sClustersReco = "\n# ClusterReco: Looping over " + ProducerUtils::str(int(MVecNHits[2].size())) + " matched clusters";
    producer->PrintInColor(sClustersReco, ProducerUtils::GetColor("green"));
    for (int i = 0; i < int(MVecNHits[2].size()); i++)
    {
      if (fClusterPreselectionSignal && MVecPur[2][i] == 0)
      {
        continue;
      }
      bool TrackMatch = false;
      bool AdjClusterMatch = false;
      std::string sFlashReco = "";
      std::string sVertexReco = "";
      std::string sClusterReco = "";
      std::string sResultColor = "white";
      std::string sAdjClusters = "";
      float OpFlashResidual = 0;
      float MatchedOpFlashPE = -1e6;
      float MatchedOpFlashResidual = 1e6;
      float MatchedOpFlashX = -1e6;

      if (MVecNHits[2][i] > fClusterPreselectionNHits)
      {
        MPrimary = true;
        MAdjClNum = 0;
        MSignalAdjClNum = 0;
        MAdjClTime = {};
        MAdjClCharge = {};
        MAdjClInd0Charge = {};
        MAdjClInd1Charge = {};
        MAdjClMaxCharge = {};
        MAdjClInd0MaxCharge = {};
        MAdjClInd1MaxCharge = {};
        MAdjClNHits = {};
        MAdjClInd0NHits = {};
        MAdjClInd1NHits = {};
        MAdjClRecoY = {};
        MAdjClRecoZ = {};
        MAdjClR = {};
        MAdjClPur = {};
        MAdjClGen = {};
        MAdjClGenPur = {};
        MAdjClMainID = {};
        MAdjClMainPDG = {};
        MAdjClMainE = {};
        MAdjClMainP = {};
        MAdjClMainK = {};
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
        MAdjFlashNHits = {};
        MAdjFlashPlane = {};
        MAdjFlashFast = {};
        MAdjFlashMaxPE = {};
        MAdjFlashRecoX = {};
        MAdjFlashRecoY = {};
        MAdjFlashRecoZ = {};
        MAdjFlashResidual = {};
        MTrackStart = {-1e6, -1e6, -1e6};
        MTrackEnd = {-1e6, -1e6, -1e6};

        for (int j = 0; j < int(MVecNHits[2].size()); j++)
        {
          if (j == i) { continue;} // Do not compare the cluster with itself

          double ClusterDistance = 0;
          producer->ComputeDistance3D(ClusterDistance, MVecTime[2][i], MVecRecoY[2][i], MVecRecoZ[2][i], MVecTime[2][j], MVecRecoY[2][j], MVecRecoZ[2][j], TPCIDdriftLength[MVecTPC[2][i]], TPCIDdriftTime[MVecTPC[2][i]]);
          if (MVecCharge[2][j] < fMinClusterCharge) { continue; } // Skip clusters that are too small
          if (ClusterDistance > fAdjClusterRad) { continue; } // Skip clusters that are too far

          sAdjClusters += "    - Cluster " + ProducerUtils::str(j) + " at distance " + ProducerUtils::str(ClusterDistance) + " with time " + ProducerUtils::str(MVecTime[2][j]) + " and charge " + ProducerUtils::str(MVecCharge[2][j]) + " in TPC " + ProducerUtils::str(MVecTPC[2][j]);
          sAdjClusters += " and hits " + ProducerUtils::str(MVecNHits[2][j]) + "\n";

          if (MVecCharge[2][j] > MVecCharge[2][i]) { MPrimary = false; }
          if (MVecGen[i] == MVecGen[j]) { MSignalAdjClNum += 1; }
          MAdjClNum += 1;
          
          // If the cluster is matched, add the information to the vectors
          AdjClusterMatch = true;
          MAdjClTime.push_back(MVecTime[2][j]);
          MAdjClInd0Charge.push_back(MVecCharge[0][j]);
          MAdjClInd1Charge.push_back(MVecCharge[1][j]);
          MAdjClCharge.push_back(MVecCharge[2][j]);
          MAdjClMaxCharge.push_back(MVecMaxCharge[2][j]);
          MAdjClInd0MaxCharge.push_back(MVecMaxCharge[0][j]);
          MAdjClInd1MaxCharge.push_back(MVecMaxCharge[1][j]);
          MAdjClNHits.push_back(MVecNHits[2][j]);
          MAdjClInd0NHits.push_back(MVecNHits[0][j]);
          MAdjClInd1NHits.push_back(MVecNHits[1][j]);
          MAdjClRecoY.push_back(MVecRecoY[2][j]);
          MAdjClRecoZ.push_back(MVecRecoZ[2][j]);
          MAdjClR.push_back(sqrt(pow(MVecRecoY[2][i] - MVecRecoY[2][j], 2) + pow(MVecRecoZ[2][i] - MVecRecoZ[2][j], 2)));
          MAdjClPur.push_back(MVecPur[2][j]);
          MAdjClGen.push_back(MVecGen[j]);
          MAdjClGenPur.push_back(MVecGenPur[j]);
          MAdjClMainID.push_back(MVecMainID[2][j]);

          // If mother exists add the mother information
          const simb::MCParticle *MAdjClTruth;
          int TerminalOutput = ProducerUtils::supress_stdout();
          MAdjClTruth = pi_serv->TrackIdToParticle_P(MVecMainID[2][j]);
          ProducerUtils::resume_stdout(TerminalOutput);
          
          if (MAdjClTruth == 0) {
            MAdjClMainPDG.push_back(0);
            MAdjClMainE.push_back(-1e6);
            MAdjClMainP.push_back(-1e6);
            MAdjClMainK.push_back(-1e6);
            MAdjClMainX.push_back(-1e6);
            MAdjClMainY.push_back(-1e6);
            MAdjClMainZ.push_back(-1e6);
            MAdjClEndX.push_back(-1e6);
            MAdjClEndY.push_back(-1e6);
            MAdjClEndZ.push_back(-1e6);
          }
          else {
            MAdjClMainPDG.push_back(MAdjClTruth->PdgCode());
            MAdjClMainE.push_back(1e3*MAdjClTruth->E());
            MAdjClMainP.push_back(1e3*MAdjClTruth->P());
            MAdjClMainK.push_back(1e3*MAdjClTruth->E() - 1e3*MAdjClTruth->Mass());
            MAdjClMainX.push_back(MAdjClTruth->Vx());
            MAdjClMainY.push_back(MAdjClTruth->Vy());
            MAdjClMainZ.push_back(MAdjClTruth->Vz());
            MAdjClEndX.push_back(MAdjClTruth->EndX());
            MAdjClEndY.push_back(MAdjClTruth->EndY());
            MAdjClEndZ.push_back(MAdjClTruth->EndZ());
          }
        }

        sResultColor = "yellow";
        if (MVecPur[2][i] > 0) {
          sResultColor = "green";
        }

        if (fClusterPreselectionPrimary && !MPrimary) { continue; }

        if (MPrimary) {
          sClusterReco += "*** Matched preselection cluster: " + ProducerUtils::str(i) + "\n";
          sClusterReco += " - MainTrackID " + ProducerUtils::str(MVecMainID[2][i]) + "\n";
          if (MVecGen[i] > 0 && int(MVecGen[i]) < (int(fLabels.size()) + 1)) {
            sClusterReco += " - Gen " + ProducerUtils::str(int(MVecGen[i])) + " -> " + fLabels[MVecGen[i] - 1];
          }
          else {
            sClusterReco += " - Gen ?? -> Unknown";
          }

          sClusterReco += " TPC " + ProducerUtils::str(MVecTPC[2][i]) + "\n";
          sClusterReco += " - Purity " + ProducerUtils::str(MVecGenPur[i]) + " Hits " + ProducerUtils::str(MVecNHits[2][i]) + "\n";
          sClusterReco += " - Charge " + ProducerUtils::str(MVecCharge[2][i]) + " ( MaxHit " + ProducerUtils::str(MVecMaxCharge[2][i]) + " )\n";
          sClusterReco += " - #AdjCl " + ProducerUtils::str(MAdjClNum) + " ( " + ProducerUtils::str(MSignalAdjClNum) + " signal ):\n";
          
          if (AdjClusterMatch) { sClusterReco += sAdjClusters; }
          sClusterReco += " - RecoCol  Time,Y,Z ( " + ProducerUtils::str(MVecTime[2][i]) + ", " + ProducerUtils::str(MVecRecoY[2][i]) + ", " + ProducerUtils::str(MVecRecoZ[2][i]) + " )\n";
          sClusterReco += " - RecoInd0 Time,Y,Z ( " + ProducerUtils::str(MVecTime[0][i]) + ", " + ProducerUtils::str(MVecRecoY[0][i]) + ", " + ProducerUtils::str(MVecRecoZ[0][i]) + " )\n";
          sClusterReco += " - RecoInd1 Time,Y,Z ( " + ProducerUtils::str(MVecTime[1][i]) + ", " + ProducerUtils::str(MVecRecoY[1][i]) + ", " + ProducerUtils::str(MVecRecoZ[1][i]) + " )\n";
          
          if (fSaveTrackInfo) {
            TVector3 ThisClVertex = {0, MVecRecoY[2][i], MVecRecoZ[2][i]};
            float MaxVertexDistance = 10; // if track is further away from ThisClVertex than
            for (int i = 0; i < TrackNum; i++)
            { // using index loop to get track idx
              recob::Track trk = *TrackList[i];
              TVector3 trk_start(0, trk.Start().Y(), trk.Start().Z());
              TVector3 trk_end(0, trk.End().Y(), trk.End().Z());
              // throw away bad tracks
              if ((trk_start - ThisClVertex).Mag() > MaxVertexDistance && (trk_end - ThisClVertex).Mag() > MaxVertexDistance) { continue; }

              MTrackNPoints = trk.NPoints();
              MTrackStart = {trk.Start().X(), trk.Start().Y(), trk.Start().Z()};
              MTrackEnd = {trk.End().X(), trk.End().Y(), trk.End().Z()};
              MTrackChi2 = trk.Chi2();
              
              sClusterReco += "*** Matched pmtrack: \n";
              sClusterReco += " - Track has start ( " + ProducerUtils::str(trk.Start().X()) + ", " + ProducerUtils::str(trk.Start().Y()) + ", " + ProducerUtils::str(trk.Start().Z()) + " )\n";
              sClusterReco += " - Track has end   ( " + ProducerUtils::str(trk.End().X()) + ", " + ProducerUtils::str(trk.End().Y()) + ", " + ProducerUtils::str(trk.End().Z()) + " )\n\n";
              TrackMatch = true;
            }; // Loop over tracks
          }; // if (fSaveTrackInfo)
        }; // if (MPrimary)
        
        if (fClusterPreselectionTrack && !TrackMatch) { continue; }
        
        std::string sFlashMatching = "";
        for (int j = 0; j < int(OpFlashPE.size()); j++)
        {
          // Skip flashes with time outside the cluster time window
          float OpFlashR = -1e6;
          double MAdjFlashX = 0;
          if ((MVecTime[2][i] - OpFlashTime[j]) < 0 || (MVecTime[2][i] - OpFlashTime[j]) > TPCIDdriftTime[MVecTPC[2][i]]) { continue; }          
          producer->ComputeDistanceX(MAdjFlashX, MVecTime[2][i], OpFlashTime[j], TPCIDdriftLength[MVecTPC[2][i]], TPCIDdriftTime[MVecTPC[2][i]]);

          if (fGeometry == "HD" && MVecTPC[2][i]%2 == 0) {
            MAdjFlashX = -MAdjFlashX;
          }
          if (fGeometry == "VD") {
            MAdjFlashX = TPCIDdriftLength[MVecTPC[2][i]] / 2 - MAdjFlashX;
          }

          // Make an eliptical cut on the flash position based on the clusters plane
          if (fGeometry == "HD") {
            OpFlashR = sqrt(pow(MVecRecoY[2][i] - OpFlashY[j], 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2));
            if (pow(MVecRecoY[2][i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1) {
              sFlashMatching += "Skipping flash " + ProducerUtils::str(j) + " at (X,Y,Z) = (" + ProducerUtils::str(OpFlashX[j]) + "," + ProducerUtils::str(OpFlashY[j]) + "," + ProducerUtils::str(OpFlashZ[j]) + ") outside cut at (X,Y,Z) = (" + ProducerUtils::str(MAdjFlashX) + "," + ProducerUtils::str(MVecRecoY[2][i]) + "," + ProducerUtils::str(MVecRecoZ[2][i]) + ") with R = " + ProducerUtils::str(OpFlashR) + " cm\n";
              continue;
            }
          }
          else if (fGeometry == "VD" && OpFlashPlane[j] == 0) { // Cathode flashes
            if (pow(MVecRecoY[2][i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1) { continue; }
            OpFlashR = sqrt(pow(MVecRecoY[2][i] - OpFlashY[j], 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2));
          }
          else if (fGeometry == "VD" && (OpFlashPlane[j] == 1 || OpFlashPlane[j] == 2)) { // Membrane flashes
            if (fAdjOpFlashMembraneProjection) {
              if (MVecRecoY[2][i] * OpFlashY[j] < 0) { continue; } // Only consider clusters and flashes on the same side of the detector
              if (pow(MAdjFlashX - OpFlashX[j], 2) / pow(fAdjOpFlashX, 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1) { continue; }
            }
            else {
              if (pow(MAdjFlashX - OpFlashX[j], 2) / pow(fAdjOpFlashX, 2) + pow(MVecRecoY[2][i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1) { continue; }
            }
            OpFlashR = sqrt(pow(MAdjFlashX - OpFlashX[j], 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2));
          } 
          else if (fGeometry == "VD" && (OpFlashPlane[j] == 3 || OpFlashPlane[j] == 4)) { // End-Cap flashes
            if (fAdjOpFlashEndCapProjection){
              if (MVecRecoZ[2][i] < fidVolZ / 2 && OpFlashPlane[j] == 3) { continue; } // Only consider clusters and flashes on the same half of the volume
              if (MVecRecoZ[2][i] > fidVolZ / 2 && OpFlashPlane[j] == 4) { continue; } // Only consider clusters and flashes on the same half of the volume
              if (pow(MAdjFlashX - OpFlashX[j], 2) / pow(fAdjOpFlashX, 2) + pow(MVecRecoY[2][i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) > 1){ continue; }
            }
            else{
              if (pow(MAdjFlashX - OpFlashX[j], 2) / pow(fAdjOpFlashX, 2) + pow(MVecRecoY[2][i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1){
                continue;
              }
            }
            OpFlashR = sqrt(pow(MAdjFlashX - OpFlashX[j], 2) + pow(MVecRecoY[2][i] - OpFlashY[j], 2));
          }
          else if (fGeometry == "VD" && OpFlashPlane[j] == -1) {
            sFlashMatching += "Skipping flash " + ProducerUtils::str(j) + " with unknown plane " + ProducerUtils::str(OpFlashPlane[j]) + "\n";            
            continue;
          }


          MAdjFlashR.push_back(OpFlashR);
          MAdjFlashPE.push_back(OpFlashPE[j]);
          MAdjFlashTime.push_back(OpFlashTime[j]);
          MAdjFlashNHits.push_back(OpFlashNHits[j]);
          MAdjFlashPlane.push_back(OpFlashPlane[j]);
          MAdjFlashMaxPE.push_back(OpFlashMaxPE[j]);
          MAdjFlashFast.push_back(OpFlashFast[j]);
          MAdjFlashRecoX.push_back(OpFlashX[j]);
          MAdjFlashRecoY.push_back(OpFlashY[j]);
          MAdjFlashRecoZ.push_back(OpFlashZ[j]);
          MAdjFlashSTD.push_back(OpFlashSTD[j]);
          MAdjFlashPur.push_back(OpFlashPur[j]);
          // Compute the residual between the predicted cluster signal and the flash
          adjophits->FlashMatchResidual(OpFlashResidual, OpHitVec[j], MAdjFlashX, double(MVecRecoY[2][i]), double(MVecRecoZ[2][i]));
          // Print the flash information for debugging
          sFlashMatching += "Matching flash " + ProducerUtils::str(j) + " with time " + ProducerUtils::str(OpFlashTime[j]) + " and PE " + ProducerUtils::str(OpFlashPE[j]) + " in plane " + ProducerUtils::str(OpFlashPlane[j]) + " at distance " + ProducerUtils::str(OpFlashR) + " with residual " + ProducerUtils::str(OpFlashResidual) + "\n";
          // Make a cut on the flash PE and MaxPE distributions
          if (OpFlashMaxPE[j] / OpFlashPE[j] > fAdjOpFlashMaxPERatioCut || OpFlashPE[j] < fAdjOpFlashMinPECut)
          {
            continue;
          }
          // Create a cut based on a sigmoid function of the opflash number of hits versus the associated drift time
          if (OpFlashNHits[j] < fAdjOpFlashMinNHitCut)
          {
            continue;
          }
          // If the residual is smaller than the minimum residual, update the minimum residual and the matched flash
          if ((fFlashMatchByResidual && OpFlashResidual < MatchedOpFlashResidual) || (!fFlashMatchByResidual && OpFlashPE[j] > MatchedOpFlashPE))
          {
            MFlashR = OpFlashR;
            MFlashPE = OpFlashPE[j];
            MFlashFast = OpFlashFast[j];
            MFlashNHits = OpFlashNHits[j];
            MFlashPlane = OpFlashPlane[j];
            MFlashMaxPE = OpFlashMaxPE[j];
            MFlashPur = OpFlashPur[j];
            MFlashSTD = OpFlashSTD[j];
            MFlashTime = OpFlashTime[j];
            MFlashRecoX = OpFlashX[j];
            MFlashRecoY = OpFlashY[j];
            MFlashRecoZ = OpFlashZ[j];
            MFlashResidual = OpFlashResidual;
            // Create an output string with the flash information
            sFlashReco = "*** Matched flash: \n - Purity " + ProducerUtils::str(OpFlashPur[j]) +
              " Plane " + ProducerUtils::str(OpFlashPlane[j]) +
              " #Hits " + ProducerUtils::str(OpFlashNHits[j]) +
              " PE " + ProducerUtils::str(OpFlashPE[j]) +
              " MaxPE " + ProducerUtils::str(OpFlashMaxPE[j]) + "\n" +
              " - Time " + ProducerUtils::str(OpFlashTime[j]) +
              " Fast " + ProducerUtils::str(OpFlashFast[j]) +
              " Residual " + ProducerUtils::str(OpFlashResidual) + "\n" +
              " - Reco Time,Y,Z ( " + ProducerUtils::str(MFlashTime) + ", " + ProducerUtils::str(OpFlashY[j]) + ", " + ProducerUtils::str(OpFlashZ[j]) + " )" + "\n";
            MatchedOpFlashX = MAdjFlashX;
            MatchedOpFlashResidual = OpFlashResidual;
            MatchedOpFlashPE = MFlashPE;
          }
          MAdjFlashResidual.push_back(OpFlashResidual);
        }
        
        producer->PrintInColor(sFlashMatching, ProducerUtils::GetColor(sResultColor), "Debug");

        if (fLabels[0] != "marley") { // If not marley, save the true signal particle information based on the mainID
          SignalParticlePDG = SignalPDGMap[SignalPrimaryMap[MVecMainID[2][i]]];
          SignalParticleE = SignalEMap[SignalPrimaryMap[MVecMainID[2][i]]];
          SignalParticleP = SignalPMap[SignalPrimaryMap[MVecMainID[2][i]]];
          SignalParticleK = SignalKMap[SignalPrimaryMap[MVecMainID[2][i]]];
          SignalParticleX = SignalStartXMap[SignalPrimaryMap[MVecMainID[2][i]]];
          SignalParticleY = SignalStartYMap[SignalPrimaryMap[MVecMainID[2][i]]];
          SignalParticleZ = SignalStartZMap[SignalPrimaryMap[MVecMainID[2][i]]];
          SignalParticleTime = SignalTimeMap[SignalPrimaryMap[MVecMainID[2][i]]];
        }

        sVertexReco += "*** Reconstructed Interaction Vertex: \n";
        sVertexReco += " - True X,Y,Z ( " + ProducerUtils::str(SignalStartXMap[SignalPrimaryMap[MVecMainID[2][i]]]) + ", " + ProducerUtils::str(SignalStartYMap[SignalPrimaryMap[MVecMainID[2][i]]]) + ", " + ProducerUtils::str(SignalStartZMap[SignalPrimaryMap[MVecMainID[2][i]]]) + " )" + "\n";
        sVertexReco += " - Main X,Y,Z ( " + ProducerUtils::str(SignalFinalXMap[SignalPrimaryMap[MVecMainID[2][i]]]) + ", " + ProducerUtils::str(SignalFinalYMap[SignalPrimaryMap[MVecMainID[2][i]]]) + ", " + ProducerUtils::str(SignalFinalZMap[SignalPrimaryMap[MVecMainID[2][i]]]) + " )" + "\n";
        sVertexReco += " - EDep X,Y,Z ( " + ProducerUtils::str(SignalMaxEDepXMap[MVecMainID[2][i]]) + ", " + ProducerUtils::str(SignalMaxEDepYMap[MVecMainID[2][i]]) + ", " + ProducerUtils::str(SignalMaxEDepZMap[MVecMainID[2][i]]) + " )" + "\n";
        sVertexReco += " - Reco X,Y,Z ( " + ProducerUtils::str(MatchedOpFlashX) + ", " + ProducerUtils::str(MVecRecoY[2][i]) + ", " + ProducerUtils::str(MVecRecoZ[2][i]) + " )";
        sClusterReco += sFlashReco;
        sClusterReco += sVertexReco;

        if (fClusterPreselectionFlashMatch && MatchedOpFlashPE < 0) {
          continue;
        } 
        // Fill the tree with the cluster information and the adjacent clusters and flashes
        MInd0Pur = MVecPur[0][i];
        MInd1Pur = MVecPur[1][i];
        MPur = MVecPur[2][i];
        MGen = MVecGen[i];
        MGenPur = MVecGenPur[i];
        MGenFrac = MVecGenFrac[i];
        MSignalFrac = {MVecFracE[i], MVecFracGa[i], MVecFracNe[i], MVecFracRest[i]};
        MInd0TPC = MVecTPC[0][i];
        MInd1TPC = MVecTPC[1][i];
        MTPC = MVecTPC[2][i];
        MInd0Charge = MVecCharge[0][i];
        MInd1Charge = MVecCharge[1][i];
        MCharge = MVecCharge[2][i];
        MInd0MaxCharge = MVecMaxCharge[0][i];
        MInd1MaxCharge = MVecMaxCharge[1][i];
        MMaxCharge = MVecMaxCharge[2][i];
        MInd0NHits = MVecNHits[0][i];
        MInd1NHits = MVecNHits[1][i];
        MNHit = MVecNHits[2][i];
        MInd0dTime = MVecdT[0][i];
        MInd1dTime = MVecdT[1][i];
        MTime = MVecTime[2][i];
        MRecX = MatchedOpFlashX;
        MInd0RecoY = MVecRecoY[0][i];
        MInd1RecoY = MVecRecoY[1][i];
        MRecY = MVecRecoY[2][i];
        MRecZ = MVecRecoZ[2][i];
        MMainID = MVecMainID[2][i];
        // If mother exists add the mother information
        const simb::MCParticle *MClTruth;
        
        int TerminalOutput = ProducerUtils::supress_stdout();
        MClTruth = pi_serv->TrackIdToParticle_P(MVecMainID[2][i]);
        ProducerUtils::resume_stdout(TerminalOutput);
        
        if (MClTruth == 0) {
          MMainVertex = {-1e6, -1e6, -1e6};
          MEndVertex = {-1e6, -1e6, -1e6};
          MMainPDG = 0;
          MMainE = -1e6;
          MMainP = -1e6;
          MMainK = -1e6;
          MMainTime = -1e6;
          MMainParentVertex = {-1e6, -1e6, -1e6};
          MMainParentPDG = 0;
          MMainParentE = -1e6;
          MMainParentP = -1e6;
          MMainParentK = -1e6;
          MMainParentTime = -1e6;
        }
        else {
          if (MFlashPur > 0) {
            MFlashCorrect = true;
          }
          MMainVertex = {MClTruth->Vx(), MClTruth->Vy(), MClTruth->Vz()};
          MEndVertex = {MClTruth->EndX(), MClTruth->EndY(), MClTruth->EndZ()};
          MMainPDG = MClTruth->PdgCode();
          MMainE = 1e3 * MClTruth->E();
          MMainP = 1e3 * MClTruth->P();
          MMainK = MMainE - 1e3 * MClTruth->Mass();
          MMainTime = MClTruth->T();
          // If exists add the parent information
          const simb::MCParticle *MClParentTruth;

          int TerminalOutput = ProducerUtils::supress_stdout();
          MClParentTruth = pi_serv->TrackIdToParticle_P(MClTruth->Mother());
          ProducerUtils::resume_stdout(TerminalOutput);
          
          if (MClParentTruth == 0) {
            MMainParentVertex = {-1e6, -1e6, -1e6};
            MMainParentPDG = 0;
            MMainParentE = -1e6;
            MMainParentP = -1e6;
            MMainParentK = -1e6;
            MMainParentTime = -1e6;
          }
          else {
            MMainParentVertex = {MClParentTruth->Vx(), MClParentTruth->Vy(), MClParentTruth->Vz()};
            MMainParentPDG = MClParentTruth->PdgCode();
            MMainParentE = 1e3 * MClParentTruth->E();
            MMainParentP = 1e3 * MClParentTruth->P();
            MMainParentK = MMainParentE - 1e3 * MClParentTruth->Mass();
            MMainParentTime = MClParentTruth->T();
          }
        }

        fSolarNuAnaTree->Fill();
        hDriftTime->Fill(MainElectronEndPointX, MTime);
        hXTruth->Fill(MRecX - SignalParticleX, SignalParticleX);
        hYTruth->Fill(MRecY - SignalParticleY, SignalParticleY);
        hZTruth->Fill(MRecZ - SignalParticleZ, SignalParticleZ);
      }
      // Check if the string sClusterReco is not empty and print it in color
      if (sClusterReco != "") { producer->PrintInColor(sClusterReco, ProducerUtils::GetColor(sResultColor)); }
    }
    fMCTruthTree->Fill();
    producer->PrintInColor("-----------------------------------------------------------------------------------------\n", ProducerUtils::GetColor("green"));
  }

  //......................................................................................................................//
  // Reset variables for each event
  void SolarNuAna::ResetVariables()
  {
    TrackNum = 0;
    OpHitNum = 0;
    OpFlashNum = 0;
    MTrackNPoints = 0;
    MTrackChi2 = 0;
    MFlashPE = -1e6;
    MFlashFast = -1e6;
    MFlashNHits = -1;
    MFlashPlane = -1;
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
    SignalParticleE = 0;
    SignalParticleP = 0;
    SignalParticleK = 0;
    SignalParticleX = 0;
    SignalParticleY = 0;
    SignalParticleZ = 0;
    SignalParticlePDG = 0;
    SignalParticleTime = 0;
    OpFlashPur.clear();
    OpFlashID.clear();
    OpFlashPE.clear();
    OpFlashSTD.clear();
    OpFlashFast.clear();
    OpFlashMaxPE.clear();
    OpFlashX.clear();
    OpFlashY.clear();
    OpFlashZ.clear();
    OpFlashTime.clear();
    OpFlashDeltaT.clear();
    OpFlashNHits.clear();
    OpFlashPlane.clear();
    HitNum = {};
    ClusterNum = {};
    SOpHitPlane = {};
    SignalPDGList = {};
    SignalPDGDepList = {};
    SignalMotherList = {};
    SignalElectronDepList = {};
    TPart = {}, GeneratorParticles = {};
    SignalIDList = {}, SignalIDDepList = {};
    SignalEDepList = {}, SignalXDepList = {}, SignalYDepList = {}, SignalZDepList = {};
    SignalMaxEDepList = {}, SignalMaxEDepXList = {}, SignalMaxEDepYList = {}, SignalMaxEDepZList = {};
    SOpHitChannel = {}, SOpHitPur = {}, SOpHitPE = {}, SOpHitX = {}, SOpHitY = {}, SOpHitZ = {}, SOpHitTime = {}, SOpHitFlashID = {};
    SignalEList = {}, SignalPList = {}, SignalKList = {}, SignalTimeList = {}, SignalEndXList = {}, SignalEndYList = {}, SignalEndZList = {};
  }
} // namespace solar
DEFINE_ART_MODULE(solar::SolarNuAna)
