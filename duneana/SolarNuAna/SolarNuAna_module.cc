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
    std::string fSignalLabel, fGeometry;
    int fDetectorSizeY, fDetectorSizeZ, fClusterAlgoAdjChannel, fClusterInd0MatchTime, fClusterInd1MatchTime, fClusterPreselectionNHits, fAdjOpFlashMinNHitCut;
    float fClusterMatchTime, fAdjClusterRad, fMinClusterCharge, fClusterMatchCharge, fAdjOpFlashX, fAdjOpFlashY, fAdjOpFlashZ, fAdjOpFlashMaxPERatioCut, fAdjOpFlashMinPECut, fClusterMatchNHit, fClusterAlgoTime;
    std::vector<std::string> fLabels, fBackgroundLabels;
    float fOpFlashAlgoMinTime, fOpFlashAlgoMaxTime, fOpFlashAlgoRad, fOpFlashAlgoPE, fOpFlashAlgoTriggerPE, fOpFlashAlgoHotVertexThld;
    bool fClusterPreselectionSignal, fClusterPreselectionPrimary, fClusterPreselectionTrack, fClusterPreselectionFlashMatch;
    bool fGenerateAdjCluster, fGenerateAdjOpFlash, fFlashMatchByResidual;
    bool fSaveSignalDaughters, fSaveSignalEDep, fSaveSignalOpHits, fSaveOpFlashInfo, fSaveTrackInfo;
    bool fAdjOpFlashMembraneProjection, fAdjOpFlashEndCapProjection; // If true, the TPC reco is projected to the membrane plane. If false, apply a 3D constraint dT, Y, Z.
    bool fTestLatestFeatures;
    // bool fOpFlashAlgoCentroid;

    // --- Our TTrees, and its associated variables.
    TTree *fConfigTree;
    TTree *fMCTruthTree;
    TTree *fSolarNuAnaTree;
    std::string TNuInteraction;
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
    std::map<unsigned int, geo::TPCID> TPCIDMap; // Key is the TPC index, value is the TPCID object
    std::map<unsigned int, float> TPCIDdriftLength; // Key is the TPC index, value is the drift length in cm
    std::map<unsigned int, float> TPCIDdriftTime; // Key is the TPC index, value is the drift time in us

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
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    std::unique_ptr<producer::ProducerUtils> producer;
    std::unique_ptr<solar::AdjOpHitsUtils> adjophits;
    std::unique_ptr<solar::LowEUtils> lowe;
  };
#endif

  //......................................................
  SolarNuAna::SolarNuAna(fhicl::ParameterSet const &p)
      : EDAnalyzer(p),
        producer(new producer::ProducerUtils(p)),
        adjophits(new solar::AdjOpHitsUtils(p)),
        lowe(new solar::LowEUtils(p))
  {
    this->reconfigure(p);
  }

  //......................................................
  void SolarNuAna::reconfigure(fhicl::ParameterSet const &p)
  {
    fSignalLabel = p.get<std::string>("SignalLabel");
    fBackgroundLabels = p.get<std::vector<std::string>>("BackgroundLabelVector");
    fHitLabel = p.get<std::string>("HitLabel");
    fOpFlashLabel = p.get<std::string>("OpFlashLabel");
    fOpHitLabel = p.get<std::string>("OpHitLabel");
    fTrackLabel = p.get<std::string>("TrackLabel");
    fGEANTLabel = p.get<std::string>("GEANT4Label");
    fGeometry = p.get<std::string>("Geometry");
    fDetectorSizeY = p.get<int>("DetectorSizeY");
    fDetectorSizeZ = p.get<int>("DetectorSizeZ");
    fClusterAlgoTime = p.get<float>("ClusterAlgoTime");
    fClusterAlgoAdjChannel = p.get<int>("ClusterAlgoAdjChannel");
    fClusterMatchNHit = p.get<float>("ClusterMatchNHit");
    fClusterMatchCharge = p.get<float>("ClusterMatchCharge");
    fClusterMatchTime = p.get<float>("ClusterMatchTime");
    fClusterInd0MatchTime = p.get<float>("ClusterInd0MatchTime");
    fClusterInd1MatchTime = p.get<float>("ClusterInd1MatchTime");
    fClusterPreselectionSignal = p.get<bool>("ClusterPreselectionSignal");
    fClusterPreselectionPrimary = p.get<bool>("ClusterPreselectionPrimary");
    fClusterPreselectionNHits = p.get<int>("ClusterPreselectionNHits");
    fClusterPreselectionTrack = p.get<bool>("ClusterPreselectionTrack");
    fClusterPreselectionFlashMatch = p.get<bool>("ClusterPreselectionFlashMatch");
    fGenerateAdjCluster = p.get<bool>("GenerateAdjCluster");
    fAdjClusterRad = p.get<float>("AdjClusterRad");
    fMinClusterCharge = p.get<float>("MinClusterCharge");
    fGenerateAdjOpFlash = p.get<bool>("GenerateAdjOpFlash");
    fOpFlashAlgoMinTime = p.get<double>("OpFlashAlgoMinTime");
    fOpFlashAlgoMaxTime = p.get<double>("OpFlashAlgoMaxTime");
    fOpFlashAlgoRad = p.get<double>("OpFlashAlgoRad");
    fOpFlashAlgoPE = p.get<float>("OpFlashAlgoPE");
    fOpFlashAlgoTriggerPE = p.get<float>("OpFlashAlgoTriggerPE");
    fOpFlashAlgoHotVertexThld = p.get<float>("OpFlashAlgoHotVertexThld");
    fAdjOpFlashMembraneProjection = p.get<bool>("AdjOpFlashMembraneProjection");
    fAdjOpFlashEndCapProjection = p.get<bool>("AdjOpFlashEndCapProjection");
    // fOpFlashAlgoCentroid = p.get<bool>("OpFlashAlgoCentroid");
    fAdjOpFlashX = p.get<float>("AdjOpFlashX");
    fAdjOpFlashY = p.get<float>("AdjOpFlashY");
    fAdjOpFlashZ = p.get<float>("AdjOpFlashZ");
    fAdjOpFlashMaxPERatioCut = p.get<float>("AdjOpFlashMaxPERatioCut");
    fAdjOpFlashMinPECut = p.get<float>("AdjOpFlashMinPECut");
    fAdjOpFlashMinNHitCut = p.get<int>("AdjOpFlashMinNHitCut");
    fFlashMatchByResidual = p.get<bool>("FlashMatchByResidual");
    fSaveSignalDaughters = p.get<bool>("SaveSignalDaughters");
    fSaveSignalEDep = p.get<bool>("SaveSignalEDep");
    fSaveSignalOpHits = p.get<bool>("SaveSignalOpHits");
    fSaveOpFlashInfo = p.get<bool>("SaveOpFlashInfo");
    fSaveTrackInfo = p.get<bool>("SaveTrackInfo");
    fTestLatestFeatures = p.get<bool>("TestLatestFeatures");
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
    fConfigTree->Branch("GEANT4Label", &fGEANTLabel);
    fConfigTree->Branch("HitLabel", &fHitLabel);
    fConfigTree->Branch("TrackLabel", &fTrackLabel);
    fConfigTree->Branch("OpHitLabel", &fOpHitLabel);
    fConfigTree->Branch("OpFlashLabel", &fOpFlashLabel);
    fConfigTree->Branch("Geometry", &fGeometry);
    fConfigTree->Branch("DetectorSizeY", &fDetectorSizeY);
    fConfigTree->Branch("DetectorSizeZ", &fDetectorSizeZ);
    fConfigTree->Branch("ClusterAlgoTime", &fClusterAlgoTime);
    fConfigTree->Branch("ClusterAlgoAdjChannel", &fClusterAlgoAdjChannel);
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
    fConfigTree->Branch("OpFlashAlgoMinTime", &fOpFlashAlgoMinTime);
    fConfigTree->Branch("OpFlashAlgoMaxTime", &fOpFlashAlgoMaxTime);
    fConfigTree->Branch("OpFlashAlgoRad", &fOpFlashAlgoRad);
    fConfigTree->Branch("OpFlashAlgoPE", &fOpFlashAlgoPE);
    fConfigTree->Branch("OpFlashAlgoTriggerPE", &fOpFlashAlgoTriggerPE);
    fConfigTree->Branch("OpFlashAlgoHotVertexThld", &fOpFlashAlgoHotVertexThld);
    fConfigTree->Branch("AdjOpFlashMembraneProjection", &fAdjOpFlashMembraneProjection);
    fConfigTree->Branch("AdjOpFlashEndCapProjection", &fAdjOpFlashEndCapProjection);
    // fConfigTree->Branch("OpFlashAlgoCentroid", &fOpFlashAlgoCentroid);
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
    fConfigTree->Branch("TestLatestFeatures", &fTestLatestFeatures);

    // MC Truth info.
    fMCTruthTree->Branch("Event", &Event, "Event/I");                                        // Event number
    fMCTruthTree->Branch("Flag", &Flag, "Flag/I");                                           // Flag used to match truth with reco tree entries
    fMCTruthTree->Branch("TruthPart", &TPart);                                               // Number particles per generator
    fMCTruthTree->Branch("Interaction", &TNuInteraction);                                    // True signal interaction process
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
    fMCTruthTree->Branch("HitNum", &HitNum);                                                 // Number of hits in each TPC plane
    fMCTruthTree->Branch("ClusterNum", &ClusterNum);                                         // Number of clusters in each TPC plane
    fMCTruthTree->Branch("TrackNum", &TrackNum, "TrackNum/I");                               // Number of PMTracks
    if (fSaveSignalDaughters)
    { // Save Signal Daughters. (Only makes sense for marley)
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
    fSolarNuAnaTree->Branch("Interaction", &TNuInteraction);                                    // True signal interaction process
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
    std::vector<recob::Hit> Ind0Hits, Ind1Hits, ColHits, GhostHits;
    std::vector<std::vector<recob::Hit>> Clusters0, Clusters1, Clusters2, Clusters3;

    // --- We want to reset all of our previous run and TTree variables ---
    ResetVariables();
    ThisGeneratorParts.clear();
    Event = evt.event();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    Flag = rand() % 10000000000;
    geo::CryostatID c(0);
    
    const geo::CryostatGeo& cryostat = geom->Cryostat(c);

    // Loop over all TPCs in the cryostat and fill the map
    std::string sTPCMap = "";
    unsigned int maxTPC = 0;
    for (auto const& tpcid : geom->Iterate<geo::TPCID>()) {
      if (tpcid.isValid) {
        // Fill the TPC map with the TPC ID
        const double driftLength = cryostat.TPC(tpcid).DriftDistance();
        const double driftTime = driftLength / art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData).DriftVelocity();
        TPCIDMap[tpcid.TPC] = tpcid;
        TPCIDdriftLength[tpcid.TPC] = driftLength;
        TPCIDdriftTime[tpcid.TPC] = driftTime;
        if (tpcid.TPC > maxTPC) {
          maxTPC = tpcid.TPC; // Keep track of the maximum TPC ID
        }
        sTPCMap += "Found TPC ID: " + std::to_string(tpcid.TPC) + " in Cryostat: " + std::to_string(c.Cryostat) + 
                   " with Drift Length: " + ProducerUtils::str(driftLength) + 
                   " cm and Drift Time: " + ProducerUtils::str(driftTime) + " us\n";
      }
    }
    // Add extra TPCID entry -1 for all clusters that are not associated with a TPC
    TPCIDMap[-1] = geo::TPCID(); // Invalid TPC ID
    // Set the drift length and time for the invalid TPC ID to the first valid TPC ID
    TPCIDdriftLength[-1] = TPCIDdriftLength.begin()->second; // Use the first valid TPC's drift length
    TPCIDdriftTime[-1] = TPCIDdriftTime.begin()->second;     // Use the first valid TPC's drift time
    producer->PrintInColor(sTPCMap, ProducerUtils::GetColor("yellow"), "Debug");



    std::string sHead = "";
    sHead = sHead + "\n#########################################";
    sHead = sHead + "\nEvent: " + ProducerUtils::str(Event) + " Flag: " + ProducerUtils::str(Flag);
    sHead = sHead + "\nTPC Map: " + ProducerUtils::str(TPCIDMap.size()) + " TPCs found";
    sHead = sHead + "\nTPC Frequency in [MHz]: " + ProducerUtils::str(clockData.TPCClock().Frequency());
    sHead = sHead + "\nTPC Tick in [us]: " + ProducerUtils::str(clockData.TPCClock().TickPeriod());
    sHead = sHead + "\nTPC DriftLength in [cm]: " + ProducerUtils::str(TPCIDdriftLength[0]);
    sHead = sHead + "\nTPC DriftTime in [us]: " + ProducerUtils::str(TPCIDdriftTime[0]);
    sHead = sHead + "\n#########################################";
    producer->PrintInColor(sHead, ProducerUtils::GetColor("magenta"));

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Create maps for ID tracking -----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // --- Fill MC Truth IDs to tracking vectors. Get a list of all of my particles in one chunk. ---
    const sim::ParticleList &PartList = pi_serv->ParticleList();
    std::string sMcTruth = "";
    sMcTruth = sMcTruth + "\nThere are a total of " + ProducerUtils::str(int(PartList.size())) + " Particles in the event\n";

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
        producer->FillMyMaps(GeneratorParticles[i], Assn, ThisValidHanlde);                          // Fill empty list with previously assigned particles
        if (GeneratorParticles[i].size() < 1000)
        {
          sMcTruth = sMcTruth + "\n# of particles " + ProducerUtils::str(int(GeneratorParticles[i].size())) + "\tfrom gen " + ProducerUtils::str(int(i) + 1) + " " + fLabels[i];
        }
        else
        {
          sMcTruth = sMcTruth + "\n# of particles " + ProducerUtils::str(int(GeneratorParticles[i].size())) + "\tfrom gen " + ProducerUtils::str(int(i) + 1) + " " + fLabels[i];
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
        sSignalTruth = sSignalTruth + "\nNumber of Signal Particles: " + ProducerUtils::str(NSignalParticles);
        if (fLabels[0] == "marley")
        {
          const simb::MCNeutrino &nue = SignalTruth.GetNeutrino();
          TNuInteraction = ProducerUtils::str(nue.InteractionType());
          SignalParticleE = 1e3 * nue.Nu().E();
          SignalParticleP = 1e3 * nue.Nu().P();
          SignalParticleK = 1e3 * nue.Nu().E() - 1e3 * nue.Nu().Mass();
          SignalParticleX = nue.Nu().Vx();
          SignalParticleY = nue.Nu().Vy();
          SignalParticleZ = nue.Nu().Vz();
          SignalParticlePDG = nue.Nu().PdgCode();
          SignalParticleTime = nue.Nu().T();
          sSignalTruth = sSignalTruth + "\nNeutrino Interaction: " + TNuInteraction;
          sSignalTruth = sSignalTruth + "\nNeutrino Energy: " + ProducerUtils::str(SignalParticleE) + " MeV";
          sSignalTruth = sSignalTruth + "\nPosition (" + ProducerUtils::str(SignalParticleX) + ", " + ProducerUtils::str(SignalParticleY) + ", " + ProducerUtils::str(SignalParticleZ) + ") cm";
        }
        if (fLabels[0] == "generator")
        {
          sSignalTruth = sSignalTruth + "\nFound generator label: " + fLabels[0] + ". Using single generator config.\n";
          if (NSignalParticles > 1)
          {
            sSignalTruth = sSignalTruth + "\n[WARNING] Multiple particles found in the Signal MCTruth. Using the first one.\n";
          }
          const simb::MCParticle &SignalParticle = SignalTruth.GetParticle(0);
          SignalParticleE = 1e3 * SignalParticle.E();
          SignalParticleP = 1e3 * SignalParticle.P();
          SignalParticleK = 1e3 * SignalParticle.E() - 1e3 * SignalParticle.Mass();
          SignalParticleX = SignalParticle.Vx();
          SignalParticleY = SignalParticle.Vy();
          SignalParticleZ = SignalParticle.Vz();
          SignalParticlePDG = SignalParticle.PdgCode();
          SignalParticleTime = SignalParticle.T();
          std::string sSignalParticle = "";
          if (abs(SignalParticle.PdgCode()) == 12)
          {
            sSignalParticle = "Neutrino";
          }
          else if (abs(SignalParticle.PdgCode()) == 11)
          {
            sSignalParticle = "Electron";
          }
          else if (abs(SignalParticle.PdgCode()) == 22)
          {
            sSignalParticle = "Photon";
          }
          else if (abs(SignalParticle.PdgCode()) == 2112)
          {
            sSignalParticle = "Neutron";
          }
          else
          {
            sSignalParticle = "Other";
          }
          TNuInteraction = "Single " + sSignalParticle;
          sSignalTruth = sSignalTruth + "\n" +  sSignalParticle + " Energy: " + ProducerUtils::str(SignalParticleE) + " MeV";
          sSignalTruth = sSignalTruth + "\nPosition (" + ProducerUtils::str(SignalParticleX) + ", " + ProducerUtils::str(SignalParticleY) + ", " + ProducerUtils::str(SignalParticleZ) + ") cm\n";
        }
      }
      art::FindManyP<simb::MCParticle> SignalAssn(Signal, evt, fGEANTLabel);
      sSignalTruth = sSignalTruth + "\nGen.\tPdgCode\t\tEnergy\t\tEndPosition\t\tMother";
      sSignalTruth = sSignalTruth + "\n------------------------------------------------------------------------";

      for (size_t i = 0; i < SignalAssn.size(); i++)
      {
        auto SignalParticles = SignalAssn.at(i);
        for (auto SignalParticle = SignalParticles.begin(); SignalParticle != SignalParticles.end(); SignalParticle++)
        {
          if (fSaveSignalDaughters)
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
          }
          std::map<int, float> SignalMaxEDepMap, SignalMaxEDepXMap, SignalMaxEDepYMap, SignalMaxEDepZMap;
          std::vector<const sim::IDE *> ides = bt_serv->TrackIdToSimIDEs_Ps((*SignalParticle)->TrackId());
          for (auto const &ide : ides)
          {
            if (ide->numElectrons < 1 || ide->energy < 1e-6 || abs(ide->x) > TPCIDdriftLength[0] || abs(ide->y) > fDetectorSizeY || abs(ide->z) > fDetectorSizeZ)
            {
              continue;
            } 
            if (ProducerUtils::InMyMap((*SignalParticle)->TrackId(), SignalMaxEDepMap) == false)
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
          SignalMaxEDepList.push_back(SignalMaxEDepMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepXList.push_back(SignalMaxEDepXMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepYList.push_back(SignalMaxEDepYMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepZList.push_back(SignalMaxEDepZMap[(*SignalParticle)->TrackId()]);
          SignalTrackIDs.emplace((*SignalParticle)->TrackId());

          if ((*SignalParticle)->PdgCode() < 1000000)
          {
            sSignalTruth = sSignalTruth + "\n" + fLabels[0] + "\t" + ProducerUtils::str((*SignalParticle)->PdgCode()) + "\t\t" + ProducerUtils::str(1e3 * (*SignalParticle)->E()) + "\t (" + ProducerUtils::str((*SignalParticle)->EndX()) + ", " + ProducerUtils::str((*SignalParticle)->EndY()) + ", " + ProducerUtils::str((*SignalParticle)->EndZ()) + ")\t" + ProducerUtils::str((*SignalParticle)->Mother());
          }
          else
          {
            sSignalTruth = sSignalTruth + "\n" + fLabels[0] + "\t" + ProducerUtils::str((*SignalParticle)->PdgCode()) + "\t" + ProducerUtils::str(1e3 * (*SignalParticle)->E()) + " (" + ProducerUtils::str((*SignalParticle)->EndX()) + ", " + ProducerUtils::str((*SignalParticle)->EndY()) + ", " + ProducerUtils::str((*SignalParticle)->EndZ()) + ")\t" + ProducerUtils::str((*SignalParticle)->Mother());
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
      sSignalTruth += ProducerUtils::str(SignalTrackID) + "; ";
    }
    producer->PrintInColor(sSignalTruth, ProducerUtils::GetColor("yellow"));

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
        AdjOpHitsUtils::FlashInfo TheFlash = FlashVec[i];
        double ThisOpFlashPur = 0;
        OpFlashPlane.push_back(TheFlash.Plane);
        OpFlashNHits.push_back(TheFlash.NHit);
        OpFlashTime.push_back(TheFlash.Time * clockData.TPCClock().TickPeriod()); // Convert to microseconds
        OpFlashDeltaT.push_back(TheFlash.TimeWidth);
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
          SOpHitTime.push_back(OpHit.PeakTime());
          SOpHitPE.push_back(OpHit.PE());
          SOpHitX.push_back(OpHitXYZ.X());
          SOpHitY.push_back(OpHitXYZ.Y());
          SOpHitZ.push_back(OpHitXYZ.Z());
          SOpHitFlashID.push_back(i);
          SOpHitPlane.push_back(TheFlash.Plane);
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
        if (abs(TheFlash.Time) < 3)
        {
          mf::LogDebug("SolarNuAna") << "Signal OpFlash PE (fast/ratio/tot/STD) " << TheFlash.FastToTotal << "/" << TheFlash.MaxPE / TheFlash.PE << "/" << TheFlash.PE << "/" << TheFlash.STD << " with purity " << ThisOpFlashPur << " time " << TheFlash.Time;
          sOpFlashTruth += "OpFlash PE " + ProducerUtils::str(TheFlash.PE) + " with purity " + ProducerUtils::str(ThisOpFlashPur) + " time " + ProducerUtils::str(TheFlash.Time) + " plane " + ProducerUtils::str(TheFlash.Plane) + "\n";
          sOpFlashTruth += " - Vertex (" + ProducerUtils::str(TheFlash.X) + ", " + ProducerUtils::str(TheFlash.Y) + ", " + ProducerUtils::str(TheFlash.Z) + ")\n";
          sOpFlashTruth += "\t*** 1st Sanity check: Ratio " + ProducerUtils::str(TheFlash.MaxPE / TheFlash.PE) + " <= " + ProducerUtils::str(fAdjOpFlashMaxPERatioCut) + "\n";
          sOpFlashTruth += "\t*** 2nd Sanity check: #OpHits " + ProducerUtils::str(int(OpHitVec[i].size())) + " >= " + ProducerUtils::str(TheFlash.NHit) + "\n";
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
        double FlashStdDev = 0.0, TotalFlashPE = 0, MaxOpHitPE = 0, FlashTime = 0;
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
            if (SignalTrackIDs.find(ThisOpHitTrackId) != SignalTrackIDs.end())
            {
              ThisOphitPurity += 1;
            }
          }
          auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(OpHit.OpChannel()).GetCenter();
          TotalFlashPE += OpHit.PE();
          varXY.push_back(sqrt(pow(TheFlash.XCenter() - OpHitXYZ.X(), 2) +  pow(TheFlash.YCenter() - OpHitXYZ.Y(), 2)) * OpHit.PE());
          varYZ.push_back(sqrt(pow(TheFlash.YCenter() - OpHitXYZ.Y(), 2) +  pow(TheFlash.ZCenter() - OpHitXYZ.Z(), 2)) * OpHit.PE());
          varXZ.push_back(sqrt(pow(TheFlash.XCenter() - OpHitXYZ.X(), 2) +  pow(TheFlash.ZCenter() - OpHitXYZ.Z(), 2)) * OpHit.PE());
          FlashTime += OpHit.PeakTime() * OpHit.PE();
          SOpHitPur.push_back(ThisOphitPurity / int(ThisOpHitTrackIds.size()));
          if (OpHit.PE() > MaxOpHitPE)
          {
            MaxOpHitPE = OpHit.PE();
          };
          SOpHitChannel.push_back(OpHit.OpChannel());
          SOpHitTime.push_back(OpHit.PeakTime());
          SOpHitPE.push_back(OpHit.PE());
          SOpHitX.push_back(OpHitXYZ.X());
          SOpHitY.push_back(OpHitXYZ.Y());
          SOpHitZ.push_back(OpHitXYZ.Z());
          SOpHitFlashID.push_back(i);
          SOpHitPlane.push_back(adjophits->GetOpHitPlane(OpHitPtr));
        } // End of OpHit loop

        // Check that all OpHits are assigned to the same plane
        int FlashPlane = -1;
        for (int j = 0; j < int(SOpHitPlane.size()); j++)
        {
          if (SOpHitPlane[j] != SOpHitPlane[0])
          {
            FlashPlane = -1;
            mf::LogError("SolarNuAna") << "OpHits are not assigned to the same plane!";
          }
          else
          {
            FlashPlane = SOpHitPlane[0];
          }
        }

        OpHitVec.push_back(MatchedHits);
        FlashTime = FlashTime / TotalFlashPE;
        
        FlashStdDev = adjophits->GetOpFlashPlaneSTD(FlashPlane, varXY, varYZ, varXZ);

        mf::LogDebug("SolarNuAna") << "Evaluating Flash purity";
        int TerminalOutput = ProducerUtils::supress_stdout();
        double ThisOpFlashPur = pbt->OpHitCollectionPurity(SignalTrackIDs, MatchedHits);
        ProducerUtils::resume_stdout(TerminalOutput);
        mf::LogDebug("SolarNuAna") << "PE of this OpFlash " << TotalFlashPE << " OpFlash time " << FlashTime;

        // Calculate the flash purity, only for the Signal events
        OpFlashID.push_back(i);
        OpFlashPlane.push_back(FlashPlane);
        OpFlashPur.push_back(ThisOpFlashPur);
        OpFlashMaxPE.push_back(MaxOpHitPE);
        OpFlashSTD.push_back(FlashStdDev);
        OpFlashTime.push_back(TheFlash.Time() * clockData.TPCClock().TickPeriod()); // Convert to microseconds
        OpFlashX.push_back(TheFlash.XCenter());
        OpFlashY.push_back(TheFlash.YCenter());
        OpFlashZ.push_back(TheFlash.ZCenter());
        OpFlashPE.push_back(TheFlash.TotalPE());
        OpFlashFast.push_back(TheFlash.FastToTotal());
        OpFlashDeltaT.push_back(TheFlash.TimeWidth());
        OpFlashNHits.push_back(MatchedHits.size());
        if (abs(TheFlash.Time()) < 3)
        {
          mf::LogDebug("SolarNuAna") << "OpFlash PE " << TheFlash.TotalPE() << " with purity " << ThisOpFlashPur << " time " << TheFlash.Time();
          sOpFlashTruth += "OpFlash PE " + ProducerUtils::str(TheFlash.TotalPE()) + " with purity " + ProducerUtils::str(ThisOpFlashPur) + " time " + ProducerUtils::str(TheFlash.Time()) + " plane " + ProducerUtils::str(FlashPlane) + "\n";
          sOpFlashTruth += " - Vertex (" + ProducerUtils::str(TheFlash.XCenter()) + ", " + ProducerUtils::str(TheFlash.YCenter()) + ", " + ProducerUtils::str(TheFlash.ZCenter()) + ")\n";
          sOpFlashTruth += "\t*** 1st Sanity check: Ratio " + ProducerUtils::str(MaxOpHitPE / TotalFlashPE) + " <= " + ProducerUtils::str(fAdjOpFlashMaxPERatioCut) + "\n";
          sOpFlashTruth += "\t*** 2nd Sanity check: #OpHits " + ProducerUtils::str(int(NMatchedHits)) + " >= " + ProducerUtils::str(int(TheFlash.PEs().size())) + "\n";
        }
      }
    }
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
      if (ThisHit.PeakTime() < 0)
        producer->PrintInColor("Negative Hit Time = " + ProducerUtils::str(ThisHit.PeakTime()), ProducerUtils::GetColor("red"));
      mf::LogDebug("SolarNuAna") << "Hit " << i << " has view " << ThisHit.View() << " and signal type " << ThisHit.SignalType();

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
    std::vector<std::vector<art::Ptr<recob::Hit>>> TestClusters;
    std::vector<std::vector<std::vector<recob::Hit>>> AllPlaneClusters;
    std::vector<std::vector<int>> ClustersIdx = {{}, {}, {}};
    std::vector<std::vector<int>> RecoHitIdx;
    // Map to associate the ClusterIdx with the position in the ClVectors
    std::map<int, std::vector<int>> ClIdxMap;
    
    if (fTestLatestFeatures)
    {
      lowe->CalcAdjHits(RecoHitsPtr, TestClusters, RecoHitIdx);
      for (int i = 0; i < int(TestClusters.size()); i++)
      {
        std::vector<recob::Hit> ThisHitVector = {}; // Convert pointer to vector
        for (int j = 0; j < int(TestClusters[i].size()); j++)
        {
          ThisHitVector.push_back(*TestClusters[i][j]);
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
    }
    else {
      // --- Now calculate the clusters ...
      lowe->CalcAdjHits(Ind0Hits, Clusters0, hAdjHits, hAdjHitsADCInt, false);
      HitNum.push_back(Ind0Hits.size());
      ClusterNum.push_back(Clusters0.size());
      lowe->CalcAdjHits(Ind1Hits, Clusters1, hAdjHits, hAdjHitsADCInt, false);
      HitNum.push_back(Ind1Hits.size());
      ClusterNum.push_back(Clusters1.size());
      lowe->CalcAdjHits(ColHits, Clusters2, hAdjHits, hAdjHitsADCInt, false);
      HitNum.push_back(ColHits.size());
      ClusterNum.push_back(Clusters2.size());
      lowe->CalcAdjHits(GhostHits, Clusters3, hAdjHits, hAdjHitsADCInt, false);
      HitNum.push_back(GhostHits.size());
      ClusterNum.push_back(Clusters3.size());
    }
    AllPlaneClusters = {Clusters0, Clusters1, Clusters2};
    fMCTruthTree->Fill();

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
          if (TPCHit.PeakTime() < 0)
            producer->PrintInColor("Negative Cluster Time = " + ProducerUtils::str(TPCHit.PeakTime()), ProducerUtils::GetColor("red"));
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
          clustT += TPCHit.Integral() * TPCHit.PeakTime();

          if (TPCHit.Integral() > maxHit)
          {
            // If clusterTPC not in TPCIDMap, set it to -1
            if (TPCIDMap.find(TPC) == TPCIDMap.end())
            {
              clustTPC = -1;
            }
            else
            {
              clustTPC = TPC;
            }
            // Look for maxHit inside cluster
            maxHit = TPCHit.Integral();            
          }
          else
          {
            clustTPC = TPC;
          }

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

          long unsigned int GeneratorType = ProducerUtils::WhichGeneratorType(GeneratorParticles, MainTrID);
          VecGenPur[int(GeneratorType)] = VecGenPur[int(GeneratorType)] + TPCHit.Integral();
          mf::LogDebug("SolarNuAna") << "\nThis particle type " << GeneratorType << "\nThis cluster's main track ID " << MainTrID;
          if (SignalTrackIDs.find(MainTrID) != SignalTrackIDs.end())
          {
            hitCharge = TPCHit.Integral();
            Pur = Pur + hitCharge;
          }
        }

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
        clustX /= ncharge;
        clustY /= ncharge;
        clustZ /= ncharge;
        clustT /= ncharge;
        Pur /= ncharge;
        mf::LogDebug("SolarNuAna") << "\ndzdy " << dzdy << " for cluster "
                                   << " (" << clustY << ", " << clustZ << ") with track ID " << MainTrID << " in plane " << idx;
        if (clustT < 0)
          producer->PrintInColor("Negative Cluster Time = " + ProducerUtils::str(clustT), ProducerUtils::GetColor("red"));
        
        if (fTestLatestFeatures)
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

        mf::LogDebug("SolarNuAna") << "\nCluster " << i << " in plane " << idx << " has #hits" << nhit << " charge, " << ncharge << " time, " << clustT;
        mf::LogDebug("SolarNuAna") << " and position (" << clustY << ", " << clustZ << ") with main track ID " << MainTrID << " and purity " << Pur;
      }
    } // Finished first cluster processing

    //-------------------------------------------------------------------- Cluster Matching -------------------------------------------------------------------------//
    std::vector<int> MVecMainID = {};
    std::vector<unsigned int> MVecGen = {};
    std::vector<std::vector<float>> MVecGenFrac = {};
    std::vector<float> MVecFracE = {}, MVecFracGa = {}, MVecFracNe = {}, MVecFracRest = {}, MVecGenPur = {};
    std::vector<std::vector<int>> MVecNHits = {{}, {}, {}}, MVecTPC = {{}, {}, {}}, MVecChannel = {{}, {}, {}};
    std::vector<std::vector<float>> MVecPur = {{}, {}, {}}, MVecMaxCharge = {{}, {}, {}}, MVecCharge = {{}, {}, {}}, MVecTime = {{}, {}, {}}, MVecRecoY = {{}, {}, {}}, MVecRecoZ = {{}, {}, {}};
    std::vector<std::vector<float>> MVecDirDir = {{}, {}, {}}, MatchedClCompleteness = {{}, {}, {}}, MVecDT = {{}, {}, {}};
    
    if (fTestLatestFeatures)
    {
      std::vector<std::vector<int>>  MatchedClustersIdx = {{}, {}, {}};
      std::vector<std::vector<std::vector<recob::Hit>>> MatchedClusters = {{}, {}, {}};
    
      std::string SolarClusterInfo = "SolarClusterInfo: ";
      SolarClusterInfo = SolarClusterInfo + "(" + ProducerUtils::str(Clusters0.size()) + "," + ProducerUtils::str(Clusters1.size()) + "," + ProducerUtils::str(Clusters2.size()) + ")";
      lowe->MatchClusters(SignalTrackIDs, MatchedClustersIdx, MatchedClusters, ClustersIdx, AllPlaneClusters, MVecNHits, MVecChannel, MVecTime, MVecRecoY, MVecRecoZ, MVecDirDir, MVecCharge, MVecPur, MatchedClCompleteness, clockData, true);
      
      int MatchedClusterNum = int(MatchedClustersIdx[2].size());
      SolarClusterInfo = SolarClusterInfo + "\nFound " + ProducerUtils::str(MatchedClusterNum) + " MatchedClusters (from col. plane loop)!";
      
      for (int ThisClIdx = 0; ThisClIdx < MatchedClusterNum; ThisClIdx++)
      {
        MVecTime[2][ThisClIdx] *= clockData.TPCClock().TickPeriod(); // Convert to microseconds
        for (int plane = 0; plane < 2; plane++) {
          int RefClIdx = ClIdxMap[MatchedClustersIdx[plane][ThisClIdx]][1]; // Get the cluster index in the plane
          if (MVecTime[plane][ThisClIdx] > -1e6)
          {
            MVecTime[plane][ThisClIdx] *= clockData.TPCClock().TickPeriod(); // Convert to microseconds
            MVecDT[plane].push_back(abs(MVecTime[2][ThisClIdx] - MVecTime[plane][ThisClIdx]));
            MVecTPC[plane].push_back(ClTPC[plane][RefClIdx]);
            MVecMaxCharge[plane].push_back(ClMaxCharge[plane][RefClIdx]);
            SolarClusterInfo = SolarClusterInfo + "\nMatched Cluster in plane " + ProducerUtils::str(plane) + " with time " + ProducerUtils::str(MVecTime[plane][ThisClIdx]) + " and charge " + ProducerUtils::str(MVecCharge[plane][ThisClIdx]) + " with TPC " + ProducerUtils::str(MVecTPC[plane][ThisClIdx]);
          }
          else
          {
            MVecDT[plane].push_back(-1e6);
            MVecTPC[plane].push_back(-1);
            MVecMaxCharge[plane].push_back(-1e6);
            SolarClusterInfo = SolarClusterInfo + "\nMatched Cluster in plane " + ProducerUtils::str(plane) + " with time -1e6 and charge -1e6 with TPC -1";
          }
        }
        int RefClIdx = ClIdxMap[MatchedClustersIdx[2][ThisClIdx]][1]; // Get the plane index of the matched cluster
        MVecTPC[2].push_back(ClTPC[2][RefClIdx]);
        MVecMaxCharge[2].push_back(ClMaxCharge[2][RefClIdx]);
        MVecGenPur.push_back(ClGenPur[2][RefClIdx]);
        MVecMainID.push_back(ClMainID[2][RefClIdx]);
        MVecGen.push_back(ClGen[2][RefClIdx]);
        MVecFracE.push_back(ClFracE[2][RefClIdx]);
        MVecFracGa.push_back(ClFracGa[2][RefClIdx]);
        MVecFracNe.push_back(ClFracNe[2][RefClIdx]);
        MVecFracRest.push_back(ClFracRest[2][RefClIdx]);  
        MVecGenFrac.push_back(ClVecGenPur[2][RefClIdx]);
        SolarClusterInfo = SolarClusterInfo + "\nMatched Cluster in plane 2 with time " + ProducerUtils::str(MVecTime[2][ThisClIdx]) + " and charge " + ProducerUtils::str(MVecCharge[2][ThisClIdx]) + " with TPC " + ProducerUtils::str(MVecTPC[2][ThisClIdx]);
      }
      producer->PrintInColor(SolarClusterInfo, ProducerUtils::GetColor("yellow"), "Debug");
    }

    else {
      for (int ii = 0; ii < int(AllPlaneClusters[2].size()); ii++)
      {
        bool match = false;
        // int ind0clustIndex = -1, ind1clustIndex = -1;
        int ind0clustNHits = 0, ind1clustNHits = 0;
        int ind0clustTPC = -1, ind1clustTPC = -1;
        double ind0clustY = -1e6, ind1clustY = -1e6, ind0clustMaxCharge = -1e6, ind1clustMaxCharge = -1e6, ind0clustCharge = -1e6, ind1clustCharge = -1e6;
        double ind0clustdT = fClusterMatchTime, ind1clustdT = fClusterMatchTime;
        if (!AllPlaneClusters[2][ii].empty())
        {
          if (!AllPlaneClusters[0].empty())
          {
            // std::cout << " - Matching cluster " << ii << " with time " << ClT[2][ii] << std::endl;
            for (int jj = 0; jj < int(AllPlaneClusters[0].size()); jj++)
            {
              if (abs(ClNHits[0][jj] - ClNHits[2][ii]) / ClNHits[2][ii] > fClusterMatchNHit || 
                  abs(ClCharge[0][jj] - ClCharge[2][ii]) / ClCharge[2][ii] > fClusterMatchCharge)
              {
                continue;
              }
              // std::cout << "    Checking Ind0 cluster with index " << jj << " and time " << ClT[0][jj] << std::endl;
              if (abs(ClT[2][ii] - ClT[0][jj]) < fClusterMatchTime && abs(fClusterInd0MatchTime - abs(ClT[2][ii] - ClT[0][jj])) < abs(fClusterInd0MatchTime - ind0clustdT))
              {
                ind0clustY = ClY[0][jj] + (ClZ[2][ii] - ClZ[0][jj]) / (Cldzdy[0][jj]);
                ind0clustdT = abs(ClT[2][ii] - ClT[0][jj]);
                ind0clustNHits = int(AllPlaneClusters[0][jj].size());
                ind0clustCharge = ClCharge[0][jj];
                ind0clustMaxCharge = ClMaxCharge[0][jj];
                ind0clustTPC = ClTPC[0][jj];
                // ind0clustIndex = jj;
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
              // std::cout << "    Checking Ind1 cluster with index " << zz << " and time " << ClT[1][zz] << std::endl;
              if (abs(ClT[2][ii] - ClT[1][zz]) < fClusterMatchTime && abs(fClusterInd1MatchTime - abs(ClT[2][ii] - ClT[1][zz])) < abs(fClusterInd1MatchTime - ind1clustdT))
              {
                ind1clustY = ClY[1][zz] + (ClZ[2][ii] - ClZ[1][zz]) / (Cldzdy[1][zz]);
                ind1clustdT = abs(ClT[2][ii] - ClT[1][zz]);
                ind1clustNHits = int(AllPlaneClusters[1][zz].size());
                ind1clustCharge = ClCharge[1][zz];
                ind1clustMaxCharge = ClMaxCharge[1][zz];
                ind1clustTPC = ClTPC[1][zz];
                // ind1clustIndex = zz;
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
          // std::cout << "    ***Matched cluster " << ii << " to Ind0: " << ind0clustIndex << ", Ind1: " << ind1clustIndex << std::endl;
          // Cluster Charge
          MVecCharge[0].push_back(ind0clustCharge);
          MVecCharge[1].push_back(ind1clustCharge);
          MVecCharge[2].push_back(ClCharge[2][ii]);
          MVecMaxCharge[0].push_back(ind0clustMaxCharge);
          MVecMaxCharge[1].push_back(ind1clustMaxCharge);
          MVecMaxCharge[2].push_back(ClMaxCharge[2][ii]);
          // Cluster Hits
          MVecNHits[0].push_back(ind0clustNHits);
          MVecNHits[1].push_back(ind1clustNHits);
          MVecNHits[2].push_back(ClNHits[2][ii]);
          // Cluster TPC
          MVecTPC[0].push_back(ind0clustTPC);
          MVecTPC[1].push_back(ind1clustTPC);
          MVecTPC[2].push_back(ClTPC[2][ii]);
          // Cluster Time
          MVecDT[0].push_back(ind0clustdT * clockData.TPCClock().TickPeriod());
          MVecDT[1].push_back(ind1clustdT * clockData.TPCClock().TickPeriod());
          MVecTime[0].push_back(ClT[0][ii] * clockData.TPCClock().TickPeriod()); // Convert to us
          MVecTime[1].push_back(ClT[1][ii] * clockData.TPCClock().TickPeriod()); // Convert to us
          MVecTime[2].push_back(ClT[2][ii] * clockData.TPCClock().TickPeriod()); // Convert to us
          // Cluster RecoY
          MVecRecoY[0].push_back(ind0clustY);
          MVecRecoY[1].push_back(ind1clustY);
          // Cluster RecoZ
          MVecRecoZ[0].push_back(ClZ[0][ii]);
          MVecRecoZ[1].push_back(ClZ[1][ii]);
          MVecRecoZ[2].push_back(ClZ[2][ii]);
          // Cluster Signal Fractions
          MVecFracE.push_back(ClFracE[2][ii]);
          MVecFracGa.push_back(ClFracGa[2][ii]);
          MVecFracNe.push_back(ClFracNe[2][ii]);
          MVecFracRest.push_back(ClFracRest[2][ii]);
          // Cluster Signal Purity
          MVecPur[0].push_back(ClPur[0][ii]);
          MVecPur[1].push_back(ClPur[1][ii]);
          MVecPur[2].push_back(ClPur[2][ii]);
          // Cluster Gen and GenFraction
          MVecMainID.push_back(ClMainID[2][ii]);
          MVecGen.push_back(ClGen[2][ii]);
          MVecGenPur.push_back(ClGenPur[2][ii]);
          MVecGenFrac.push_back(ClVecGenPur[2][ii]);

          float buffer = 1;
          if ((ind0clustY > -buffer * fDetectorSizeY && ind0clustY < buffer * fDetectorSizeY) && (ind1clustY > -buffer * fDetectorSizeY && ind1clustY < buffer * fDetectorSizeY))
          {
            mf::LogDebug("SolarNuAna") << "BOTH IND RECO INSIDE OF DETECTOR";
            MVecRecoY[2].push_back((ind0clustY + ind1clustY) / 2);
          }
          else if (ind0clustY > -buffer * fDetectorSizeY && ind0clustY < buffer * fDetectorSizeY)
          {
            mf::LogDebug("SolarNuAna") << "IND1 OUTSIDE OF DETECTOR";
            MVecRecoY[2].push_back(ind0clustY);
          }
          else if (ind1clustY > -buffer * fDetectorSizeY && ind1clustY < buffer * fDetectorSizeY)
          {
            mf::LogDebug("SolarNuAna") << "IND0 OUTSIDE OF DETECTOR";
            MVecRecoY[2].push_back(ind1clustY);
          }
          else
          {
            mf::LogDebug("SolarNuAna") << "RECO OUTSIDE OF DETECTOR";
            MVecRecoY[2].push_back((ind0clustY + ind1clustY) / 2);
            if (ClGen[2][ii] == 1)
            {
              mf::LogWarning("SolarNuAna") << "Signal cluster reconstructed outside of detector volume! RecoY = " << ProducerUtils::str((ind0clustY + ind1clustY) / 2);
            }
          }

          // Print in color if the cluster is matched
          mf::LogDebug("SolarNuAna") << "¡¡¡ Matched cluster !!! ";
          mf::LogDebug("SolarNuAna") << " - Cluster " << ProducerUtils::str(ClMainID[2][ii]) << " Gen " << ProducerUtils::str(ClGen[2][ii]) << " Purity " << ProducerUtils::str(ClGenPur[2][ii]) << " Hits " << ProducerUtils::str(ClNHits[2][ii]);
          mf::LogDebug("SolarNuAna") << " - Hits(ind0, ind1, col) " << ProducerUtils::str(ind0clustNHits) << ", " << ProducerUtils::str(ind1clustNHits) << ", " << ProducerUtils::str(ClNHits[2][ii]);
          mf::LogDebug("SolarNuAna") << " - Positions y(ind0, ind1) = " << ProducerUtils::str(ind0clustY) << ", " << ProducerUtils::str(ind1clustY) << ", z = " << ProducerUtils::str(ClZ[2][ii]) << "\n";
        } // if (match == true)
        else
        {
          // std::cout << "    No match found for index " << ii << ". Skipping..." << std::endl;
        }
      } // Loop over collection plane clusters
    } // Finished cluster matching
    //-------------------------------------------------------------------- Cluster Tree Export -------------------------------------------------------------------------//
    // Loop over matched clusters and export to tree if all conditions are satisfied
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
          if (j == i)
          {
            continue;
          }

          double ClusterDistance = 0;
          producer->ComputeDistance3D(ClusterDistance, MVecTime[2][i], MVecRecoY[2][i], MVecRecoZ[2][i], MVecTime[2][j], MVecRecoY[2][j], MVecRecoZ[2][j], TPCIDdriftLength[MVecTPC[2][i]], TPCIDdriftTime[MVecTPC[2][i]]);
          if (ClusterDistance > fAdjClusterRad)
          {
            continue;
          }
          sAdjClusters += "    - Cluster " + ProducerUtils::str(j) + " at distance " + ProducerUtils::str(ClusterDistance) + " with time " + ProducerUtils::str(MVecTime[2][j]) + " and charge " + ProducerUtils::str(MVecCharge[2][j]) + " in TPC " + ProducerUtils::str(MVecTPC[2][j]);
          if (MVecCharge[2][j] < fMinClusterCharge)
          {
            continue;
          }
          sAdjClusters += " and hits " + ProducerUtils::str(MVecNHits[2][j]) + "\n";
          if (MVecCharge[2][j] > MVecCharge[2][i])
          {
            MPrimary = false;
          }
          MAdjClNum += 1;
          if (MVecGen[i] == MVecGen[j])
          {
            MSignalAdjClNum += 1;
          }
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
          MAdjClMainID.push_back(MVecMainID[j]);

          // If mother exists add the mother information
          const simb::MCParticle *MAdjClTruth;
          int TerminalOutput = ProducerUtils::supress_stdout();
          MAdjClTruth = pi_serv->TrackIdToParticle_P(MVecMainID[j]);
          ProducerUtils::resume_stdout(TerminalOutput);
          if (MAdjClTruth == 0)
          {
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
          else
          {
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
        if (MVecPur[2][i] > 0)
        {
          sResultColor = "green";
        }
        if (fClusterPreselectionPrimary && !MPrimary)
        {
          continue;
        }
        if (MPrimary)
        {
          sClusterReco += "*** Matched preselection cluster: " + ProducerUtils::str(i) + "\n";
          sClusterReco += " - MainTrackID " + ProducerUtils::str(MVecMainID[i]) + "\n";
          if (MVecGen[i] > 0 && int(MVecGen[i]) < (int(fLabels.size()) + 1))
          {
            sClusterReco += " - Gen " + ProducerUtils::str(int(MVecGen[i])) + " -> " + fLabels[MVecGen[i] - 1];
          }
          else
          {
            sClusterReco += " - Gen ?? -> Unknown";
          }
          sClusterReco += " TPC " + ProducerUtils::str(MVecTPC[2][i]) + "\n";
          // sClusterReco += " - Truth X,Y,Z ( " + ProducerUtils::str(MVecMainX[i]) + ", " + ProducerUtils::str(MVecMainY[i]) + ", " + ProducerUtils::str(MVecMainZ[i]) + " )\n";
          sClusterReco += " - Purity " + ProducerUtils::str(MVecGenPur[i]) + " Hits " + ProducerUtils::str(MVecNHits[2][i]) + "\n";
          sClusterReco += " - #AdjCl " + ProducerUtils::str(MAdjClNum) + " ( " + ProducerUtils::str(MSignalAdjClNum) + " signal ):\n";
          if (AdjClusterMatch) {sClusterReco += sAdjClusters;}
          sClusterReco += " - RecoCol  Time,Y,Z ( " + ProducerUtils::str(MVecTime[2][i]) + ", " + ProducerUtils::str(MVecRecoY[2][i]) + ", " + ProducerUtils::str(MVecRecoZ[2][i]) + " )\n";
          sClusterReco += " - RecoInd0 Time,Y,Z ( " + ProducerUtils::str(MVecTime[0][i]) + ", " + ProducerUtils::str(MVecRecoY[0][i]) + ", " + ProducerUtils::str(MVecRecoZ[0][i]) + " )\n";
          sClusterReco += " - RecoInd1 Time,Y,Z ( " + ProducerUtils::str(MVecTime[1][i]) + ", " + ProducerUtils::str(MVecRecoY[1][i]) + ", " + ProducerUtils::str(MVecRecoZ[1][i]) + " )\n";
          
          if (fSaveTrackInfo){
            TVector3 ThisClVertex = {0, MVecRecoY[2][i], MVecRecoZ[2][i]};
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
              sClusterReco += " - Track has start ( " + ProducerUtils::str(trk.Start().X()) + ", " + ProducerUtils::str(trk.Start().Y()) + ", " + ProducerUtils::str(trk.Start().Z()) + " )\n";
              sClusterReco += " - Track has end   ( " + ProducerUtils::str(trk.End().X()) + ", " + ProducerUtils::str(trk.End().Y()) + ", " + ProducerUtils::str(trk.End().Z()) + " )\n\n";
              TrackMatch = true;
            }; // Loop over tracks
          }; // if (fSaveTrackInfo)
        }; // if (MPrimary)
        if (fClusterPreselectionTrack && !TrackMatch)
        {
          continue;
        }

        for (int j = 0; j < int(OpFlashPE.size()); j++)
        {
          float OpFlashR = -1e6;
          // Skip flashes with time outside the cluster time window
          if ((MVecTime[2][i] - OpFlashTime[j]) < 0 || (MVecTime[2][i] - OpFlashTime[j]) > TPCIDdriftTime[MVecTPC[2][i]])
          {
            continue;
          }
          
          // Compute the time distance between the cluster and the flash. Use factor 2 to convert us to TPC tics
          double MAdjFlashX = 0;
          
          // For HD 1x2x6 (only 1 APA) the x coordinate is determined by the TPC number
          if (fGeometry == "HD" && MVecTPC[2][i]%2 != 0)
          {
            producer->ComputeDistanceX(MAdjFlashX, MVecTime[2][i], OpFlashTime[j], TPCIDdriftLength[MVecTPC[2][i]], TPCIDdriftTime[MVecTPC[2][i]]);
          }
          else if (fGeometry == "HD" && MVecTPC[2][i]%2 == 0)
          {
            producer->ComputeDistanceX(MAdjFlashX, MVecTime[2][i], OpFlashTime[j], TPCIDdriftLength[MVecTPC[2][i]], TPCIDdriftTime[MVecTPC[2][i]]);
            MAdjFlashX = -MAdjFlashX;
          }
          else if (fGeometry == "VD")
          {
            producer->ComputeDistanceX(MAdjFlashX, MVecTime[2][i], OpFlashTime[j], TPCIDdriftLength[MVecTPC[2][i]], TPCIDdriftTime[MVecTPC[2][i]]);
            // Change the x distance to the x coordinate of the VD geometry [-driftLength/2, driftLength/2]
            MAdjFlashX = TPCIDdriftLength[MVecTPC[2][i]] / 2 - MAdjFlashX;
          }
          else
          {
            producer->PrintInColor("Unknown geometry " + fGeometry, ProducerUtils::GetColor("red"));
            producer->ComputeDistanceX(MAdjFlashX, MVecTime[2][i], OpFlashTime[j], TPCIDdriftLength[MVecTPC[2][i]], TPCIDdriftTime[MVecTPC[2][i]]);
          }

          // Make an eliptical cut on the flash position based on the clusters plane
          if (fGeometry == "VD" && OpFlashPlane[j] == 0) // Cathode flashes
          {
            if (pow(MVecRecoY[2][i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1)
            {
              continue;
            }
            OpFlashR = sqrt(pow(MVecRecoY[2][i] - OpFlashY[j], 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2));
          }
          else if (fGeometry == "VD" && (OpFlashPlane[j] == 1 || OpFlashPlane[j] == 2)) // Membrane flashes
          {
            if (fAdjOpFlashMembraneProjection){
              if (pow(MAdjFlashX - OpFlashX[j], 2) / pow(fAdjOpFlashX, 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1){
                continue;
              }
            }
            else{
              if (pow(MAdjFlashX - OpFlashX[j], 2) / pow(fAdjOpFlashX, 2) + pow(MVecRecoY[2][i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1){
                continue;
              }
            }
            OpFlashR = sqrt(pow(MAdjFlashX - OpFlashX[j], 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2));
          } 
          else if (fGeometry == "VD" && (OpFlashPlane[j] == 3 || OpFlashPlane[j] == 4)) // End-Cap flashes
          {
            if (fAdjOpFlashEndCapProjection){
              if (pow(MAdjFlashX - OpFlashX[j], 2) / pow(fAdjOpFlashX, 2) + pow(MVecRecoY[2][i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) > 1){
                continue;
              }
            }
            else{
              if (pow(MAdjFlashX - OpFlashX[j], 2) / pow(fAdjOpFlashX, 2) + pow(MVecRecoY[2][i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1){
                continue;
              }
            }
            OpFlashR = sqrt(pow(MAdjFlashX - OpFlashX[j], 2) + pow(MVecRecoY[2][i] - OpFlashY[j], 2));
          }
          else if (fGeometry == "VD" && OpFlashPlane[j] == -1)
          {
            continue;
          }

          if (fGeometry == "HD")
          {
            if (pow(MVecRecoY[2][i] - OpFlashY[j], 2) / pow(fAdjOpFlashY, 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2) / pow(fAdjOpFlashZ, 2) > 1)
            {
              continue;
            }
            OpFlashR = sqrt(pow(MVecRecoY[2][i] - OpFlashY[j], 2) + pow(MVecRecoZ[2][i] - OpFlashZ[j], 2));
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
          std::string sFlashMatching = "Testing flash " + ProducerUtils::str(j) + " with time " + ProducerUtils::str(OpFlashTime[j]) + " and PE " + ProducerUtils::str(OpFlashPE[j]);
          producer->PrintInColor(sFlashMatching, ProducerUtils::GetColor(sResultColor), "Debug");
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
        sVertexReco += "*** Reconstructed Interaction Vertex: \n";
        sVertexReco += " - True X,Y,Z ( " + ProducerUtils::str(SignalParticleX) + ", " + ProducerUtils::str(SignalParticleY) + ", " + ProducerUtils::str(SignalParticleZ) + " )" + "\n";
        sVertexReco += " - Reco X,Y,Z ( " + ProducerUtils::str(MatchedOpFlashX) + ", " + ProducerUtils::str(MVecRecoY[2][i]) + ", " + ProducerUtils::str(MVecRecoZ[2][i]) + " )";
        sClusterReco += sFlashReco;
        sClusterReco += sVertexReco;
        if (fClusterPreselectionFlashMatch && MatchedOpFlashPE < 0)
        {
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
        // Cluster TPC
        MInd0TPC = MVecTPC[0][i];
        MInd1TPC = MVecTPC[1][i];
        MTPC = MVecTPC[2][i];
        // Cluster Charge
        MInd0Charge = MVecCharge[0][i];
        MInd1Charge = MVecCharge[1][i];
        MCharge = MVecCharge[2][i];
        MInd0MaxCharge = MVecMaxCharge[0][i];
        MInd1MaxCharge = MVecMaxCharge[1][i];
        MMaxCharge = MVecMaxCharge[2][i];
        // Cluster Hits
        MInd0NHits = MVecNHits[0][i];
        MInd1NHits = MVecNHits[1][i];
        MNHit = MVecNHits[2][i];
        // Cluster Time
        MInd0dTime = MVecDT[0][i];
        MInd1dTime = MVecDT[1][i];
        MTime = MVecTime[2][i];
        // Cluster RecoX
        MRecX = MatchedOpFlashX;
        // Cluster RecoY
        MInd0RecoY = MVecRecoY[0][i];
        MInd1RecoY = MVecRecoY[1][i];
        MRecY = MVecRecoY[2][i];
        // Cluster RecoZ
        MRecZ = MVecRecoZ[2][i];
        // Cluster MainID
        MMainID = MVecMainID[i];

        // If mother exists add the mother information
        const simb::MCParticle *MClTruth;
        int TerminalOutput = ProducerUtils::supress_stdout();
        MClTruth = pi_serv->TrackIdToParticle_P(MVecMainID[i]);
        ProducerUtils::resume_stdout(TerminalOutput);
        if (MClTruth == 0)
        {
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
        else
        {
          if (MFlashPur > 0)
          {
            MFlashCorrect = true;
          };
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
          if (MClParentTruth == 0)
          {
            MMainParentVertex = {-1e6, -1e6, -1e6};
            MMainParentPDG = 0;
            MMainParentE = -1e6;
            MMainParentP = -1e6;
            MMainParentK = -1e6;
            MMainParentTime = -1e6;
          }
          else
          {
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
      if (sClusterReco != "")
      {
        producer->PrintInColor(sClusterReco, ProducerUtils::GetColor(sResultColor));
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
