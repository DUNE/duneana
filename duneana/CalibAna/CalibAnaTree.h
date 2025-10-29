#ifndef DUNEANA_CalibAnaTree
#define DUNEANA_CalibAnaTree

////////////////////////////////////////////////////////////////////////
// Class:       CalibAnaTree
// Plugin Type: analyzer
//
// Based on TrackCaloSkimmer.h by Gray Putnam,
// https://github.com/SBNSoftware/sbncode/blob/develop/sbncode/Calibration/TrackCaloSkimmer.h.
//
////////////////////////////////////////////////////////////////////////


#include <cstdlib>
#include <iostream>

#include <TTree.h>
#include "TFitter.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/ServicePack.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"

#include "larcorealg/GeoAlgo/GeoAlgo.h"
#include "larcorealg/Geometry/fwd.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "CalibAnaTreeObj.h"
#include "ICATSelectionTool.h"

namespace dune {
  class CalibAnaTree;
  enum EDet {kNOTDEFINED, kDUNED, kICARUS};
}

class dune::CalibAnaTree : public art::EDAnalyzer {
public:
  explicit CalibAnaTree(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CalibAnaTree(CalibAnaTree const&) = delete;
  CalibAnaTree(CalibAnaTree&&) = delete;
  CalibAnaTree& operator=(CalibAnaTree const&) = delete;
  CalibAnaTree& operator=(CalibAnaTree&&) = delete;

  ~CalibAnaTree();

  void analyze(art::Event const& e) override;

  void respondToOpenInputFile(const art::FileBlock& fb) override {
    (void) fb;
    fMeta.ifile ++;
  }

private:
  // Internal data struct
  struct GlobalTrackInfo {
    geo::Point_t start;
    geo::Point_t end;
    geo::Vector_t dir;
    geo::Vector_t enddir;
    int ID;
  };


  // Represents a "Snippet" of ADCs shared by a set of hits on a wire
  struct Snippet {
    geo::WireID wire;
    int start_tick;
    int end_tick;

    bool operator==(const Snippet &rhs) const {
      return wire == rhs.wire && start_tick == rhs.start_tick && end_tick == rhs.end_tick;
    }

    bool operator<(const Snippet &rhs) const {
      return wire < rhs.wire || (wire == rhs.wire && start_tick < rhs.start_tick);
    }
  };

  /*
  typedef struct //!< 2D point for clustering : - group (number of the associated cluster)
                 //!<                           - index (index to retreive the info like energy to the associated hit)
  {
    float y, z; //!<  [cm]
    int group, index;
  } point_t, *point;

  typedef struct //!< Temporary Clusters used to construct the Clusters (ClusterInfo)
  {
    float Sumy , Sumz , ECol , EInd1 , EInd2 , PeakTime;
    int Npoint, NCol , NInd1 , NInd2;
    std::list<int> lChannelCol , lChannelInd1 , lChannelInd2;
    std::vector<int> vMCPDG , vMCMOMpdg , vNOF;
    std::vector<float> vMCWEI;
    std::vector<float> vMCX , vMCY , vMCZ;
    std::vector<std::string> vMCGenTag;
    std::vector<float> vMCE;
    std::vector<int> vMCNe;
  } TempCluster ;

  //!< stucture used for association between a grid coordinate and a int
  //!< this structure is used for voxelisation and isolation search
  struct GridHasher
  {
      std::size_t operator()(const std::tuple<int, int, int>& key) const {
          // Use a combination of prime numbers to create a unique hash from 3D coordinates
          return std::get<0>(key) * 73856093 ^ std::get<1>(key) * 19349663 ^ std::get<2>(key) * 83492791;
      }
  };
*/

  // Fill vars
  void FillTrack(const recob::Track &track,
    const recob::PFParticle &pfp, float t0,
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<const recob::TrackHitMeta*> &thms,
    const std::vector<art::Ptr<recob::SpacePoint>> &sps,
    const std::vector<art::Ptr<anab::Calorimetry>> &calo,
    const std::map<geo::WireID, art::Ptr<raw::RawDigit>> &rawdigits,
    const std::vector<GlobalTrackInfo> &tracks,
    const geo::WireReadoutGeom *wireReadout,
    const detinfo::DetectorClocksData &clock_data,
    const cheat::BackTrackerService *bt_serv,
    const dune::EDet det);

  void FillTrackDaughterRays(const recob::Track &trk,
    const recob::PFParticle &pfp,
    const std::vector<art::Ptr<recob::PFParticle>> &PFParticleList,
    const art::FindManyP<recob::SpacePoint> &PFParticleSPs);

  void FillTrackEndHits(const geo::GeometryCore *geometry,
    const geo::WireReadoutGeom *wireReadout,
    const detinfo::DetectorPropertiesData &dprop,
    const recob::Track &track,
    const std::vector<art::Ptr<recob::Hit>> &allHits,
    const art::FindManyP<recob::SpacePoint> &allHitSPs);

  void FillTrackTruth(const detinfo::DetectorClocksData &clock_data,
    const std::vector<art::Ptr<recob::Hit>> &trkHits,
    const std::vector<art::Ptr<simb::MCParticle>> &mcparticles,
    const std::vector<geo::BoxBoundedGeo> &active_volumes,
    const std::vector<std::vector<geo::BoxBoundedGeo>> &tpc_volumes,
    const std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> id_to_ide_map,
    const std::map<int, std::vector<art::Ptr<recob::Hit>>> id_to_truehit_map,
    const detinfo::DetectorPropertiesData &dprop,
    const geo::GeometryCore *geo,
    const geo::WireReadoutGeom *wireReadout);

  TrackHitInfo MakeHit(const recob::Hit &hit,
    unsigned hkey,
    const recob::TrackHitMeta &thm,
    const recob::Track &trk,
    const art::Ptr<recob::SpacePoint> &sp,
    const std::vector<art::Ptr<anab::Calorimetry>> &calo,
    const geo::WireReadoutGeom *wireReadout,
    const detinfo::DetectorClocksData &dclock,
    const cheat::BackTrackerService *bt_serv);

  void DoTailFit();

  // declare truth utils, ported from CAFana in SBNCode
  std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>>
  PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels,
	    const geo::WireReadoutGeom &geo);
  std::map<int, std::vector<art::Ptr<recob::Hit>>>
  PrepTrueHits(const std::vector<art::Ptr<recob::Hit>> &allHits,
		 const detinfo::DetectorClocksData &clockData,
		 const cheat::BackTrackerService &backtracker);
  std::vector<std::pair<int, float>>
  AllTrueParticleIDEnergyMatches(const detinfo::DetectorClocksData &clockData,
				   const std::vector<art::Ptr<recob::Hit> >& hits,
				   bool rollup_unsaved_ids);
  float TotalHitEnergy(const detinfo::DetectorClocksData &clockData,
			 const std::vector<art::Ptr<recob::Hit> >& hits);
  int GetShowerPrimary(const int g4ID);




  // config

  // tags
 
  bool fLowEnergyClusterAnalysis;
  bool fTrackAnalysis;
  art::InputTag fRDTLabel;

  float fRadiusInt;
  float fRadiusExt;
  //float fElectronVelocity; //You have a local variable of the same name, and it's only used in analyze().  So commenting this out to make Clang happy.
  float fCoincidenceWd1_left;
  float fCoincidenceWd1_right;
  float fCoincidenceWd2_left;
  float fCoincidenceWd2_right;

  float fPitch;
  float fPitchMultiplier;

  bool  bIs3ViewsCoincidence;
  bool  bIsPD;
  bool  bIsHD;

  float fNumberInitClusters;
  float fMaxSizeCluster;
  float fMinSizeCluster;
  float fClusterSizeMulti;
  float fNumberConvStep;
  float fCovering;

  art::InputTag fPFPproducer;
  std::vector<art::InputTag> fT0producers;
  art::InputTag fCALOproducer;
  art::InputTag fTRKproducer;
  art::InputTag fTRKHMproducer;
  art::InputTag fHITproducer;
  std::vector<art::InputTag> fRawDigitproducers;
  std::string fG4producer;
  std::string fSimChannelproducer;
  bool fRequireT0;
  bool fDoTailFit;
  bool fVerbose;
  bool fSilenceMissingDataProducts;
  double fHitRawDigitsTickCollectWidth;
  double fTailFitResidualRange;
  int fHitRawDigitsWireCollectWidth;
  bool fFillTrackEndHits;
  float fTrackEndHitWireBox;
  float fTrackEndHitTimeBox;
  std::string fDetector;

  // tools
  std::vector<std::unique_ptr<dune::ICATSelectionTool>> fSelectionTools;

  // persistent info
  MetaInfo fMeta;

  std::map<Snippet, int> fSnippetCount;
  std::map<geo::WireID, std::pair<int, int>> fWiresToSave;

  const geo::WireReadoutGeom& fWireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  const geo::Geometry* fGeom;
  float fgeoXmin = 1e6;
  float fgeoXmax =-1e6;
  float fgeoYmin = 1e6;
  float fgeoYmax =-1e6;
  float fgeoZmin = 1e6;
  float fgeoZmax =-1e6;

  // Output
  TTree *fTree_track;
  TTree *fTree_LE;
  TrackInfo *fTrack;
  ClusterInfo *fCluster;

  // Fitting info
  TFitter fFitExp;
  TFitter fFitConst;
};
#endif
