////////////////////////////////////////////////////////////////////////
//
// \file CAFMaker_module.cc
//
// Chris Marshall's version
// Largely based on historical FDSensOpt/CAFMaker_module.cc
// Overhauled by Pierre Granger to adapt it to the new CAF format
//
///////////////////////////////////////////////////////////////////////

#ifndef CAFMaker_H
#define CAFMaker_H

// Generic C++ includes
#include <iostream>
#include <array>
#include <cmath>
#include <deque>
#include <limits>
#include <tuple>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"

#include "duneanaobj/StandardRecord/Flat/FlatRecord.h"

//#include "Utils/AppInit.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dunereco/CVN/func/InteractionType.h"
#include "dunereco/CVN/func/Result.h"
#include "dunereco/RegCNN/func/RegCNNResult.h"
#include "dunereco/FDSensOpt/FDSensOptData/AngularRecoOutput.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larcore/Geometry/Geometry.h"
#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/AnaUtils/DUNEAnaSliceUtils.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"


// root
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

// pdg
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"

// genie
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/GHEP/GHepParticle.h"


namespace caf {

  class CAFMaker : public art::EDAnalyzer {

    public:

      explicit CAFMaker(fhicl::ParameterSet const& pset);
      virtual ~CAFMaker();
      void beginJob() override;
      void endJob() override;
      void beginSubRun(const art::SubRun& sr) override;
      void endSubRun(const art::SubRun& sr) override;
      void analyze(art::Event const & evt) override;


    private:
      /// Truth category of the ionisation contributing to a hit.
      /// Mutually exclusive and exhaustive, mirroring the members of caf::SRHitSummary.
      enum class HitOrigin { kBeam = 0, kBeamDaughter, kBeamContam, kCosmic, kOther, kUnmatched, kN };

      /// Which truth category each GEANT4 track ID belongs to, for one event.
      struct TruthClassification {
        std::map<int, HitOrigin> g4ToOrigin;  ///< abs(G4 TrackId) -> category
        int mainBeamTid = -1;                 ///< G4 TrackId of the main (trigger) beam particle
      };

      /// Per-hit truth breakdown.
      struct HitTruth {
        bool computed = false;
        HitOrigin dominant = HitOrigin::kUnmatched;  ///< category of the largest contributor
        std::array<float, static_cast<std::size_t>(HitOrigin::kN)> frac{};  ///< energyFrac per category; sums to <= 1
      };

      /// Memoisation of HitTruth, keyed by (art product, hit index within it).
      /// Back-tracking a hit is expensive and every hit is visited at least twice (once
      /// event-wide, once within its slice), so the cache must be shared across both passes.
      /// Keying on the ProductID as well as the index means hits from a different product
      /// simply get their own table rather than colliding.
      class HitTruthCache {
        public:
          const HitTruth& Get(const art::Ptr<recob::Hit> &hit,
                              detinfo::DetectorClocksData const& clockData,
                              cheat::BackTrackerService const& bt,
                              const TruthClassification &cls);
        private:
          std::map<art::ProductID, std::vector<HitTruth>> fByProduct;
      };

      /// The SRHitCategory member of `summary` corresponding to `origin`.
      static caf::SRHitCategory& CategoryRef(caf::SRHitSummary &summary, HitOrigin origin);

      void PreLoadMCParticlesInfo(art::Event const& evt);
      void FillTruthInfo(caf::SRTruthBranch& sr,
                         art::Event const& evt);
      void FillMetaInfo(caf::SRDetectorMeta &meta, art::Event const& evt) const;
      void FillBeamInfo(caf::SRBeamBranch &beam, const art::Event &evt) const;

      /// @name Ported from protoana::ProtoDUNEBeamlineUtils
      ///
      /// Reimplemented here so the module keeps no protoduneana dependency. They read only
      /// `beam::ProtoDUNEBeamEvent`, which comes from dunecore::DuneObj. Kept faithful to
      /// protoduneana v10_20_09d00; if the upstream cuts are retuned these must be revisited.
      ///@{
      /// Beamline trigger is the beam trigger, and is matched to the DAQ trigger.
      static bool IsGoodBeamlineTrigger(const beam::ProtoDUNEBeamEvent &beamEvent);
      /// Exactly one non-glitching fiber in each momentum monitor, so the momentum is unambiguous.
      static bool HasPerfectBeamMomentum(const beam::ProtoDUNEBeamEvent &beamEvent);
      /// CERN-calibration TOF + Cherenkov PID selection; returns PDG codes.
      std::vector<int> GetBeamPDGCandidates(const beam::ProtoDUNEBeamEvent &beamEvent,
                                            double nominal_momentum) const;
      ///@}
      void FillRecoInfo(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, caf::SRTruthBranch &truthBranch, const art::Event &evt, const art::Ptr<recob::Slice> &slicePtr, const art::FindManyP<recob::PFParticle> &sliceToPFP, const std::vector<caf::SRHitSummary> &sliceSummaries, const caf::SRBeamInstrumentation &beamInst) const;
      void FillRecoInfoSliceLoop(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, caf::SRTruthBranch &truthBranch, const art::Event &evt, const caf::SRBeamInstrumentation &beamInst) const;
      void FillCVNInfo(caf::SRCVNScoreBranch &cvnBranch, const art::Event &evt) const;
      void FillEnergyInfo(caf::SRNeutrinoEnergyBranch &ErecBranch, const art::Event &evt) const;
      void FillRecoParticlesInfo(caf::SRRecoParticlesBranch &recoParticlesBranch, caf::SRFD &fdBranch, const art::Event &evt, const art::Ptr<recob::Slice> &slicePtr, const art::FindManyP<recob::PFParticle> &sliceToPFP) const;
      void FillDirectionInfo(caf::SRDirectionBranch &dirBranch, const art::Event &evt) const;
      int FillGENIERecord(simb::MCTruth const& mctruth, simb::GTruth const& gtruth);
      double GetVisibleEnergy(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const;
      void FillTruthMatchingAndOverlap(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt, std::vector<TrueParticleID> &truth, std::vector<float> &truthOverlap) const;
      void FillPFPMetadata(caf::SRPFP &pfpBranch, art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const;
      void GetMVAResults(caf::SRPFP & output_pfp, art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt, anab::MVAReader<recob::Hit,4> * hitResults, int planeid, bool charge_weighted) const;
      double GetWallDistance(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const;
      double GetWallDistance(recob::SpacePoint const& sp) const;
      void ComputeActiveBounds();
      double GetSingleHitsEnergy(art::Event const& evt, const art::Ptr<recob::Slice> &slicePtr, int plane) const;
      std::map<int, std::vector<const sim::IDE*>> slice_IDEs(
          std::vector<const sim::IDE*> ides,
          double the_z0, double the_pitch, double true_endZ) const;
      double ComputeTrueInteractingEnergy(
          const art::Event& evt,
          const simb::MCParticle* true_beam_particle) const;
      void ComputeRecoInteractingEnergy(
          const art::Event& evt,
          const recob::Track* thisTrack,
          detinfo::DetectorClocksData const& clockData,
          caf::SRInteractionBranch& interaction,
          caf::SRTruthBranch& truthBranch,
          const caf::SRBeamInstrumentation& beamInst) const;
      bool IsVertexContained(caf::SRVector3D const& vtx) const;

      /// Collection-plane index. Charge-based energies are only meaningful on this plane,
      /// and summing ionisation across all three planes would triple-count it.
      int CollectionPlane() const { return fVPlaneAsCollector ? 1 : 2; }
      /// Lifetime-corrected ADC area -> energy [GeV].
      double ChargeToEnergyGeV(double charge, int plane) const;
      /// Map every GEANT4 track ID in the event onto a truth category, by generator module label.
      TruthClassification BuildTruthClassification(const art::Event &evt) const;
      /// Accumulate the hit composition of `hits` (one slice, or the whole event) into `out`.
      void AccumulateHitSummary(const std::vector<art::Ptr<recob::Hit>> &hits,
                                detinfo::DetectorClocksData const& clockData,
                                detinfo::DetectorPropertiesData const& detProp,
                                const TruthClassification &cls,
                                HitTruthCache &cache,
                                caf::SRHitSummary &out) const;
      /// Fill the event-wide hit summary, the per-slice summaries and the slice bookkeeping.
      void ComputeHitSummaries(const art::Event &evt,
                               const std::vector<art::Ptr<recob::Slice>> &slicePtrs,
                               caf::SRInteractionBranch &ixn,
                               std::vector<caf::SRHitSummary> &sliceSummaries) const;
      /// True ionisation energy of the beam particle, its daughters and the beam halo,
      /// straight from the SimChannels (i.e. independent of hit finding).
      void FillSimChannelBeamEnergy(const TruthClassification &cls,
                                    caf::SRInteractionBranch &ixn) const;

      std::string fCVNLabel;
      bool fIsAtmoCVN;
      std::string fRegCNNLabel;

      // std::string fMCTruthLabel;
      std::vector<std::string> fMCTruthLabel; //To be able to read both the neutrino and cosmic MC truth
      std::string fGTruthLabel;
      std::string fMCFluxLabel;
      std::string fPOTSummaryLabel;
      std::string fEnergyRecoCaloLabel;
      std::string fEnergyRecoLepCaloLabel;
      std::string fEnergyRecoMuRangeLabel;
      std::string fEnergyRecoMuMcsLabel;
      std::string fEnergyRecoMuMcsLLHDLabel;
      std::string fEnergyRecoECaloLabel;
      std::string fDirectionRecoLabelNue;
      std::string fDirectionRecoLabelNumu;
      std::string fDirectionRecoLabelCalo;
      std::string fPandoraLabel;
      std::string fParticleIDLabel;
      std::string fTrackLabel;
      std::string fShowerLabel;
      std::string fSpacePointLabel;
      double fContainedDistThreshold;
      std::string fHitLabel;
      std::string fG4Label;
      std::string fMVALabel;


      std::map<int, std::tuple<art::Ptr<simb::MCParticle>, int, bool, int>> fMCParticlesMap; //[tid] = (MCParticle, interaction, isPrimary, SRParticle ID)

      TTree* fTree = nullptr;
      TTree* fMetaTree = nullptr;
      TTree* fGENIETree = nullptr;

      std::unique_ptr<TFile> fFlatFile;
      TTree* fFlatTree = nullptr; //Ownership will be managed directly by ROOT
      std::unique_ptr<flat::Flat<caf::StandardRecord>> fFlatRecord;

      genie::NtpMCEventRecord *fEventRecord = nullptr;

      double fMetaPOT;
      int fMetaRun, fMetaSubRun, fMetaVersion;

      const geo::Geometry* fGeom;
      std::vector<double> fActiveBounds;
      std::vector<double> fVertexFiducialVolumeCut;

      calo::CalorimetryAlg fCalorimetryAlg;                    ///< the calorimetry algorithm
      double fRecombFactor; ///< recombination factor for the isolated hits
      std::string fCalorimetryLabelSCE;    ///< SCE-corrected calorimetry label for beam track
      bool fVPlaneAsCollector;             ///< Use V-plane (1) instead of W-plane (2) as collection plane
      std::string fPFParticleLabel;        ///< PFParticle label used to retrieve beam particle track
      bool fMCHasBI;                       ///< MC sample has beam instrumentation info (enables data-path KE in ComputeRecoInteractingEnergy)
      std::string fBeamModuleLabel;        ///< Label for beam instrumentation data product (real data)
      double fBeamInstPFix;               ///< Momentum correction factor for MC beam instrumentation
      double fBeamPIDMomentum;            ///< Nominal beam momentum [GeV/c] the PID cuts are defined at
      bool   fUseCERNCalibSelection;      ///< Use the CERN-calibrated TOF cuts (reco after ~v08_07_00) rather than the older values
      geo::WireReadoutGeom const* fWireReadout = nullptr; ///< Wire readout geometry for wire pitch

      std::vector<std::string> fBeamTruthLabels;   ///< Generator module labels producing the beam particle(s)
      std::vector<std::string> fCosmicTruthLabels; ///< Generator module labels producing cosmic rays
      bool fFillHitSummaries;                      ///< Fill the per-slice / event-wide hit truth composition

      const std::map<simb::Generator_t, caf::Generator> fgenMap = {
        {simb::Generator_t::kUnknown, caf::Generator::kUnknownGenerator},
        {simb::Generator_t::kGENIE,   caf::Generator::kGENIE},
        {simb::Generator_t::kCRY,     caf::Generator::kCRY},
        {simb::Generator_t::kGIBUU,   caf::Generator::kGIBUU},
        {simb::Generator_t::kNuWro,   caf::Generator::kNuWro},
        {simb::Generator_t::kMARLEY,  caf::Generator::kMARLEY},
        {simb::Generator_t::kNEUT,    caf::Generator::kNEUT},
        {simb::Generator_t::kCORSIKA, caf::Generator::kCORSIKA},
        {simb::Generator_t::kGEANT,   caf::Generator::kGEANT}
      };


  }; // class CAFMaker


  //------------------------------------------------------------------------------
  CAFMaker::CAFMaker(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset),
      fCVNLabel(pset.get<std::string>("CVNLabel")),
      fIsAtmoCVN(pset.get<bool>("IsAtmoCVN")),
      fRegCNNLabel(pset.get<std::string>("RegCNNLabel")),
      // fMCTruthLabel(pset.get<std::string>("MCTruthLabel")),
      fMCTruthLabel(pset.get<std::vector<std::string>>("MCTruthLabel")),
      fGTruthLabel(pset.get<std::string>("GTruthLabel")),
      fMCFluxLabel(pset.get<std::string>("MCFluxLabel")),
      fPOTSummaryLabel(pset.get<std::string>("POTSummaryLabel")),
      fEnergyRecoCaloLabel(pset.get<std::string>("EnergyRecoCaloLabel")),
      fEnergyRecoLepCaloLabel(pset.get<std::string>("EnergyRecoLepCaloLabel")),
      fEnergyRecoMuRangeLabel(pset.get<std::string>("EnergyRecoMuRangeLabel")),
      fEnergyRecoMuMcsLabel(pset.get<std::string>("EnergyRecoMuMcsLabel")),
      fEnergyRecoMuMcsLLHDLabel(pset.get<std::string>("EnergyRecoMuMcsLLHDLabel")),
      fEnergyRecoECaloLabel(pset.get<std::string>("EnergyRecoECaloLabel")),
      fDirectionRecoLabelNue(pset.get<std::string>("DirectionRecoLabelNue")),
      fDirectionRecoLabelNumu(pset.get<std::string>("DirectionRecoLabelNumu")),
      fDirectionRecoLabelCalo(pset.get<std::string>("DirectionRecoLabelCalo")),
      fPandoraLabel(pset.get< std::string >("PandoraLabel")),
      fParticleIDLabel(pset.get< std::string >("ParticleIDLabel")),
      fTrackLabel(pset.get< std::string >("TrackLabel")),
      fShowerLabel(pset.get< std::string >("ShowerLabel")),
      fSpacePointLabel(pset.get< std::string >("SpacePointLabel")),
      fContainedDistThreshold(pset.get< double >("ContainedDistThreshold")),
      fHitLabel(pset.get< std::string >("HitLabel")),
      fG4Label(pset.get< std::string >("G4Label")),
      fMVALabel(pset.get<std::string>("MVALabel")),
      fEventRecord(new genie::NtpMCEventRecord),
      fGeom(&*art::ServiceHandle<geo::Geometry>()),
      fVertexFiducialVolumeCut(pset.get<std::vector<double>>("VertexFiducialVolumeCut")),
      fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
      fRecombFactor(pset.get<double>("RecombFactor")),
      fCalorimetryLabelSCE(pset.get<std::string>("CalorimetryLabelSCE", "")),
      fVPlaneAsCollector(pset.get<bool>("VPlaneAsCollector", false)),
      fPFParticleLabel(pset.get<std::string>("PFParticleLabel", "pandora")),
      fMCHasBI(pset.get<bool>("MCHasBI", false)),
      fBeamModuleLabel(pset.get<std::string>("BeamModuleLabel", "beamevent")),
      fBeamInstPFix(pset.get<double>("BeamInstPFix", 1.)),
      // Defaults match protoduneana's ProtoDUNEBeamlineUtils.fcl, whose selection these reproduce.
      fBeamPIDMomentum(pset.get<double>("BeamPIDMomentum", 1.)),
      fUseCERNCalibSelection(pset.get<bool>("UseCERNCalibSelection", true)),
      fBeamTruthLabels(pset.get<std::vector<std::string>>("BeamTruthLabels", {"generator"})),
      fCosmicTruthLabels(pset.get<std::vector<std::string>>("CosmicTruthLabels", {"cosmicgenerator"})),
      fFillHitSummaries(pset.get<bool>("FillHitSummaries", true))
  {

    if(pset.get<bool>("CreateFlatCAF")){
      // LZ4 is the fastest format to decompress. I get 3x faster loading with
      // this compared to the default, and the files are only slightly larger.
      fFlatFile = std::make_unique<TFile>("flatcaf.root", "RECREATE", "",
                            ROOT::CompressionSettings(ROOT::kLZ4, 1));
    }

    ComputeActiveBounds();
    fWireReadout = &art::ServiceHandle<geo::WireReadout>()->Get();

    if(fVertexFiducialVolumeCut.size() != 6){
      throw cet::exception("CAFMaker") << "VertexFiducialVolumeCut must be a vector of 6 elements";
    }
  }

  //------------------------------------------------------------------------------
  caf::CAFMaker::~CAFMaker()
  {
  }

  //------------------------------------------------------------------------------
  void CAFMaker::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("cafTree", "cafTree");

    // Create the branch. We will update the address before we write the tree
    caf::StandardRecord* rec = 0;
    fTree->Branch("rec", "caf::StandardRecord", &rec);

    fMetaTree = tfs->make<TTree>("meta", "meta");

    fMetaTree->Branch("pot", &fMetaPOT, "pot/D");
    fMetaTree->Branch("run", &fMetaRun, "run/I");
    fMetaTree->Branch("subrun", &fMetaSubRun, "subrun/I");
    fMetaTree->Branch("version", &fMetaVersion, "version/I");

    fMetaPOT = 0.;
    fMetaVersion = 1;

    fGENIETree = tfs->make<TTree>("genieEvt", "genieEvt");

    fGENIETree->Branch("genie_record", "genie::NtpMCEventRecord", &fEventRecord);

    if(fFlatFile){
      fFlatFile->cd();
      fFlatTree = new TTree("cafTree", "cafTree");

      fFlatRecord = std::make_unique<flat::Flat<caf::StandardRecord>>(fFlatTree, "rec", "", nullptr);
    }

  }

  //------------------------------------------------------------------------------

  void CAFMaker::PreLoadMCParticlesInfo(art::Event const& evt)
  {
    //Preloading the MCParticles info to be able to access them later
    std::vector<art::Ptr< simb::MCParticle>> mcparticles = dune_ana::DUNEAnaEventUtils::GetMCParticles(evt, fG4Label);

    //Creating a map of the MC particles to easily access them by their TrackId
    fMCParticlesMap.clear();
    //Not the most efficient way, but I prefer good readibility over performance here
    for(art::Ptr<simb::MCParticle> const& mcpart: mcparticles) {
      if(mcpart->TrackId() == 0) continue; //Skip the neutrino particle (TrackId 0)
      fMCParticlesMap[mcpart->TrackId()] = {mcpart, -1, false, -1}; // The only useful information here is the pointer to the MCParticle, the other values will be filled in FillTruthInfo
    }
  }

  //------------------------------------------------------------------------------

  caf::SRHitCategory& CAFMaker::CategoryRef(caf::SRHitSummary &summary, HitOrigin origin)
  {
    switch(origin){
      case HitOrigin::kBeam:         return summary.beam;
      case HitOrigin::kBeamDaughter: return summary.beam_daughters;
      case HitOrigin::kBeamContam:   return summary.beam_contamination;
      case HitOrigin::kCosmic:       return summary.cosmic;
      case HitOrigin::kOther:        return summary.other;
      default:                       return summary.unmatched;
    }
  }

  //------------------------------------------------------------------------------

  double CAFMaker::ChargeToEnergyGeV(double charge, int plane) const
  {
    return fCalorimetryAlg.ElectronsFromADCArea(charge, plane)*1./fRecombFactor/util::kGeVToElectrons;
  }

  //------------------------------------------------------------------------------

  CAFMaker::TruthClassification CAFMaker::BuildTruthClassification(const art::Event &evt) const
  {
    TruthClassification cls;

    // Which generator produced each MCTruth tells us what its particles are. Using the module
    // label rather than simb::Origin_t matters here: the ProtoDUNE beam gun and the radiological
    // generators both report kSingleParticle, so Origin_t alone cannot separate beam from
    // radioactivity. Same approach as GetGeneratorTag() in duneana/CalibAna/LEClustersFunctions.h.
    std::vector<art::Handle<std::vector<simb::MCTruth>>> mctruthHandles = evt.getMany<std::vector<simb::MCTruth>>();

    // Track IDs from the beam generator(s), kept aside so the beam family can be split into
    // the main particle, its descendants, and unrelated beam-halo activity.
    std::set<int> beamTids;

    for(auto const& mctruthHandle : mctruthHandles){
      if(!mctruthHandle.isValid()) continue;
      const std::string label = mctruthHandle.provenance()->moduleLabel();

      HitOrigin category = HitOrigin::kOther;
      const bool isBeamLabel = std::find(fBeamTruthLabels.begin(), fBeamTruthLabels.end(), label) != fBeamTruthLabels.end();
      if(isBeamLabel){
        category = HitOrigin::kBeamContam; //Refined below once the main beam particle is known
      }
      else if(std::find(fCosmicTruthLabels.begin(), fCosmicTruthLabels.end(), label) != fCosmicTruthLabels.end()){
        category = HitOrigin::kCosmic;
      }
      else if(!mctruthHandle->empty()){
        // Not a configured label: fall back on the generator's own claim about its origin.
        const simb::Origin_t origin = mctruthHandle->front().Origin();
        if(origin == simb::kCosmicRay) category = HitOrigin::kCosmic;
      }

      // largeant associates *every* particle it tracked to its parent MCTruth, secondaries
      // included, so this classifies each full ancestry tree without walking mothers.
      art::FindManyP<simb::MCParticle> truthToParticles(mctruthHandle, evt, fG4Label);
      if(!truthToParticles.isValid()){
        mf::LogWarning("CAFMaker") << "No MCTruth->MCParticle associations for generator '" << label
                                   << "' under G4 label '" << fG4Label << "'. Its hits will be counted as unmatched.";
        continue;
      }

      for(size_t i = 0; i < mctruthHandle->size(); ++i){
        for(art::Ptr<simb::MCParticle> const& part : truthToParticles.at(i)){
          const int tid = std::abs(part->TrackId());
          cls.g4ToOrigin[tid] = category;
          if(isBeamLabel) beamTids.insert(tid);
        }
      }
    }

    if(beamTids.empty()) return cls;

    // The ProtoDUNE beam generator injects the whole spill: the particle that fired the trigger
    // *plus* the beam halo, which can be far more energetic (tens of GeV muons that never enter
    // the TPC). The trigger particle is by convention the first generator primary the beam
    // MCTruth wrote, which is what protoana::ProtoDUNETruthUtils::GetGeantGoodParticle picks by
    // scanning GetParticle(t) in order. So order by (MCTruth index, generated particle index)
    // and take the lowest -- NOT the most energetic.
    std::tuple<size_t, size_t, int> bestKey{std::numeric_limits<size_t>::max(),
                                            std::numeric_limits<size_t>::max(),
                                            std::numeric_limits<int>::max()};
    for(auto const& mctruthHandle : mctruthHandles){
      if(!mctruthHandle.isValid()) continue;
      const std::string label = mctruthHandle.provenance()->moduleLabel();
      if(std::find(fBeamTruthLabels.begin(), fBeamTruthLabels.end(), label) == fBeamTruthLabels.end()) continue;

      art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> truthToParticles(mctruthHandle, evt, fG4Label);
      const bool haveGenIndex = truthToParticles.isValid();

      // Without the GeneratedParticleInfo metadata we cannot tell generator primaries from G4
      // secondaries via the association, so fall back on the particle's own provenance and on
      // the G4 track ID ordering (G4 numbers primaries in generator order).
      std::unique_ptr<art::FindManyP<simb::MCParticle>> plainToParticles;
      if(!haveGenIndex) plainToParticles = std::make_unique<art::FindManyP<simb::MCParticle>>(mctruthHandle, evt, fG4Label);

      for(size_t i = 0; i < mctruthHandle->size(); ++i){
        std::vector<art::Ptr<simb::MCParticle>> parts;
        std::vector<const sim::GeneratedParticleInfo*> infos;
        if(haveGenIndex){
          parts = truthToParticles.at(i);
          infos = truthToParticles.data(i);
        }
        else{
          if(!plainToParticles->isValid()) continue;
          parts = plainToParticles->at(i);
        }

        for(size_t j = 0; j < parts.size(); ++j){
          const bool isGenPrimary = haveGenIndex ? infos[j]->hasGeneratedParticleIndex()
                                                 : (parts[j]->Process() == "primary");
          if(!isGenPrimary || parts[j]->Mother() != 0) continue;

          const int tid = std::abs(parts[j]->TrackId());
          const size_t genIndex = haveGenIndex ? static_cast<size_t>(infos[j]->generatedParticleIndex())
                                               : static_cast<size_t>(tid);
          const std::tuple<size_t, size_t, int> key{i, genIndex, tid};
          if(key < bestKey){
            bestKey = key;
            cls.mainBeamTid = tid;
          }
        }
      }
    }

    if(cls.mainBeamTid < 0){
      mf::LogWarning("CAFMaker") << "Found " << beamTids.size() << " beam-generator track IDs but no generator-level "
                                 << "primary among them. All beam activity will be recorded as beam contamination.";
      return cls;
    }

    cls.g4ToOrigin[cls.mainBeamTid] = HitOrigin::kBeam;

    // Everything descending from the main beam particle is beam daughter activity; anything else
    // the beam generator produced (the halo) stays as beam contamination. Walking down from the
    // main particle visits each descendant once, unlike walking every particle's mother chain up.
    std::deque<int> toVisit{cls.mainBeamTid};
    while(!toVisit.empty()){
      const int tid = toVisit.front();
      toVisit.pop_front();

      auto it = fMCParticlesMap.find(tid);
      if(it == fMCParticlesMap.end()) continue;
      art::Ptr<simb::MCParticle> const& part = std::get<0>(it->second);

      for(int d = 0; d < part->NumberDaughters(); ++d){
        const int daughter = std::abs(part->Daughter(d));
        auto dIt = cls.g4ToOrigin.find(daughter);
        if(dIt == cls.g4ToOrigin.end()) continue;                    //Not tracked by the beam generator
        if(dIt->second == HitOrigin::kBeamDaughter) continue;        //Already visited
        dIt->second = HitOrigin::kBeamDaughter;
        toVisit.push_back(daughter);
      }
    }

    return cls;
  }

  //------------------------------------------------------------------------------

  const CAFMaker::HitTruth& CAFMaker::HitTruthCache::Get(const art::Ptr<recob::Hit> &hit,
                                                         detinfo::DetectorClocksData const& clockData,
                                                         cheat::BackTrackerService const& bt,
                                                         const TruthClassification &cls)
  {
    std::vector<HitTruth> &table = fByProduct[hit.id()];
    if(table.size() <= hit.key()) table.resize(hit.key() + 1);

    HitTruth &out = table[hit.key()];
    if(out.computed) return out;
    out.computed = true;

    float sumFrac = 0;
    for(sim::TrackIDE const& ide : bt.HitToTrackIDEs(clockData, hit)){
      // G4 rolls EM daughters up into their parent and flags them with a negative track ID.
      auto it = cls.g4ToOrigin.find(std::abs(ide.trackID));
      const HitOrigin origin = (it != cls.g4ToOrigin.end()) ? it->second : HitOrigin::kUnmatched;
      out.frac[static_cast<std::size_t>(origin)] += ide.energyFrac;
      sumFrac += ide.energyFrac;
    }

    // BackTracker normalises energyFrac against *all* the ionisation in the hit's time window,
    // including deposits it cannot attribute to a particle, so the fractions sum to <= 1 and the
    // remainder is genuinely unmatched charge rather than a rounding artefact.
    out.frac[static_cast<std::size_t>(HitOrigin::kUnmatched)] += std::max(0.f, 1.f - sumFrac);

    float best = -1;
    for(std::size_t c = 0; c < static_cast<std::size_t>(HitOrigin::kN); ++c){
      if(out.frac[c] > best){
        best = out.frac[c];
        out.dominant = static_cast<HitOrigin>(c);
      }
    }

    return out;
  }

  //------------------------------------------------------------------------------

  void CAFMaker::AccumulateHitSummary(const std::vector<art::Ptr<recob::Hit>> &hits,
                                      detinfo::DetectorClocksData const& clockData,
                                      detinfo::DetectorPropertiesData const& detProp,
                                      const TruthClassification &cls,
                                      HitTruthCache &cache,
                                      caf::SRHitSummary &out) const
  {
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    const int coll = CollectionPlane();

    // Switch from the "not computed" defaults to real counters.
    out.nhits = 0;
    out.nhits_coll = 0;
    out.charge_coll = 0;
    for(std::size_t c = 0; c < static_cast<std::size_t>(HitOrigin::kN); ++c){
      caf::SRHitCategory &cat = CategoryRef(out, static_cast<HitOrigin>(c));
      cat.nhits = 0;
      cat.nhits_coll = 0;
      cat.nhits_frac_coll = 0;
      cat.charge_coll = 0;
    }

    for(art::Ptr<recob::Hit> const& hit : hits){
      const HitTruth &truth = cache.Get(hit, clockData, *bt_serv, cls);

      ++out.nhits;
      ++CategoryRef(out, truth.dominant).nhits;

      // Charge only means something on the collection plane, and summing ionisation over all
      // three planes would count each deposit three times.
      if(static_cast<int>(hit->WireID().Plane) != coll) continue;

      ++out.nhits_coll;
      ++CategoryRef(out, truth.dominant).nhits_coll;

      const double charge = dune_ana::DUNEAnaHitUtils::LifetimeCorrection(clockData, detProp, hit)*hit->Integral();
      out.charge_coll += charge;

      // Split mixed hits between categories by true energy fraction rather than assigning the
      // whole hit to its dominant contributor: that is the weighting energy resolution cares about.
      for(std::size_t c = 0; c < static_cast<std::size_t>(HitOrigin::kN); ++c){
        caf::SRHitCategory &cat = CategoryRef(out, static_cast<HitOrigin>(c));
        cat.nhits_frac_coll += truth.frac[c];
        cat.charge_coll     += charge*truth.frac[c];
      }
    }

    out.E_reco = ChargeToEnergyGeV(out.charge_coll, coll);
    for(std::size_t c = 0; c < static_cast<std::size_t>(HitOrigin::kN); ++c){
      caf::SRHitCategory &cat = CategoryRef(out, static_cast<HitOrigin>(c));
      cat.E_reco = ChargeToEnergyGeV(cat.charge_coll, coll);
    }
  }

  //------------------------------------------------------------------------------

  void CAFMaker::FillSimChannelBeamEnergy(const TruthClassification &cls,
                                          caf::SRInteractionBranch &ixn) const
  {
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    const geo::View_t collView = static_cast<geo::View_t>(CollectionPlane());

    double eMain = 0, eDaughters = 0, eContam = 0;

    // One sweep over the SimChannels. Asking BackTracker for the IDEs of each beam track ID
    // instead would rescan every SimChannel once per track ID, which is minutes per event.
    for(art::Ptr<sim::SimChannel> const& simchannel : bt_serv->SimChannels()){
      if(fWireReadout->View(simchannel->Channel()) != collView) continue; //Avoid triple counting

      for(auto const& tdcide : simchannel->TDCIDEMap()){
        for(sim::IDE const& ide : tdcide.second){
          auto it = cls.g4ToOrigin.find(std::abs(ide.trackID));
          if(it == cls.g4ToOrigin.end()) continue;

          if     (it->second == HitOrigin::kBeam)         eMain      += ide.energy;
          else if(it->second == HitOrigin::kBeamDaughter) eDaughters += ide.energy;
          else if(it->second == HitOrigin::kBeamContam)   eContam    += ide.energy;
        }
      }
    }

    ixn.beam_E_true_simch               = eMain*1e-3;      //MeV -> GeV
    ixn.beam_daughters_E_true_simch     = eDaughters*1e-3;
    ixn.beam_contamination_E_true_simch = eContam*1e-3;
  }

  //------------------------------------------------------------------------------

  void CAFMaker::ComputeHitSummaries(const art::Event &evt,
                                     const std::vector<art::Ptr<recob::Slice>> &slicePtrs,
                                     caf::SRInteractionBranch &ixn,
                                     std::vector<caf::SRHitSummary> &sliceSummaries) const
  {
    sliceSummaries.clear();
    if(!fFillHitSummaries) return;
    if(evt.isRealData()) return; //No SimChannels to back-track against

    auto hitHandle = evt.getHandle<std::vector<recob::Hit>>(fHitLabel);
    if(!hitHandle){
      mf::LogWarning("CAFMaker") << "No Hit collection found with label " << fHitLabel
                                 << ". Hit summaries will be left unfilled.";
      return;
    }

    const TruthClassification cls = BuildTruthClassification(evt);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    // Local to this call, so no mutable member is needed despite this being a const method.
    HitTruthCache cache;

    std::vector<art::Ptr<recob::Hit>> allHits;
    art::fill_ptr_vector(allHits, hitHandle);
    AccumulateHitSummary(allHits, clockData, detProp, cls, cache, ixn.allhits);

    sliceSummaries.resize(slicePtrs.size());
    bool productsConsistent = true;

    int bestBeamHits = -1;
    int nSlicesWithBeam = 0;

    for(art::Ptr<recob::Slice> const& slicePtr : slicePtrs){
      std::vector<art::Ptr<recob::Hit>> sliceHits = dune_ana::DUNEAnaSliceUtils::GetHits(slicePtr, evt, fPandoraLabel);

      // If Pandora clustered a different hit collection than HitLabel, the event-wide numbers
      // are not a valid denominator for the per-slice ones. Flag it rather than silently
      // producing completeness values above 1.
      if(!sliceHits.empty() && sliceHits.front().id() != hitHandle.id()) productsConsistent = false;

      caf::SRHitSummary &summary = sliceSummaries[slicePtr.key()];
      AccumulateHitSummary(sliceHits, clockData, detProp, cls, cache, summary);

      const int beamTreeHits = summary.beam.nhits + summary.beam_daughters.nhits;
      if(beamTreeHits > 0) ++nSlicesWithBeam;
      if(beamTreeHits > bestBeamHits){
        bestBeamHits = beamTreeHits;
        ixn.best_beam_slice_id = static_cast<int>(slicePtr.key());
      }
    }

    if(!productsConsistent){
      mf::LogWarning("CAFMaker") << "Slice hits come from a different art product than HitLabel ('" << fHitLabel
                                 << "'). Per-slice hit counts are not bounded by the event-wide ones; "
                                 << "SRHitSummary::consistent_products is set to false.";
      ixn.allhits.consistent_products = false;
      for(caf::SRHitSummary &summary : sliceSummaries) summary.consistent_products = false;
    }

    ixn.nslices = static_cast<int>(slicePtrs.size());
    ixn.nslices_with_beam_hits = nSlicesWithBeam;

    const int eventBeamTreeHits = ixn.allhits.beam.nhits + ixn.allhits.beam_daughters.nhits;
    if(eventBeamTreeHits > 0 && bestBeamHits >= 0){
      ixn.beam_hit_completeness_best = static_cast<float>(bestBeamHits)/eventBeamTreeHits;
    }

    FillSimChannelBeamEnergy(cls, ixn);
  }


  //------------------------------------------------------------------------------
  void CAFMaker::FillTruthInfo(caf::SRTruthBranch& truthBranch, art::Event const& evt)
  {
    truthBranch.nu.clear();
    truthBranch.nnu = 0;

    size_t cumulativeInteractionIndex = 0;

    art::Handle<std::vector<simb::GTruth>> gtruthHandle;
    art::Handle<std::vector<simb::MCFlux>> fluxHandle;

    const bool haveGTruthHandle = evt.getByLabel(fGTruthLabel, gtruthHandle) && gtruthHandle.isValid();
    const bool haveFluxHandle   = evt.getByLabel(fMCFluxLabel, fluxHandle)   && fluxHandle.isValid();

    if (!haveGTruthHandle) {
      mf::LogWarning("CAFMaker") << "GTruth product not valid. GENIE-derived truth fields will be left unfilled.";
    }
    if (!haveFluxHandle) {
      mf::LogWarning("CAFMaker") << "MCFlux product not valid. Flux-derived truth fields will be left unfilled.";
    }

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();

    for (const auto& MCTruthLabel : fMCTruthLabel) {
      art::Handle<std::vector<simb::MCTruth>> mctruthHandle;
      evt.getByLabel(MCTruthLabel, mctruthHandle);

      if (!mctruthHandle.isValid()) {
        mf::LogWarning("CAFMaker")
          << "MCTruth vector for label '" << MCTruthLabel
          << "' is not valid. Skipping this label.";
        continue;
      }

      const auto& mctruthVec = *mctruthHandle;

      art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> fmParticles(mctruthHandle, evt, fG4Label);
      if (!fmParticles.isValid()) {
        mf::LogWarning("CAFMaker") << "MCTruth->MCParticle associations for label '"
          << MCTruthLabel << "' not found. Prim/sec filling will be skipped for this label.";
      }

      for (size_t i = 0; i < mctruthVec.size(); ++i) {
        const simb::MCTruth& mct = mctruthVec[i];
        const bool hasNu = mct.NeutrinoSet();

        const bool haveThisGTruth = haveGTruthHandle && i < gtruthHandle->size();
        const bool haveThisFlux   = haveFluxHandle   && i < fluxHandle->size();

        if (hasNu) {
          caf::SRTrueInteraction inter;
          inter.id = cumulativeInteractionIndex++;

          // Optional GENIE record
          if (haveThisGTruth) {
            inter.genieIdx = FillGENIERecord(mct, (*gtruthHandle)[i]);
          }

          // Generator info is meaningful for both neutrino and non-neutrino MCTruth
          const simb::MCGeneratorInfo& genInfo = mct.GeneratorInfo();
          auto it = fgenMap.find(genInfo.generator);
          inter.generator = (it != fgenMap.end()) ? it->second : caf::Generator::kUnknownGenerator;

          inter.genVersion.clear();
          if (!genInfo.generatorVersion.empty()) {
            size_t last = 0;
            size_t next = 0;
            std::string s(genInfo.generatorVersion);
            while ((next = s.find('.', last)) != std::string::npos) {
              inter.genVersion.push_back(std::stoi(s.substr(last, next - last)));
              last = next + 1;
            }
            inter.genVersion.push_back(std::stoi(s.substr(last)));
          }

          // Neutrino-specific block
          const simb::MCNeutrino& neutrino = mct.GetNeutrino();

          inter.pdg = neutrino.Nu().PdgCode();
          inter.iscc = !(neutrino.CCNC());
          inter.mode = static_cast<caf::ScatteringMode>(neutrino.Mode());
          inter.E = neutrino.Nu().E();

          inter.vtx.SetX(neutrino.Lepton().Vx());
          inter.vtx.SetY(neutrino.Lepton().Vy());
          inter.vtx.SetZ(neutrino.Lepton().Vz());
          inter.time = neutrino.Lepton().T();

          inter.momentum.SetX(neutrino.Nu().Momentum().X());
          inter.momentum.SetY(neutrino.Nu().Momentum().Y());
          inter.momentum.SetZ(neutrino.Nu().Momentum().Z());

          inter.W = neutrino.W();
          inter.Q2 = neutrino.QSqr();
          inter.bjorkenX = neutrino.X();
          inter.inelasticity = neutrino.Y();

          {
            TLorentzVector q = neutrino.Nu().Momentum() - neutrino.Lepton().Momentum();
            inter.q0 = q.E();
            inter.modq = q.Vect().Mag();
          }

          inter.isvtxcont = IsVertexContained(inter.vtx);

          // Optional GTruth-dependent neutrino fields
          if (haveThisGTruth) {
            const simb::GTruth& gt = (*gtruthHandle)[i];
            inter.targetPDG = gt.ftgtPDG;
            inter.ischarm = gt.fIsCharm;
            inter.isseaquark = gt.fIsSeaQuark;
            inter.resnum = gt.fResNum;
            inter.xsec = gt.fXsec;
            inter.genweight = gt.fweight;
            inter.t = gt.fgT;

            inter.hitnuc = neutrino.HitNuc();
            auto* pdef = genie::PDGLibrary::Instance()->Find(inter.hitnuc);
            if (pdef) {
              const double nucMass = pdef->Mass();
              inter.removalE = nucMass - gt.fHitNucP4.E();
            }
          }
          else {
            inter.hitnuc = neutrino.HitNuc();
          }

          // Optional flux-dependent fields
          if (haveThisFlux) {
            inter.pdgorig = (*fluxHandle)[i].fntype;
          }

          // TODO fields still unavailable in the current chain:
          // inter.baseline
          // inter.prod_vtx
          // inter.parent_dcy_mom
          // inter.parent_dcy_mode
          // inter.parent_pdg
          // inter.parent_dcy_E
          // inter.imp_weight

          inter.nproton = 0;
          inter.nneutron = 0;
          inter.npip = 0;
          inter.npim = 0;
          inter.npi0 = 0;
          inter.nprim = 0;
          inter.nprefsi = 0;
          inter.nsec = 0;

          // Pre-FSI particles: only meaningful for neutrino generator record
          for (int p = 0; p < mct.NParticles(); ++p) {
            const simb::MCParticle& mcpart = mct.GetParticle(p);
            if (mcpart.StatusCode() != genie::EGHepStatus::kIStHadronInTheNucleus) continue;

            caf::SRTrueParticle part;
            part.pdg = mcpart.PdgCode();
            part.G4ID = mcpart.TrackId();
            part.interaction_id = inter.id;
            part.time = mcpart.T();
            part.p = caf::SRLorentzVector(mcpart.Momentum());
            part.end_p = caf::SRLorentzVector(mcpart.EndMomentum());
            part.start_pos = caf::SRVector3D(mcpart.Position().Vect());
            part.end_pos = caf::SRVector3D(mcpart.EndPosition().Vect());
            part.parent = mcpart.Mother();

            inter.prefsi.push_back(std::move(part));
            ++inter.nprefsi;
          }

          if (fmParticles.isValid()) {
            for (const art::Ptr<simb::MCParticle>& mcpart : fmParticles.at(i)) {
              const int tid = mcpart->TrackId();
              if (tid == 0) continue;

              bool isPrimary = (mcpart->Mother() == 0);
              fMCParticlesMap[tid] = std::make_tuple(mcpart, inter.id, isPrimary,
                  isPrimary ? (int)inter.prim.size() : (int)inter.sec.size());

              caf::SRTrueParticle part;
              part.pdg = mcpart->PdgCode();
              part.G4ID = tid;
              part.interaction_id = inter.id;
              part.time = mcpart->T();
              part.p = caf::SRLorentzVector(mcpart->Momentum());
              part.end_p = caf::SRLorentzVector(mcpart->EndMomentum());
              part.start_pos = caf::SRVector3D(mcpart->Position().Vect());
              part.end_pos = caf::SRVector3D(mcpart->EndPosition().Vect());
              part.parent = mcpart->Mother();

              for (int d = 0; d < mcpart->NumberDaughters(); ++d) {
                part.daughters.push_back(mcpart->Daughter(d));
              }

              if (fMCParticlesMap.count(mcpart->Mother()) > 0) {
                art::Ptr<simb::MCParticle> ancestor_ptr;
                bool ancestor_is_primary = false;
                int ancestor_id = -1;
                int ancestor_ixn = -1;
                std::tie(ancestor_ptr, ancestor_ixn, ancestor_is_primary, ancestor_id) =
                  fMCParticlesMap.at(mcpart->Mother());

                caf::TrueParticleID ancestor;
                ancestor.type = ancestor_is_primary ? caf::TrueParticleID::kPrimary
                                                    : caf::TrueParticleID::kSecondary;
                ancestor.ixn = ancestor_ixn;
                ancestor.part = ancestor_id;
                part.ancestor_id = ancestor;
              }

              if (isPrimary) {
                inter.prim.push_back(std::move(part));
                switch (mcpart->PdgCode()) {
                  case 2212: ++inter.nproton; break;
                  case 2112: ++inter.nneutron; break;
                  case 211:  ++inter.npip;    break;
                  case -211: ++inter.npim;    break;
                  case 111:  ++inter.npi0;    break;
                  default: break;
                }
              }
              else {
                inter.sec.push_back(std::move(part));
              }
            }
          }

          inter.nprim = inter.prim.size();
          inter.nsec  = inter.sec.size();

          // Fill parentID and daughtersID now that prim/sec vectors are complete
          // and fMCParticlesMap holds the final SR indices for all particles in this interaction.
          auto fill_id_refs = [&](std::vector<caf::SRTrueParticle>& particles) {
            for (auto& p : particles) {
              if (fMCParticlesMap.count(p.parent) > 0) {
                bool par_isPrimary; int par_ixn, par_idx;
                std::tie(std::ignore, par_ixn, par_isPrimary, par_idx) = fMCParticlesMap.at(p.parent);
                p.parentID.ixn  = par_ixn;
                p.parentID.type = par_isPrimary ? caf::TrueParticleID::kPrimary : caf::TrueParticleID::kSecondary;
                p.parentID.part = par_idx;
              }
              for (unsigned int dTid : p.daughters) {
                caf::TrueParticleID dID;
                if (fMCParticlesMap.count(dTid) > 0) {
                  bool d_isPrimary; int d_ixn, d_idx;
                  std::tie(std::ignore, d_ixn, d_isPrimary, d_idx) = fMCParticlesMap.at(dTid);
                  dID.ixn  = d_ixn;
                  dID.type = d_isPrimary ? caf::TrueParticleID::kPrimary : caf::TrueParticleID::kSecondary;
                  dID.part = d_idx;
                }
                p.daughtersID.push_back(dID);
              }
            }
          };
          fill_id_refs(inter.prim);
          fill_id_refs(inter.sec);

          truthBranch.nu.push_back(std::move(inter));
        }
        else {
          // Non-neutrino MCTruth:
          // each visible source particle becomes its own truth interaction entry.
          auto const& assoc_parts = fmParticles.at(i);
          auto const& assoc_infos = fmParticles.data(i);
          for (size_t j = assoc_parts.size(); j-- > 0; ) {
            if (!assoc_infos[j]->hasGeneratedParticleIndex()) continue;
            const simb::MCParticle& mcpart = *assoc_parts[j];

            caf::SRTrueInteraction inter;
            inter.id = cumulativeInteractionIndex++;
            const int trackID = std::abs(mcpart.TrackId());
            auto mcPartIt = fMCParticlesMap.find(trackID);
            if (mcPartIt != fMCParticlesMap.end()) {
              mcPartIt->second = std::make_tuple(std::get<0>(mcPartIt->second), inter.id, true, inter.prim.size());
            }
            else {
             // create warning
              mf::LogWarning("CAFMaker") << "MCParticle with TrackId " << trackID << " not found in preloaded map. This particle will not be included in the truth record.";
              continue;
            }
            
            // Generator info can still be filled
            const simb::MCGeneratorInfo& genInfo = mct.GeneratorInfo();
            auto it = fgenMap.find(genInfo.generator);
            inter.generator = (it != fgenMap.end()) ? it->second : caf::Generator::kUnknownGenerator;

            inter.genVersion.clear();
            if (!genInfo.generatorVersion.empty()) {
              size_t last = 0;
              size_t next = 0;
              std::string s(genInfo.generatorVersion);
              while ((next = s.find('.', last)) != std::string::npos) {
                inter.genVersion.push_back(std::stoi(s.substr(last, next - last)));
                last = next + 1;
              }
              inter.genVersion.push_back(std::stoi(s.substr(last)));
            }

            caf::SRTrueParticle part;
            part.pdg = mcpart.PdgCode();
            part.G4ID = std::abs(mcpart.TrackId());
            part.interaction_id = inter.id;
            part.time = mcpart.T();
            part.p = caf::SRLorentzVector(mcpart.Momentum());
            part.end_p = caf::SRLorentzVector(mcpart.EndMomentum());
            part.start_pos = caf::SRVector3D(mcpart.Position().Vect());
            part.end_pos = caf::SRVector3D(mcpart.EndPosition().Vect());
            part.parent = mcpart.Mother();

            if (fMCParticlesMap.count(mcpart.Mother()) > 0) {
              art::Ptr<simb::MCParticle> ancestor_ptr;
              bool ancestor_is_primary = false;
              int ancestor_id = -1;
              int ancestor_ixn = -1;
              std::tie(ancestor_ptr, ancestor_ixn, ancestor_is_primary, ancestor_id) =
                fMCParticlesMap.at(mcpart.Mother());

              caf::TrueParticleID ancestor;
              ancestor.type = ancestor_is_primary ? caf::TrueParticleID::kPrimary
                                                  : caf::TrueParticleID::kSecondary;
              ancestor.ixn = ancestor_ixn;
              ancestor.part = ancestor_id;
              part.ancestor_id = ancestor;
            }
            std::deque<int> secondaries_to_add;

            const simb::MCParticle* the_g4_part = nullptr;
            
            const int thisTrackId = std::abs(mcpart.TrackId());
            auto g4It = plist.find(thisTrackId);
            if (g4It != plist.end()) {
              the_g4_part = g4It->second;
            }

            if (the_g4_part) {
              for (int d = 0; d < the_g4_part->NumberDaughters(); ++d) {
                int daughter = the_g4_part->Daughter(d);
                part.daughters.push_back(daughter);
                secondaries_to_add.push_back(daughter);
              }
            }

            part.true_beam_interactingEnergy = ComputeTrueInteractingEnergy(evt, &mcpart);
            inter.prim.push_back(std::move(part));
            inter.nprim = 1;

            while (!secondaries_to_add.empty()) {
              int this_sec_ID = secondaries_to_add.front();
              secondaries_to_add.pop_front();

              if (fMCParticlesMap.count(this_sec_ID) == 0) continue;
              fMCParticlesMap[this_sec_ID] = std::make_tuple(std::get<0>(fMCParticlesMap.at(this_sec_ID)), inter.id, false, inter.sec.size());

              auto sec_mcpart = std::get<0>(fMCParticlesMap.at(this_sec_ID));

              caf::SRTrueParticle sec;
              sec.pdg = sec_mcpart->PdgCode();
              sec.G4ID = std::abs(sec_mcpart->TrackId());
              sec.interaction_id = inter.id;
              sec.time = sec_mcpart->T();
              sec.p = caf::SRLorentzVector(sec_mcpart->Momentum());
              sec.end_p = caf::SRLorentzVector(sec_mcpart->EndMomentum());
              sec.start_pos = caf::SRVector3D(sec_mcpart->Position().Vect());
              sec.end_pos = caf::SRVector3D(sec_mcpart->EndPosition().Vect());
              sec.parent = sec_mcpart->Mother();

              if (fMCParticlesMap.count(sec_mcpart->Mother()) > 0) {
                art::Ptr<simb::MCParticle> ancestor_ptr;
                bool ancestor_is_primary = false;
                int ancestor_id = -1;
                int ancestor_ixn = -1;
                std::tie(ancestor_ptr, ancestor_ixn, ancestor_is_primary, ancestor_id) =
                  fMCParticlesMap.at(sec_mcpart->Mother());

                caf::TrueParticleID ancestor;
                ancestor.type = ancestor_is_primary ? caf::TrueParticleID::kPrimary
                                                    : caf::TrueParticleID::kSecondary;
                ancestor.ixn = ancestor_ixn;
                ancestor.part = ancestor_id;
                sec.ancestor_id = ancestor;
              }


              for (int d = 0; d < sec_mcpart->NumberDaughters(); ++d) {
                int daughter = sec_mcpart->Daughter(d);
                sec.daughters.push_back(daughter);
                secondaries_to_add.push_back(daughter);
              }

              sec.true_beam_interactingEnergy = ComputeTrueInteractingEnergy(evt, sec_mcpart.get());
              inter.sec.push_back(std::move(sec));
            }

            inter.nsec = inter.sec.size();

            // Fill parentID and daughtersID now that prim/sec are complete
            auto fill_id_refs_nonu = [&](std::vector<caf::SRTrueParticle>& particles) {
              for (auto& p : particles) {
                if (fMCParticlesMap.count(p.parent) > 0) {
                  bool par_isPrimary; int par_ixn, par_idx;
                  std::tie(std::ignore, par_ixn, par_isPrimary, par_idx) = fMCParticlesMap.at(p.parent);
                  p.parentID.ixn  = par_ixn;
                  p.parentID.type = par_isPrimary ? caf::TrueParticleID::kPrimary : caf::TrueParticleID::kSecondary;
                  p.parentID.part = par_idx;
                }
                for (unsigned int dTid : p.daughters) {
                  caf::TrueParticleID dID;
                  if (fMCParticlesMap.count(dTid) > 0) {
                    bool d_isPrimary; int d_ixn, d_idx;
                    std::tie(std::ignore, d_ixn, d_isPrimary, d_idx) = fMCParticlesMap.at(dTid);
                    dID.ixn  = d_ixn;
                    dID.type = d_isPrimary ? caf::TrueParticleID::kPrimary : caf::TrueParticleID::kSecondary;
                    dID.part = d_idx;
                  }
                  p.daughtersID.push_back(dID);
                }
              }
            };
            fill_id_refs_nonu(inter.prim);
            fill_id_refs_nonu(inter.sec);

            // Optional visible-particle vertex proxy from the particle start point
            inter.vtx.SetX(mcpart.Vx());
            inter.vtx.SetY(mcpart.Vy());
            inter.vtx.SetZ(mcpart.Vz());
            inter.time = mcpart.T();
            inter.momentum.SetX(mcpart.Px());
            inter.momentum.SetY(mcpart.Py());
            inter.momentum.SetZ(mcpart.Pz());
            inter.isvtxcont = IsVertexContained(inter.vtx);

            truthBranch.nu.push_back(std::move(inter));
          }
        }
      }
    }

    truthBranch.nnu = truthBranch.nu.size();
  }

  //------------------------------------------------------------------------------

  void CAFMaker::FillPFPMetadata(caf::SRPFP &pfpBranch, art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const {
    art::Ptr<larpandoraobj::PFParticleMetadata> metadata = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfp, evt, fPandoraLabel);
      if(metadata.isNull()){
        mf::LogWarning("CAFMaker") << "No metadata found for PFP with ID " << pfp->Self();
        return;
      }

      std::map<std::string, float> properties = metadata->GetPropertiesMap();

      std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fPandoraLabel);
      std::vector<art::Ptr<recob::Hit>> hits_U = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, 0);
      std::vector<art::Ptr<recob::Hit>> hits_V = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, 1);
      std::vector<art::Ptr<recob::Hit>> hits_W = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, 2);
      std::vector<art::Ptr<recob::SpacePoint>> sps = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfp, evt, fSpacePointLabel);

      pfpBranch.nhits_U = hits_U.size();
      pfpBranch.nhits_V = hits_V.size();
      pfpBranch.nhits_W = hits_W.size();
      pfpBranch.nhits_3D = sps.size();

      std::map<std::string, float*> property_mapping = {
        {"LArPfoHierarchyFeatureTool_DaughterParentHitRatio", &pfpBranch.daughter_parent_hit_ratio},
        {"LArPfoHierarchyFeatureTool_NDaughterHits3D", &pfpBranch.ndaughters_hit_3d},
        {"LArThreeDChargeFeatureTool_EndFraction", &pfpBranch.charge_end_fraction},
        {"LArThreeDChargeFeatureTool_FractionalSpread", &pfpBranch.charge_fractional_spread},
        {"LArThreeDLinearFitFeatureTool_DiffStraightLineMean", &pfpBranch.diff_straight_line_mean},
        {"LArThreeDLinearFitFeatureTool_Length", &pfpBranch.line_length},
        {"LArThreeDLinearFitFeatureTool_MaxFitGapLength", &pfpBranch.max_fit_gap_length},
        {"LArThreeDLinearFitFeatureTool_SlidingLinearFitRMS", &pfpBranch.sliding_linear_fit_rms},
        {"LArThreeDOpeningAngleFeatureTool_AngleDiff", &pfpBranch.angle_diff_3d},
        {"LArThreeDPCAFeatureTool_SecondaryPCARatio", &pfpBranch.secondary_pca_ratio},
        {"LArThreeDPCAFeatureTool_TertiaryPCARatio", &pfpBranch.tertiary_pca_ratio},
        {"LArThreeDVertexDistanceFeatureTool_VertexDistance", &pfpBranch.vertex_distance},
        {"TrackScore", &pfpBranch.track_score},
      };

      for(const auto& [property_name, property_value] : property_mapping) {
        if(properties.find(property_name) != properties.end()) {
          *(property_value) = properties[property_name];
        } 
      }

  }

  //------------------------------------------------------------------------------

  void CAFMaker::GetMVAResults(caf::SRPFP & output_pfp, art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt, anab::MVAReader<recob::Hit,4> * hitResults, int planeid, bool charge_weighted) const
  {
    if (!hitResults) {
      mf::LogWarning("CAFMaker") << "Null pointer provided for hitResults in GetMVAResults";
      return;
    }

    if (planeid < -1 || planeid > 2) {
      std::stringstream ss;
      ss << "Unknown planeid (" << planeid << ") provided to GetMVAResults";
      throw std::runtime_error(
        ss.str()
      );
    }

    //First getting all the hits belonging to that PFP
    std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fPandoraLabel);
    if (planeid != -1) {
      hits = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, planeid);
    }

    float denom = 0.;
    output_pfp.cnn_stem_scores.charge_weighted = charge_weighted;
    output_pfp.cnn_stem_scores.plane_ID = planeid;
    for (const auto & hit : hits){
      auto output = hitResults->getOutput(hit);

      float scale = (charge_weighted ? hit->Integral() : 1.);

      output_pfp.cnn_stem_scores.track_score  += scale*output[hitResults->getIndex("track")];
      output_pfp.cnn_stem_scores.shower_score += scale*output[hitResults->getIndex("em")];
      output_pfp.cnn_stem_scores.empty_score  += scale*output[hitResults->getIndex("none")];
      output_pfp.cnn_stem_scores.michel_score += scale*output[hitResults->getIndex("michel")];
      denom += scale;
    }

    if (denom > 0.) {
      output_pfp.cnn_stem_scores /= denom;
    }
    else {
      output_pfp.cnn_stem_scores.track_score  = output_pfp.NaN;
      output_pfp.cnn_stem_scores.shower_score = output_pfp.NaN;
      output_pfp.cnn_stem_scores.empty_score  = output_pfp.NaN;
      output_pfp.cnn_stem_scores.michel_score = output_pfp.NaN;
    }
  }

  //------------------------------------------------------------------------------
  void CAFMaker::FillRecoInfoSliceLoop(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, caf::SRTruthBranch &truthBranch, const art::Event &evt, const caf::SRBeamInstrumentation &beamInst) const
  {
    // get handle to slices
    auto sliceHandle = evt.getHandle<std::vector<recob::Slice>>(fPandoraLabel);
    if (!sliceHandle) {
      mf::LogWarning("CAFMaker") << "No Slice collection found with label " << fPandoraLabel;
      return;
    }
    std::vector<art::Ptr<recob::Slice>> slicePtrs;
    art::fill_ptr_vector(slicePtrs, sliceHandle);

    art::FindManyP<recob::PFParticle> sliceToPFP(sliceHandle, evt, fPandoraLabel);
    if (!sliceToPFP.isValid()) {
       mf::LogWarning("CAFMaker") << "Slice->PFParticle associations not found for label " << fPandoraLabel;
       return;
    }
  
    // Done once for the whole event: back-tracking every hit is expensive, and the per-slice
    // summaries and the event-wide denominators share the same cache.
    std::vector<caf::SRHitSummary> sliceSummaries;
    ComputeHitSummaries(evt, slicePtrs, recoBranch.ixn, sliceSummaries);

    for (const auto& slicePtr : slicePtrs) {
      FillRecoInfo(recoBranch, fdBranch, truthBranch, evt, slicePtr, sliceToPFP, sliceSummaries, beamInst);
    }

  }

  void CAFMaker::FillRecoInfo(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, caf::SRTruthBranch &truthBranch, const art::Event &evt, const art::Ptr<recob::Slice> &slicePtr, const art::FindManyP<recob::PFParticle> &sliceToPFP, const std::vector<caf::SRHitSummary> &sliceSummaries, const caf::SRBeamInstrumentation &beamInst) const {
    SRInteractionBranch &ixn = recoBranch.ixn;

    //Only filling with Pandora Reco for the moment
    std::vector<SRInteraction> &pandora = ixn.pandora;

    lar_pandora::PFParticleVector particleVector = sliceToPFP.at(slicePtr.key());

    lar_pandora::VertexVector vertexVector;
    lar_pandora::PFParticlesToVertices particlesToVertices;
    lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraLabel, vertexVector, particlesToVertices);

    auto pfParticles_beam = evt.getValidHandle<std::vector<recob::PFParticle>>(fPandoraLabel);
    const art::FindManyP<larpandoraobj::PFParticleMetadata> findMetaData_beam(pfParticles_beam, evt, fPandoraLabel);

    bool found_beam = false;
    for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      bool is_test_beam = false;
      if (findMetaData_beam.isValid()) {
        const auto& mdVec = findMetaData_beam.at(particle->Self());
        if (!mdVec.empty() && mdVec.at(0).isNonnull()) {
          const auto& mdMap = mdVec.at(0)->GetPropertiesMap();
          is_test_beam = (mdMap.find("IsTestBeam") != mdMap.end());
        }
      }
      if (is_test_beam) {
        found_beam = true;
        break;
      }
    }

    bool fill_info_condition = false;
    for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      bool is_test_beam = false;
      if (findMetaData_beam.isValid()) {
        const auto& mdVec = findMetaData_beam.at(particle->Self());
        if (!mdVec.empty() && mdVec.at(0).isNonnull()) {
          const auto& mdMap = mdVec.at(0)->GetPropertiesMap();
          is_test_beam = (mdMap.find("IsTestBeam") != mdMap.end());
        }
      }

      fill_info_condition = false;
      if (found_beam) {
        fill_info_condition = is_test_beam && particle->IsPrimary();
      } else {
        fill_info_condition = particle->IsPrimary();
      }

      if(fill_info_condition){
      // if(particle->IsPrimary()){

        SRInteraction reco;

        reco.vtx = SRVector3D(-999, -999, -999); //Setting an unambiguous default value if no vertex is found

        //Retrieving the reco vertex
        lar_pandora::PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(particle);
        if (particlesToVertices.end() != vIter) {
          const lar_pandora::VertexVector &vertexVector = vIter->second;
          if (vertexVector.size() == 1) {
            const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
            double xyz[3] = {0.0, 0.0, 0.0} ;
            vertex->XYZ(xyz);
            reco.vtx = SRVector3D(xyz[0], xyz[1], xyz[2]);
          }
        }

        SRDirectionBranch &dir = reco.dir;
        FillDirectionInfo(dir, evt);


        //Neutrino flavours hypotheses
        SRNeutrinoHypothesisBranch &nuhyp = reco.nuhyp;
        //Filling only CVN at the moment.
        FillCVNInfo(nuhyp.cvn, evt);

        //Neutrino energy hypothese
        SRNeutrinoEnergyBranch &Enu = reco.Enu;
        FillEnergyInfo(Enu, evt);

        //List of reconstructed particles
        SRRecoParticlesBranch &part = reco.part;
        FillRecoParticlesInfo(part, fdBranch, evt, slicePtr, sliceToPFP);

        //Assuming a single TrueInteraction for now. TODO: Change this if several interactions end up being simulated in the same event
        reco.truth = {0}; 
        reco.truthOverlap = {1.};
        reco.isFromTrigger = is_test_beam;

        if (is_test_beam) {
          art::Ptr<recob::Track> beamTrack;
          if (dune_ana::DUNEAnaPFParticleUtils::HasTrack(particle, evt, fPandoraLabel, fTrackLabel)) {
            beamTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(particle, evt, fPandoraLabel, fTrackLabel);
          }
          if (beamTrack.isNonnull()) {
            auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
            ComputeRecoInteractingEnergy(evt, beamTrack.get(), clockData, ixn, truthBranch, beamInst);
          } else {
            mf::LogWarning("CAFMaker") << "No track associated to beam particle — skipping reco interacting energy.";
          }
        }

        //Record which slice this interaction came from, so it can be matched up offline with
        //the per-slice hit composition and with the other slices in the event.
        reco.id = static_cast<long int>(slicePtr.key());
        if (slicePtr.key() < sliceSummaries.size()) {
          reco.hits = sliceSummaries[slicePtr.key()];

          if (is_test_beam) {
            ixn.beam_slice_id = static_cast<int>(slicePtr.key());
            ixn.beam_slice_index = static_cast<int>(pandora.size()); //Index this record is about to take

            const int sliceBeamTreeHits = reco.hits.beam.nhits + reco.hits.beam_daughters.nhits;
            const int eventBeamTreeHits = ixn.allhits.beam.nhits + ixn.allhits.beam_daughters.nhits;
            if (eventBeamTreeHits > 0) {
              ixn.beam_hit_completeness_beam = static_cast<float>(sliceBeamTreeHits)/eventBeamTreeHits;
            }
            if (reco.hits.nhits > 0) {
              ixn.beam_hit_purity_beam = static_cast<float>(sliceBeamTreeHits)/reco.hits.nhits;
            }
          }
        }

        pandora.emplace_back(reco);
        break; //We do this work only the first primary particle in the slice to avoid repetitions
      }
    }

    ixn.npandora = pandora.size();
    ixn.ndlp = ixn.dlp.size();

  }

  //------------------------------------------------------------------------------


  void CAFMaker::FillTruthMatchingAndOverlap(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt, std::vector<TrueParticleID> &truth, std::vector<float> &truthOverlap) const{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    //First getting all the hits belonging to that PFP
    std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fPandoraLabel);

    TruthMatchUtils::IDToEDepositMap idToEDepositMap;
    for (const art::Ptr<recob::Hit>& pHit : hits){
      TruthMatchUtils::FillG4IDToEnergyDepositMap(idToEDepositMap, clockData, pHit, true);
    }

    float totalEDeposit = 0;
    for (const auto& [id, eDeposit] : idToEDepositMap) {
      totalEDeposit += eDeposit;
    }

    if (totalEDeposit <= 0) {
      mf::LogWarning("CAFMaker") << "No energy deposit found for PFP with ID " << pfp->Self() << ". Skipping truth matching.";
      return;
    }

    for (const auto& [id, eDeposit] : idToEDepositMap) {
      if (fMCParticlesMap.count(id) == 0) {
        mf::LogWarning("CAFMaker") << "No MCParticle found with ID " << id << " for PFP with ID " << pfp->Self() << ". Skipping.";
        continue;
      }

      art::Ptr<simb::MCParticle> mcpart;
      bool isPrimary;
      int srID;
      int ixnID=-1;
      std::tie(mcpart, ixnID, isPrimary, srID) = fMCParticlesMap.at(id);
      caf::TrueParticleID::PartType type = isPrimary ? TrueParticleID::kPrimary : TrueParticleID::kSecondary;

      TrueParticleID truePart;
      truePart.type = type;
      // truePart.ixn = 0; //This assumes only 1 interaction
      truePart.ixn = ixnID;
      truePart.part = srID;

      truth.push_back(truePart);
      truthOverlap.push_back(eDeposit / totalEDeposit);
    }


  }

  //------------------------------------------------------------------------------


  double CAFMaker::GetWallDistance(recob::SpacePoint const& sp) const{
    //Get the position of the space point in world coordinates
    double x = sp.XYZ()[0];
    double y = sp.XYZ()[1];
    double z = sp.XYZ()[2];

    //Get the distance to the wall
    double dist = std::numeric_limits<double>::max();

    dist = std::min(dist, x - fActiveBounds[0]);
    dist = std::min(dist, fActiveBounds[1] - x);
    dist = std::min(dist, y - fActiveBounds[2]);
    dist = std::min(dist, fActiveBounds[3] - y);
    dist = std::min(dist, z - fActiveBounds[4]);
    dist = std::min(dist, fActiveBounds[5] - z);

    return dist;
  }

  //------------------------------------------------------------------------------

  double CAFMaker::GetWallDistance(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const{
    double minDist = std::numeric_limits<double>::max();
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfp, evt, fSpacePointLabel);

    if(spacePoints.empty()){
      mf::LogWarning("CAFMaker") << "No space points found with label '" << fSpacePointLabel << "'";
      return minDist;
    }

    //Returning the minimum distance to the wall
    
    for(auto const& sp: spacePoints){
      double dist = GetWallDistance(*sp);
      if(dist < minDist){
        minDist = dist;
      }
    }
    return minDist;
  }

  //------------------------------------------------------------------------------


  void CAFMaker::ComputeActiveBounds(){
    double minx = 99999;
    double maxx = -99999;
    double miny = 99999;
    double maxy = -99999;
    double minz = 99999;
    double maxz = -99999;
  
    fActiveBounds = {minx, maxx, miny, maxy, minz, maxz};
  
     for (geo::TPCGeo const& TPC: fGeom->Iterate<geo::TPCGeo>()) {
      // get center in world coordinates
      auto const center = TPC.GetCenter();
      double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
  
      if( center.X() - tpcDim[0] < fActiveBounds[0] ) fActiveBounds[0] = center.X() - tpcDim[0];
      if( center.X() + tpcDim[0] > fActiveBounds[1] ) fActiveBounds[1] = center.X() + tpcDim[0];
      if( center.Y() - tpcDim[1] < fActiveBounds[2] ) fActiveBounds[2] = center.Y() - tpcDim[1];
      if( center.Y() + tpcDim[1] > fActiveBounds[3] ) fActiveBounds[3] = center.Y() + tpcDim[1];
      if( center.Z() - tpcDim[2] < fActiveBounds[4] ) fActiveBounds[4] = center.Z() - tpcDim[2];
      if( center.Z() + tpcDim[2] > fActiveBounds[5] ) fActiveBounds[5] = center.Z() + tpcDim[2];
    } // for all TPC
  
    //Note that on the y axis there is an extra on non-instrumented 8cm on each side for the HD workspace geom. Not tweaking it here as this would be too hacky
  }

  //------------------------------------------------------------------------------

  double CAFMaker::GetSingleHitsEnergy(art::Event const& evt, const art::Ptr<recob::Slice> &slicePtr, int plane) const{
    std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaSliceUtils::GetHits(slicePtr, evt, fPandoraLabel);

    std::vector<art::Ptr<recob::Hit>> collection_plane_hits = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, plane);

    const art::FindManyP<recob::SpacePoint> sp_assoc(collection_plane_hits, evt, fSpacePointLabel);

    if (!sp_assoc.isValid()) {
      mf::LogWarning("CAFMaker") << "No space points found with label '" << fSpacePointLabel << "'";
      return 0;
    }

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

    double charge = 0;

    for (uint i = 0; i < collection_plane_hits.size(); i++){
      std::vector<art::Ptr<recob::SpacePoint>> matching_sps = sp_assoc.at(i);
      if(!matching_sps.empty()){ // If there are some matched spacepoints, this means that the hit is associated to some PFP and we don't want it here
        continue;
      }
      // Get the energy deposited in the hit
      charge += dune_ana::DUNEAnaHitUtils::LifetimeCorrection(clockData, detProp, collection_plane_hits[i])*collection_plane_hits[i]->Integral();
    }

    return ChargeToEnergyGeV(charge, 2);

  }

  //------------------------------------------------------------------------------

  std::map<int, std::vector<const sim::IDE*>> CAFMaker::slice_IDEs(
      std::vector<const sim::IDE*> ides,
      double the_z0, double the_pitch, double true_endZ) const {

    std::map<int, std::vector<const sim::IDE*>> results;

    for (size_t i = 0; i < ides.size(); ++i) {
      int slice_num = std::floor(
          (ides[i]->z - (the_z0 - the_pitch/2.)) / the_pitch);
      results[slice_num].push_back(ides[i]);
    }

    return results;
  }

  //------------------------------------------------------------------------------

  double CAFMaker::ComputeTrueInteractingEnergy(
      const art::Event& evt,
      const simb::MCParticle* true_beam_particle) const {

    const simb::MCTrajectory& true_beam_trajectory = true_beam_particle->Trajectory();
    double true_beam_mass   = true_beam_particle->Mass() * 1.e3;
    double true_beam_startP = true_beam_particle->P();
    int    true_beam_ID     = true_beam_particle->TrackId();
    double true_beam_endZ   = true_beam_particle->EndZ();

    constexpr geo::PlaneID planeID{0, 1, 2};
    double fZ0    = fWireReadout->Wire(geo::WireID(planeID, 0)).GetCenter().Z();
    double fPitch = fWireReadout->Plane(planeID).WirePitch();

    double init_KE = std::sqrt(1.e6 * true_beam_startP * true_beam_startP +
                               true_beam_mass * true_beam_mass) - true_beam_mass;

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    auto view2_IDEs = bt_serv->TrackIdToSimIDEs_Ps(true_beam_ID, geo::View_t(2));

    std::sort(view2_IDEs.begin(), view2_IDEs.end(),
              [](const sim::IDE* i1, const sim::IDE* i2){ return i1->z < i2->z; });

    size_t remove_index = 0;
    bool   do_remove    = false;
    if (view2_IDEs.size()) {
      for (size_t i = 1; i < view2_IDEs.size() - 1; ++i) {
        const sim::IDE* prev_IDE = view2_IDEs[i-1];
        const sim::IDE* this_IDE = view2_IDEs[i];
        if (this_IDE->trackID < 0 && (this_IDE->z - prev_IDE->z) > 5) {
          remove_index = i;
          do_remove    = true;
          break;
        }
      }
    }
    if (do_remove)
      view2_IDEs.erase(view2_IDEs.begin() + remove_index, view2_IDEs.end());

    auto sliced_ides = slice_IDEs(view2_IDEs, fZ0, fPitch, true_beam_endZ);
    if (sliced_ides.size()) {
      auto  first_slice = sliced_ides.begin();
      auto  theIDEs     = first_slice->second;
      if (theIDEs.empty()) return -1.;

      double ide_z = theIDEs[0]->z;
      for (size_t i = 1; i < true_beam_trajectory.size(); ++i) {
        double z0 = true_beam_trajectory.Z(i-1);
        double z1 = true_beam_trajectory.Z(i);
        if (z0 < ide_z && z1 > ide_z) {
          init_KE = 1.e3 * true_beam_trajectory.E(i-1) - true_beam_mass;
          break;
        }
      }

      double true_beam_interactingEnergy = init_KE;
      for (auto it = sliced_ides.begin(); it != sliced_ides.end(); ++it) {
        auto sliceIDEs = it->second;
        if (sliceIDEs.size()) {
          double deltaE = 0.;
          for (size_t i = 0; i < sliceIDEs.size(); ++i)
            deltaE += sliceIDEs[i]->energy;
          true_beam_interactingEnergy -= deltaE;
        }
      }
      return true_beam_interactingEnergy;
    } else {
      return -1.;
    }
  }

  //------------------------------------------------------------------------------

  void CAFMaker::ComputeRecoInteractingEnergy(
      const art::Event& evt,
      const recob::Track* thisTrack,
      detinfo::DetectorClocksData const& clockData,
      caf::SRInteractionBranch& interaction,
      caf::SRTruthBranch& truthBranch,
      const caf::SRBeamInstrumentation& beamInst) const {

    // Retrieve calorimetry — direct art association, no protoduneana dependency
    auto tracksHandle = evt.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
    art::FindManyP<anab::Calorimetry> fmCalo(tracksHandle, evt, fCalorimetryLabelSCE);
    std::vector<anab::Calorimetry> calo;
    if (fmCalo.isValid()) {
      for (auto caloPtr : fmCalo.at(thisTrack->ID()))
        calo.push_back(*caloPtr);
    }

    bool found_calo = false;
    size_t index = 0;
    for (index = 0; index < calo.size(); ++index) {
      auto this_plane = calo[index].PlaneID().Plane;
      if ((this_plane == 2 && !fVPlaneAsCollector) ||
          (this_plane == 1 &&  fVPlaneAsCollector)) {
        found_calo = true;
        break;
      }
    }
    if (!found_calo) return;

    auto calo_dEdX  = calo[index].dEdx();
    auto calo_range = calo[index].ResidualRange();
    auto TpIndices  = calo[index].TpIndices();
    auto theXYZPoints = calo[index].XYZ();

    auto allHits = evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel);

    struct calo_point {
      calo_point() = default;
      calo_point(size_t w, double in_tick, double p, double dqdx, double dedx, double dq,
                 double cali_dqdx, double cali_dedx, double r, size_t idx, double wire_z_in,
                 int t, double efield, double input_x, double input_y, double input_z)
          : wire(w), tick(in_tick), pitch(p), dQdX(dqdx), dEdX(dedx), dQ(dq),
            calibrated_dQdX(cali_dqdx), calibrated_dEdX(cali_dedx),
            res_range(r), hit_index(idx), wire_z(wire_z_in), tpc(t),
            EField(efield), x(input_x), y(input_y), z(input_z) {}
      size_t wire; double tick, pitch, dQdX, dEdX, dQ;
      double calibrated_dQdX, calibrated_dEdX, res_range;
      size_t hit_index; double wire_z; int tpc; double EField, x, y, z;
    };

    std::vector<calo_point> reco_beam_calo_points;
    for (size_t i = 0; i < calo_dEdX.size(); ++i) {
      const recob::Hit& theHit = (*allHits)[TpIndices[i]];
      reco_beam_calo_points.push_back(calo_point(
        theHit.WireID().Wire, theHit.PeakTime(),
        calo[index].TrkPitchVec()[i],
        calo[index].dQdx()[i], calo_dEdX[i],
        theHit.Integral(),
        calo[index].dQdx()[i], calo_dEdX[i],
        calo_range[i], TpIndices[i],
        0., 0., 0.,
        theXYZPoints[i].X(), theXYZPoints[i].Y(), theXYZPoints[i].Z()
      ));
    }

    std::sort(reco_beam_calo_points.begin(), reco_beam_calo_points.end(),
              [](const calo_point& a, const calo_point& b){ return a.z < b.z; });

    // Beam instrumentation momentum, already retrieved and stored by FillBeamInfo. Unset (NaN)
    // when no beam event was found, which stands in for the 0. the local retrieval used to yield.
    const double beam_inst_P = std::isfinite(beamInst.P) ? beamInst.P : 0.;

    double init_KE = 0.;
    if (evt.isRealData() || fMCHasBI) {
      double mass = 139.57;
      init_KE = std::sqrt(1.e6 * beam_inst_P * beam_inst_P + mass * mass) - mass;
    } else {
      if (truthBranch.nu.empty()) {
        mf::LogWarning("CAFMaker") << "ComputeRecoInteractingEnergy: truthBranch.nu is empty.";
        return;
      }
      auto& nu = truthBranch.nu[0];
      if (nu.prim.empty()) {
        mf::LogWarning("CAFMaker") << "ComputeRecoInteractingEnergy: nu.prim is empty.";
        return;
      }
      auto* pdef = genie::PDGLibrary::Instance()->Find(nu.pdg);
      if (!pdef) {
        mf::LogWarning("CAFMaker") << "ComputeRecoInteractingEnergy: PDGLibrary returned nullptr for PDG " << nu.pdg << ".";
        return;
      }
      double true_beam_mass = pdef->Mass();
      init_KE = nu.prim[0].p.E - true_beam_mass;
      std::cout << "Not real data, the true mass is " << true_beam_mass << " MeV/c^2 and the true energy is " << nu.prim[0].p.E << " MeV" << std::endl;
    }
    
    std::cout << "Initial KE for reco beam: " << init_KE << " MeV" << std::endl;
    if (reco_beam_calo_points.size() > 0) {
      std::vector<double> reco_beam_incidentEnergies;
      reco_beam_incidentEnergies.push_back(init_KE);
      for (size_t i = 0; i < reco_beam_calo_points.size() - 1; ++i) {
        if (reco_beam_calo_points[i].dEdX < 0.) continue;
        double this_energy = reco_beam_incidentEnergies.back() - (reco_beam_calo_points[i].dEdX * reco_beam_calo_points[i].pitch);
        reco_beam_incidentEnergies.push_back(this_energy);
      }
      if (!reco_beam_incidentEnergies.empty())
        interaction.reco_beam_interactingEnergy = reco_beam_incidentEnergies.back();
    } else {
      mf::LogWarning("CAFMaker") << "ComputeRecoInteractingEnergy: no calorimetry points found.";
    }
  }

  //------------------------------------------------------------------------------


  void CAFMaker::beginSubRun(const art::SubRun& sr)
  {
    art::Handle<sumdata::POTSummary> pots = sr.getHandle<sumdata::POTSummary>(fPOTSummaryLabel);
    if( pots ) fMetaPOT += pots->totpot;

    fMetaSubRun = sr.id().subRun();
    fMetaRun = sr.id().run();

  }

  //------------------------------------------------------------------------------

  void CAFMaker::FillMetaInfo(caf::SRDetectorMeta &meta, const art::Event &evt) const
  {
    meta.enabled = true;
    meta.run = evt.id().run();
    meta.subrun = evt.id().subRun();
    meta.event = evt.id().event();
    meta.subevt = 0; //Hardcoded to 0, only makes sense in ND where multiple interactions can occur in the same event

    //Nothing is filled about the trigger for the moment
  }

  //------------------------------------------------------------------------------

  // The three functions below are ports of protoana::ProtoDUNEBeamlineUtils
  // (protoduneana v10_20_09d00). They are reproduced rather than called so that CAFMaker keeps no
  // dependency on protoduneana; every input they touch lives on beam::ProtoDUNEBeamEvent, which
  // dunecore::DuneObj already provides.

  bool CAFMaker::IsGoodBeamlineTrigger(const beam::ProtoDUNEBeamEvent &beamEvent)
  {
    // Timing trigger 12 is the beam trigger; CheckIsMatched() says the beamline record was
    // matched to this DAQ trigger.
    return (beamEvent.GetTimingTrigger() == 12 && beamEvent.CheckIsMatched());
  }

  //------------------------------------------------------------------------------

  bool CAFMaker::HasPerfectBeamMomentum(const beam::ProtoDUNEBeamEvent &beamEvent)
  {
    // One active fiber per momentum monitor means an unambiguous momentum. Fibers flagged in the
    // monitor's glitch mask are dropped first, so this is stricter than the raw fiber counts.
    auto countGood = [&beamEvent](const std::string &monitor) {
      const std::vector<short> &fibers = beamEvent.GetActiveFibers(monitor);
      const std::array<short, 192> &glitch_mask = beamEvent.GetFBM(monitor).glitch_mask;
      int n = 0;
      for (size_t i = 0; i < fibers.size(); ++i) {
        if (!glitch_mask[fibers[i]]) ++n;
      }
      return n;
    };

    return (countGood("XBPF022697") == 1 &&
            countGood("XBPF022701") == 1 &&
            countGood("XBPF022702") == 1);
  }

  //------------------------------------------------------------------------------

  std::vector<int> CAFMaker::GetBeamPDGCandidates(const beam::ProtoDUNEBeamEvent &beamEvent,
                                                  double nominal_momentum) const
  {
    std::vector<int> pdgs;

    // The cuts are only defined at the momenta the beam line was calibrated at.
    const std::vector<double> valid_momenta = {1., 2., 3., 6., 7.};
    if (std::find(valid_momenta.begin(), valid_momenta.end(), nominal_momentum) == valid_momenta.end()) {
      mf::LogWarning("CAFMaker") << "Beam PID: reference momentum " << nominal_momentum
                                 << " GeV/c is not one of 1, 2, 3, 6, 7; no PID assigned.";
      return pdgs;
    }

    // Naming follows the reference: CKov0 is the high-pressure counter, CKov1 the low-pressure one.
    const int high_pressure_status = beamEvent.GetCKov0Status();
    const int low_pressure_status  = beamEvent.GetCKov1Status();

    if (nominal_momentum == 1. || nominal_momentum == 2.) {
      // At 1 and 2 GeV/c the separation is driven by time of flight, with the low-pressure
      // Cherenkov only tagging electrons.
      if (beamEvent.GetTOFChan() == -1) return pdgs;   // no valid TOF
      if (low_pressure_status == -1)    return pdgs;   // no valid Cherenkov

      const double tof = beamEvent.GetTOF();

      // The CERN-calibrated cuts differ between the two momenta; the pre-calibration fallback
      // uses a single boundary.
      const double e_cut    = 105.;   // same at both momenta
      const double mip_cut  = (nominal_momentum == 1. ? 110. : 103.);
      const double p_cut_hi = 160.;
      const double old_cut  = (nominal_momentum == 1. ? 170. : 160.);

      if (((fUseCERNCalibSelection && tof < e_cut) || (!fUseCERNCalibSelection && tof < old_cut))
          && low_pressure_status == 1) {
        pdgs.push_back(11);
      }
      else if (((fUseCERNCalibSelection && tof < mip_cut) || (!fUseCERNCalibSelection && tof < old_cut))
               && low_pressure_status == 0) {
        pdgs.push_back(13);
        pdgs.push_back(211);
      }
      else if (((fUseCERNCalibSelection && tof > mip_cut && tof < p_cut_hi)
                || (!fUseCERNCalibSelection && tof > old_cut))
               && low_pressure_status == 0) {
        pdgs.push_back(2212);
      }
    }
    else if (nominal_momentum == 3.) {
      // At 3 GeV/c both Cherenkovs separate, and TOF is no longer used.
      if (high_pressure_status == -1 || low_pressure_status == -1) return pdgs;

      if (low_pressure_status == 1 && high_pressure_status == 1) {
        pdgs.push_back(11);
      }
      else if (low_pressure_status == 0 && high_pressure_status == 1) {
        pdgs.push_back(13);
        pdgs.push_back(211);
      }
      else { // 0, 0
        pdgs.push_back(321);
        pdgs.push_back(2212);
      }
    }
    else { // 6 or 7 GeV/c
      if (high_pressure_status == -1 || low_pressure_status == -1) return pdgs;

      if (low_pressure_status == 1 && high_pressure_status == 1) {
        pdgs.push_back(11);
        pdgs.push_back(13);
        pdgs.push_back(211);
      }
      else if (low_pressure_status == 0 && high_pressure_status == 1) {
        pdgs.push_back(321);
      }
      else { // 0, 0
        pdgs.push_back(2212);
      }
    }

    return pdgs;
  }

  //------------------------------------------------------------------------------

  void CAFMaker::FillBeamInfo(caf::SRBeamBranch &beam, const art::Event &evt) const
  {
    beam.ismc = !evt.isRealData();

    // The NuMI fields of SRBeamBranch (toroids, horn current, batch positions) have no
    // test-beam analogue and stay unfilled. Everything below is the H4-VLE beam line, read
    // straight out of beam::ProtoDUNEBeamEvent.
    //
    // No protoana::ProtoDUNEBeamlineUtils dependency: the selections that would need it are
    // ported above (IsGoodBeamlineTrigger, HasPerfectBeamMomentum, GetBeamPDGCandidates).
    caf::SRBeamInstrumentation &inst = beam.inst;

    std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
    if (evt.isRealData()) {
      auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamModuleLabel);
      if (beamHandle.isValid()) art::fill_ptr_vector(beamVec, beamHandle);
    }
    else {
      try {
        auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("generator");
        if (beamHandle.isValid()) art::fill_ptr_vector(beamVec, beamHandle);
      }
      catch (const cet::exception &) {
        mf::LogWarning("CAFMaker") << "BeamEvent generator object not found; beam instrumentation left empty.";
      }
    }

    if (beamVec.empty()) return;

    const beam::ProtoDUNEBeamEvent &beamEvent = *(beamVec.at(0));

    inst.trigger = beamEvent.GetTimingTrigger();
    // Simulation has no DAQ trigger to match against, so the beamline-trigger requirement only
    // makes sense on data.
    inst.valid   = evt.isRealData() ? IsGoodBeamlineTrigger(beamEvent) : true;

    const std::vector<double> momenta = beamEvent.GetRecoBeamMomenta();
    inst.nmomenta = static_cast<int>(momenta.size());
    inst.momenta.reserve(momenta.size());
    for (double p : momenta) inst.momenta.push_back(static_cast<float>(p));
    if (!momenta.empty()) {
      inst.P_raw = static_cast<float>(momenta[0]);
      // fBeamInstPFix puts simulated beam instrumentation on the data momentum scale. inst.P is
      // the corrected value, which is what ComputeRecoInteractingEnergy consumes; inst.P_raw is
      // kept so the correction can be undone or varied offline.
      inst.P = static_cast<float>(evt.isRealData() ? momenta[0] : momenta[0] * fBeamInstPFix);
    }

    // The singular TOF is the one the PID cuts act on; the vectors are all the candidates.
    inst.TOF      = static_cast<float>(beamEvent.GetTOF());
    inst.TOF_chan = beamEvent.GetTOFChan();

    const std::vector<double> tofs  = beamEvent.GetTOFs();
    const std::vector<int>    chans = beamEvent.GetTOFChans();
    inst.TOFs.reserve(tofs.size());
    inst.TOF_chans.reserve(tofs.size());
    for (size_t i = 0; i < tofs.size(); ++i) {
      inst.TOFs.push_back(static_cast<float>(tofs[i]));
      if (i < chans.size()) inst.TOF_chans.push_back(chans[i]);
    }

    inst.C0          = beamEvent.GetCKov0Status();
    inst.C1          = beamEvent.GetCKov1Status();
    inst.C0_pressure = beamEvent.GetCKov0Pressure();
    inst.C1_pressure = beamEvent.GetCKov1Pressure();

    // Position/direction of the beamline track projected towards the TPC face, leading track only.
    const auto &beamTracks = beamEvent.GetBeamTracks();
    inst.ntracks = static_cast<int>(beamTracks.size());
    if (!beamTracks.empty()) {
      const auto &traj = beamTracks[0].Trajectory();
      inst.pos = caf::SRVector3D(traj.End().X(), traj.End().Y(), traj.End().Z());
      inst.dir = caf::SRVector3D(traj.EndDirection().X(),
                                 traj.EndDirection().Y(),
                                 traj.EndDirection().Z());
    }

    inst.nfibers_p1 = static_cast<int>(beamEvent.GetActiveFibers("XBPF022697").size());
    inst.nfibers_p2 = static_cast<int>(beamEvent.GetActiveFibers("XBPF022701").size());
    inst.nfibers_p3 = static_cast<int>(beamEvent.GetActiveFibers("XBPF022702").size());
    inst.perfect_momentum = HasPerfectBeamMomentum(beamEvent);

    inst.PDG_candidates = GetBeamPDGCandidates(beamEvent, fBeamPIDMomentum);
  }

  //------------------------------------------------------------------------------

  void CAFMaker::FillCVNInfo(caf::SRCVNScoreBranch &cvnBranch, const art::Event &evt) const
  {
    art::Handle<std::vector<cvn::Result>> cvnin = evt.getHandle<std::vector<cvn::Result>>(fCVNLabel);

    if( !cvnin.failedToGet() && !cvnin->empty()) {
      if(fIsAtmoCVN){ //Hotfix to take care of the fact that the CVN for atmospherics is storing results in a weird way...
        const std::vector<std::vector<float>> &scores = (*cvnin)[0].fOutput;
        cvnBranch.nc = scores[0][2];
        cvnBranch.nue = scores[0][1];
        cvnBranch.numu = scores[0][0];
      }

      else{ //Normal code
        cvnBranch.isnubar = (*cvnin)[0].GetIsAntineutrinoProbability();
        cvnBranch.nue = (*cvnin)[0].GetNueProbability();
        cvnBranch.numu = (*cvnin)[0].GetNumuProbability();
        cvnBranch.nutau = (*cvnin)[0].GetNutauProbability();
        cvnBranch.nc = (*cvnin)[0].GetNCProbability();

        cvnBranch.protons0 = (*cvnin)[0].Get0protonsProbability();
        cvnBranch.protons1 = (*cvnin)[0].Get1protonsProbability();
        cvnBranch.protons2 = (*cvnin)[0].Get2protonsProbability();
        cvnBranch.protonsN = (*cvnin)[0].GetNprotonsProbability();

        cvnBranch.chgpi0 = (*cvnin)[0].Get0pionsProbability();
        cvnBranch.chgpi1 = (*cvnin)[0].Get1pionsProbability();
        cvnBranch.chgpi2 = (*cvnin)[0].Get2pionsProbability();
        cvnBranch.chgpiN = (*cvnin)[0].GetNpionsProbability();

        cvnBranch.pizero0 = (*cvnin)[0].Get0pizerosProbability();
        cvnBranch.pizero1 = (*cvnin)[0].Get1pizerosProbability();
        cvnBranch.pizero2 = (*cvnin)[0].Get2pizerosProbability();
        cvnBranch.pizeroN = (*cvnin)[0].GetNpizerosProbability();

        cvnBranch.neutron0 = (*cvnin)[0].Get0neutronsProbability();
        cvnBranch.neutron1 = (*cvnin)[0].Get1neutronsProbability();
        cvnBranch.neutron2 = (*cvnin)[0].Get2neutronsProbability();
        cvnBranch.neutronN = (*cvnin)[0].GetNneutronsProbability();
      }
      
    }
  }

  //------------------------------------------------------------------------------

  void CAFMaker::FillDirectionInfo(caf::SRDirectionBranch &dirBranch, const art::Event &evt) const
  {
    art::Handle<dune::AngularRecoOutput> dirReco = evt.getHandle<dune::AngularRecoOutput>(fDirectionRecoLabelNumu);
    if(!dirReco.failedToGet()){
      dirBranch.lngtrk.SetX(dirReco->fRecoDirection.X());
      dirBranch.lngtrk.SetY(dirReco->fRecoDirection.Y());
      dirBranch.lngtrk.SetZ(dirReco->fRecoDirection.Z());
    }
    else{
      mf::LogWarning("CAFMaker") << "No AngularRecoOutput found with label '" << fDirectionRecoLabelNumu << "'";
    }

    dirReco = evt.getHandle<dune::AngularRecoOutput>(fDirectionRecoLabelNue);
    if(!dirReco.failedToGet()){
      dirBranch.heshw.SetX(dirReco->fRecoDirection.X());
      dirBranch.heshw.SetY(dirReco->fRecoDirection.Y());
      dirBranch.heshw.SetZ(dirReco->fRecoDirection.Z());
    }
    else{
      mf::LogWarning("CAFMaker") << "No AngularRecoOutput found with label '" << fDirectionRecoLabelNue << "'";
    }

    dirReco = evt.getHandle<dune::AngularRecoOutput>(fDirectionRecoLabelCalo);
    if(!dirReco.failedToGet()){
      dirBranch.calo.SetX(dirReco->fRecoDirection.X());
      dirBranch.calo.SetY(dirReco->fRecoDirection.Y());
      dirBranch.calo.SetZ(dirReco->fRecoDirection.Z());
    }
    else{
      mf::LogWarning("CAFMaker") << "No AngularRecoOutput found with label '" << fDirectionRecoLabelCalo << "'";
    }

  }

  //------------------------------------------------------------------------------

  void CAFMaker::FillEnergyInfo(caf::SRNeutrinoEnergyBranch &ErecBranch, const art::Event &evt) const
  {
    //Filling the reg CNN results
    art::InputTag itag(fRegCNNLabel, "regcnnresult");
    art::Handle<std::vector<cnn::RegCNNResult>> regcnn = evt.getHandle<std::vector<cnn::RegCNNResult>>(itag);
    if(!regcnn.failedToGet() && !regcnn->empty()){
        const std::vector<float>& cnnResults = (*regcnn)[0].fOutput;
        ErecBranch.regcnn = cnnResults[0];
    }
    else{
      mf::LogWarning("CAFMaker") << itag << " does not correspond to a valid RegCNNResult product";
    }

    std::map<std::string, float*> ereco_map = {
      {fEnergyRecoCaloLabel, &(ErecBranch.calo)},
      {fEnergyRecoLepCaloLabel, &(ErecBranch.lep_calo)},
      {fEnergyRecoMuRangeLabel, &(ErecBranch.mu_range)},
      {fEnergyRecoMuMcsLabel, &(ErecBranch.mu_mcs)},
      {fEnergyRecoMuMcsLLHDLabel, &(ErecBranch.mu_mcs_llhd)},
      {fEnergyRecoECaloLabel, &(ErecBranch.e_calo)}
    };

    for(auto [label, record] : ereco_map){
       art::Handle<dune::EnergyRecoOutput> ereco = evt.getHandle<dune::EnergyRecoOutput>(label);
       if(ereco.failedToGet()){
        mf::LogWarning("CAFMaker") << label << " does not correspond to a valid EnergyRecoOutput product";
       }
       else{
        *record = ereco->fNuLorentzVector.E();

        if(label == fEnergyRecoECaloLabel){
          //Adding the had information
          ErecBranch.e_had = ereco->fHadLorentzVector.E();
        }
        else if(label == fEnergyRecoMuRangeLabel){
          //Adding the had information
          ErecBranch.mu_had = ereco->fHadLorentzVector.E();
        }
       }
    }

  }

  //------------------------------------------------------------------------------

  void CAFMaker::FillRecoParticlesInfo(caf::SRRecoParticlesBranch &recoParticlesBranch, caf::SRFD &fdBranch, const art::Event &evt, const art::Ptr<recob::Slice> &slicePtr, const art::FindManyP<recob::PFParticle> &sliceToPFP) const
  {
    //Doing quite a lot of things here related to saving the reco particles
    //Will try to be pedagogical in the comments

    //Getting Ar density in g/cm3
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
    double lar_density = detProp.Density();

    //Getting all the PFParticles from this slice
    lar_pandora::PFParticleVector particleVector = sliceToPFP.at(slicePtr.key());

    // if MVALabel is "", the hitResults pointer will be a null pointer and GetMVAResults
    // will return NaN for all the scores without crashing
    anab::MVAReader<recob::Hit,4> * hitResults;
    if(fMVALabel == ""){
      mf::LogWarning("CAFMaker") << "MVALabel is empty, the MVA scores will not be filled";
      hitResults = nullptr;
    }
    else{
      hitResults = new anab::MVAReader<recob::Hit, 4>(evt, fMVALabel);
    }

    unsigned int nuID = std::numeric_limits<unsigned int>::max();
    for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      if(particle->IsPrimary() && (std::abs(particle->PdgCode()) == 12 || std::abs(particle->PdgCode()) == 14 || std::abs(particle->PdgCode()) == 16)){
        nuID = particle->Self(); //Finding the ID of neutrino's particle
        break;
      }
    }

    //Getting the PID information (only used for the tracks)
    art::Handle<std::vector<recob::Track>> tracks_handle = evt.getHandle<std::vector<recob::Track>>(fTrackLabel);
    if(!tracks_handle.isValid()){
      mf::LogWarning("CAFMaker") << "No Track found with label '" << fTrackLabel << "'";
    }

    //Creating a FindManyP object to link the tracks to the PIDs
    const art::FindManyP<anab::ParticleID> fmPID(tracks_handle, evt, fParticleIDLabel);

    //Creating the FD interaction record where we are going to save the tracks/showers in parallel to the PFPs objects
    caf::SRFDInt fdIxn;

    //Mapping particle idx in the vector to their Pandora ID for later use
    std::map<unsigned int, int> pandoraIDToPFPIdx;

    //Iterating on all the PFParticles to fill the reco particles
    for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      if(particle->Self() == nuID){ //Skipping the neutrino that is not a "real" reco particle
        continue;
      }

      //Creating the particle record for this PFP
      caf::SRRecoParticle particle_record;

      particle_record.primary = (particle->Parent() == nuID); //Is primary if the parent is the neutrino
      //TODO: For now PDG is taken from the PFP only which does not include PID info beyond track/shower. Would be good to have the user specifying a specific PID module that takes care of that.
      particle_record.pdg = particle->PdgCode();
      particle_record.tgtA = 40; //Interaction on Ar40. TODO: Maybe to improve if we want to consider interactions outside the detector active volume
      //TODO: Not filling the momentum information as it requires some specific PID to be made.
      // particle_record.p;

      //Pre-fill the parent and daughter fields with the Pandora IDs. They will be converted to SR indices later.
      particle_record.parent = (particle_record.primary) ? -1 : particle->Parent(); //Setting to -1 if primary as there is no saved record for the neutrino

      for(auto const& daughter_id : particle->Daughters()){
        particle_record.daughters.push_back(daughter_id);
      }
      pandoraIDToPFPIdx[particle->Self()] = fdIxn.npfps;

      FillTruthMatchingAndOverlap(particle, evt, particle_record.truth, particle_record.truthOverlap);

      particle_record.walldist = GetWallDistance(particle, evt); //Getting the distance to the wall for this PFP
      particle_record.contained = (particle_record.walldist < fContainedDistThreshold); //Setting the contained flag based on the distance to the wall

      //Getting the track and shower objects associated to the PFP
      art::Ptr<recob::Track> track;
      if(dune_ana::DUNEAnaPFParticleUtils::HasTrack(particle, evt, fPandoraLabel, fTrackLabel)){ 
        track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(particle, evt, fPandoraLabel, fTrackLabel);
      }
      art::Ptr<recob::Shower> shower;
      if(dune_ana::DUNEAnaPFParticleUtils::HasShower(particle, evt, fPandoraLabel, fShowerLabel)){
        shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(particle, evt, fPandoraLabel, fShowerLabel);
      }
      //Seeing which option Pandora prefers
      bool isTrack = lar_pandora::LArPandoraHelper::IsTrack(particle);
      bool isShower = lar_pandora::LArPandoraHelper::IsShower(particle);
      //For every PFP we create a track and a shower object and save it, independently of the existence of a track object to keep the PFP/Track/Shower parallel indexing
      SRTrack srtrack;
      SRShower srshower;

      //We define Evis at the PFP level as the hits are the same for shower and track
      double Evis = GetVisibleEnergy(particle, evt);

      //This variable will be updated correctly during the recob::Track processing and will be used to fill the reco particle energy method if isTrack.
      caf::PartEMethod trackErecoMethod = caf::PartEMethod::kUnknownMethod;

      if(track.isNonnull()){
        srtrack.start.SetX(track->Start().X());
        srtrack.start.SetY(track->Start().Y());
        srtrack.start.SetZ(track->Start().Z());

        srtrack.end.SetX(track->End().X());
        srtrack.end.SetY(track->End().Y());
        srtrack.end.SetZ(track->End().Z());

        srtrack.dir.SetX(track->StartDirection().X());
        srtrack.dir.SetY(track->StartDirection().Y());
        srtrack.dir.SetZ(track->StartDirection().Z());

        srtrack.enddir.SetX(track->EndDirection().X());
        srtrack.enddir.SetY(track->EndDirection().Y());
        srtrack.enddir.SetZ(track->EndDirection().Z());

        srtrack.Evis = Evis; //Using the visible energy of the PFP
        //srtrack.qual TODO: Not sure we have anything relevant to put on the FD side for this at the moment

        srtrack.len_gcm2 = track->Length() * lar_density; //Length in g/cm2
        srtrack.len_cm = track->Length();

        //TODO: I would prefer to use some unified module that the user can setup and that will decide how to compute the energy rather than making a specific choice here
        //Putting Evis as placeholder to not confuse the user too much
        srtrack.E = Evis;
        trackErecoMethod = caf::PartEMethod::kCalorimetry; //Using the visible energy of the PFP

        //Truth matching already filled at the PFP level, no need to do it again here
        //srtrack.truth
        //srtrack.truthOverlap

        //Filling the PFP score with the PIDA score computed for the associated track
        if(fmPID.isValid()){
          std::vector<art::Ptr<anab::ParticleID>> pid_vec = fmPID.at(track.key());
          if(pid_vec.empty()){
            mf::LogWarning("CAFMaker") << "No ParticleID found for track with key " << track.key();
          }
          else{
            art::Ptr<anab::ParticleID> pid = pid_vec[0];
            const std::vector<anab::sParticleIDAlgScores> pScores = pid->ParticleIDAlgScores();
            for(const anab::sParticleIDAlgScores &pScore : pScores){ //Several scores are saved for different assumptions
              if(pScore.fAssumedPdg == 0){ //PIDA score is when there is no assumed pdf
                particle_record.score = pScore.fValue;
                break;
              }
            }
          }
          
        }
      }

      if(shower.isNonnull()){
        //Filling the shower information
        srshower.start.SetX(shower->ShowerStart().X());
        srshower.start.SetY(shower->ShowerStart().Y());
        srshower.start.SetZ(shower->ShowerStart().Z());

        srshower.direction.SetX(shower->Direction().X());
        srshower.direction.SetY(shower->Direction().Y());
        srshower.direction.SetZ(shower->Direction().Z());
        
        srshower.Evis = Evis; //Using the visible energy of the PFP
        //Truth matching already filled at the PFP level, no need to do it again here
        //srshower.truth
        //srshower.truthOverlap
      }

      if (!isTrack && !isShower){
        mf::LogWarning("CAFMaker") << "PFP with ID " << particle->Self() << " is not associated to either a track or a shower according to Pandora. This particle will be saved without kinematic information.";
      }

      if(isTrack){
        if(track.isNonnull()){ //I hope this condition is always fullfilled is the particle is tagged at track, but who knows...
          particle_record.start = SRVector3D(track->Start().X(), track->Start().Y(), track->Start().Z());
          particle_record.end = SRVector3D(track->End().X(), track->End().Y(), track->End().Z());
          particle_record.E = srtrack.E;
          particle_record.E_method = trackErecoMethod;
        }
        particle_record.origRecoObjType = caf::RecoObjType::kTrack;
      }
      else{
        if(shower.isNonnull()){ //I hope this condition is always fullfilled is the particle is tagged at shower, but who knows...
          particle_record.start = SRVector3D(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
          //Only filling the start, no defined end for a shower

          particle_record.E = srshower.Evis; //Using the visible energy of the PFP
          particle_record.E_method = caf::PartEMethod::kCalorimetry;
          
        }
        particle_record.origRecoObjType = caf::RecoObjType::kShower;
      }

      //Saving the track record
      fdIxn.tracks.push_back(std::move(srtrack));
      fdIxn.ntracks++;

      //Saving the shower record
      fdIxn.showers.push_back(std::move(srshower));
      fdIxn.nshowers++;

      //Also saving PFP metadata
      caf::SRPFP pfp_metarecord;
      GetMVAResults(pfp_metarecord, particle, evt, hitResults, 2, true ); //TODO -- make configurable
      FillPFPMetadata(pfp_metarecord, particle, evt);
      fdIxn.pfps.push_back(std::move(pfp_metarecord));
      fdIxn.npfps++;
        
      //Saving the particle record for this PFP
      recoParticlesBranch.pandora.push_back(std::move(particle_record));
      recoParticlesBranch.npandora++;

    }


    //Now that all particles are saved, we can convert the parent/daughter fields from Pandora IDs to SR indices
    for(auto &particle : recoParticlesBranch.pandora){
    //Daughters
      // for(auto &daughter_idx : particle.daughters){
      for(size_t i = 0; i < particle.daughters.size(); i++){
        unsigned int daughter_pfpID = particle.daughters[i];
        //Finding the SR index in the map
        if(pandoraIDToPFPIdx.count(daughter_pfpID) == 0){
          mf::LogWarning("CAFMaker") << "No SR index found for daughter PFP ID " << daughter_pfpID;
          particle.daughters[i] = -1; //Setting to -1 to avoid confusion
        }
        else{
          particle.daughters[i] = pandoraIDToPFPIdx.at(daughter_pfpID);
        }
      }

      //Parent
      if(particle.parent == -1) continue; //Skipping if primary
      unsigned int parent_pfpID = particle.parent;
      //Finding the SR index in the map
      if(pandoraIDToPFPIdx.count(parent_pfpID) == 0){
        mf::LogWarning("CAFMaker") << "No SR index found for parent PFP ID " << parent_pfpID;
        particle.parent = -1; //Setting to -1 to avoid confusion
      }
      else{
        particle.parent = pandoraIDToPFPIdx.at(parent_pfpID);
      }
    }

    //Saving the FD interaction record
    fdBranch.pandora.push_back(std::move(fdIxn));
    fdBranch.npandora++;

    //Adding some extra record with all the leftover hits not associated to any particle
    caf::SRRecoParticle single_hits;
    single_hits.primary = false;
    single_hits.pdg = 0; //Not a real particle
    single_hits.tgtA = 40; //Interaction on Ar40.
    single_hits.E = GetSingleHitsEnergy(evt, slicePtr, 2); //Using the collection plane for now
    single_hits.origRecoObjType = caf::RecoObjType::kHitCollection;

    recoParticlesBranch.pandora.push_back(std::move(single_hits));
    recoParticlesBranch.npandora++;

    // particle_record.origRecoObjType

    delete hitResults;
  }


  //------------------------------------------------------------------------------


  double CAFMaker::GetVisibleEnergy(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const
  {
    if(!dune_ana::DUNEAnaPFParticleUtils::HasShower(pfp, evt, fPandoraLabel, fShowerLabel)){
      return 0;
    }
    //Using the shower version of the PFP to compute the visible energy for the particle
    art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, evt, fPandoraLabel, fShowerLabel);
    if(!shower){ //Should always exist in theory but who knows...
      return 0;
    }

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
    // Get the hits on the collection plane
    const std::vector<art::Ptr<recob::Hit> > showerHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(shower,evt,fShowerLabel),2));
    // Compute the charge
    const double showerCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, showerHits));

    return ChargeToEnergyGeV(showerCharge, 2);

  }

  //------------------------------------------------------------------------------

  bool CAFMaker::IsVertexContained(caf::SRVector3D const& vtx) const
  {
    //Checking if the vertex is contained in the fiducial volume
    return (vtx.X() > fActiveBounds[0] + fVertexFiducialVolumeCut[0] && vtx.X() < fActiveBounds[1] - fVertexFiducialVolumeCut[1] &&
            vtx.Y() > fActiveBounds[2] + fVertexFiducialVolumeCut[2] && vtx.Y() < fActiveBounds[3] - fVertexFiducialVolumeCut[3] &&
            vtx.Z() > fActiveBounds[4] + fVertexFiducialVolumeCut[4] && vtx.Z() < fActiveBounds[5] - fVertexFiducialVolumeCut[5]);
  }

  //------------------------------------------------------------------------------

  int CAFMaker::FillGENIERecord(simb::MCTruth const& mctruth, simb::GTruth const& gtruth)
  {
    std::unique_ptr<const genie::EventRecord> record(evgb::RetrieveGHEP(mctruth, gtruth));
    int cur_idx = fGENIETree->GetEntries();
    fEventRecord->Fill(cur_idx, record.get());
    fGENIETree->Fill();

    return cur_idx;
  }


  //------------------------------------------------------------------------------
  
  void CAFMaker::analyze(art::Event const & evt)
  {
    caf::StandardRecord sr;
    caf::StandardRecord* psr = &sr;

    PreLoadMCParticlesInfo(evt);
    

    if(fTree){
      fTree->SetBranchAddress("rec", &psr);
    }

    std::string geoName = fGeom->DetectorName();

    mf::LogInfo("CAFMaker") << "Geo name is: " << geoName;

    SRDetectorMeta *detector;
    SRFD *fdBranch;

    if(geoName.find("dunevd10kt") != std::string::npos){
      detector = &(sr.meta.fd_vd);
      fdBranch = &(sr.fd.vd);
      mf::LogInfo("CAFMaker") << "Assuming the FD VD detector";
    }
    else if (geoName.find("dune10kt") != std::string::npos)
    {
      detector = &(sr.meta.fd_hd);
      fdBranch = &(sr.fd.hd);
      mf::LogInfo("CAFMaker") << "Assuming the FD HD detector";
    }
    else if (geoName.find("protodune") != std::string::npos)
    {
      detector = &(sr.meta.pd_hd);
      fdBranch = &(sr.fd.pd_hd);
      mf::LogInfo("CAFMaker") << "Assuming the PDUNE detector";
    }
    else {
      mf::LogWarning("CAFMaker") << "Didn't detect a know geometry. Defaulting to FD HD!";
      detector = &(sr.meta.fd_hd);
      fdBranch = &(sr.fd.hd);
    }
    



    FillMetaInfo(*detector, evt);

    FillBeamInfo(sr.beam, evt);

    FillTruthInfo(sr.mc, evt);

    FillRecoInfoSliceLoop(sr.common, *fdBranch, sr.mc, evt, sr.beam.inst);

    if(fTree){
      fTree->Fill();
    }

    if(fFlatTree){
      fFlatRecord->Clear();
      fFlatRecord->Fill(sr);
      fFlatTree->Fill();
    }
  }

  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  void CAFMaker::endSubRun(const art::SubRun& sr){
  }

  void CAFMaker::endJob()
  {
    fMetaTree->Fill();

    if(fFlatFile){
      fFlatFile->cd();
      fFlatTree->Write();
      fMetaTree->CloneTree()->Write();
      fGENIETree->CloneTree()->Write();
      fFlatFile->Close();
    }

    delete fEventRecord; //Making this a unique_pointer requires too many circonvolutions because of TTree->Branch requiring a pointer to a pointer

  }

  DEFINE_ART_MODULE(CAFMaker)

} // namespace caf

#endif // CAFMaker_H
