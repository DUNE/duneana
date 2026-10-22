////////////////////////////////////////////////////////////////////////
/// \file    PDSPSelection.h
/// \brief   ProtoDUNE-SP pi+ absorption / charge-exchange selection for the CAFMaker
/// \author  D. Pullia <dariopullia@gmail.com>
///
/// Port of protoduneana's PDSPAnalyzer (the inputs) plus runEventSelection.cc /
/// SelectionDefinitions.h / EventSelection.fcl (the cuts), as of protoduneana
/// v10_20_09d00. Line references in comments ("PDSPA:L", "SD:L", "rES:L",
/// "PTU:L") point into that version; old_analysis_reference/*.md documents them.
///
/// The chain runs twice per event with the same inputs:
///   - legacy:   all Variant switches off, i.e. the old code, quirks included;
///   - improved: the switches set in the Improved fcl table.
/// Both share every line of code, so any difference between them comes from a
/// switch. No protoduneana dependency: the helpers are re-implemented here.
////////////////////////////////////////////////////////////////////////

#ifndef CAFMAKER_PDSPSELECTION_H
#define CAFMAKER_PDSPSELECTION_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "cetlib/filesystem.h"
#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TFile.h"
#include "TMath.h"
#include "TProfile.h"

#include "duneanaobj/StandardRecord/SRBeamInstrumentation.h"
#include "duneanaobj/StandardRecord/SRPDSPSelection.h"

namespace caf {
namespace pdspsel {

  /// Beam-spot constants of the TPC beam cut (old beam_cut_TPC). One set for MC, one for data.
  struct BeamSpot {
    double meanX, meanY, meanZ;    ///< Mean SCE-corrected calo start [cm]
    double sigmaX, sigmaY, sigmaZ; ///< Its spread [cm]
    double thetaX, thetaY, thetaZ; ///< Mean beam direction: angles to the x, y, z axes [deg]
  };

  /// What separates the improved implementation from the legacy one. All false = the old chain.
  struct Variant {
    bool beamPFPByScore      = false; ///< Beam = the tagged primary with the best TestBeamScore (legacy: first primary of the first tagged slice)
    bool michelAtEnd         = false; ///< Michel window at the downstream end of the beam track (legacy: hit farthest from the origin)
    bool requireTwoCaloPts   = false; ///< Fewer than 2 distinct calo points fails the beam cut (legacy: NaN angle passes)
    bool requireShowerEnergy = false; ///< A shower without collection calorimetry cannot tag charge exchange
    bool biQualityInMC       = false; ///< One beam momentum required in MC as well as in data
  };

  class PDSPSelection {
    public:
      PDSPSelection(fhicl::ParameterSet const& pset,
                    std::string const& pfpLabel, std::string const& trackLabel,
                    std::string const& hitLabel, std::string const& mvaLabel,
                    int collectionPlane);

      /// Fill both implementations and the shared MC truth for one event.
      void Fill(art::Event const& evt, caf::SRBeamInstrumentation const& inst,
                caf::SRPDSPSelectionBranch& out) const;

    private:
      // The old chain's sentinels, kept so the stored values match the old ntuples.
      static constexpr double kMissing = -999.;
      static constexpr double kEmptydEdx = -9999.;
      static constexpr double kPionMass = 139.57; // [MeV], hardcoded in the old code (PDSPA:L3330)

      /// Products and associations, looked up once per event and shared by both implementations.
      struct EventData;

      struct BeamChoice {
        int   pfpIndex = -1;      ///< Index into the PFParticle collection
        int   slice    = -1;
        int   nTaggedSlices = 0;
        bool  tagged   = false;
        float score    = std::numeric_limits<float>::signaling_NaN();
      };

      void FillTruth(art::Event const& evt, caf::SRPDSPSelectionBranch& out) const;
      void RunVariant(art::Event const& evt, EventData const& ed,
                      caf::SRBeamInstrumentation const& inst, int truePDG,
                      Variant const& var, caf::SRPDSPSelection& sel) const;

      BeamChoice ChooseBeamPFP(EventData const& ed, bool byScore) const;
      std::vector<anab::Calorimetry> TrackCalo(EventData const& ed, bool forced, recob::Track const& track) const;
      std::pair<double, int> MichelScore(EventData const& ed, recob::Track const& track,
                                         geo::Point_t const& ref, caf::SRVector3D& vtx) const;
      static std::pair<double, int> Chi2PID(std::vector<float> const& dedx,
                                            std::vector<float> const& range, TProfile* profile);
      double TruncatedMean(std::vector<float> const& dedx) const;
      static bool ShowerDistEnergyCheck(double dist, double energy);

      std::string fPFPLabel, fTrackLabel, fHitLabel, fMVALabel;
      std::string fForcedTrackLabel, fForcedShowerLabel;
      std::string fBeamCaloLabel, fDaughterCaloLabel;
      std::string fGeneratorLabel;
      int  fCollectionPlane;

      bool   fCheckCalo;
      double fEndZLow, fEndZHigh, fThreshold;
      bool   fCexNPi0, fSignalPastFV;
      double fTrackScoreCut, fLimit, fChi2Cut, fdEdXLow, fdEdXMed, fdEdXHigh;
      bool   fDoMichel;
      double fMichelCut;
      bool   fDoNTracks;
      std::vector<int> fExcludeRuns;
      bool   fDoAngle;
      double fXYZCut, fCosCut;
      BeamSpot fMCSpot, fDataSpot;

      Variant fImproved;

      std::unique_ptr<TFile> fTemplateFile;
      TProfile* fProtonTemplate = nullptr; ///< Owned by fTemplateFile
  };

  //------------------------------------------------------------------------------
  // Implementation
  //------------------------------------------------------------------------------

  namespace detail {
    /// Build an association finder, or nullptr when the association product is absent.
    template <typename F, typename... Args>
    std::unique_ptr<F> TryMake(Args&&... args) {
      try {
        auto p = std::make_unique<F>(std::forward<Args>(args)...);
        if (!p->isValid()) return nullptr;
        return p;
      }
      catch (cet::exception const&) {
        return nullptr;
      }
    }

    inline BeamSpot ReadSpot(fhicl::ParameterSet const& p, BeamSpot const& def) {
      return BeamSpot{p.get<double>("MeanX", def.meanX),   p.get<double>("MeanY", def.meanY),   p.get<double>("MeanZ", def.meanZ),
                      p.get<double>("SigmaX", def.sigmaX), p.get<double>("SigmaY", def.sigmaY), p.get<double>("SigmaZ", def.sigmaZ),
                      p.get<double>("ThetaX", def.thetaX), p.get<double>("ThetaY", def.thetaY), p.get<double>("ThetaZ", def.thetaZ)};
    }
  } // detail

  struct PDSPSelection::EventData {
    art::Handle<std::vector<recob::PFParticle>> pfps;
    art::Handle<std::vector<recob::Track>>      tracks;        // pandoraTrack
    art::Handle<std::vector<recob::Track>>      forcedTracks;  // pandora2Track
    art::Handle<std::vector<recob::Shower>>     forcedShowers; // pandora2Shower
    art::Handle<std::vector<recob::Cluster>>    clusters;
    art::Handle<std::vector<recob::Hit>>        hits;
    std::vector<art::Ptr<recob::Hit>>           hitPtrs;

    std::unique_ptr<art::FindManyP<larpandoraobj::PFParticleMetadata>> pfpMeta;
    std::unique_ptr<art::FindOneP<recob::Slice>>      pfpSlice;
    std::unique_ptr<art::FindManyP<recob::Track>>     pfpTrack;
    std::unique_ptr<art::FindManyP<recob::Track>>     pfpForcedTrack;
    std::unique_ptr<art::FindManyP<recob::Shower>>    pfpForcedShower;
    std::unique_ptr<art::FindManyP<recob::Cluster>>   pfpClusters;
    std::unique_ptr<art::FindManyP<recob::Hit>>       clusterHits;
    std::unique_ptr<art::FindManyP<anab::Calorimetry>> beamCalo;
    std::unique_ptr<art::FindManyP<anab::Calorimetry>> forcedTrackCalo;
    std::unique_ptr<art::FindManyP<anab::Calorimetry>> forcedShowerCalo;
    std::unique_ptr<art::FindManyP<recob::Hit>>       trackHits;
    std::unique_ptr<art::FindManyP<recob::Hit, recob::TrackHitMeta>> trackHitMetas;
    std::unique_ptr<art::FindManyP<recob::Track>>     hitTracks;
    std::unique_ptr<anab::MVAReader<recob::Hit, 4>>   cnn;

    /// Pandora metadata of one PFParticle; empty if it has none.
    std::map<std::string, float> Meta(size_t self) const {
      if (!pfpMeta) return {};
      auto const& v = pfpMeta->at(self);
      if (v.empty() || v.at(0).isNull()) return {};
      return v.at(0)->GetPropertiesMap();
    }
  };

  inline PDSPSelection::PDSPSelection(fhicl::ParameterSet const& pset,
                                      std::string const& pfpLabel, std::string const& trackLabel,
                                      std::string const& hitLabel, std::string const& mvaLabel,
                                      int collectionPlane)
    : fPFPLabel(pfpLabel), fTrackLabel(trackLabel), fHitLabel(hitLabel), fMVALabel(mvaLabel),
      fForcedTrackLabel(pset.get<std::string>("ForcedTrackLabel", "pandora2Track")),
      fForcedShowerLabel(pset.get<std::string>("ForcedShowerLabel", "pandora2Shower")),
      fBeamCaloLabel(pset.get<std::string>("BeamCaloLabel", "pandoracali")),
      fDaughterCaloLabel(pset.get<std::string>("DaughterCaloLabel", "pandora2cali")),
      fGeneratorLabel(pset.get<std::string>("BeamGeneratorLabel", "generator")),
      fCollectionPlane(collectionPlane),
      // Defaults are the values of protoduneana's EventSelection.fcl.
      fCheckCalo(pset.get<bool>("CheckCalo", true)),
      fEndZLow(pset.get<double>("EndZLow", -0.49375)),
      fEndZHigh(pset.get<double>("EndZHigh", 222.10561)),
      fThreshold(pset.get<double>("Threshold", 0.150)),
      fCexNPi0(pset.get<bool>("CexNPi0", true)),
      fSignalPastFV(pset.get<bool>("SignalPastFV", true)),
      fTrackScoreCut(pset.get<double>("TrackScoreCut", 0.3)),
      fLimit(pset.get<double>("Limit", 0.16)),
      fChi2Cut(pset.get<double>("Chi2Cut", 70.)),
      fdEdXLow(pset.get<double>("dEdXLow", 0.5)),
      fdEdXMed(pset.get<double>("dEdXMed", 2.8)),
      fdEdXHigh(pset.get<double>("dEdXHigh", 3.4)),
      fDoMichel(pset.get<bool>("DoMichel", true)),
      fMichelCut(pset.get<double>("MichelCut", 0.6)),
      fDoNTracks(pset.get<bool>("DoNTracks", false)),
      fExcludeRuns(pset.get<std::vector<int>>("ExcludeRuns", {})),
      fDoAngle(pset.get<bool>("DoAngle", true)),
      fXYZCut(pset.get<double>("XYZCut", 3.)),
      fCosCut(pset.get<double>("CosCut", 0.95)),
      fMCSpot(detail::ReadSpot(pset.get<fhicl::ParameterSet>("MCBeamSpot", {}),
                               {-30.81, 422.41, 0.11, 5.02, 4.51, 0.22, 101.58, 101.19, 16.59})),
      fDataSpot(detail::ReadSpot(pset.get<fhicl::ParameterSet>("DataBeamSpot", {}),
                                 {-27.91, 424.36, 3.78, 4.71, 5.16, 1.10, 100.95, 104.17, 18.05}))
  {
    const auto imp = pset.get<fhicl::ParameterSet>("Improved", {});
    fImproved.beamPFPByScore      = imp.get<bool>("BeamPFPByScore", true);
    fImproved.michelAtEnd         = imp.get<bool>("MichelAtEnd", true);
    fImproved.requireTwoCaloPts   = imp.get<bool>("RequireTwoCaloPoints", true);
    fImproved.requireShowerEnergy = imp.get<bool>("RequireShowerEnergy", true);
    fImproved.biQualityInMC       = imp.get<bool>("BIQualityInMC", true);

    // Same lookup as the old analyzer's OpenFile (PDSPA:L359-L389): a path, else FW_SEARCH_PATH.
    const std::string name = pset.get<std::string>("dEdXTemplateFile", "dEdxrestemplates.root");
    std::string path = name;
    if (!cet::file_exists(path)) {
      cet::search_path sp{"FW_SEARCH_PATH"};
      if (!sp.find_file(name, path)) {
        throw cet::exception("PDSPSelection") << "Could not find the dE/dx template file " << name
                                              << " directly or in FW_SEARCH_PATH";
      }
    }
    fTemplateFile.reset(TFile::Open(path.c_str(), "READ"));
    if (!fTemplateFile || fTemplateFile->IsZombie()) {
      throw cet::exception("PDSPSelection") << "Could not open " << path;
    }
    fProtonTemplate = dynamic_cast<TProfile*>(fTemplateFile->Get("dedx_range_pro"));
    if (!fProtonTemplate) {
      throw cet::exception("PDSPSelection") << "No TProfile dedx_range_pro in " << path;
    }
    mf::LogInfo("PDSPSelection") << "Proton dE/dx template from " << path;
  }

  //------------------------------------------------------------------------------

  inline void PDSPSelection::Fill(art::Event const& evt, caf::SRBeamInstrumentation const& inst,
                                  caf::SRPDSPSelectionBranch& out) const
  {
    try {
      FillTruth(evt, out);

      EventData ed;
      ed.pfps = evt.getHandle<std::vector<recob::PFParticle>>(fPFPLabel);
      ed.tracks = evt.getHandle<std::vector<recob::Track>>(fTrackLabel);
      ed.hits = evt.getHandle<std::vector<recob::Hit>>(fHitLabel);
      if (!ed.pfps.isValid() || !ed.tracks.isValid() || !ed.hits.isValid()) {
        mf::LogWarning("PDSPSelection") << "Missing PFParticles, tracks or hits; selection not evaluated.";
        return;
      }
      ed.forcedTracks  = evt.getHandle<std::vector<recob::Track>>(fForcedTrackLabel);
      ed.forcedShowers = evt.getHandle<std::vector<recob::Shower>>(fForcedShowerLabel);
      ed.clusters      = evt.getHandle<std::vector<recob::Cluster>>(fPFPLabel);
      art::fill_ptr_vector(ed.hitPtrs, ed.hits);

      using namespace detail;
      ed.pfpMeta         = TryMake<art::FindManyP<larpandoraobj::PFParticleMetadata>>(ed.pfps, evt, fPFPLabel);
      ed.pfpSlice        = TryMake<art::FindOneP<recob::Slice>>(ed.pfps, evt, fPFPLabel);
      ed.pfpTrack        = TryMake<art::FindManyP<recob::Track>>(ed.pfps, evt, fTrackLabel);
      ed.pfpForcedTrack  = TryMake<art::FindManyP<recob::Track>>(ed.pfps, evt, fForcedTrackLabel);
      ed.pfpForcedShower = TryMake<art::FindManyP<recob::Shower>>(ed.pfps, evt, fForcedShowerLabel);
      ed.pfpClusters     = TryMake<art::FindManyP<recob::Cluster>>(ed.pfps, evt, fPFPLabel);
      if (ed.clusters.isValid())
        ed.clusterHits   = TryMake<art::FindManyP<recob::Hit>>(ed.clusters, evt, fPFPLabel);
      ed.beamCalo        = TryMake<art::FindManyP<anab::Calorimetry>>(ed.tracks, evt, fBeamCaloLabel);
      ed.trackHits       = TryMake<art::FindManyP<recob::Hit>>(ed.tracks, evt, fTrackLabel);
      ed.trackHitMetas   = TryMake<art::FindManyP<recob::Hit, recob::TrackHitMeta>>(ed.tracks, evt, fTrackLabel);
      ed.hitTracks       = TryMake<art::FindManyP<recob::Track>>(ed.hits, evt, fTrackLabel);
      if (ed.forcedTracks.isValid())
        ed.forcedTrackCalo  = TryMake<art::FindManyP<anab::Calorimetry>>(ed.forcedTracks, evt, fDaughterCaloLabel);
      if (ed.forcedShowers.isValid())
        ed.forcedShowerCalo = TryMake<art::FindManyP<anab::Calorimetry>>(ed.forcedShowers, evt, fDaughterCaloLabel);
      ed.cnn = std::make_unique<anab::MVAReader<recob::Hit, 4>>(evt, fMVALabel);

      if (!ed.pfpMeta || !ed.pfpTrack || !ed.beamCalo || !ed.pfpForcedTrack || !ed.pfpForcedShower ||
          !ed.forcedTrackCalo || !ed.forcedShowerCalo || !ed.pfpClusters || !ed.clusterHits ||
          !ed.trackHits || !ed.trackHitMetas || !ed.hitTracks) {
        mf::LogWarning("PDSPSelection")
          << "Missing association(s) needed by the selection (check that " << fBeamCaloLabel << " and "
          << fDaughterCaloLabel << " are produced in this job); selection not evaluated.";
        return;
      }

      const int truePDG = out.true_beam_pdg;
      RunVariant(evt, ed, inst, truePDG, Variant{}, out.legacy);
      RunVariant(evt, ed, inst, truePDG, fImproved, out.improved);
    }
    catch (cet::exception const& e) {
      mf::LogWarning("PDSPSelection") << "Selection not evaluated for this event: " << e.what();
      out.legacy.event_class = out.improved.event_class = caf::kPDSPNotEvaluated;
      out.legacy.event_class_inclusive = out.improved.event_class_inclusive = caf::kPDSPNotEvaluated;
    }
  }

  //------------------------------------------------------------------------------

  inline void PDSPSelection::FillTruth(art::Event const& evt, caf::SRPDSPSelectionBranch& out) const
  {
    if (evt.isRealData()) return;

    // True beam particle, as protoana::ProtoDUNETruthUtils::GetGeantGoodParticle: the first
    // "primary" of MCTruth[0], then the first G4 particle (TrackId order) with the same PDG
    // and |dE| < 1e-5. The old analyzer dropped the event when this fails (PDSPA:L1335-L1338).
    auto truthHandle = evt.getHandle<std::vector<simb::MCTruth>>(fGeneratorLabel);
    if (!truthHandle.isValid() || truthHandle->empty()) return;
    simb::MCTruth const& genTruth = truthHandle->front();

    const simb::MCParticle* goodGen = nullptr;
    for (int t = 0; t < genTruth.NParticles(); ++t) {
      if (genTruth.GetParticle(t).Process() == "primary") { goodGen = &genTruth.GetParticle(t); break; }
    }
    if (!goodGen) return;

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    sim::ParticleList const& plist = pi_serv->ParticleList();
    const simb::MCParticle* beam = nullptr;
    for (auto const& entry : plist) {
      if (entry.second && entry.second->PdgCode() == goodGen->PdgCode() &&
          std::fabs(entry.second->E() - goodGen->E()) < 1e-5) {
        beam = entry.second;
        break;
      }
    }
    if (!beam) return;

    out.true_beam_pdg  = beam->PdgCode();
    out.true_beam_g4id = beam->TrackId();
    out.true_beam_end  = caf::SRVector3D(beam->EndX(), beam->EndY(), beam->EndZ());
    const std::string endProcess = beam->EndProcess();
    out.true_beam_pip_inelastic = (endProcess == "pi+Inelastic");

    // Direct daughters, skipping muIoni / hIoni (PDSPA:L3827-L3875).
    int nPi0 = 0;
    int nChPiAbove = 0;
    for (int i = 0; i < beam->NumberDaughters(); ++i) {
      auto it = plist.find(beam->Daughter(i));
      if (it == plist.end() || !it->second) continue; // the old code would have crashed here
      simb::MCParticle const* d = it->second;
      const std::string process = d->Process();
      if (process == "muIoni" || process == "hIoni") continue;
      if (std::abs(d->PdgCode()) == 211 && d->P() > fThreshold) ++nChPiAbove;
      if (process.find("Inelastic") != std::string::npos && d->PdgCode() == 111) ++nPi0;
    }
    out.true_npi0 = nPi0;
    out.true_nchpi_above_thr = nChPiAbove;

    // new_interaction_topology (SD:L13-L67).
    const int pdg = out.true_beam_pdg;
    const double endZ = beam->EndZ();
    int topology = -1;
    if (pdg == 211) {
      if (endZ < fEndZLow) topology = 4;
      else if (endZ > fEndZHigh && (!fSignalPastFV || endProcess != "pi+Inelastic")) topology = 6;
      else if (endProcess == "pi+Inelastic") {
        if (nChPiAbove > 0 || (!fCexNPi0 && nPi0 > 1)) topology = 3;
        else if (nPi0 == 0) topology = 1;
        else if ((fCexNPi0 && nPi0 > 0) || (!fCexNPi0 && nPi0 == 1)) topology = 2;
      }
      else topology = 7;
    }
    else if (pdg == -13) topology = 5;
    else topology = 7;
    out.true_event_class = topology;
  }

  //------------------------------------------------------------------------------

  inline PDSPSelection::BeamChoice PDSPSelection::ChooseBeamPFP(EventData const& ed, bool byScore) const
  {
    // Slice map of the primary PFParticles, ordered by slice index (ProtoDUNEPFParticleUtils
    // SliceMapHelper): slice from the Slice association, else metadata "SliceIndex", 9999 dropped.
    std::map<unsigned int, std::vector<size_t>> sliceMap;
    for (size_t p = 0; p < ed.pfps->size(); ++p) {
      recob::PFParticle const& pfp = ed.pfps->at(p);
      if (!pfp.IsPrimary()) continue;
      unsigned int slice = 9999;
      art::Ptr<recob::Slice> slicePtr;
      if (ed.pfpSlice) slicePtr = ed.pfpSlice->at(pfp.Self());
      if (slicePtr.isNonnull()) slice = static_cast<unsigned short>(slicePtr->ID());
      else {
        const auto md = ed.Meta(pfp.Self());
        auto it = md.find("SliceIndex");
        if (it != md.end()) slice = static_cast<unsigned short>(it->second);
      }
      if (slice != 9999) sliceMap[slice].push_back(p);
    }

    BeamChoice choice;
    for (auto const& [slice, members] : sliceMap) {
      bool sliceTagged = false;
      for (size_t p : members) {
        const auto md = ed.Meta(ed.pfps->at(p).Self());
        if (md.find("IsTestBeam") == md.end()) continue;
        sliceTagged = true;
        if (!byScore) break;
        auto it = md.find("TestBeamScore");
        const float score = (it != md.end()) ? it->second : kMissing;
        if (choice.pfpIndex < 0 || score > choice.score) {
          choice.pfpIndex = static_cast<int>(p);
          choice.slice = static_cast<int>(slice);
          choice.tagged = true;
          choice.score = score;
        }
      }
      if (!sliceTagged) continue;
      ++choice.nTaggedSlices;
      // Legacy (PDSPA:L1427, L1453-L1456): the first primary of the lowest-index tagged slice,
      // which need not be the tagged PFParticle itself.
      if (!byScore && choice.pfpIndex < 0) {
        choice.pfpIndex = static_cast<int>(members.front());
        choice.slice = static_cast<int>(slice);
        const auto md = ed.Meta(ed.pfps->at(members.front()).Self());
        choice.tagged = (md.find("IsTestBeam") != md.end());
        auto it = md.find("TestBeamScore");
        choice.score = (it != md.end()) ? it->second : kMissing;
      }
    }
    return choice;
  }

  //------------------------------------------------------------------------------

  inline std::vector<anab::Calorimetry> PDSPSelection::TrackCalo(EventData const& ed, bool forced,
                                                                 recob::Track const& track) const
  {
    // As ProtoDUNETrackUtils::GetRecoTrackCalorimetry: the association is indexed by track ID,
    // and any failure gives an empty result.
    std::vector<anab::Calorimetry> result;
    auto const& finder = forced ? ed.forcedTrackCalo : ed.beamCalo;
    if (!finder || track.ID() < 0 || static_cast<size_t>(track.ID()) >= finder->size()) return result;
    for (auto const& c : finder->at(track.ID())) result.push_back(*c);
    return result;
  }

  //------------------------------------------------------------------------------

  inline std::pair<double, int> PDSPSelection::MichelScore(EventData const& ed, recob::Track const& track,
                                                           geo::Point_t const& ref, caf::SRVector3D& vtx) const
  {
    // Port of ProtoDUNETrackUtils::GetVertexMichelScore with its default arguments (PTU:L653-L776).
    auto start = track.Vertex();
    auto end = track.End();
    if (std::max(start.X(), end.X()) > 0. || std::min(start.X(), end.X()) < -200. ||
        std::max(start.Y(), end.Y()) > 500. || std::min(start.Y(), end.Y()) < 200. ||
        std::min(start.Z(), end.Z()) < 25. || track.Length() < 5.) {
      return {-1., 0};
    }

    // Collection-plane hits of the track, keyed by trajectory point (GetRecoHitsFromTrajPoints).
    size_t trackIndex = 0;
    for (size_t i = 0; i < ed.tracks->size(); ++i) {
      if ((*ed.tracks)[i].ID() == track.ID()) { trackIndex = i; break; }
    }
    std::map<size_t, recob::Hit const*> trajHits;
    if (track.ID() >= 0 && static_cast<size_t>(track.ID()) < ed.trackHits->size()) {
      auto const& hitsOfTrack = ed.trackHits->at(track.ID());
      auto const& metaHits = ed.trackHitMetas->at(trackIndex);
      auto const& metas = ed.trackHitMetas->data(trackIndex);
      for (size_t i = 0; i < metaHits.size(); ++i) {
        const unsigned int idx = metas[i]->Index();
        if (idx == static_cast<unsigned int>(std::numeric_limits<int>::max())) continue;
        if (!track.HasValidPoint(idx)) continue;
        for (auto const& h : hitsOfTrack) {
          if (metaHits[i].key() == h.key() && idx < track.NumberTrajectoryPoints()) trajHits[idx] = h.get();
        }
      }
    }
    std::vector<recob::Hit const*> view2Hits;
    std::vector<size_t> view2Index;
    for (auto const& [idx, hit] : trajHits) {
      if (hit->View() != 2) continue;
      view2Hits.push_back(hit);
      view2Index.push_back(idx);
    }

    // Vertex hit: the trajectory hit farthest from `ref`.
    int vertexTPC = -1, vertexWire = -1;
    float vertexTime = -1.;
    double highest = -1.;
    for (size_t i = 0; i < view2Hits.size(); ++i) {
      auto const loc = track.Trajectory().LocationAtPoint(view2Index[i]);
      const double d = std::sqrt(std::pow(loc.X() - ref.X(), 2) + std::pow(loc.Y() - ref.Y(), 2) +
                                 std::pow(loc.Z() - ref.Z(), 2));
      if (d > highest) {
        highest = d;
        vertexTPC = view2Hits[i]->WireID().TPC;
        vertexTime = view2Hits[i]->PeakTime();
        vertexWire = view2Hits[i]->WireID().Wire;
        vtx = caf::SRVector3D(loc.X(), loc.Y(), loc.Z());
      }
    }

    // Window of +-15 wires, +-100 ticks on the same TPC, collection plane, excluding the track's
    // own trajectory hits and hits belonging to another track longer than 25 cm.
    std::pair<double, int> result{0., 0};
    for (size_t i = 0; i < ed.hitPtrs.size(); ++i) {
      recob::Hit const* hit = ed.hitPtrs[i].get();
      if (std::find(view2Hits.begin(), view2Hits.end(), hit) != view2Hits.end()) continue;
      const int wire = hit->WireID().Wire;
      if (std::abs(wire - vertexWire) > 15 || std::abs(hit->PeakTime() - vertexTime) > 100. ||
          static_cast<int>(hit->WireID().TPC) != vertexTPC || hit->View() != 2) continue;
      auto const& owners = ed.hitTracks->at(ed.hitPtrs[i].key());
      if (!owners.empty() && owners[0]->ID() != track.ID() && owners[0]->Length() > 25.) continue;
      const auto cnnOut = ed.cnn->getOutput(ed.hitPtrs[i]);
      result.first += cnnOut[ed.cnn->getIndex("michel")];
      result.second += 1;
    }
    return result;
  }

  //------------------------------------------------------------------------------

  inline std::pair<double, int> PDSPSelection::Chi2PID(std::vector<float> const& dedx,
                                                       std::vector<float> const& range,
                                                       TProfile* profile)
  {
    // Port of ProtoDUNETrackUtils::Chi2PID (PTU:L424-L470).
    if (dedx.empty() || range.empty()) return {9999., -1};
    double chi2 = 0.;
    int npt = 0;
    for (size_t i = 1; i + 1 < dedx.size(); ++i) {  // first and last point ignored
      if (dedx[i] > 1000. || dedx[i] < 0.) continue;
      const int bin = profile->FindBin(range[i]);
      if (bin < 1 || bin > profile->GetNbinsX()) continue;
      double t = profile->GetBinContent(bin);
      if (t < 1.e-6) t = (profile->GetBinContent(bin - 1) + profile->GetBinContent(bin + 1)) / 2.;
      double terr = profile->GetBinError(bin);
      if (terr < 1.e-6) terr = (profile->GetBinError(bin - 1) + profile->GetBinError(bin + 1)) / 2.;
      const double res = (0.04231 + 0.0001783 * dedx[i] * dedx[i]) * dedx[i];
      chi2 += std::pow(dedx[i] - t, 2) / (terr * terr + res * res);
      ++npt;
    }
    if (npt == 0) return {9999., -1};
    return {chi2, npt};
  }

  inline double PDSPSelection::TruncatedMean(std::vector<float> const& dedx) const
  {
    // Port of truncatedMean_pos (SD:L367-L446).
    std::vector<double> v;
    for (float d : dedx) if (d >= 0.) v.push_back(d);
    if (v.empty()) return kEmptydEdx;
    std::sort(v.begin(), v.end());
    const int lo = std::rint(v.size() * fLimit);
    const int hi = std::rint(v.size() * (1. - fLimit));
    double sum = 0.;
    int n = 0;
    for (int i = lo; i < hi; ++i) { sum += v[i]; ++n; }
    return sum / n;
  }

  inline bool PDSPSelection::ShowerDistEnergyCheck(double dist, double energy)
  {
    // Port of shower_dist_energy_check (SD:L479-L508); all values hardcoded there too.
    const double x1 = 7.5, y1 = 880.;
    const double x2 = 90., y2 = 305.;
    if (energy < 140. && dist < 2.5) return false;
    if (energy < 60. && dist < 5.) return false;
    if (energy < 20. && dist < 15.) return false;
    if (energy > 880.) return false;
    if (dist > 90.) return false;
    if (((dist - x1) * (y2 - y1) - (energy - y1) * (x2 - x1)) < 0) return false;
    return true;
  }

  //------------------------------------------------------------------------------

  inline void PDSPSelection::RunVariant(art::Event const& evt, EventData const& ed,
                                        caf::SRBeamInstrumentation const& inst, int truePDG,
                                        Variant const& var, caf::SRPDSPSelection& sel) const
  {
    const bool isData = evt.isRealData();
    BeamSpot const& spot = isData ? fDataSpot : fMCSpot;

    // ---- Prefilter: the old Filter() calls (rES:L187, L309, L394) ----
    const bool biQuality = inst.valid && inst.nmomenta == 1 && (!fDoNTracks || inst.ntracks == 1);
    if (isData) {
      const bool pid = std::find(inst.PDG_candidates.begin(), inst.PDG_candidates.end(), 211) != inst.PDG_candidates.end();
      const bool goodRun = std::find(fExcludeRuns.begin(), fExcludeRuns.end(), static_cast<int>(evt.run())) == fExcludeRuns.end();
      sel.pass_prefilter = pid && biQuality && goodRun;
    }
    else {
      sel.pass_prefilter = (truePDG == 211 || truePDG == -13) && (!var.biQualityInMC || biQuality);
    }

    // ---- Beam particle ----
    const BeamChoice beam = ChooseBeamPFP(ed, var.beamPFPByScore);
    sel.nbeam_slices = beam.nTaggedSlices;
    recob::PFParticle const* beamPFP = nullptr;
    art::Ptr<recob::Track> beamTrack;
    if (beam.pfpIndex >= 0) {
      beamPFP = &ed.pfps->at(beam.pfpIndex);
      sel.beam_slice = beam.slice;
      sel.beam_pfp_id = static_cast<int>(beamPFP->Self());
      sel.beam_pfp_is_tagged = beam.tagged;
      sel.beam_score = beam.score;
      auto const& trks = ed.pfpTrack->at(beamPFP->Self());
      if (!trks.empty()) beamTrack = trks.front();
    }

    // Old reco_beam_end*: -999 without a beam track (PDSPA:L2274-L2276).
    double endX = kMissing, endY = kMissing, endZ = kMissing;
    bool hasCalo = false;
    if (beamTrack) {
      sel.beam_track_id = beamTrack->ID();
      auto s = beamTrack->Trajectory().Start();
      auto e = beamTrack->Trajectory().End();
      if (s.Z() > e.Z()) std::swap(s, e); // PDSPA:L2906-L2914
      sel.beam_start = caf::SRVector3D(s.X(), s.Y(), s.Z());
      sel.beam_end = caf::SRVector3D(e.X(), e.Y(), e.Z());
      sel.beam_len = beamTrack->Length();
      endX = e.X(); endY = e.Y(); endZ = e.Z();

      // Collection-plane calorimetry of the beam track.
      const auto calos = TrackCalo(ed, false, *beamTrack);
      for (auto const& c : calos) {
        if (static_cast<int>(c.PlaneID().Plane) != fCollectionPlane) continue;
        hasCalo = true;
        auto pts = c.XYZ();
        sel.beam_calo_npts = static_cast<int>(pts.size());

        // Thin-slice interacting energy (PDSPA:L3324-L3349), from the beam-line momentum.
        struct P { double dedx, pitch, z; };
        std::vector<P> cp;
        for (size_t i = 0; i < c.dEdx().size(); ++i) cp.push_back({c.dEdx()[i], c.TrkPitchVec()[i], pts[i].Z()});
        std::sort(cp.begin(), cp.end(), [](P const& a, P const& b) { return a.z < b.z; });
        if (std::isfinite(inst.P)) {
          double E = std::sqrt(1.e6 * inst.P * inst.P + kPionMass * kPionMass) - kPionMass;
          if (cp.size() >= 2) {
            for (size_t i = 0; i + 1 < cp.size(); ++i) {
              if (cp[i].dedx < 0.) continue;
              E -= cp[i].dedx * cp[i].pitch;
            }
          }
          sel.beam_interacting_energy = E;
        }

        std::sort(pts.begin(), pts.end(), [](auto const& a, auto const& b) { return a.Z() < b.Z(); });
        if (!pts.empty()) {
          sel.beam_calo_start = caf::SRVector3D(pts.front().X(), pts.front().Y(), pts.front().Z());
          sel.beam_calo_end = caf::SRVector3D(pts.back().X(), pts.back().Y(), pts.back().Z());
        }
        break;
      }
    }

    // ---- S1 primary_isBeamType (SD:L539-L552) ----
    sel.beam_is_track = beamTrack && (!fCheckCalo || hasCalo);

    // ---- S2 passBeamCut, TPC version (SD:L671-L723) ----
    if (sel.beam_calo_npts > 0) {
      const double sx = sel.beam_calo_start.X(), sy = sel.beam_calo_start.Y(), sz = sel.beam_calo_start.Z();
      const double dx = sel.beam_calo_end.X() - sx, dy = sel.beam_calo_end.Y() - sy, dz = sel.beam_calo_end.Z() - sz;
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      const double deg = TMath::Pi() / 180.;
      const double cosTheta = std::cos(spot.thetaX * deg) * dx / r + std::cos(spot.thetaY * deg) * dy / r +
                              std::cos(spot.thetaZ * deg) * dz / r;  // NaN when r == 0
      sel.beam_dx_sigma = (sx - spot.meanX) / spot.sigmaX;
      sel.beam_dy_sigma = (sy - spot.meanY) / spot.sigmaY;
      sel.beam_dz_sigma = (sz - spot.meanZ) / spot.sigmaZ;
      sel.beam_costheta = cosTheta;
      bool pass = !(std::fabs(sel.beam_dx_sigma) > fXYZCut) && !(std::fabs(sel.beam_dy_sigma) > fXYZCut) &&
                  !(std::fabs(sel.beam_dz_sigma) > fXYZCut);
      if (fDoAngle && cosTheta < fCosCut) pass = false;       // a NaN cosTheta does not fail here
      if (var.requireTwoCaloPts && !(r > 0.)) pass = false;   // improved: a point-like track fails
      sel.pass_beam_cut = pass;
    }
    // Without calorimetry points the old start was -999, which always fails the position cut.

    // ---- S3 primary_ends_inAPA3 (SD:L554-L564), raw end z ----
    sel.ends_in_fv = (endZ < fEndZHigh);

    // ---- S4 vertex_cut (SD:L725-L734) ----
    if (beamTrack) {
      // Legacy: farthest from the origin (check_x/y/z = 0). Improved: farthest from the upstream
      // end, i.e. the downstream (interaction) end.
      const geo::Point_t ref = var.michelAtEnd
        ? geo::Point_t(sel.beam_start.X(), sel.beam_start.Y(), sel.beam_start.Z())
        : geo::Point_t(0., 0., 0.);
      const auto michel = MichelScore(ed, *beamTrack, ref, sel.michel_vtx);
      sel.michel_score = michel.first;
      sel.michel_nhits = michel.second;
      sel.michel_tagged = (michel.second > 0 && (michel.first / michel.second) > fMichelCut);
    }

    // ---- Daughters: S5 has_noPion_daughter, S6 has_shower_dist_energy ----
    if (beamPFP) {
      for (size_t daughterID : beamPFP->Daughters()) {
        recob::PFParticle const& dPFP = ed.pfps->at(daughterID);
        caf::SRPDSPSelDaughter d;
        d.pfp_id = static_cast<int>(dPFP.Self());

        // Collection-plane, charge-weighted CNN track score (PDSPA:L493-L531, plane 2 hardcoded there).
        double qTrack = 0., qTot = 0.;
        int nColl = 0;
        for (auto const& cl : ed.pfpClusters->at(dPFP.Self())) {
          if (cl->ID() < 0 || static_cast<size_t>(cl->ID()) >= ed.clusterHits->size()) continue;
          for (auto const& h : ed.clusterHits->at(cl->ID())) {
            if (h->WireID().Plane != 2) continue;
            const auto out = ed.cnn->getOutput(h);
            qTrack += h->Integral() * out[ed.cnn->getIndex("track")];
            qTot += h->Integral();
            ++nColl;
          }
        }
        d.nhits_coll = nColl;
        const double score = (qTot != 0.) ? qTrack / qTot : kMissing;
        d.track_score = score;

        // Forced track fit (PDSPA:L4428-L4530).
        std::vector<float> dedx, range;
        double chi2 = kMissing;
        int ndof = static_cast<int>(kMissing);
        auto const& ftrk = ed.pfpForcedTrack->at(dPFP.Self());
        if (!ftrk.empty()) {
          d.track_id = ftrk.front()->ID();
          for (auto const& c : TrackCalo(ed, true, *ftrk.front())) {
            if (c.PlaneID().Plane != 2) continue;
            dedx = c.dEdx();
            range = c.ResidualRange();
            d.ncalo_pts = static_cast<int>(dedx.size());
            std::tie(chi2, ndof) = Chi2PID(dedx, range, fProtonTemplate);
            break;
          }
        }
        const double trunc = TruncatedMean(dedx);
        d.trunc_dEdx = trunc;
        d.chi2_proton = chi2;
        d.chi2_ndof = ndof;

        // secondary_noPion (SD:L585-L596).
        if (d.track_id != -999 && score > fTrackScoreCut) {
          if (trunc < fdEdXMed && trunc > fdEdXLow) d.is_pion_like = true;
          else if (trunc < fdEdXHigh && ndof > 0 && chi2 / ndof > fChi2Cut) d.is_pion_like = true;
        }

        // Forced shower fit (PDSPA:L4734-L4846) and has_shower_dist_energy (SD:L510-L537).
        double shX = kMissing, shY = kMissing, shZ = kMissing, shE = kMissing;
        bool shCalo = false;
        auto const& fshw = ed.pfpForcedShower->at(dPFP.Self());
        if (!fshw.empty()) {
          auto const& shw = fshw.front();
          d.shower_id = shw->ID();
          shX = shw->ShowerStart().X(); shY = shw->ShowerStart().Y(); shZ = shw->ShowerStart().Z();
          if (shw->ID() >= 0 && static_cast<size_t>(shw->ID()) < ed.forcedShowerCalo->size()) {
            for (auto const& c : ed.forcedShowerCalo->at(shw->ID())) {
              if (c->PlaneID().Plane != 2) continue;
              shE = c->KineticEnergy();
              shCalo = true;
              break;
            }
          }
        }
        d.shower_start = caf::SRVector3D(shX, shY, shZ);
        const double dist = std::sqrt(std::pow(shX - endX, 2) + std::pow(shY - endY, 2) + std::pow(shZ - endZ, 2));
        d.shower_dist = dist;
        d.shower_energy = shE;
        d.is_pi0_shower = score < fTrackScoreCut && score > 0. && ShowerDistEnergyCheck(dist, shE) &&
                          (!var.requireShowerEnergy || shCalo);

        if (d.is_pion_like) sel.has_pion_daughter = true;
        if (d.is_pi0_shower) sel.has_pi0_shower = true;
        sel.daughters.push_back(d);
      }
    }
    sel.ndaughters = static_cast<int>(sel.daughters.size());

    // ---- selection_ID / selection_ID_inclusive (SD:L112-L178) ----
    if (!sel.pass_prefilter) {
      sel.event_class = sel.event_class_inclusive = caf::kPDSPFailPrefilter;
      return;
    }
    const bool michel = fDoMichel && sel.michel_tagged;
    if (!sel.beam_is_track)         sel.event_class = caf::kPDSPNoBeamTrack;
    else if (!sel.pass_beam_cut)    sel.event_class = caf::kPDSPFailBeamCut;
    else if (!sel.ends_in_fv)       sel.event_class = caf::kPDSPPastFV;
    else if (michel)                sel.event_class = caf::kPDSPMichel;
    else if (sel.has_pion_daughter) sel.event_class = caf::kPDSPRejectedInt;
    else if (sel.has_pi0_shower)    sel.event_class = caf::kPDSPCex;
    else                            sel.event_class = caf::kPDSPAbs;

    if (!sel.beam_is_track)         sel.event_class_inclusive = 4;
    else if (!sel.pass_beam_cut)    sel.event_class_inclusive = 3;
    else if (!sel.ends_in_fv)       sel.event_class_inclusive = 2;
    else if (michel)                sel.event_class_inclusive = 5;
    else                            sel.event_class_inclusive = 1;
  }

} // pdspsel
} // caf

#endif // CAFMAKER_PDSPSELECTION_H
