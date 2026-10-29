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
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
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
      void PreLoadMCParticlesInfo(art::Event const& evt);
      void FillTruthInfo(caf::SRTruthBranch& sr,
                         art::Event const& evt);
      void FillParticleIDRefs(std::vector<caf::SRTrueParticle>& particles) const;
      void FillAncestorIDs(caf::SRTrueInteraction& inter) const;
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
      void FillRecoInfo(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, const caf::SRTruthBranch &truthBranch, const caf::SRBeamInstrumentation &beamInst, const art::Event &evt, const art::Ptr<recob::Slice> &slicePtr, const art::FindManyP<recob::PFParticle> &sliceToPFP) const;
      void FillRecoInfoSliceLoop(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, const caf::SRTruthBranch &truthBranch, const caf::SRBeamInstrumentation &beamInst, const art::Event &evt) const;
      void FillCVNInfo(caf::SRCVNScoreBranch &cvnBranch, const art::Event &evt) const;
      void FillEnergyInfo(caf::SRNeutrinoEnergyBranch &ErecBranch, const art::Event &evt) const;
      void FillRecoParticlesInfo(caf::SRRecoParticlesBranch &recoParticlesBranch, caf::SRFD &fdBranch, const art::Event &evt, const art::Ptr<recob::Slice> &slicePtr, const art::FindManyP<recob::PFParticle> &sliceToPFP) const;
      void FillDirectionInfo(caf::SRDirectionBranch &dirBranch, const art::Event &evt) const;
      int FillGENIERecord(simb::MCTruth const& mctruth, simb::GTruth const& gtruth);
      void FillGeneratorInfo(simb::MCGeneratorInfo const& genInfo, caf::SRTrueInteraction &inter) const;
      double GetVisibleEnergy(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const;
      void FillTruthMatchingAndOverlap(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt, std::vector<TrueParticleID> &truth, std::vector<float> &truthOverlap) const;
      void FillInteractionTruthMatchingAndOverlap(lar_pandora::PFParticleVector const& particles, const art::Event &evt, std::vector<std::size_t> &truth, std::vector<float> &truthOverlap) const;
      void FillPFPMetadata(caf::SRPFP &pfpBranch, art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const;
      double GetWallDistance(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const;
      double GetWallDistance(recob::SpacePoint const& sp) const;
      void ComputeActiveBounds();
      double GetSingleHitsEnergy(art::Event const& evt, int plane) const;
      bool IsVertexContained(caf::SRVector3D const& vtx) const;

      /// @name ProtoDUNE thin-target interacting energies, ported from protoduneana's PDSPAnalyzer
      ///@{
      /// Collection-plane index (also the view index used for the true energy depositions).
      int CollectionPlane() const { return fVPlaneAsCollector ? 1 : 2; }
      /// Group IDEs into wire-pitch slices along z, starting from the first wire at z0.
      std::map<int, std::vector<const sim::IDE*>> SliceIDEs(std::vector<const sim::IDE*> const& ides, double z0, double pitch) const;
      /// True kinetic energy [MeV] of a G4 particle at its interaction point, from its collection-view IDEs.
      double ComputeTrueInteractingEnergy(const simb::MCParticle &trueParticle) const;
      /// Reco kinetic energy [MeV] of the beam track at its interaction point, from its SCE-corrected calorimetry.
      double ComputeRecoInteractingEnergy(const art::Event &evt, const art::Ptr<recob::Track> &beamTrack, const caf::SRTruthBranch &truthBranch, const caf::SRBeamInstrumentation &beamInst) const;
      ///@}
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

      std::string fBeamModuleLabel;       ///< Label of the beam::ProtoDUNEBeamEvent product in real data
      std::string fMCBeamModuleLabel;     ///< Label of the beam::ProtoDUNEBeamEvent product in simulation
      double fBeamInstPFix;               ///< Momentum correction factor for MC beam instrumentation
      double fBeamPIDMomentum;            ///< Nominal beam momentum [GeV/c] the PID cuts are defined at
      bool   fUseCERNCalibSelection;      ///< Use the CERN-calibrated TOF cuts (reco after ~v08_07_00) rather than the older values

      std::string fBeamTrackLabel;        ///< Track label for the beam particle track (its calorimetry is fCalorimetryLabelSCE)
      std::string fCalorimetryLabelSCE;   ///< SCE-corrected calorimetry label for the beam track
      bool fVPlaneAsCollector;            ///< Use the V plane (1) instead of the W plane (2) as collection plane
      bool fMCHasBI;                      ///< MC has beam instrumentation: start the reco interacting energy from inst.P instead of truth
      geo::WireReadoutGeom const* fWireReadout; ///< Wire readout geometry, for the wire pitch

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
      fEventRecord(new genie::NtpMCEventRecord),
      fGeom(&*art::ServiceHandle<geo::Geometry>()),
      fVertexFiducialVolumeCut(pset.get<std::vector<double>>("VertexFiducialVolumeCut")),
      fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
      fRecombFactor(pset.get<double>("RecombFactor")),
      fBeamModuleLabel(pset.get<std::string>("BeamModuleLabel")),
      fMCBeamModuleLabel(pset.get<std::string>("MCBeamModuleLabel")),
      fBeamInstPFix(pset.get<double>("BeamInstPFix")),
      fBeamPIDMomentum(pset.get<double>("BeamPIDMomentum")),
      fUseCERNCalibSelection(pset.get<bool>("UseCERNCalibSelection")),
      fBeamTrackLabel(pset.get<std::string>("BeamTrackLabel")),
      fCalorimetryLabelSCE(pset.get<std::string>("CalorimetryLabelSCE")),
      fVPlaneAsCollector(pset.get<bool>("VPlaneAsCollector")),
      fMCHasBI(pset.get<bool>("MCHasBI")),
      fWireReadout(&art::ServiceHandle<geo::WireReadout>()->Get())
  {

    if(pset.get<bool>("CreateFlatCAF")){
      // LZ4 is the fastest format to decompress. I get 3x faster loading with
      // this compared to the default, and the files are only slightly larger.
      fFlatFile = std::make_unique<TFile>("flatcaf.root", "RECREATE", "",
                            ROOT::CompressionSettings(ROOT::kLZ4, 1));
    }

    ComputeActiveBounds();

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
  // Fill parentID and daughtersID for a set of particles once fMCParticlesMap
  // holds the final SR indices for all particles in their interaction.
  void CAFMaker::FillParticleIDRefs(std::vector<caf::SRTrueParticle>& particles) const
  {
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
  }

  //------------------------------------------------------------------------------
  // Fill ancestor_id: the root primary particle each particle ultimately descended
  // from, found by walking up parentID links. Must run after FillParticleIDRefs
  // has populated parentID for every particle in the interaction.
  void CAFMaker::FillAncestorIDs(caf::SRTrueInteraction& inter) const
  {
    auto getParticle = [&](const caf::TrueParticleID& id) -> const caf::SRTrueParticle* {
      if (id.type == caf::TrueParticleID::kPrimary && id.part >= 0 && (size_t)id.part < inter.prim.size())
        return &inter.prim[id.part];
      if (id.type == caf::TrueParticleID::kSecondary && id.part >= 0 && (size_t)id.part < inter.sec.size())
        return &inter.sec[id.part];
      return nullptr;
    };

    auto fillAncestor = [&](std::vector<caf::SRTrueParticle>& particles) {
      for (auto& p : particles) {
        caf::TrueParticleID current = p.parentID;
        while (current.type == caf::TrueParticleID::kSecondary) {
          const caf::SRTrueParticle* parent = getParticle(current);
          if (!parent) {
            current = caf::TrueParticleID();
            break;
          }
          current = parent->parentID;
        }
        if (current.type == caf::TrueParticleID::kPrimary) {
          p.ancestor_id = current;
        }
      }
    };

    fillAncestor(inter.prim);
    fillAncestor(inter.sec);
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

      art::FindManyP<simb::MCParticle> fmParticles(mctruthHandle, evt, fG4Label);
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
          FillGeneratorInfo(genInfo, inter);

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
              part.start_pos = caf::SRVector3D(mcpart->Position().Vect());
              part.end_pos = caf::SRVector3D(mcpart->EndPosition().Vect());
              part.parent = mcpart->Mother();

              for (int d = 0; d < mcpart->NumberDaughters(); ++d) {
                part.daughters.push_back(mcpart->Daughter(d));
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
          FillParticleIDRefs(inter.prim);
          FillParticleIDRefs(inter.sec);
          FillAncestorIDs(inter);

          truthBranch.nu.push_back(std::move(inter));
        }
        else {
          // Non-neutrino MCTruth:
          // each visible source particle becomes its own truth interaction entry.
          for (int p = 0; p < mct.NParticles(); ++p) {
            const simb::MCParticle& mcpart = mct.GetParticle(p);

            caf::SRTrueInteraction inter;
            inter.id = cumulativeInteractionIndex++;
            // abs is needed because if you look at the genenerator truth instead of the G4 truth, the TrackId is negative for the primaries. It should be always positive, the abs is just a safety measure
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
            FillGeneratorInfo(genInfo, inter);

            caf::SRTrueParticle part;
            part.pdg = mcpart.PdgCode();
            // abs is needed because if you look at the genenerator truth instead of the G4 truth, the TrackId is negative for the primaries. It should be always positive, the abs is just a safety measure
            part.G4ID = std::abs(mcpart.TrackId());
            part.interaction_id = inter.id;
            part.time = mcpart.T();
            part.p = caf::SRLorentzVector(mcpart.Momentum());
            part.start_pos = caf::SRVector3D(mcpart.Position().Vect());
            part.end_pos = caf::SRVector3D(mcpart.EndPosition().Vect());
            part.parent = mcpart.Mother();

            std::deque<int> secondaries_to_add;

            const simb::MCParticle* the_g4_part = nullptr;
            // abs is needed because if you look at the genenerator truth instead of the G4 truth, the TrackId is negative for the primaries. It should be always positive, the abs is just a safety measure            
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
              //Only for the primary: each call scans the SimChannels through the BackTracker
              part.true_beam_interactingEnergy = ComputeTrueInteractingEnergy(*the_g4_part);
            }

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
              // abs is needed because if you look at the genenerator truth instead of the G4 truth, the TrackId is negative for the primaries. It should be always positive, the abs is just a safety measure
              sec.G4ID = std::abs(sec_mcpart->TrackId());
              sec.interaction_id = inter.id;
              sec.time = sec_mcpart->T();
              sec.p = caf::SRLorentzVector(sec_mcpart->Momentum());
              sec.start_pos = caf::SRVector3D(sec_mcpart->Position().Vect());
              sec.end_pos = caf::SRVector3D(sec_mcpart->EndPosition().Vect());
              sec.parent = sec_mcpart->Mother();

              for (int d = 0; d < sec_mcpart->NumberDaughters(); ++d) {
                int daughter = sec_mcpart->Daughter(d);
                sec.daughters.push_back(daughter);
                secondaries_to_add.push_back(daughter);
              }

              inter.sec.push_back(std::move(sec));
            }

            inter.nsec = inter.sec.size();

            // Fill parentID and daughtersID now that prim/sec are complete
            FillParticleIDRefs(inter.prim);
            FillParticleIDRefs(inter.sec);
            FillAncestorIDs(inter);

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
  void CAFMaker::FillRecoInfoSliceLoop(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, const caf::SRTruthBranch &truthBranch, const caf::SRBeamInstrumentation &beamInst, const art::Event &evt) const
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
  
    for (const auto& slicePtr : slicePtrs) {
      FillRecoInfo(recoBranch, fdBranch, truthBranch, beamInst, evt, slicePtr, sliceToPFP);
    }

  }

  void CAFMaker::FillRecoInfo(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, const caf::SRTruthBranch &truthBranch, const caf::SRBeamInstrumentation &beamInst, const art::Event &evt, const art::Ptr<recob::Slice> &slicePtr, const art::FindManyP<recob::PFParticle> &sliceToPFP) const {
    SRInteractionBranch &ixn = recoBranch.ixn;

    //Only filling with Pandora Reco for the moment
    std::vector<SRInteraction> &pandora = ixn.pandora;

    lar_pandora::PFParticleVector particleVector = sliceToPFP.at(slicePtr.key());

    lar_pandora::VertexVector vertexVector;
    lar_pandora::PFParticlesToVertices particlesToVertices;
    lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraLabel, vertexVector, particlesToVertices);

    auto pfParticles_beam = evt.getValidHandle<std::vector<recob::PFParticle>>(fPandoraLabel);
    const art::FindManyP<larpandoraobj::PFParticleMetadata> findMetaData_beam(pfParticles_beam, evt, fPandoraLabel);

    std::vector<bool> isTestBeamVec(particleVector.size(), false);
    bool found_beam = false;
    for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      if (findMetaData_beam.isValid()) {
        const auto& mdVec = findMetaData_beam.at(particle.key());
        if (!mdVec.empty() && mdVec.at(0).isNonnull()) {
          const auto& mdMap = mdVec.at(0)->GetPropertiesMap();
          isTestBeamVec[n] = (mdMap.find("IsTestBeam") != mdMap.end());
        }
      }
      found_beam = found_beam || isTestBeamVec[n];
    }

    bool fill_info_condition = false;
    for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      bool is_test_beam = isTestBeamVec[n];

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

        FillInteractionTruthMatchingAndOverlap(particleVector, evt, reco.truth, reco.truthOverlap);
        reco.isFromTrigger = is_test_beam;

        if (is_test_beam) {
          if (dune_ana::DUNEAnaPFParticleUtils::HasTrack(particle, evt, fPandoraLabel, fBeamTrackLabel)) {
            art::Ptr<recob::Track> beamTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(particle, evt, fPandoraLabel, fBeamTrackLabel);
            ixn.reco_beam_interactingEnergy = ComputeRecoInteractingEnergy(evt, beamTrack, truthBranch, beamInst);
          }
          else {
            mf::LogWarning("CAFMaker") << "No track associated to the beam particle; reco interacting energy left unfilled.";
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
      truePart.ixn = ixnID;
      truePart.part = srID;

      truth.push_back(truePart);
      truthOverlap.push_back(eDeposit / totalEDeposit);
    }


  }

  //------------------------------------------------------------------------------


  void CAFMaker::FillInteractionTruthMatchingAndOverlap(lar_pandora::PFParticleVector const& particles, const art::Event &evt, std::vector<std::size_t> &truth, std::vector<float> &truthOverlap) const{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    //Collecting the hits of every PFParticle in the slice, since a reco interaction can be matched to more than one true interaction
    std::vector<art::Ptr<recob::Hit>> hits;
    for (const art::Ptr<recob::PFParticle>& particle : particles) {
      std::vector<art::Ptr<recob::Hit>> particleHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(particle, evt, fPandoraLabel);
      hits.insert(hits.end(), particleHits.begin(), particleHits.end());
    }

    TruthMatchUtils::IDToEDepositMap idToEDepositMap;
    for (const art::Ptr<recob::Hit>& pHit : hits){
      TruthMatchUtils::FillG4IDToEnergyDepositMap(idToEDepositMap, clockData, pHit, true);
    }

    float totalEDeposit = 0;
    for (const auto& [id, eDeposit] : idToEDepositMap) {
      totalEDeposit += eDeposit;
    }

    if (totalEDeposit <= 0) {
      mf::LogWarning("CAFMaker") << "No energy deposit found for this reco interaction. Skipping interaction-level truth matching.";
      return;
    }

    //Aggregating the energy deposit by true interaction, since several G4 tracks can belong to the same true interaction
    std::map<std::size_t, float> ixnEDeposit;
    for (const auto& [id, eDeposit] : idToEDepositMap) {
      if (fMCParticlesMap.count(id) == 0) {
        continue;
      }

      art::Ptr<simb::MCParticle> mcpart;
      bool isPrimary;
      int srID;
      int ixnID = -1;
      std::tie(mcpart, ixnID, isPrimary, srID) = fMCParticlesMap.at(id);
      if (ixnID < 0) {
        continue;
      }

      ixnEDeposit[static_cast<std::size_t>(ixnID)] += eDeposit;
    }

    for (const auto& [ixnID, eDeposit] : ixnEDeposit) {
      truth.push_back(ixnID);
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

  double CAFMaker::GetSingleHitsEnergy(art::Event const& evt, int plane) const{
    std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitLabel);

    std::vector<art::Ptr<recob::Hit>> collection_plane_hits = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, plane);

    const art::FindManyP<recob::SpacePoint> sp_assoc(hits, evt, fSpacePointLabel);

    if (!sp_assoc.isValid()) {
      mf::LogWarning("CAFMaker") << "No space points found with label '" << fSpacePointLabel << "'";
      return 0;
    }

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

    double charge = 0;

    for (uint i = 0; i < collection_plane_hits.size(); i++){
      std::vector<art::Ptr<recob::SpacePoint>> matching_sps = sp_assoc.at(collection_plane_hits[i].key());
      if(!matching_sps.empty()){ // If there are some matched spacepoints, this means that the hit is associated to some PFP and we don't want it here
        continue;
      }
      // Get the energy deposited in the hit
      charge += dune_ana::DUNEAnaHitUtils::LifetimeCorrection(clockData, detProp, collection_plane_hits[i])*collection_plane_hits[i]->Integral();
    }

    return fCalorimetryAlg.ElectronsFromADCArea(charge,2)*1./fRecombFactor/util::kGeVToElectrons;

  }

  //------------------------------------------------------------------------------

  // The three functions below port the thin-target interacting energy of protoduneana's
  // PDSPAnalyzer (true_beam_interactingEnergy and reco_beam_interactingEnergy), without the
  // protoduneana dependency. All energies are kinetic energies in MeV.

  std::map<int, std::vector<const sim::IDE*>> CAFMaker::SliceIDEs(std::vector<const sim::IDE*> const& ides, double z0, double pitch) const
  {
    std::map<int, std::vector<const sim::IDE*>> results;
    for (const sim::IDE* ide : ides) {
      int slice_num = std::floor((ide->z - (z0 - pitch/2.)) / pitch);
      results[slice_num].push_back(ide);
    }
    return results;
  }

  //------------------------------------------------------------------------------

  double CAFMaker::ComputeTrueInteractingEnergy(const simb::MCParticle &trueParticle) const
  {
    const simb::MCTrajectory &trajectory = trueParticle.Trajectory();
    const double mass = trueParticle.Mass()*1.e3; //GeV -> MeV

    //Wire positions of TPC 1, the beam-side TPC the beam enters through, as in PDSPAnalyzer
    const geo::PlaneID planeID{0, 1, static_cast<unsigned int>(CollectionPlane())};
    const double z0    = fWireReadout->Wire(geo::WireID(planeID, 0)).GetCenter().Z();
    const double pitch = fWireReadout->Plane(planeID).WirePitch();

    double init_KE = std::sqrt(1.e6*trueParticle.P()*trueParticle.P() + mass*mass) - mass;

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    auto ides = bt_serv->TrackIdToSimIDEs_Ps(trueParticle.TrackId(), geo::View_t(CollectionPlane()));
    std::sort(ides.begin(), ides.end(), [](const sim::IDE* a, const sim::IDE* b){ return a->z < b->z; });

    //Drop everything after the first gap of more than 5 cm made by an EM daughter (negative track ID)
    for (size_t i = 1; i + 1 < ides.size(); ++i) {
      if (ides[i]->trackID < 0 && (ides[i]->z - ides[i-1]->z) > 5) {
        ides.erase(ides.begin() + i, ides.end());
        break;
      }
    }

    auto sliced_ides = SliceIDEs(ides, z0, pitch);
    if (sliced_ides.empty()) return -999.;

    //Start from the trajectory energy just before the first slice with depositions
    const double first_z = sliced_ides.begin()->second.front()->z;
    for (size_t i = 1; i < trajectory.size(); ++i) {
      if (trajectory.Z(i-1) < first_z && trajectory.Z(i) > first_z) {
        init_KE = 1.e3*trajectory.E(i-1) - mass;
        break;
      }
    }

    double interactingEnergy = init_KE;
    for (auto const& [slice, sliceIDEs] : sliced_ides) {
      for (const sim::IDE* ide : sliceIDEs) interactingEnergy -= ide->energy;
    }
    return interactingEnergy;
  }

  //------------------------------------------------------------------------------

  double CAFMaker::ComputeRecoInteractingEnergy(const art::Event &evt, const art::Ptr<recob::Track> &beamTrack, const caf::SRTruthBranch &truthBranch, const caf::SRBeamInstrumentation &beamInst) const
  {
    const art::FindManyP<anab::Calorimetry> fmCalo(std::vector<art::Ptr<recob::Track>>{beamTrack}, evt, fCalorimetryLabelSCE);
    if (!fmCalo.isValid()) {
      mf::LogWarning("CAFMaker") << "No calorimetry found with label '" << fCalorimetryLabelSCE << "' for the beam track.";
      return -999.;
    }

    art::Ptr<anab::Calorimetry> calo;
    for (art::Ptr<anab::Calorimetry> const& c : fmCalo.at(0)) {
      if (static_cast<int>(c->PlaneID().Plane) == CollectionPlane()) {
        calo = c;
        break;
      }
    }
    if (calo.isNull() || calo->dEdx().empty()) return -999.;

    //(z, dE/dx, pitch) of each calorimetry point, ordered along the beam direction
    std::vector<std::tuple<double, double, double>> points;
    for (size_t i = 0; i < calo->dEdx().size(); ++i) {
      points.emplace_back(calo->XYZ()[i].Z(), calo->dEdx()[i], calo->TrkPitchVec()[i]);
    }
    std::sort(points.begin(), points.end(), [](auto const& a, auto const& b){ return std::get<0>(a) < std::get<0>(b); });

    double init_KE = 0.;
    if (evt.isRealData() || fMCHasBI) {
      //Beam instrumentation momentum, assuming a pion as in PDSPAnalyzer. NaN if no beam event was found.
      const double mass = 139.57;
      const double p = std::isfinite(beamInst.P) ? 1.e3*beamInst.P : 0.;
      init_KE = std::sqrt(p*p + mass*mass) - mass;
    }
    else {
      //Assumes the first true interaction is the beam particle, as with MCTruthLabel: [ "generator" ]
      if (truthBranch.nu.empty() || truthBranch.nu[0].prim.empty()) {
        mf::LogWarning("CAFMaker") << "No true beam particle found; reco interacting energy left unfilled.";
        return -999.;
      }
      const caf::SRLorentzVector &p = truthBranch.nu[0].prim[0].p;
      const double p2 = p.px*p.px + p.py*p.py + p.pz*p.pz;
      init_KE = (p.E - std::sqrt(std::max(0., p.E*p.E - p2)))*1.e3; //GeV -> MeV
    }

    //The last point is not subtracted, as in PDSPAnalyzer's reco_beam_incidentEnergies
    double interactingEnergy = init_KE;
    for (size_t i = 0; i + 1 < points.size(); ++i) {
      const double dEdx = std::get<1>(points[i]);
      if (dEdx < 0.) continue;
      interactingEnergy -= dEdx*std::get<2>(points[i]);
    }
    return interactingEnergy;
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
    caf::SRBeamInstrumentation &inst = beam.inst;

    const std::string &beamLabel = evt.isRealData() ? fBeamModuleLabel : fMCBeamModuleLabel;
    auto beamHandle = evt.getHandle<std::vector<beam::ProtoDUNEBeamEvent>>(beamLabel);
    if (!beamHandle || beamHandle->empty()) {
      mf::LogWarning("CAFMaker") << "No beam::ProtoDUNEBeamEvent found with label '" << beamLabel
                                 << "'; beam instrumentation left empty.";
      return;
    }

    const beam::ProtoDUNEBeamEvent &beamEvent = beamHandle->front();

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
      // fBeamInstPFix puts simulated beam instrumentation on the data momentum scale; inst.P_raw
      // is kept so the correction can be undone or varied offline.
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
      particle_record.contained = (particle_record.walldist > fContainedDistThreshold); //Setting the contained flag based on the distance to the wall

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
    single_hits.E = GetSingleHitsEnergy(evt, 2); //Using the collection plane for now
    single_hits.origRecoObjType = caf::RecoObjType::kHitCollection;

    recoParticlesBranch.pandora.push_back(std::move(single_hits));
    recoParticlesBranch.npandora++;

    // particle_record.origRecoObjType
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

    return fCalorimetryAlg.ElectronsFromADCArea(showerCharge,2)*1./fRecombFactor/util::kGeVToElectrons;

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

  void CAFMaker::FillGeneratorInfo(simb::MCGeneratorInfo const& genInfo, caf::SRTrueInteraction &inter) const
  {
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

    FillRecoInfoSliceLoop(sr.common, *fdBranch, sr.mc, sr.beam.inst, evt);

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
