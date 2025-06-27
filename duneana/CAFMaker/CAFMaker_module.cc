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

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

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
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"


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
                         std::vector<simb::MCTruth> const& mctruth,
                         std::vector<simb::GTruth> const& gtruth,
                         std::vector<simb::MCFlux> const& flux,
                         art::Event const& evt);

      void FillMetaInfo(caf::SRDetectorMeta &meta, art::Event const& evt) const;
      void FillBeamInfo(caf::SRBeamBranch &beam, const art::Event &evt) const;
      void FillRecoInfo(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, const art::Event &evt) const;
      void FillCVNInfo(caf::SRCVNScoreBranch &cvnBranch, const art::Event &evt) const;
      void FillEnergyInfo(caf::SRNeutrinoEnergyBranch &ErecBranch, const art::Event &evt) const;
      void FillRecoParticlesInfo(caf::SRRecoParticlesBranch &recoParticlesBranch, caf::SRFD &fdBranch, const art::Event &evt) const;
      void FillDirectionInfo(caf::SRDirectionBranch &dirBranch, const art::Event &evt) const;
      int FillGENIERecord(simb::MCTruth const& mctruth, simb::GTruth const& gtruth);
      double GetVisibleEnergy(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const;
      void FillTruthMatchingAndOverlap(art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt, std::vector<TrueParticleID> &truth, std::vector<float> &truthOverlap) const;
      void FillPFPMetadata(caf::SRPFP &pfpBranch, art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const;
      double GetWallDistance(recob::PFParticle const& pfp, const art::Event &evt) const;
      double GetWallDistance(recob::SpacePoint const& sp) const;
      void ComputeActiveBounds();
      double GetSingleHitsEnergy(art::Event const& evt, int plane) const;
      bool IsVertexContained(caf::SRVector3D const& vtx) const;

      std::string fCVNLabel;
      bool fIsAtmoCVN;
      std::string fRegCNNLabel;

      std::string fMCTruthLabel;
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
      std::string fContainedDistThreshold;
      std::string fHitLabel;
      std::string fG4Label;


      std::map<int, std::tuple<art::Ptr<simb::MCParticle>, bool, int>> fMCParticlesMap; //[tid] = (MCParticle, isPrimary, SRParticle ID)
      uint fNprimaries = 0;
      uint fNsecondaries = 0;

      TTree* fTree = nullptr;
      TTree* fMetaTree = nullptr;
      TTree* fGENIETree = nullptr;

      std::unique_ptr<TFile> fFlatFile;
      TTree* fFlatTree; //Ownership will be managed directly by ROOT
      std::unique_ptr<flat::Flat<caf::StandardRecord>> fFlatRecord;

      genie::NtpMCEventRecord *fEventRecord = nullptr;

      double fMetaPOT;
      int fMetaRun, fMetaSubRun, fMetaVersion;

      const geo::Geometry* fGeom;
      std::vector<double> fActiveBounds;
      std::vector<double> fVertexFiducialVolumeCut;

      calo::CalorimetryAlg fCalorimetryAlg;                    ///< the calorimetry algorithm
      double fRecombFactor; ///< recombination factor for the isolated hits
      
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
      fMCTruthLabel(pset.get<std::string>("MCTruthLabel")),
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
      fContainedDistThreshold(pset.get< std::string >("ContainedDistThreshold")),
      fHitLabel(pset.get< std::string >("HitLabel")),
      fG4Label(pset.get< std::string >("G4Label")),
      fEventRecord(new genie::NtpMCEventRecord),
      fGeom(&*art::ServiceHandle<geo::Geometry>()),
      fVertexFiducialVolumeCut(pset.get<std::vector<double>>("VertexFiducialVolumeCut")),
      fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
      fRecombFactor(pset.get<double>("RecombFactor"))
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
    fNprimaries = 0;
    fNsecondaries = 0;
    for(art::Ptr<simb::MCParticle> const& mcpart: mcparticles) {
      if(mcpart->TrackId() == 0) continue; //Skip the neutrino particle (TrackId 0)
      bool isPrimary = (mcpart->Mother() == 0);
      if(isPrimary) {
        fMCParticlesMap[mcpart->TrackId()] = {mcpart, true, fNprimaries};
        fNprimaries++;
      } else {
        fMCParticlesMap[mcpart->TrackId()] = {mcpart, false, fNsecondaries};
        fNsecondaries++;
      } 
    }
  }

  //------------------------------------------------------------------------------

  void CAFMaker::FillTruthInfo(caf::SRTruthBranch& truthBranch,
                               std::vector<simb::MCTruth> const& mctruth,
                               std::vector<simb::GTruth> const& gtruth,
                               std::vector<simb::MCFlux> const& flux,
                               art::Event const& evt)
  {
    for(size_t i=0; i<mctruth.size(); i++){
      caf::SRTrueInteraction inter;

      inter.id = i;
      inter.genieIdx = FillGENIERecord(mctruth[i], gtruth[i]); //Filling the GENIE EventRecord tree and associating the right index here

      const simb::MCNeutrino &neutrino = mctruth[i].GetNeutrino();

      inter.pdg = neutrino.Nu().PdgCode();
      inter.pdgorig = flux[i].fntype;
      inter.iscc = !(neutrino.CCNC()); // ccnc is 0=CC 1=NC
      inter.mode = static_cast<caf::ScatteringMode>(neutrino.Mode());
      inter.targetPDG = gtruth[i].ftgtPDG;
      inter.hitnuc = neutrino.HitNuc(); //gtruth[i].fHitNucPDG seems not to be set properly so using the MCNeutrino info
      double nucMass = genie::PDGLibrary::Instance()->Find(inter.hitnuc)->Mass();
      inter.removalE = nucMass - gtruth[i].fHitNucP4.E(); //Estimating the removal energy as hitNucleonMass - hitNucleonEnergy as the actual value is not stored...
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

      TLorentzVector q = neutrino.Nu().Momentum()-neutrino.Lepton().Momentum();
      inter.q0 = q.E();
      inter.modq = q.Vect().Mag();
      inter.t = gtruth[i].fgT;
      inter.isvtxcont = IsVertexContained(inter.vtx);

      inter.ischarm = gtruth[i].fIsCharm;
      inter.isseaquark = gtruth[i].fIsSeaQuark;
      inter.resnum = gtruth[i].fResNum;
      inter.xsec = gtruth[i].fXsec;
      inter.genweight = gtruth[i].fweight;

      //TODO: To be done later when the info will be propagated/available
      // inter.baseline ///< Distance from decay to interaction [m]
      // inter.prod_vtx ///< Neutrino production vertex [cm; beam coordinates]
      // inter.parent_dcy_mom ///< Neutrino parent momentum at decay [GeV; beam coordinates]
      // inter.parent_dcy_mode ///< Parent hadron/muon decay mode
      // inter.parent_pdg ///< PDG Code of parent particle ID
      // inter.parent_dcy_E ///< Neutrino parent energy at decay [GeV]
      // inter.imp_weight ///< Importance weight from flux file

      const simb::MCGeneratorInfo &genInfo = mctruth[i].GeneratorInfo();

      std::map<simb::Generator_t, caf::Generator>::const_iterator it = fgenMap.find(genInfo.generator);
      if (it != fgenMap.end())
      {
        inter.generator = it->second;
      }
      else{
        inter.generator = caf::Generator::kUnknownGenerator;
      }

      //Parsing the GENIE version because it is stored as a vector of uint.
      size_t last = 0;
      size_t next = 0;
      std::string s(genInfo.generatorVersion);
      char delimiter = '.';
      while ((next = s.find(delimiter, last)) != string::npos){
        inter.genVersion.push_back(std::stoi(s.substr(last, next - last)));  
        last = next + 1;
      }
      inter.genVersion.push_back(std::stoi(s.substr(last)));

      //TODO: Ask to implement a map in the StandardRecord to put everything there
      // CURRENTLY DISABLED genConfigString FIELD
      // if(genInfo.generatorConfig.find("tune") != genInfo.generatorConfig.end()){
      //   inter.genConfigString = genInfo.generatorConfig.at("tune");
      // }

      inter.nproton = 0;
      inter.nneutron = 0;
      inter.npip = 0;
      inter.npim = 0;
      inter.npi0 = 0;
      inter.nprim = 0;
      inter.nprefsi = 0;
      inter.nsec = 0;

      //Looping on the particles in the MC truth to get the prefsi particles only
      for( int p = 0; p < mctruth[i].NParticles(); p++ ) {
        const simb::MCParticle &mcpart = mctruth[i].GetParticle(p);
        if( mcpart.StatusCode() != genie::EGHepStatus::kIStHadronInTheNucleus) continue; //We only fill the prefsi particles here, the primaries will be filled later

        caf::SRTrueParticle part;
        int pdg = mcpart.PdgCode();
        part.pdg = pdg;
        part.G4ID = mcpart.TrackId();
        part.interaction_id = inter.id;
        part.time = mcpart.T();
        part.p = caf::SRLorentzVector(mcpart.Momentum());
        part.start_pos = caf::SRVector3D(mcpart.Position().Vect());
        part.end_pos = caf::SRVector3D(mcpart.EndPosition().Vect());
        part.parent = mcpart.Mother();

        inter.prefsi.push_back(std::move(part));
        inter.nprefsi++;

        //start and end processes/subprocesses are currently not filled as they are strings and it's tricky to convert them back to uint
      }

      //Now filling the primaries and secondaries infos
      inter.prim.resize(fNprimaries);
      inter.nprim = fNprimaries;
      inter.sec.resize(fNsecondaries);
      inter.nsec = fNsecondaries;

      for (auto const& [tid, mcpart_tuple] : fMCParticlesMap) {
        art::Ptr<simb::MCParticle> mcpart = std::get<0>(mcpart_tuple);
        uint srID = std::get<2>(mcpart_tuple);
        bool isPrimary = (mcpart->Mother() == 0);
        SRTrueParticle part;
        part.pdg = mcpart->PdgCode();
        part.G4ID = mcpart->TrackId();
        part.interaction_id = 0; //TODO: Only considering a single interaction in the FD currently
        part.time = mcpart->T();
        part.p = SRLorentzVector(mcpart->Momentum());
        part.start_pos = SRVector3D(mcpart->Position().Vect());
        part.end_pos = SRVector3D(mcpart->EndPosition().Vect());
        part.parent = mcpart->Mother();

        for(int i = 0; i < mcpart->NumberDaughters(); i++){
          int daughter = mcpart->Daughter(i);
          part.daughters.push_back(daughter);
        }

        caf::TrueParticleID ancestor;
        if(fMCParticlesMap.count(mcpart->Mother()) > 0){
          art::Ptr<simb::MCParticle> ancestor_ptr;
          bool ancestor_is_primary;
          int ancestor_id;
          std::tie(ancestor_ptr, ancestor_is_primary, ancestor_id) = fMCParticlesMap[mcpart->Mother()];

          ancestor.type = ancestor_is_primary ? caf::TrueParticleID::kPrimary : caf::TrueParticleID::kSecondary;
          ancestor.ixn = part.interaction_id; //Same interaction as the secondary particle
          ancestor.part = ancestor_id; //The ID of the ancestor particle in the StandardRecord

          part.ancestor_id = ancestor;
        }

        //TODO: Not filling the start and end processes/subprocesses as they are strings and it's tricky to convert them back to uint
        if(isPrimary){
          inter.prim[srID] = std::move(part);

          switch(mcpart->PdgCode()){
            case 2212: // Proton
              inter.nproton++;
              break;
            case 2112: // Neutron
              inter.nneutron++;
              break;
            case 211: // Pi+
              inter.npip++;
              break;
            case -211: // Pi-
              inter.npim++;
              break;
            case 111: // Pi0
              inter.npi0++;
              break;
          }
        }
        else{
          inter.sec[srID] = std::move(part);
        }

      }

      truthBranch.nu.push_back(std::move(inter));
    } // loop through MC truth i

    truthBranch.nnu = mctruth.size();
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
        } else {
          mf::LogWarning("CAFMaker") << "Unknown property '" << property_name << "' for PFP with ID " << pfp->Self();
        }
      }

  }


  //------------------------------------------------------------------------------

  void CAFMaker::FillRecoInfo(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, const art::Event &evt) const {
    SRInteractionBranch &ixn = recoBranch.ixn;

    //Only filling with Pandora Reco for the moment
    std::vector<SRInteraction> &pandora = ixn.pandora;
    
    lar_pandora::PFParticleVector particleVector;
    lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraLabel, particleVector);
    lar_pandora::VertexVector vertexVector;
    lar_pandora::PFParticlesToVertices particlesToVertices;
    lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraLabel, vertexVector, particlesToVertices);

    for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      if(particle->IsPrimary()){
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
        FillRecoParticlesInfo(part, fdBranch, evt);

        //Assuming a single TrueInteraction for now. TODO: Change this if several interactions end up being simulated in the same event
        reco.truth = {0}; 
        reco.truthOverlap = {1.};

        pandora.emplace_back(reco);
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
      std::tie(mcpart, isPrimary, srID) = fMCParticlesMap.at(id);

      TrueParticleID truePart;
      truePart.type = isPrimary ? TrueParticleID::kPrimary : TrueParticleID::kSecondary;
      truePart.ixn = 0; //Assuming a single interaction for now
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

  double CAFMaker::GetWallDistance(recob::PFParticle const& pfp, const art::Event &evt) const{
    double minDist = std::numeric_limits<double>::max();
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaEventUtils::GetSpacePoints(evt, fSpacePointLabel);

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
 
 
  void CAFMaker::beginSubRun(const art::SubRun& sr)
  {
    art::Handle<sumdata::POTSummary> pots = sr.getHandle<sumdata::POTSummary>(fPOTSummaryLabel);
    if( pots ) fMetaPOT += pots->totpot;

    fMetaRun = sr.id().subRun();
    fMetaSubRun = sr.id().run();

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

  void CAFMaker::FillBeamInfo(caf::SRBeamBranch &beam, const art::Event &evt) const
  {
    //This part will only be relevant when working on real data with real beam.
    beam.ismc = true; //Hardcoded to true at the moment.
  }

  //------------------------------------------------------------------------------

  void CAFMaker::FillCVNInfo(caf::SRCVNScoreBranch &cvnBranch, const art::Event &evt) const
  {
    art::Handle<std::vector<cvn::Result>> cvnin = evt.getHandle<std::vector<cvn::Result>>(fCVNLabel);

    if( !cvnin.failedToGet() && !cvnin->empty()) {
      if(fIsAtmoCVN){ //Hotfix to take care of the fact that the CVN for atmospherics is storing results in a weird way...
        const std::vector<std::vector<float>> &scores = (*cvnin)[0].fOutput;
        cvnBranch.nc = scores[0][0];
        cvnBranch.nue = scores[0][1];
        cvnBranch.numu = scores[0][2];
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

  void CAFMaker::FillRecoParticlesInfo(caf::SRRecoParticlesBranch &recoParticlesBranch, caf::SRFD &fdBranch, const art::Event &evt) const
  {
    //Doing quite a lot of things here related to saving the reco particles
    //Will try to be pedagogical in the comments

    //Getting Ar density in g/cm3
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
    double lar_density = detProp.Density();

    //Getting all the PFParticles from the event
    lar_pandora::PFParticleVector particleVector;
    lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraLabel, particleVector);
    unsigned int nuID = std::numeric_limits<unsigned int>::max();
    for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      if(particle->IsPrimary()){
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

      FillTruthMatchingAndOverlap(particle, evt, particle_record.truth, particle_record.truthOverlap);
      particle_record.walldist = GetWallDistance(*particle, evt); //Getting the distance to the wall for this PFP
      particle_record.contained = (particle_record.walldist < fContainedDistThreshold); //Setting the contained flag based on the distance to the wall

      //Getting the track and shower objects associated to the PFP
      art::Ptr<recob::Track> track;
      if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(particle, evt, fPandoraLabel, fTrackLabel)){ //Unlike what its name suggets, this function only checks if an associated track exists
        track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(particle, evt, fPandoraLabel, fTrackLabel);
      }
      
      art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(particle, evt, fPandoraLabel, fShowerLabel);

      //Seeing which option Pandora prefers
      bool isTrack = lar_pandora::LArPandoraHelper::IsTrack(particle);

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
    else {
      mf::LogWarning("CAFMaker") << "Didn't detect a know geometry. Defaulting to FD HD!";
      detector = &(sr.meta.fd_hd);
      fdBranch = &(sr.fd.hd);
    }
    



    FillMetaInfo(*detector, evt);

    FillBeamInfo(sr.beam, evt);
    art::Handle<std::vector<simb::MCTruth>> mct = evt.getHandle< std::vector<simb::MCTruth> >(fMCTruthLabel);
    art::Handle<std::vector<simb::GTruth>> gt = evt.getHandle< std::vector<simb::GTruth> >(fGTruthLabel);
    art::Handle<std::vector<simb::MCFlux>> mcft = evt.getHandle< std::vector<simb::MCFlux> >(fMCFluxLabel);
    if ( !mct ) {
      mf::LogWarning("CAFMaker") << "No MCTruth. SRTruthBranch will be empty!";
    }
    else if ( !gt ) {
      mf::LogWarning("CAFMaker") << "No GTruth. SRTruthBranch will be empty!";
    }
    else if ( !mcft ) {
      mf::LogWarning("CAFMaker") << "No MCFlux. SRTruthBranch will be empty!";
    }
    else {
      FillTruthInfo(sr.mc, *mct, *gt, *mcft, evt);
    }

    FillRecoInfo(sr.common, *fdBranch, evt);

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
