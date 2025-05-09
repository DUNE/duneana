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
#include "larcore/Geometry/Geometry.h"
#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

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
      float GetVisibleEnergy(recob::PFParticle const& pfp, const art::Event &evt) const;
      std::tuple<int, float> GetTruthMatchingAndOverlap(recob::PFParticle const& pfp, const art::Event &evt) const;
      float GetWallDistance(recob::PFParticle const& pfp, const art::Event &evt) const;
      float GetWallDistance(recob::SpacePoint const& sp) const;
      void ComputeActiveBounds();
      float GetSingleHitsEnergy(art::Event const& evt) const;

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
      std::string fEnergyRecoECaloLabel;
      std::string fDirectionRecoLabelNue;
      std::string fDirectionRecoLabelNumu;
      std::string fPandoraLabel;
      std::string fParticleIDLabel;
      std::string fPIDACut;
      std::string fTrackLabel;
      std::string fShowerLabel;
      std::string fSpacePointLabel;
      std::string fContainedDistThreshold;
      std::string fHitLabel;

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

      calo::CalorimetryAlg fCalorimetryAlg;                    ///< the calorimetry algorithm
      float fRecombFactor; ///< recombination factor for the isolated hits

      const std::map<simb::Generator_t, caf::Generator> fgenMap = {
        {simb::Generator_t::kUnknown, caf::Generator::kUnknownGenerator},
        {simb::Generator_t::kGENIE,   caf::Generator::kGENIE},
        {simb::Generator_t::kCRY,     caf::Generator::kUnknownGenerator},
        {simb::Generator_t::kGIBUU,   caf::Generator::kGIBUU},
        {simb::Generator_t::kNuWro,   caf::Generator::kUnknownGenerator},
        {simb::Generator_t::kMARLEY,  caf::Generator::kUnknownGenerator},
        {simb::Generator_t::kNEUT,    caf::Generator::kNEUT},
        {simb::Generator_t::kCORSIKA, caf::Generator::kUnknownGenerator},
        {simb::Generator_t::kGEANT,   caf::Generator::kUnknownGenerator}
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
      fEnergyRecoECaloLabel(pset.get<std::string>("EnergyRecoECaloLabel")),
      fDirectionRecoLabelNue(pset.get<std::string>("DirectionRecoLabelNue")),
      fDirectionRecoLabelNumu(pset.get<std::string>("DirectionRecoLabelNumu")),
      fPandoraLabel(pset.get< std::string >("PandoraLabel")),
      fParticleIDLabel(pset.get< std::string >("ParticleIDLabel")),
      fPIDACut(pset.get< std::string >("PIDACut")),
      fTrackLabel(pset.get< std::string >("TrackLabel")),
      fShowerLabel(pset.get< std::string >("ShowerLabel")),
      fSpacePointLabel(pset.get< std::string >("SpacePointLabel")),
      fContainedDistThreshold(pset.get< std::string >("ContainedDistThreshold")),
      fHitLabel(pset.get< std::string >("HitLabel")),
      fRecombFactor(pset.get<float>("RecombFactor")),
      fEventRecord(new genie::NtpMCEventRecord),
      fGeom(&*art::ServiceHandle<geo::Geometry>()),
      fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
  {

    if(pset.get<bool>("CreateFlatCAF")){
      // LZ4 is the fastest format to decompress. I get 3x faster loading with
      // this compared to the default, and the files are only slightly larger.
      fFlatFile = std::make_unique<TFile>("flatcaf.root", "RECREATE", "",
                            ROOT::CompressionSettings(ROOT::kLZ4, 1));
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
      inter.hitnuc = gtruth[i].fHitNucPDG;
      //TODO inter.removalE ; Not sure the info can be retrieved from Gtruth and MCTruth (at least not trivially)
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
      //TODO: inter.isvtxcont ; Not sure this is the best place to define containment

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

      //TODO: Ask to add all the generators in the StandardRecord

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

      //Filling the same fields as ND-CAFMaker. Some fields are not filled at this stage
      for( int p = 0; p < mctruth[i].NParticles(); p++ ) {
        const simb::MCParticle &mcpart = mctruth[i].GetParticle(p);
        if( mcpart.StatusCode() != genie::EGHepStatus::kIStStableFinalState
          && mcpart.StatusCode() != genie::EGHepStatus::kIStHadronInTheNucleus) continue;

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

        for(int daughterID = 0; daughterID < mcpart.NumberDaughters(); daughterID++){
          int daughter = mcpart.Daughter(daughterID);
          part.daughters.push_back(daughter);
        }

        //TODO: start and end processes/subprocesses are currently not filled

        if( mcpart.StatusCode() == genie::EGHepStatus::kIStStableFinalState )
        {
          inter.prim.push_back(std::move(part));
          inter.nprim++;

          if( pdg == 2212 ) inter.nproton++;
            else if( pdg == 2112 ) inter.nneutron++;
            else if( pdg ==  211 ) inter.npip++;
            else if( pdg == -211 ) inter.npim++;
            else if( pdg ==  111 ) inter.npi0++;
          }
          else // kIStHadronInTheNucleus
          {
            inter.prefsi.push_back(std::move(part));
            inter.nprefsi++;
        }

      }

      truthBranch.nu.push_back(std::move(inter));
    } // loop through MC truth i

    truthBranch.nnu = mctruth.size();
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
            // fData->nuvtxpdg[iv] = particle->PdgCode(); TODO: Reuse this elsewhere, add a branch in SRNeutrinoHypothesisBranch for it
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

        reco.truth = {0}; //Assuming a single TrueInteraction for now
        //TODO: not sure of what to put there..
        // std::vector<float>   truthOverlap;              ///< Fractional overlap between this reco interaction and each true interaction

        pandora.emplace_back(reco);
      }
    }

    ixn.npandora = pandora.size();
    ixn.ndlp = ixn.dlp.size();
  }


  //------------------------------------------------------------------------------

  std::tuple<int, float> CAFMaker::GetTruthMatchingAndOverlap(recob::PFParticle const& pfp, const art::Event &evt) const{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    //First getting all the hits belonging to that PFP
    std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fPandoraLabel);

    TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, true));
  }

  //------------------------------------------------------------------------------


  float CAFMaker::GetWallDistance(recob::SpacePoint const& sp) const{
    //Get the position of the space point in world coordinates
    double x = sp.XYZ()[0];
    double y = sp.XYZ()[1];
    double z = sp.XYZ()[2];

    //Get the distance to the wall
    double dist = 99999;
    for (size_t i=0; i<6; i++){
      dist = std::min(dist, std::abs(fActiveBounds[i] - x));
    }

    return dist;
  }

  //------------------------------------------------------------------------------

  float CAFMaker::GetWallDistance(recob::PFParticle const& pfp, const art::Event &evt) const{
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaEventUtils::GetSpacePoints(evt, fSpacePointLabel);

    if(spacePoints.empty()){
      mf::LogWarning("CAFMaker") << "No space points found with label '" << fSpacePointLabel << "'";
      return 99999;
    }

    //Returning the minimum distance to the wall
    return std::min_element(spacePoints.begin(), spacePoints.end(),
      [this](const art::Ptr<recob::SpacePoint> &sp1, const art::Ptr<recob::SpacePoint> &sp2){
        return GetWallDistance(*sp1) < GetWallDistance(*sp2);
      })->get();

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
  
    //TODO: Add this extra margin to the CAFMaker config for FD
    //Tweak the bounds on the y axis as there is an extra on non-instrumented 8cm on each side...
    // fActiveBounds[2] += 8;
    // fActiveBounds[3] -= 8;
  }

  //------------------------------------------------------------------------------

  float GetSingleHitsEnergy(art::Event const& evt) const{
    std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitLabel);

    std::vector<art::Ptr<recob::Hit>> collection_plane_hits = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, 2);

    std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaHitUtils::GetSpacePoints(collection_plane_hits, evt, fSpacePointLabel);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

    float charge = 0;

    for (uint i = 0; i < collection_plane_hits.size(); i++){
      if(spacePoints[i]->IsValid()){ // If the spacepoint is valid, this means that the hit is associated to some PFP and we don't want it here
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
      {fEnergyRecoECaloLabel, &(ErecBranch.e_calo)}
    };

    for(auto [label, record] : ereco_map){
       art::Handle<dune::EnergyRecoOutput> ereco = evt.getHandle<dune::EnergyRecoOutput>(label);
       if(ereco.failedToGet()){
        mf::LogWarning("CAFMaker") << label << " does not correspond to a valid EnergyRecoOutput product";
       }
       else{
        *record = ereco->fNuLorentzVector.E();
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
    art::Handle<std::vector<anab::ParticleID>> PIDs = evt.getHandle<std::vector<anab::ParticleID>>(fParticleIDLabel);
    if(!PIDs.isValid()){
      mf::LogWarning("CAFMaker") << "No ParticleID found with label '" << fParticleIDLabel << "'";
    }

    //Creating a FindManyP object to link the tracks to the PIDs
    const art::FindManyP<recob::Track, anab::ParticleID> fmPID(PIDs, evt, fParticleIDLabel);

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

      //TODO: Pain/Usefulness ratio too bad to implement truth matching for now
      // particle_record.truth;
      // particle_record.truthOverlap;
      particle_record.walldist = GetWallDistance(*particle, evt); //Getting the distance to the wall for this PFP
      particle_record.contained = (particle_record.walldist < fContainedDistThreshold); //Setting the contained flag based on the distance to the wall

      //Getting the track and shower objects associated to the PFP
      art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(particle, evt, fPandoraLabel, fTrackLabel);
      art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(particle, evt, fPandoraLabel, fShowerLabel);

      //Seeing which option Pandora prefers
      bool isTrack = dune_ana::DUNEAnaPFParticleUtils::IsTrack(particle, evt, fPandoraLabel, fTrackLabel);

      //For every PFP we create a track and a shower object and save it, independently of the existence of a track object to keep the PFP/Track/Shower parallel indexing
      SRTrack srtrack;
      SRShower srshower;

      //We define Evis at the PFP level as the hits are the same for shower and track
      float Evis = GetVisibleEnergy(*particle, evt);

      //This variable will be updated correctly during the recob::Track processing and will be used to fill the reco particle energy method if isTrack.
      caf::PartEMethod  = caf::PartEMethod::kUnknownMethod;

      if(track){
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

        if(particle_record.contained){ //If the particle is contained we apply the range method to determine its energy

        }

        //TODO: I would prefer to use some unified module that the user can setup and that will decide how to compute the energy rather than making a specific choice here
        //Putting Evis as placeholder to not confuse the user too much
        srtrack.E = Evis;
        trackErecoMethod = caf::PartEMethod::kCalorimetry; //Using the visible energy of the PFP

        //TODO: Pain/Usefulness ratio too bad to implement truth matching for now
        //srtrack.truth
        //srtrack.truthOverlap

        //Filling the PFP score with the PIDA score computed for the associated track
        if(fmPID.isValid()){
          art::Ptr<anab::ParticleID> pid = fmPID.at(track.key());
          if(pid.isValid()){
            const std::vector<anab::sParticleIDAlgScores> pScores = pid->ParticleIDAlgScores();
            for(const anab::sParticleIDAlgScores &pScore : pScores){ //Several scores are saved for different assumptions
              if(pScore.fAssumedPdg == 0){ //PIDA score is when there is no assumed pdf
                particle_record.score = pScore.fScore;
                break;
              }
            }
          }
        }
      }

      if(shower){
        srshower.start.SetX(shower->ShowerStart().X());
        srshower.start.SetY(shower->ShowerStart().Y());
        srshower.start.SetZ(shower->ShowerStart().Z());

        srshower.direction.SetX(shower->Direction().X());
        srshower.direction.SetY(shower->Direction().Y());
        srshower.direction.SetZ(shower->Direction().Z());
        
        srshower.Evis = Evis; //Using the visible energy of the PFP

        //TODO: Pain/Usefulness ratio too bad to implement truth matching for now
        //srshower.truth
        //srshower.truthOverlap
      }

      if(isTrack){
        if(track){ //I hope this condition is always fullfilled is the particle is tagged at track, but who knows...
          particle_record.start = SRVector3D(track->Start().X(), track->Start().Y(), track->Start().Z());
          particle_record.end = SRVector3D(track->End().X(), track->End().Y(), track->End().Z());
          particle_record.E = srtrack.E;
          particle_record.E_method = trackErecoMethod;
        }
        particle_record.origRecoObjType = caf::RecoObjectType::kTrack;
      }
      else{
        if(shower){ //I hope this condition is always fullfilled is the particle is tagged at shower, but who knows...
          particle_record.start = SRVector3D(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
          //Only filling the start, no defined end for a shower

          particle_record.E = srshower.Evis; //Using the visible energy of the PFP
          particle_record.E_method = caf::PartEMethod::kCalorimetry;
          
        }
        particle_record.origRecoObjType = caf::RecoObjectType::kShower;
      }

      //Saving the track record
      fdIxn.tracks.push_back(std::move(srtrack));
      fdIxn.ntracks++;

      //Saving the shower record
      fdIxn.showers.push_back(std::move(srshower));
      fdIxn.nshowers++;
        
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
    single_hits.E = GetSingleHitsEnergy(evt);
    single_hits.origRecoObjType = caf::RecoObjectType::kHitCollection;

    recoParticlesBranch.pandora.push_back(std::move(single_hits));
    recoParticlesBranch.npandora++;

    // particle_record.origRecoObjType
  }


  //------------------------------------------------------------------------------


  float CAFMaker::GetVisibleEnergy(recob::PFParticle const& pfp, const art::Event &evt) const
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

    return showerCharge;

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
