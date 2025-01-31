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
      void FillRecoInfo(caf::SRCommonRecoBranch &recoBranch, const art::Event &evt) const;
      void FillCVNInfo(caf::SRCVNScoreBranch &cvnBranch, const art::Event &evt) const;
      void FillEnergyInfo(caf::SRNeutrinoEnergyBranch &ErecBranch, const art::Event &evt) const;
      void FillRecoParticlesInfo(caf::SRRecoParticlesBranch &recoParticlesBranch, const art::Event &evt) const;
      void FillDirectionInfo(caf::SRDirectionBranch &dirBranch, const art::Event &evt) const;
      int FillGENIERecord(simb::MCTruth const& mctruth, simb::GTruth const& gtruth);

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
      std::string fPandoraNuVertexModuleLabel;

      TTree* fTree = nullptr;
      TTree* fMetaTree = nullptr;
      TTree* fGENIETree = nullptr;

      std::unique_ptr<TFile> fFlatFile;
      TTree* fFlatTree; //Ownership will be managed directly by ROOT
      std::unique_ptr<flat::Flat<caf::StandardRecord>> fFlatRecord;

      genie::NtpMCEventRecord *fEventRecord = nullptr;

      double fMetaPOT;
      int fMetaRun, fMetaSubRun, fMetaVersion;

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
      fPandoraNuVertexModuleLabel(pset.get< std::string >("PandoraNuVertexModuleLabel")),
      fEventRecord(new genie::NtpMCEventRecord)
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

  void CAFMaker::FillRecoInfo(caf::SRCommonRecoBranch &recoBranch, const art::Event &evt) const {
    SRInteractionBranch &ixn = recoBranch.ixn;

    //Only filling with Pandora Reco for the moment
    std::vector<SRInteraction> &pandora = ixn.pandora;
    
    lar_pandora::PFParticleVector particleVector;
    lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraNuVertexModuleLabel, particleVector);
    lar_pandora::VertexVector vertexVector;
    lar_pandora::PFParticlesToVertices particlesToVertices;
    lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraNuVertexModuleLabel, vertexVector, particlesToVertices);

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
        FillRecoParticlesInfo(part, evt);

        reco.truth = {0}; //Assuming a single TrueInteraction for now
        //TODO, not sure of what to put there..
        // std::vector<float>   truthOverlap;              ///< Fractional overlap between this reco interaction and each true interaction

        pandora.emplace_back(reco);
      }
    }

    ixn.npandora = pandora.size();
    ixn.ndlp = ixn.dlp.size();
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

  void CAFMaker::FillRecoParticlesInfo(caf::SRRecoParticlesBranch &recoParticlesBranch, const art::Event &evt) const
  {
    lar_pandora::PFParticleVector particleVector;
    lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraNuVertexModuleLabel, particleVector);
    unsigned int nuID = std::numeric_limits<unsigned int>::max();
    for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      if(particle->IsPrimary()){
        nuID = particle->Self(); //Finding the ID of neutrino's particle
        break;
      }
    }

    for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      if(particle->Self() == nuID){ //Discarding the neutrino
        continue;
      }

      caf::SRRecoParticle particle_record;
      particle_record.primary = (particle->Parent() == nuID);
      //For now PDG is taken from the PFP only which does not include PID info beyond track/shower
      particle_record.pdg = particle->PdgCode();

      //TODO: Add the energy reco information. Requires some uniform framework and specific choice about the method to use.

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
    

    if(fTree){
      fTree->SetBranchAddress("rec", &psr);
    }

    art::ServiceHandle<geo::Geometry const> fGeometry;
    std::string geoName = fGeometry->DetectorName();

    mf::LogInfo("CAFMaker") << "Geo name is: " << geoName;

    SRDetectorMeta *detector;

    if(geoName.find("dunevd10kt") != std::string::npos){
      detector = &(sr.meta.fd_vd);
      mf::LogInfo("CAFMaker") << "Assuming the FD VD detector";
    }
    else if (geoName.find("dune10kt") != std::string::npos)
    {
      detector = &(sr.meta.fd_hd);
      mf::LogInfo("CAFMaker") << "Assuming the FD HD detector";
    }
    else {
      mf::LogWarning("CAFMaker") << "Didn't detect a know geometry. Defaulting to FD HD!";
      detector = &(sr.meta.fd_hd);
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

    FillRecoInfo(sr.common, evt);

    //TODO -> Coordinate with sim/reco to see what to put there
    //SRFDBranch

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
