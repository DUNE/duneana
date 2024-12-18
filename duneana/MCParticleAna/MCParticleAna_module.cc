////////////////////////////////////////////////////////////////////////////////////
// Class:       MCParticleAna                                                     //
// Module Type: analyzer                                                          //
// File:        MCParticleAna_module.cc                                           //
//                                                                                //
// Written by Sergio Manthey Corchado                                             //
// This module can be used to dump truth particle info for MCParticle analysis    //
////////////////////////////////////////////////////////////////////////////////////

// C++ includes
#ifndef MCParticleAna_h
#define MCParticleAna_h

// ROOT includes
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TTree.h"
#include "TRandom.h"
#include <fcntl.h>

// Larsoft includes (not all might be necessary)
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
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

namespace solar
{
  class MCParticleAna : public art::EDAnalyzer
  {
  public:
    // --- Standard constructor and destructor for an ART module.
    explicit MCParticleAna(fhicl::ParameterSet const &p);
    MCParticleAna(MCParticleAna const &) = delete;
    MCParticleAna(MCParticleAna &&) = delete;
    MCParticleAna &operator=(MCParticleAna const &) = delete;
    MCParticleAna &operator=(MCParticleAna &&) = delete;
    void analyze(art::Event const &evt) override;
    void reconfigure(fhicl::ParameterSet const &p);
    void beginJob() override;

  private:
    // --- Some of our own functions.
    void ResetVariables();

    // --- Input settings imported from the fcl
    std::string fGeometry;
    int fDetectorSizeX, fDetectorSizeY, fDetectorSizeZ, fDetectorDriftTime;
    std::vector<std::string> fParticleLabels;

    // --- Our TTrees, and its associated variables.
    TTree *fMCTruthTree;
    int Event, Flag, ParticlePDG, ParticleLabelID;
    float ParticleE, ParticleP, ParticleK, ParticleX, ParticleY, ParticleZ, ParticleEndX, ParticleEndY, ParticleEndZ, ParticleTime;
    std::vector<int> ParticleCount;
    std::vector<std::map<int, simb::MCParticle>> GeneratorParticles = {};

    // --- Declare our services
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    std::unique_ptr<solar::SolarAuxUtils> solaraux;
  };
#endif

  //......................................................
  MCParticleAna::MCParticleAna(fhicl::ParameterSet const &p)
      : EDAnalyzer(p),
        solaraux(new solar::SolarAuxUtils(p))
  {
    this->reconfigure(p);
  }

  //......................................................
  void MCParticleAna::reconfigure(fhicl::ParameterSet const &p)
  {
    fParticleLabels = p.get<std::vector<std::string>>("ParticleLabelVector");
    fGeometry = p.get<std::string>("Geometry");
    fDetectorSizeX = p.get<int>("DetectorSizeX");
    fDetectorSizeY = p.get<int>("DetectorSizeY");
    fDetectorSizeZ = p.get<int>("DetectorSizeZ");
    fDetectorDriftTime = p.get<float>("DetectorDriftTime");
  } // Reconfigure

  //......................................................
  void MCParticleAna::beginJob()
  {
    // --- Make our handle to the TFileService
    art::ServiceHandle<art::TFileService> tfs;
    fMCTruthTree = tfs->make<TTree>("MCTruthTree", "MCTruthTree");

    // Repeated Truth info.
    fMCTruthTree->Branch("Event", &Event, "Event/I");                      // Event number
    fMCTruthTree->Branch("Flag", &Flag, "Flag/I");                         // Flag used to match truth with reco tree entries
    fMCTruthTree->Branch("ParticleLabelID", &ParticleLabelID, "ParticleLabelID/I"); // Label ID for the particle. 0 = not generated, 1 = first generator, 2 = second generator, ...
    fMCTruthTree->Branch("ParticleE", &ParticleE, "ParticleE/F");          // True signal energy [MeV]
    fMCTruthTree->Branch("ParticleP", &ParticleP, "ParticleP/F");          // True signal momentum [MeV]
    fMCTruthTree->Branch("ParticleK", &ParticleK, "ParticleK/F");          // True signal K.E. [MeV] 
    fMCTruthTree->Branch("ParticleX", &ParticleX, "ParticleX/F");          // True signal X [cm]
    fMCTruthTree->Branch("ParticleY", &ParticleY, "ParticleY/F");          // True signal Y [cm]
    fMCTruthTree->Branch("ParticleZ", &ParticleZ, "ParticleZ/F");          // True signal Z [cm]
    fMCTruthTree->Branch("ParticleEndX", &ParticleEndX, "ParticleEndX/F"); // True signal EndX [cm]
    fMCTruthTree->Branch("ParticleEndY", &ParticleEndY, "ParticleEndY/F"); // True signal EndY [cm]
    fMCTruthTree->Branch("ParticleEndZ", &ParticleEndZ, "ParticleEndZ/F"); // True signal EndZ [cm]
    fMCTruthTree->Branch("ParticlePDG", &ParticlePDG, "ParticlePDG/I");    // True signal PDG
    fMCTruthTree->Branch("ParticleTime", &ParticleTime, "ParticleTime/F"); // True signal time [tick]
  } // BeginJob

  //......................................................
  void MCParticleAna::analyze(art::Event const &evt)
  {
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------- Prepare everything for new event ----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    std::vector<std::set<int>> trackids = {};
    std::map<int, simb::MCParticle> ThisGeneratorParts;

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
    std::vector<int> ParticelTypeCount = {0, 0, 0, 0, 0}; // {"alpha", "electron", "gamma", "neutron", "other"}

    std::string sMcTruth = "";
    sMcTruth = sMcTruth + "\nThere are a total of " + SolarAuxUtils::str(int(PartList.size())) + " Particles in the event\n";

    // Loop over all bkg handles and collect track IDs
    for (size_t i = 0; i < fParticleLabels.size(); i++)
    {
      GeneratorParticles.push_back(ThisGeneratorParts); // For each label insert empty list

      art::Handle<std::vector<simb::MCTruth>> ThisHandle;
      evt.getByLabel(fParticleLabels[i], ThisHandle);

      if (ThisHandle)
      {
        for (auto const &ParticleTruth : *ThisHandle)
        {
          int NParticles = ParticleTruth.NParticles();
          sMcTruth = sMcTruth + "\n# of particles " + SolarAuxUtils::str(NParticles) + "\tfrom gen " + SolarAuxUtils::str(int(i) + 1) + " " + fParticleLabels[i];
          ParticleLabelID = i+1;
          for (int j = 0; j < NParticles; j++)
          {
            const simb::MCParticle &Particle = ParticleTruth.GetParticle(j);
            ParticleE = 1e3 * Particle.E();
            ParticleP = 1e3 * Particle.P();
            ParticleK = 1e3 * Particle.E() - 1e3 * Particle.Mass();
            ParticleX = Particle.Vx();
            ParticleY = Particle.Vy();
            ParticleZ = Particle.Vz();
            ParticleEndX = Particle.EndX();
            ParticleEndY = Particle.EndY();
            ParticleEndZ = Particle.EndZ();
            ParticlePDG = Particle.PdgCode();
            ParticleTime = Particle.T();
            std::string sParticle = "";
            if (abs(Particle.PdgCode()) == 1000020040)
            {
              ParticelTypeCount[0]++;
              sParticle = "Alpha";
            }
            else if (abs(Particle.PdgCode()) == 11)
            {
              ParticelTypeCount[1]++;
              sParticle = "Electron";
            }
            else if (abs(Particle.PdgCode()) == 22)
            {
              ParticelTypeCount[2]++;
              sParticle = "Gamma";
            }
            else if (abs(Particle.PdgCode()) == 2112)
            {
              ParticelTypeCount[3]++;
              sParticle = "Neutron";
            }
            else
            {
              ParticelTypeCount[4]++;
              sParticle = "Other";
            }
            fMCTruthTree->Fill();
          }
        }
      }
      else
      {
        sMcTruth = sMcTruth + "\n# of particles " + SolarAuxUtils::str(int(GeneratorParticles[i].size())) + "\tfrom gen " + SolarAuxUtils::str(int(i) + 1) + " " + fParticleLabels[i] + " *not generated!";
        ParticleCount.push_back(0);
        std::set<int> ThisGeneratorIDs = {};
        trackids.push_back(ThisGeneratorIDs);
      }
    }
    sMcTruth = sMcTruth + "\nParticle Type Count: " + SolarAuxUtils::str(ParticelTypeCount[0]) + " Alphas, " + SolarAuxUtils::str(ParticelTypeCount[1]) + " Electrons, " + SolarAuxUtils::str(ParticelTypeCount[2]) + " Gammas, " + SolarAuxUtils::str(ParticelTypeCount[3]) + " Neutrons, " + SolarAuxUtils::str(ParticelTypeCount[4]) + " Others";
    solaraux->PrintInColor(sMcTruth, SolarAuxUtils::GetColor("bright_red"));
  } // Analyze

  // ########################################################################################################################################//
  // _FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_//
  // ########################################################################################################################################//

  //......................................................
  // Reset variables for each event
  void MCParticleAna::ResetVariables()
  {
    ParticleLabelID = 0; // Label ID for the particle. 0 = not generated, 1 = first generator, 2 = second generator, ...
    ParticleE = 0;
    ParticleK = 0;
    ParticleX = 0;
    ParticleY = 0;
    ParticleZ = 0;
    ParticleEndX = 0;
    ParticleEndY = 0;
    ParticleEndZ = 0;
    ParticlePDG = 0;
    ParticleTime = 0;
    ParticleCount.clear();
    GeneratorParticles.clear();
  }
} // namespace solar
DEFINE_ART_MODULE(solar::MCParticleAna)
