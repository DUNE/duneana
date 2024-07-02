//Analyzer for flash<--->track matching (testing against MC)

#include <string>
#include <vector>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larana/OpticalDetector/OpDigiProperties.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "TH1.h"
#include "TTree.h"

namespace opreco {

  //bookkeeping on all the matches
  struct flash_match {
    std::vector<size_t> particle_indices;
    std::vector<size_t> track_indices;
  };
  struct track_match {
    std::vector<size_t> particle_indices;
    std::vector<size_t> flash_indices;
  };
  struct particle_match {
    std::vector<size_t> flash_indices;
    std::vector<size_t> track_indices;
  };

  class OpticalRecoAna : public art::EDAnalyzer {
  public:
    OpticalRecoAna(const fhicl::ParameterSet&);

    void beginJob();
    void analyze(const art::Event&);

  private:
    //fcl parameters here
    std::string fFlashModuleLabel;
    std::string fTrackModuleLabel;
    float fKineticEnergyMin;
    float fPEMin;
    float fTimeMatchMax;

    std::vector<flash_match> fFlash_match_vector;
    std::vector<track_match> fTrack_match_vector;
    std::vector<particle_match> fParticle_match_vector;

    void get_MC_particle_list(sim::ParticleList const&, std::vector<simb::MCParticle>&);
    simb::Origin_t get_MC_particle_origin(simb::MCParticle const&);
    float update_MC_particle_time(simb::MCParticle const&, bool&, geo::Geometry const&);

    //matching functions
    void match_flashes_to_tracks(std::vector<recob::OpFlash> const&,
                                 std::vector<recob::Track> const&);

    void match_flashes_to_particles(std::vector<recob::OpFlash> const&,
                                    std::vector<simb::MCParticle> const&,
                                    float const&,
                                    geo::Geometry const&);

    TH1F* fTimeDiff;
    TH1F* fTimeDiff_fine;
    TH1I* fNBeamFlashes;

    TTree* fMatchTree_PF;
    TTree* fMatchTree_PF_NotNu;
    int fRun;
    int fEvent;
    int fParticleID;
    int fParticleMother;
    int fParticleTrackID;
    float fParticleTime;
    float fParticleVx;
    float fParticleVy;
    float fParticleVz;
    int fFlashID;
    float fFlashTime;
    float fFlashTimeWidth;
    float fFlashY;
    float fFlashZ;
    float fFlashYWidth;
    float fFlashZWidth;
    float fFlashPE;
    int fFlashOnBeamTime;
    int fMatchIndex;
    int fMatchIndex_NotNu;
  };

}

namespace opreco {
  DEFINE_ART_MODULE(OpticalRecoAna)
}

// Constructor
opreco::OpticalRecoAna::OpticalRecoAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
{

  // Indicate that the Input Module comes from .fcl
  fFlashModuleLabel = pset.get<std::string>("FlashModule");
  fTrackModuleLabel = pset.get<std::string>("TrackModule");
  fKineticEnergyMin = pset.get<float>("KineticEnergyMin");
  fPEMin = pset.get<float>("PEMin");
  fTimeMatchMax = pset.get<float>("TimeMatchMax");

  fFlash_match_vector.clear();
  fTrack_match_vector.clear();
  fParticle_match_vector.clear();
}

// Do something here to setup the file (like make a TTree)
void opreco::OpticalRecoAna::beginJob()
{
  art::ServiceHandle<art::TFileService const> tfs;
  fTimeDiff = tfs->make<TH1F>(
    "htdiff",
    "Time difference between particles and flashes; t_diff (ns); flash/particle pairs",
    1e3,
    -5e6,
    5e6);
  fTimeDiff_fine = tfs->make<TH1F>(
    "htdiff_fine",
    "Time difference between particles and flashes; t_diff (ns); flash/particle pairs",
    100,
    -1000,
    1000);
  fNBeamFlashes = tfs->make<TH1I>(
    "hNBeamFlashes", "Number of flashes OnBeamTime per event; N_{Flashes}; Events", 5, 0, 5);

  fMatchTree_PF = tfs->make<TTree>("MatchTree_PF", "MatchTree_PF");
  fMatchTree_PF->Branch("Run", &fRun, "Run/I");
  fMatchTree_PF->Branch("Event", &fEvent, "Event/I");
  fMatchTree_PF->Branch("p_id", &fParticleID, "p_id/I");
  fMatchTree_PF->Branch("p_time", &fParticleTime, "p_time/F");
  fMatchTree_PF->Branch("p_vx", &fParticleVx, "p_vx/F");
  fMatchTree_PF->Branch("p_vy", &fParticleVy, "p_vy/F");
  fMatchTree_PF->Branch("p_vz", &fParticleVz, "p_vz/F");
  fMatchTree_PF->Branch("p_mid", &fParticleMother, "p_mid/I");
  fMatchTree_PF->Branch("p_trackid", &fParticleTrackID, "p_trackid/I");
  fMatchTree_PF->Branch("FlashID", &fFlashID, "flash_id/I");
  fMatchTree_PF->Branch("f_time", &fFlashTime, "f_time/F");
  fMatchTree_PF->Branch("f_timewidth", &fFlashTimeWidth, "f_timewidth/F");
  fMatchTree_PF->Branch("f_y", &fFlashY, "f_y/F");
  fMatchTree_PF->Branch("f_z", &fFlashZ, "f_z/F");
  fMatchTree_PF->Branch("f_ywidth", &fFlashYWidth, "f_ywidth/F");
  fMatchTree_PF->Branch("f_zwidth", &fFlashYWidth, "f_zwidth/F");
  fMatchTree_PF->Branch("f_onbeamtime", &fFlashOnBeamTime, "f_onbeamtime/I");
  fMatchTree_PF->Branch("f_pe", &fFlashPE, "f_pe/F");
  fMatchTree_PF->Branch("matchIndex", &fMatchIndex, "matchIndex/I");

  fMatchTree_PF_NotNu = tfs->make<TTree>("MatchTree_PF_NotNu", "MatchTree_PF_NotNu");
  fMatchTree_PF_NotNu->Branch("Run", &fRun, "Run/I");
  fMatchTree_PF_NotNu->Branch("Event", &fEvent, "Event/I");
  fMatchTree_PF_NotNu->Branch("p_id", &fParticleID, "p_id/I");
  fMatchTree_PF_NotNu->Branch("p_time", &fParticleTime, "p_time/F");
  fMatchTree_PF_NotNu->Branch("p_vx", &fParticleVx, "p_vx/F");
  fMatchTree_PF_NotNu->Branch("p_vy", &fParticleVy, "p_vy/F");
  fMatchTree_PF_NotNu->Branch("p_vz", &fParticleVz, "p_vz/F");
  fMatchTree_PF_NotNu->Branch("p_mid", &fParticleMother, "p_mid/I");
  fMatchTree_PF_NotNu->Branch("p_trackid", &fParticleTrackID, "p_trackid/I");
  fMatchTree_PF_NotNu->Branch("FlashID", &fFlashID, "flash_id/I");
  fMatchTree_PF_NotNu->Branch("f_time", &fFlashTime, "f_time/F");
  fMatchTree_PF_NotNu->Branch("f_timewidth", &fFlashTimeWidth, "f_timewidth/F");
  fMatchTree_PF_NotNu->Branch("f_y", &fFlashY, "f_y/F");
  fMatchTree_PF_NotNu->Branch("f_z", &fFlashZ, "f_z/F");
  fMatchTree_PF_NotNu->Branch("f_ywidth", &fFlashYWidth, "f_ywidth/F");
  fMatchTree_PF_NotNu->Branch("f_zwidth", &fFlashYWidth, "f_zwidth/F");
  fMatchTree_PF_NotNu->Branch("f_onbeamtime", &fFlashOnBeamTime, "f_onbeamtime/I");
  fMatchTree_PF_NotNu->Branch("f_pe", &fFlashPE, "f_pe/F");
  fMatchTree_PF_NotNu->Branch("matchIndex", &fMatchIndex_NotNu, "matchIndex/I");
}

// The analyzer itself
void opreco::OpticalRecoAna::analyze(const art::Event& evt)
{

  //get flag for MC
  const bool is_MC = !(evt.isRealData());

  fRun = evt.id().run();
  fEvent = evt.id().event();

  //first get out track, flash, and particle list handles
  art::Handle<std::vector<recob::OpFlash>> flash_handle;
  evt.getByLabel(fFlashModuleLabel, flash_handle);
  std::vector<recob::OpFlash> const& flash_vector(*flash_handle);
  fFlash_match_vector.resize(flash_vector.size());

  MF_LOG_INFO("OpticalRecoAna") << "Number of flashes is " << flash_vector.size() << std::flush;

  //all this for the MC matching
  if (is_MC) {

    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    std::vector<simb::MCParticle> particle_list;
    get_MC_particle_list(pi_serv->ParticleList(), particle_list);
    fParticle_match_vector.resize(particle_list.size());

    mf::LogInfo("OpticalRecoAna") << "Number of MC particles is " << particle_list.size();

    art::ServiceHandle<opdet::OpDigiProperties const> odp;
    const float ns_per_PMT_tick =
      1e3; // already converted to microseconds//( 1e3 / odp->SampleFreq()) ; //SampleFreq is in MHz
    art::ServiceHandle<geo::Geometry const> geometry_handle;

    match_flashes_to_particles(flash_vector, particle_list, ns_per_PMT_tick, *geometry_handle);
    int n_flashes = 0;
    for (auto const& my_flash : flash_vector) {
      if (my_flash.OnBeamTime() == 1) n_flashes++;
    }
    fNBeamFlashes->Fill(n_flashes);
  }
}

simb::Origin_t opreco::OpticalRecoAna::get_MC_particle_origin(simb::MCParticle const& particle)
{
  art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
  return (pi_serv->TrackIdToMCTruth_P(particle.TrackId()))->Origin();
}

void opreco::OpticalRecoAna::get_MC_particle_list(sim::ParticleList const& plist,
                                                  std::vector<simb::MCParticle>& particle_vector)
{

  for (sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {

    const simb::MCParticle* particle = (*ipart).second;
    //do not store if it's below our energy cut
    if (particle->E() < 0.001 * particle->Mass() + fKineticEnergyMin) continue;

    //check to see if it's a charged particle we expect to leave ionization
    const int pdg = std::abs(particle->PdgCode());
    if (pdg == 11      // electron
        || pdg == 13   //muon
        || pdg == 15   //tau
        || pdg == 211  //charged pion
        || pdg == 321  //charged kaon
        || pdg == 2212 //proton
        || pdg == 2214 //delta
        || pdg == 2224 //delta
        || pdg == 3222 //sigma
        || pdg == 3112 //sigma
        || pdg == 3312 //xi
        || pdg == 3334 //omega
    ) {
      particle_vector.emplace_back(*particle);
    }
  }
}

float opreco::OpticalRecoAna::update_MC_particle_time(simb::MCParticle const& particle,
                                                      bool& pass_check,
                                                      geo::Geometry const& geometry)
{

  pass_check = false;

  const size_t numtraj = particle.NumberTrajectoryPoints();
  size_t t_iter = 0;
  while (t_iter < numtraj) {
    try {
      // check if the particle is inside a TPC
      geo::Point_t const pos{particle.Vx(t_iter), particle.Vy(t_iter), particle.Vz(t_iter)};
      geometry.PositionToTPC(pos);
    }
    catch (cet::exception& e) {
      ++t_iter;
      continue;
    }
    break;
  }

  if (t_iter == numtraj) return 0;

  pass_check = true;
  return particle.T(t_iter);
}

void opreco::OpticalRecoAna::match_flashes_to_particles(
  std::vector<recob::OpFlash> const& flash_vector,
  std::vector<simb::MCParticle> const& particle_vector,
  float const& ns_per_PMT_tick,
  geo::Geometry const& geometry)
{
  fMatchIndex = 0;
  fMatchIndex_NotNu = 0;
  int lastFlashID = -1;

  for (size_t i_particle = 0; i_particle < particle_vector.size(); i_particle++) {

    simb::MCParticle const& my_particle(particle_vector.at(i_particle));
    bool pass_check = false;
    const float particle_time = update_MC_particle_time(my_particle, pass_check, geometry);
    if (!pass_check) continue;

    for (size_t i_flash = 0; i_flash < flash_vector.size(); i_flash++) {

      recob::OpFlash const& my_flash(flash_vector.at(i_flash));
      if (my_flash.TotalPE() < fPEMin) continue;

      const float flash_time = my_flash.Time() * ns_per_PMT_tick;

      fTimeDiff->Fill(particle_time - flash_time);
      fTimeDiff_fine->Fill(particle_time - flash_time);

      fParticleID = i_particle;
      fParticleTime = particle_time;
      fParticleVx = my_particle.Vx(0);
      fParticleVy = my_particle.Vy(0);
      fParticleVz = my_particle.Vz(0);
      fParticleMother = my_particle.Mother();
      fParticleTrackID = my_particle.TrackId();
      fFlashID = i_flash;
      fFlashTime = flash_time;
      fFlashTimeWidth = my_flash.TimeWidth() * ns_per_PMT_tick;
      fFlashY = my_flash.YCenter();
      fFlashZ = my_flash.ZCenter();
      fFlashYWidth = my_flash.YWidth();
      fFlashZWidth = my_flash.ZWidth();
      fFlashPE = my_flash.TotalPE();
      fFlashOnBeamTime = my_flash.OnBeamTime();

      bool beam_match = ((std::abs(particle_time - flash_time) / fFlashTimeWidth) < fTimeMatchMax &&
                         get_MC_particle_origin(my_particle) == simb::kBeamNeutrino);
      bool nonbeam_match =
        ((std::abs(particle_time - flash_time) / fFlashTimeWidth) < fTimeMatchMax &&
         get_MC_particle_origin(my_particle) != simb::kBeamNeutrino);

      if (beam_match) {
        fMatchTree_PF->Fill();
        fMatchIndex++;

        fFlash_match_vector.at(i_flash).particle_indices.push_back(i_particle);
        fParticle_match_vector.at(i_particle).flash_indices.push_back(i_flash);
      }
      else if (nonbeam_match) {
        if (lastFlashID != fFlashID) {
          fMatchTree_PF_NotNu->Fill();
          fMatchIndex_NotNu++;
          lastFlashID = fFlashID;
        }
      }

    } //end inner loop over particles

  } //end loop over flashes

} //end match_flashes_to_particles
