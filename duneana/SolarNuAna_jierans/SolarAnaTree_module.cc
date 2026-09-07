////////////////////////////////////////////////////////////////////////
// Class:       SolarAnaTree
// Plugin Type: analyzer (Unknown Unknown)
// File:        SolarAnaTree_module.cc
//
// Generated at Fri Aug 30 14:50:19 2024 by jierans using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>

#include <map>
#include <vector>

#define INVALID -99999

namespace duneana {
class SolarAnaTree;

struct EventDataBuffer {
  int event;
  int run;
  int subrun;

  void branch_on(TTree *tree) {
    tree->Branch("event", &event, "event/i");
    tree->Branch("run", &run, "run/i");
    tree->Branch("subrun", &subrun, "subrun/i");
  }

  void clear() {
    event = -1;
    run = -1;
    subrun = -1;
  }
};

// One row per recob::Track. trackID is a synthetic key (unique within the
// event, across all track producers) shared with HitBuffer::trackID so that
// hits associated to this track can be joined back to it.
struct TrackBuffer {
  int trackID;
  int track_id;
  std::string producer;

  int n_trajectory_points;
  double length;
  double chi2;
  int ndof;
  int particle_id;

  double start_x, start_y, start_z;
  double end_x, end_y, end_z;
  double start_dir_x, start_dir_y, start_dir_z;
  double end_dir_x, end_dir_y, end_dir_z;

  // Per-trajectory-point information, one entry per index in
  // [0, n_trajectory_points). point_valid flags points excluded from the fit
  // (e.g. by outlier rejection); their positions are still stored.
  std::vector<double> point_x, point_y, point_z;
  std::vector<int> point_valid;

  // Sums of charge (summed ADC and integral) over all recob::Hits
  // associated to this track, split out per view.
  double sadc_U, sadc_V, sadc_Z;
  double integral_U, integral_V, integral_Z;

  void branch_on(TTree *tree) {
    tree->Branch("trackID", &trackID);
    tree->Branch("track_id", &track_id);
    tree->Branch("producer", &producer);
    tree->Branch("n_trajectory_points", &n_trajectory_points);
    tree->Branch("length", &length);
    tree->Branch("chi2", &chi2);
    tree->Branch("ndof", &ndof);
    tree->Branch("particle_id", &particle_id);
    tree->Branch("start_x", &start_x);
    tree->Branch("start_y", &start_y);
    tree->Branch("start_z", &start_z);
    tree->Branch("end_x", &end_x);
    tree->Branch("end_y", &end_y);
    tree->Branch("end_z", &end_z);
    tree->Branch("start_dir_x", &start_dir_x);
    tree->Branch("start_dir_y", &start_dir_y);
    tree->Branch("start_dir_z", &start_dir_z);
    tree->Branch("end_dir_x", &end_dir_x);
    tree->Branch("end_dir_y", &end_dir_y);
    tree->Branch("end_dir_z", &end_dir_z);
    tree->Branch("point_x", &point_x);
    tree->Branch("point_y", &point_y);
    tree->Branch("point_z", &point_z);
    tree->Branch("point_valid", &point_valid);
    tree->Branch("sadc_U", &sadc_U);
    tree->Branch("sadc_V", &sadc_V);
    tree->Branch("sadc_Z", &sadc_Z);
    tree->Branch("integral_U", &integral_U);
    tree->Branch("integral_V", &integral_V);
    tree->Branch("integral_Z", &integral_Z);
  }

  void from_track(const recob::Track &track, int key,
                   const std::string &prod) {
    trackID = key;
    track_id = track.ID();
    producer = prod;
    n_trajectory_points = track.NumberTrajectoryPoints();
    length = track.Length();
    chi2 = track.Chi2();
    ndof = track.Ndof();
    particle_id = track.ParticleId();
    start_x = track.Start().X();
    start_y = track.Start().Y();
    start_z = track.Start().Z();
    end_x = track.End().X();
    end_y = track.End().Y();
    end_z = track.End().Z();
    start_dir_x = track.StartDirection().X();
    start_dir_y = track.StartDirection().Y();
    start_dir_z = track.StartDirection().Z();
    end_dir_x = track.EndDirection().X();
    end_dir_y = track.EndDirection().Y();
    end_dir_z = track.EndDirection().Z();

    point_x.clear();
    point_y.clear();
    point_z.clear();
    point_valid.clear();
    point_x.reserve(n_trajectory_points);
    point_y.reserve(n_trajectory_points);
    point_z.reserve(n_trajectory_points);
    point_valid.reserve(n_trajectory_points);
    for (size_t i = 0; i < track.NumberTrajectoryPoints(); i++) {
      recob::Track::Point_t const &pos = track.LocationAtPoint(i);
      point_x.push_back(pos.X());
      point_y.push_back(pos.Y());
      point_z.push_back(pos.Z());
      point_valid.push_back(track.HasValidPoint(i) ? 1 : 0);
    }
  }

  // Sums HitSummedADC()/Integral() over the given hits (typically the hits
  // associated to this track), split out per view (U, V, Z).
  void sum_hit_charges(const std::vector<art::Ptr<recob::Hit>> &hits) {
    sadc_U = sadc_V = sadc_Z = 0;
    integral_U = integral_V = integral_Z = 0;
    for (art::Ptr<recob::Hit> const &hit : hits) {
      switch (hit->View()) {
      case geo::kU:
        sadc_U += hit->HitSummedADC();
        integral_U += hit->Integral();
        break;
      case geo::kV:
        sadc_V += hit->HitSummedADC();
        integral_V += hit->Integral();
        break;
      case geo::kZ: // == geo::kW
        sadc_Z += hit->HitSummedADC();
        integral_Z += hit->Integral();
        break;
      default:
        break;
      }
    }
  }
};

// Geometric info about a channel, looked up via geo::WireReadout the same
// way TriggerAnaTree does for its trigger primitives.
struct ChannelInfo {
  unsigned int rop_id;
  unsigned int tpcset_id;
  int view;
  double wire_x, wire_y, wire_z;
};

// One row per recob::Hit. trackID matches TrackBuffer::trackID for the track
// this hit is associated with, or -1 if the hit is not associated to any
// track.
struct HitBuffer {
  int trackID;
  std::string producer;

  unsigned int channel;
  int view;
  int cryostat, tpc, plane, wire;

  unsigned int readout_plane_id;
  unsigned int tpcset_id;
  double wire_x, wire_y, wire_z;

  int start_tick, end_tick;
  double peak_time, sigma_peak_time;
  double rms;
  double peak_amplitude, sigma_peak_amplitude;
  double summed_adc, integral, sigma_integral;
  short multiplicity, local_index;
  double goodness_of_fit;
  int ndf;

  void branch_on(TTree *tree) {
    tree->Branch("trackID", &trackID);
    tree->Branch("producer", &producer);
    tree->Branch("channel", &channel);
    tree->Branch("view", &view);
    tree->Branch("cryostat", &cryostat);
    tree->Branch("tpc", &tpc);
    tree->Branch("plane", &plane);
    tree->Branch("wire", &wire);
    tree->Branch("readout_plane_id", &readout_plane_id);
    tree->Branch("tpcset_id", &tpcset_id);
    tree->Branch("wire_x", &wire_x);
    tree->Branch("wire_y", &wire_y);
    tree->Branch("wire_z", &wire_z);
    tree->Branch("start_tick", &start_tick);
    tree->Branch("end_tick", &end_tick);
    tree->Branch("peak_time", &peak_time);
    tree->Branch("sigma_peak_time", &sigma_peak_time);
    tree->Branch("rms", &rms);
    tree->Branch("peak_amplitude", &peak_amplitude);
    tree->Branch("sigma_peak_amplitude", &sigma_peak_amplitude);
    tree->Branch("summed_adc", &summed_adc);
    tree->Branch("integral", &integral);
    tree->Branch("sigma_integral", &sigma_integral);
    tree->Branch("multiplicity", &multiplicity);
    tree->Branch("local_index", &local_index);
    tree->Branch("goodness_of_fit", &goodness_of_fit);
    tree->Branch("ndf", &ndf);
  }

  void from_hit(const recob::Hit &hit, int track_key,
                const std::string &prod, const ChannelInfo &chinfo) {
    trackID = track_key;
    producer = prod;
    channel = hit.Channel();
    view = hit.View();
    cryostat = hit.WireID().Cryostat;
    tpc = hit.WireID().TPC;
    plane = hit.WireID().Plane;
    wire = hit.WireID().Wire;
    readout_plane_id = chinfo.rop_id;
    tpcset_id = chinfo.tpcset_id;
    wire_x = chinfo.wire_x;
    wire_y = chinfo.wire_y;
    wire_z = chinfo.wire_z;
    start_tick = hit.StartTick();
    end_tick = hit.EndTick();
    peak_time = hit.PeakTime();
    sigma_peak_time = hit.SigmaPeakTime();
    rms = hit.RMS();
    peak_amplitude = hit.PeakAmplitude();
    sigma_peak_amplitude = hit.SigmaPeakAmplitude();
    summed_adc = hit.HitSummedADC();
    integral = hit.Integral();
    sigma_integral = hit.SigmaIntegral();
    multiplicity = hit.Multiplicity();
    local_index = hit.LocalIndex();
    goodness_of_fit = hit.GoodnessOfFit();
    ndf = hit.DegreesOfFreedom();
  }
};

} // namespace duneana

class duneana::SolarAnaTree : public art::EDAnalyzer {
public:
  explicit SolarAnaTree(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SolarAnaTree(SolarAnaTree const &) = delete;
  SolarAnaTree(SolarAnaTree &&) = delete;
  SolarAnaTree &operator=(SolarAnaTree const &) = delete;
  SolarAnaTree &operator=(SolarAnaTree &&) = delete;

  // Required functions.
  void beginJob() override;
  void analyze(art::Event const &e) override;

private:
  art::ServiceHandle<art::TFileService> tfs;
  EventDataBuffer ev_buf;

  std::unordered_map<int, int> trkId_to_truthBlockId;

  bool dump_mctruths;
  TTree *mctruth_tree;
  int mctruth_pdg;
  std::string mctruth_process;
  int mctruth_status, mctruth_id, mctruth_trackid;
  std::string mctruth_gen_name;
  double mctruth_x, mctruth_y, mctruth_z, mctruth_t;
  double mctruth_Px, mctruth_Py, mctruth_Pz, mctruth_P;
  double mctruth_en, mctruth_ek;

  TTree *mcneutrino_tree;
  int mcneutrino_nupdg, mcneutrino_leptonpdg, mcneutrino_ccnc, mcneutrino_mode,
      mcneutrino_iteractionType, mcneutrino_target, mcneutrino_hitnuc,
      mcneutrino_hitquark;
  double mcneutrino_w, mcneutrino_x, mcneutrino_y, mcneutrino_qsqr,
      mcneutrino_pt, mcneutrino_theta;

  bool dump_mcparticles;
  TTree *mcparticle_tree;
  int mcparticle_pdg;
  std::string mcparticle_process;
  int mcparticle_status, mcparticle_trackid, mcparticle_truthid,
      mcparticle_mother;
  std::string mcparticle_gen_name;
  double mcparticle_x, mcparticle_y, mcparticle_z, mcparticle_t;
  double mcparticle_end_x, mcparticle_end_y, mcparticle_end_z, mcparticle_end_t;
  double mcparticle_Px, mcparticle_Py, mcparticle_Pz;
  double mcparticle_en, mcparticle_ek;

  bool dump_recotracks;
  TTree *track_tree;
  TrackBuffer track_buf;

  bool dump_recohits;
  TTree *hit_tree;
  HitBuffer hit_buf;

  const geo::WireReadoutGeom *fWireReadout;
  ChannelInfo get_channel_info_for_channel(raw::ChannelID_t channel,
                                            const geo::WireID &wireid);
};

duneana::SolarAnaTree::SolarAnaTree(fhicl::ParameterSet const &p)
    : EDAnalyzer{p}, dump_mctruths(p.get<bool>("dump_mctruths", true)),
      dump_mcparticles(p.get<bool>("dump_mcparticles", true)),
      dump_recotracks(p.get<bool>("dump_recotracks", true)),
      dump_recohits(p.get<bool>("dump_recohits", true)) {}

void duneana::SolarAnaTree::beginJob() {
  if (dump_mctruths) {
    mctruth_tree = tfs->make<TTree>("mctruths", "mctruths");
    ev_buf.branch_on(mctruth_tree);
    mctruth_tree->Branch("block_id", &mctruth_id);
    mctruth_tree->Branch("truth_track_id", &mctruth_trackid);
    mctruth_tree->Branch("pdg", &mctruth_pdg);
    mctruth_tree->Branch("generator_name", &mctruth_gen_name);
    mctruth_tree->Branch("status_code", &mctruth_status);
    mctruth_tree->Branch("x", &mctruth_x);
    mctruth_tree->Branch("y", &mctruth_y);
    mctruth_tree->Branch("z", &mctruth_z);
    mctruth_tree->Branch("t", &mctruth_t);
    mctruth_tree->Branch("px", &mctruth_Px);
    mctruth_tree->Branch("py", &mctruth_Py);
    mctruth_tree->Branch("pz", &mctruth_Pz);
    mctruth_tree->Branch("p", &mctruth_P);
    mctruth_tree->Branch("energy", &mctruth_en);
    mctruth_tree->Branch("kinetic_energy", &mctruth_ek);
    mctruth_tree->Branch("process", &mctruth_process);

    mcneutrino_tree = tfs->make<TTree>("mcneutrinos", "mcneutrinos");
    ev_buf.branch_on(mcneutrino_tree);

    mcneutrino_tree->Branch("block_id", &mctruth_id);
    mcneutrino_tree->Branch("generator_name", &mctruth_gen_name);
    mcneutrino_tree->Branch("nupdg", &mcneutrino_nupdg);
    mcneutrino_tree->Branch("leptonpdg", &mcneutrino_leptonpdg);
    mcneutrino_tree->Branch("ccnc", &mcneutrino_ccnc);
    mcneutrino_tree->Branch("mode", &mcneutrino_mode);
    mcneutrino_tree->Branch("interactionType", &mcneutrino_iteractionType);
    mcneutrino_tree->Branch("target", &mcneutrino_target);
    mcneutrino_tree->Branch("hitnuc", &mcneutrino_hitnuc);
    mcneutrino_tree->Branch("hitquark", &mcneutrino_hitquark);
    mcneutrino_tree->Branch("w", &mcneutrino_w);
    mcneutrino_tree->Branch("x", &mcneutrino_x);
    mcneutrino_tree->Branch("y", &mcneutrino_y);
    mcneutrino_tree->Branch("qsqr", &mcneutrino_qsqr);
    mcneutrino_tree->Branch("pt", &mcneutrino_pt);
    mcneutrino_tree->Branch("theta", &mcneutrino_theta);
  }
  if (dump_mcparticles) {
    mcparticle_tree = tfs->make<TTree>("mcparticles", "mcparticles");
    ev_buf.branch_on(mcparticle_tree);
    mcparticle_tree->Branch("pdg", &mcparticle_pdg);
    mcparticle_tree->Branch("generator_name", &mcparticle_gen_name);
    mcparticle_tree->Branch("status_code", &mcparticle_status);
    mcparticle_tree->Branch("g4_track_id", &mcparticle_trackid);
    mcparticle_tree->Branch("mother", &mcparticle_mother);
    mcparticle_tree->Branch("truth_block_id", &mcparticle_truthid);
    mcparticle_tree->Branch("x", &mcparticle_x);
    mcparticle_tree->Branch("y", &mcparticle_y);
    mcparticle_tree->Branch("z", &mcparticle_z);
    mcparticle_tree->Branch("t", &mcparticle_t);
    mcparticle_tree->Branch("end_x", &mcparticle_end_x);
    mcparticle_tree->Branch("end_y", &mcparticle_end_y);
    mcparticle_tree->Branch("end_z", &mcparticle_end_z);
    mcparticle_tree->Branch("end_t", &mcparticle_end_t);
    mcparticle_tree->Branch("px", &mcparticle_Px);
    mcparticle_tree->Branch("py", &mcparticle_Py);
    mcparticle_tree->Branch("pz", &mcparticle_Pz);
    mcparticle_tree->Branch("energy", &mcparticle_en);
    mcparticle_tree->Branch("kinetic_energy", &mcparticle_ek);
    mcparticle_tree->Branch("process", &mcparticle_process);
  }
  if (dump_recotracks) {
    track_tree = tfs->make<TTree>("tracks", "tracks");
    ev_buf.branch_on(track_tree);
    track_buf.branch_on(track_tree);
  }
  if (dump_recohits) {
    hit_tree = tfs->make<TTree>("hits", "hits");
    ev_buf.branch_on(hit_tree);
    hit_buf.branch_on(hit_tree);
  }
}

void duneana::SolarAnaTree::analyze(art::Event const &e) {
  ev_buf.run = e.run();
  ev_buf.subrun = e.subRun();
  ev_buf.event = e.event();

  fWireReadout = &art::ServiceHandle<geo::WireReadout>()->Get();

  if (dump_mctruths) {
    std::vector<art::Handle<std::vector<simb::MCTruth>>> mctruthHandles =
        e.getMany<std::vector<simb::MCTruth>>();

    int truth_block_counter = 0;
    trkId_to_truthBlockId.clear();

    for (auto const &mctruthHandle : mctruthHandles) {
      // Extract the generator name from the truth handle input label
      std::string generator_name =
          mctruthHandle.provenance()->inputTag().label();

      mctruth_id = truth_block_counter;
      // NOTE: here we are making an assumption that the geant4 stage's process
      // name is largeant. This should be safe mostly.
      art::FindManyP<simb::MCParticle> assns(mctruthHandle, e, "largeant");
      for (size_t i = 0; i < mctruthHandle->size(); i++) {
        const simb::MCTruth &truthblock =
            *art::Ptr<simb::MCTruth>(mctruthHandle, i);

        std::vector<art::Ptr<simb::MCParticle>> matched_mcparts = assns.at(i);
        for (art::Ptr<simb::MCParticle> mcpart : matched_mcparts) {
          trkId_to_truthBlockId[mcpart->TrackId()] = truth_block_counter;
        }

        if (truthblock.NeutrinoSet()) {
          const simb::MCNeutrino &mcneutrino = truthblock.GetNeutrino();
          mcneutrino_nupdg = mcneutrino.Nu().PdgCode();
          mcneutrino_leptonpdg = mcneutrino.Lepton().PdgCode();
          mcneutrino_ccnc = mcneutrino.CCNC();
          mcneutrino_mode = mcneutrino.Mode();
          mcneutrino_iteractionType = mcneutrino.InteractionType();
          mcneutrino_target = mcneutrino.Target();
          mcneutrino_hitnuc = mcneutrino.HitNuc();
          mcneutrino_hitquark = mcneutrino.HitQuark();
          mcneutrino_w = mcneutrino.W();
          mcneutrino_x = mcneutrino.X();
          mcneutrino_y = mcneutrino.Y();
          mcneutrino_qsqr = mcneutrino.QSqr();
          mcneutrino_pt = mcneutrino.Pt();
          mcneutrino_theta = mcneutrino.Theta();
          mcneutrino_tree->Fill();
        }

        int nparticles = truthblock.NParticles();

        for (int ipart = 0; ipart < nparticles; ipart++) {
          const simb::MCParticle &part = truthblock.GetParticle(ipart);
          mctruth_pdg = part.PdgCode();
          mctruth_gen_name = generator_name;
          mctruth_status = part.StatusCode();
          mctruth_process = part.Process();
          mctruth_trackid = part.TrackId();
          mctruth_x = part.Vx();
          mctruth_y = part.Vy();
          mctruth_z = part.Vz();
          mctruth_t = part.T();
          mctruth_Px = part.Px();
          mctruth_Py = part.Py();
          mctruth_Pz = part.Pz();
          mctruth_P = part.P();
          mctruth_en = part.E();
          mctruth_ek = part.E() - part.Mass();
          mctruth_tree->Fill();
        }
        truth_block_counter++;
      }
    }
  }

  if (dump_mcparticles) {
    std::vector<art::Handle<std::vector<simb::MCParticle>>> mcparticleHandles =
        e.getMany<std::vector<simb::MCParticle>>();

    for (auto const &mcparticleHandle : mcparticleHandles) {
      std::string generator_name =
          mcparticleHandle.provenance()->inputTag().label();

      for (const simb::MCParticle &part : *mcparticleHandle) {
        mcparticle_pdg = part.PdgCode();
        mcparticle_gen_name = generator_name;
        mcparticle_status = part.StatusCode();
        mcparticle_trackid = part.TrackId();
        mcparticle_mother = part.Mother();
        mcparticle_truthid =
            dump_mctruths ? trkId_to_truthBlockId.at(part.TrackId()) : -1;
        mcparticle_process = part.Process();
        mcparticle_x = part.Vx();
        mcparticle_y = part.Vy();
        mcparticle_z = part.Vz();
        mcparticle_t = part.T();
        mcparticle_end_x = part.EndX();
        mcparticle_end_y = part.EndY();
        mcparticle_end_z = part.EndZ();
        mcparticle_end_t = part.EndT();
        mcparticle_Px = part.Px();
        mcparticle_Py = part.Py();
        mcparticle_Pz = part.Pz();
        mcparticle_en = part.E();
        mcparticle_ek = part.E() - part.Mass();
        mcparticle_tree->Fill();
      }
    }
  }

  // Maps each associated recob::Hit to the synthetic trackID of the track it
  // belongs to. Populated below when dumping tracks, and consulted when
  // dumping hits so that hits carry a shared trackID pointing back to their
  // track (or -1 if unassociated).
  std::map<art::Ptr<recob::Hit>, int> hit_to_trackID;

  if (dump_recotracks) {
    std::vector<art::Handle<std::vector<recob::Track>>> trackHandles =
        e.getMany<std::vector<recob::Track>>();

    int track_key = 0;
    for (auto const &trackHandle : trackHandles) {
      std::string producer = trackHandle.provenance()->inputTag().encode();
      art::FindManyP<recob::Hit> assns(trackHandle, e,
                                        trackHandle.provenance()->moduleLabel());

      for (size_t i = 0; i < trackHandle->size(); i++) {
        const recob::Track &track = trackHandle->at(i);
        track_buf.from_track(track, track_key, producer);

        std::vector<art::Ptr<recob::Hit>> matched_hits;
        if (assns.isValid()) {
          matched_hits = assns.at(i);
          for (art::Ptr<recob::Hit> const &hit_ptr : matched_hits) {
            hit_to_trackID[hit_ptr] = track_key;
          }
        }
        track_buf.sum_hit_charges(matched_hits);
        track_tree->Fill();

        track_key++;
      }
    }
  }

  if (dump_recohits) {
    std::vector<art::Handle<std::vector<recob::Hit>>> hitHandles =
        e.getMany<std::vector<recob::Hit>>();

    for (auto const &hitHandle : hitHandles) {
      std::string producer = hitHandle.provenance()->inputTag().encode();

      for (size_t i = 0; i < hitHandle->size(); i++) {
        art::Ptr<recob::Hit> hit_ptr(hitHandle, i);
        auto it = hit_to_trackID.find(hit_ptr);
        int track_key = it != hit_to_trackID.end() ? it->second : -1;
        ChannelInfo chinfo = get_channel_info_for_channel(hit_ptr->Channel(),
                                                           hit_ptr->WireID());
        hit_buf.from_hit(*hit_ptr, track_key, producer, chinfo);
        hit_tree->Fill();
      }
    }
  }
}

duneana::ChannelInfo duneana::SolarAnaTree::get_channel_info_for_channel(
    raw::ChannelID_t channel, const geo::WireID &wireid) {
  readout::ROPID rop = fWireReadout->ChannelToROP(channel);
  ChannelInfo result;
  result.rop_id = rop.ROP;
  result.tpcset_id = rop.asTPCsetID().TPCset;
  result.view = fWireReadout->View(rop);

  // Hits already know which wire they were reconstructed on (GausHitFinder
  // et al. set this from the same channel via ChannelToWire()[0]), so reuse
  // it here instead of re-deriving it from the channel.
  if (fWireReadout->HasWire(wireid)) {
    geo::WireGeo const &wire = fWireReadout->Wire(wireid);
    result.wire_x = wire.GetCenter().X();
    result.wire_y = wire.GetCenter().Y();
    result.wire_z = wire.GetCenter().Z();
  } else {
    result.wire_x = INVALID;
    result.wire_y = INVALID;
    result.wire_z = INVALID;
  }
  return result;
}

DEFINE_ART_MODULE(duneana::SolarAnaTree)
