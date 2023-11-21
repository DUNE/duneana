////////////////////////////////////////////////////////////////////////
// Class:       EMPi0Energy
// Module Type: analyzer
// File:        EMPi0Energy_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), August 2015
//
// Analyser module to produce useful information for characterising
// pi0 showers.
////////////////////////////////////////////////////////////////////////

// Framework includes:
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDAnalyzer.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/Filters/ChannelFilter.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "lardataobj/Simulation/sim.h"

// ROOT & C++ includes
#include <string>
#include <vector>
#include <map>
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TVector3.h"

const int kMaxClusters = 100;

namespace emshower {
  class EMPi0Energy;
}

class emshower::EMPi0Energy : public art::EDAnalyzer {
public:

  EMPi0Energy(fhicl::ParameterSet const& pset);
  void analyze(art::Event const& evt);
  void reconfigure(fhicl::ParameterSet const& p);
  void reset();
  double ConvertChargeToEnergy(double charge, int plane);
  double FindDepositedEnergy(int trackID);
  int FindTrackID(detinfo::DetectorClocksData const& clockData,
                  art::Ptr<recob::Hit> const& hit);
  int FindTrueTrack(detinfo::DetectorClocksData const& clockData,
                    std::vector<art::Ptr<recob::Hit> > const& clusterHits);
  double FindVertexDetectorDistance(const simb::MCParticle* particle);

private:

  // Variables for the tree
  TTree* fTree;
  double trueEnergyPi0;
  double trueEnergyHighEPhoton;
  double trueEnergyLowEPhoton;
  int    trueTrackIDPi0;
  int    trueTrackIDHighEPhoton;
  int    trueTrackIDLowEPhoton;
  double depositHighEPhoton;
  double depositLowEPhoton;
  double vertexDetectorDistHighEPhoton;
  double vertexDetectorDistLowEPhoton;
  int    nclusters;
  int    cluster_plane      [kMaxClusters];
  int    cluster_size       [kMaxClusters];
  int    cluster_truetrackid[kMaxClusters];
  double cluster_charge     [kMaxClusters];
  double cluster_energy     [kMaxClusters];

  std::string fHitsModuleLabel;
  std::string fClusterModuleLabel;
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  art::ServiceHandle<geo::Geometry> geom;

  // For converting charge to energy
  double Uintercept, Ugradient;
  double Vintercept, Vgradient;
  double Zintercept, Zgradient;

};

emshower::EMPi0Energy::EMPi0Energy(fhicl::ParameterSet const& pset) : EDAnalyzer(pset) {
  this->reconfigure(pset);
  fTree = tfs->make<TTree>("EMPi0Energy","EMPi0Energy");
  fTree->Branch("TrueEnergyPi0",                &trueEnergyPi0);
  fTree->Branch("TrueEnergyHighEPhoton",        &trueEnergyHighEPhoton);
  fTree->Branch("TrueEnergyLowEPhoton",         &trueEnergyLowEPhoton);
  fTree->Branch("TrueTrackIDPi0",               &trueTrackIDPi0);
  fTree->Branch("TrueTrackIDHighEPhoton",       &trueTrackIDHighEPhoton);
  fTree->Branch("TrueTrackIDLowEPhoton",        &trueTrackIDLowEPhoton);
  fTree->Branch("DepositHighEPhoton",           &depositHighEPhoton);
  fTree->Branch("DepositLowEPhoton",            &depositLowEPhoton);
  fTree->Branch("VertexDetectorDistHighEPhoton",&vertexDetectorDistHighEPhoton);
  fTree->Branch("VertexDetectorDistLowEPhoton", &vertexDetectorDistLowEPhoton);
  fTree->Branch("NClusters",                    &nclusters);
  fTree->Branch("Cluster_Plane",                cluster_plane,      "cluster_plane[NClusters]/I");
  fTree->Branch("Cluster_Size",                 cluster_size,       "cluster_size[NClusters]/I");
  fTree->Branch("Cluster_TrueTrackID",          cluster_truetrackid,"cluster_truetrackid[NClusters]/I");
  fTree->Branch("Cluster_Charge",               cluster_charge,     "cluster_charge[NClusters]/D");
  fTree->Branch("Cluster_Energy",               cluster_energy,     "cluster_energy[NClusters]/D");

  Uintercept = -1519.33; Ugradient = 148867;
  Vintercept = -1234.91; Vgradient = 149458;
  Zintercept = -1089.73; Zgradient = 145372;

}

void emshower::EMPi0Energy::reconfigure(fhicl::ParameterSet const& pset) {
  fHitsModuleLabel = pset.get<std::string>("HitsModuleLabel");
  fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
}

void emshower::EMPi0Energy::analyze(art::Event const& evt) {

  /// Analyse function to save information for calibrating shower energies
  /// This is written assuming single pi0 per event

  this->reset();

  // Get the hits out of the event record
  std::vector<art::Ptr<recob::Hit> > hits;
  auto hitHandle = evt.getHandle<std::vector<recob::Hit> >(fHitsModuleLabel);
  if (hitHandle)
    art::fill_ptr_vector(hits, hitHandle);

  // Get the clusters out of the event record
  std::vector<art::Ptr<recob::Cluster> > clusters;
  auto clusterHandle = evt.getHandle<std::vector<recob::Cluster> >(fClusterModuleLabel);
  if (clusterHandle)
    art::fill_ptr_vector(clusters, clusterHandle);

  art::FindManyP<recob::Hit> fmh(clusterHandle, evt, fClusterModuleLabel);

  // Look at the clusters
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  for (unsigned int clus = 0; clus < clusters.size(); ++clus) {

    art::Ptr<recob::Cluster> cluster = clusters.at(clus);
    std::vector<art::Ptr<recob::Hit> > hits = fmh.at(clus);

    // Find the charge deposited by hits in this cluster in this plane
    double charge = 0;
    for (std::vector<art::Ptr<recob::Hit> >::iterator hit = hits.begin(); hit != hits.end(); ++hit)
      charge += ((*hit)->Integral() * TMath::Exp((500 * (*hit)->PeakTime())/3e6));

    cluster_plane      [clus] = cluster->Plane().Plane;
    cluster_size       [clus] = hits.size();
    cluster_truetrackid[clus] = this->FindTrueTrack(clockData, hits);
    cluster_charge     [clus] = charge;
    cluster_energy     [clus] = this->ConvertChargeToEnergy(charge, cluster->Plane().Plane);

  }

  // Event level information

  nclusters = clusters.size();

  // Get the pi0 and the decay photons
  const sim::ParticleList& trueParticles = particleinventory->ParticleList();
  const simb::MCParticle* truePi0 = trueParticles.Primary(0);
  if (truePi0->NumberDaughters() != 2) return;
  const simb::MCParticle* truePhoton1 = particleinventory->TrackIdToParticle_P(truePi0->Daughter(0));
  const simb::MCParticle* truePhoton2 = particleinventory->TrackIdToParticle_P(truePi0->Daughter(1));

  // Make sure photon 1 energy > photon 2 energy
  if (truePhoton1->Momentum().E() < truePhoton2->Momentum().E()) {
    const simb::MCParticle* tmp = truePhoton2;
    truePhoton2 = truePhoton1;
    truePhoton1 = tmp;
  }

  trueEnergyPi0  = truePi0->Momentum().E();
  trueTrackIDPi0 = truePi0->TrackId();

  trueEnergyHighEPhoton = truePhoton1->Momentum().E();
  trueEnergyLowEPhoton = truePhoton2->Momentum().E();
  trueTrackIDHighEPhoton = truePhoton1->TrackId();
  trueTrackIDLowEPhoton = truePhoton2->TrackId();

  // Find the energy deposited on each plane in the TPC
  depositHighEPhoton = this->FindDepositedEnergy(trueTrackIDHighEPhoton)/3;
  depositLowEPhoton = this->FindDepositedEnergy(trueTrackIDLowEPhoton)/3;

  // Find the distance between the particle vertex and the edge of the detector
  vertexDetectorDistHighEPhoton = this->FindVertexDetectorDistance(truePhoton1);
  vertexDetectorDistLowEPhoton = this->FindVertexDetectorDistance(truePhoton2);

  fTree->Fill();

  return;

}

double emshower::EMPi0Energy::ConvertChargeToEnergy(double charge, int plane) {

  /// Converts charge to energy

  double energy = 0;

  switch (plane) {
  case 0:
    energy = (double)(charge - Uintercept)/(double)Ugradient;
    break;
  case 1:
    energy = (double)(charge - Vintercept)/(double)Vgradient;
    break;
  case 2:
    energy = (double)(charge - Zintercept)/(double)Zgradient;
    break;
  }

  return energy;

}

double emshower::EMPi0Energy::FindDepositedEnergy(int trackID) {

      std::vector<sim::IDE> ides;
    for( auto ide_P : backtracker->TrackIdToSimIDEs_Ps(trackID) ){
      ides.push_back(*ide_P);
    }

  double deposit = 0;

  for (std::vector<sim::IDE>::iterator ideIt = ides.begin(); ideIt != ides.end(); ++ideIt)
    deposit += ideIt->energy;

  // Put energies in GeV units
  deposit /= 1000;

  return deposit;

}

int emshower::EMPi0Energy::FindTrackID(detinfo::DetectorClocksData const& clockData,
                                       art::Ptr<recob::Hit> const& hit) {

  /// Find the true track ID this hit is associated with

  double particleEnergy = 0;
  int likelyTrackID = 0;
  std::vector<sim::TrackIDE> trackIDs = backtracker->HitToTrackIDEs(clockData, hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = TMath::Abs(trackIDs.at(idIt).trackID);
    }
  }
  return likelyTrackID;
}

int emshower::EMPi0Energy::FindTrueTrack(detinfo::DetectorClocksData const& clockData,
                                         std::vector<art::Ptr<recob::Hit> > const& clusterHits) {

  /// Find the true track which is most associated to this cluster

  std::map<int,double> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator clusHitIt = clusterHits.begin(); clusHitIt != clusterHits.end(); ++clusHitIt) {
    art::Ptr<recob::Hit> hit = *clusHitIt;
    int trackID = FindTrackID(clockData, hit);
    trackMap[trackID] += hit->Integral();
  }
  //return std::max_element(trackMap.begin(), trackMap.end(), [](const std::pair<int,double>& p1, const std::pair<int,double>& p2) {return p1.second < p2.second;} )->first;
  double highestCharge = 0;
  int clusterTrack = 0;
  for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt)
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      clusterTrack  = trackIt->first;
    }
  return clusterTrack;
}

double emshower::EMPi0Energy::FindVertexDetectorDistance(const simb::MCParticle* particle) {

  /// Finds the rough distance from the end of a particle track to the edge of the detector, along the direction of the particle

  TVector3 end = TVector3(particle->EndX(),particle->EndY(),particle->EndZ());
  TVector3 direction = TVector3(particle->Px(),particle->Py(),particle->Pz()).Unit();

  int distanceStep = 1, steps = 0;
  TVector3 pos;
  bool inTPC = true;

  while (inTPC) {
    pos = end + ( (steps*distanceStep) * direction );
    if (!geom->FindTPCAtPosition(geo::vect::toPoint(pos)).isValid)
      inTPC = false;
    ++steps;
  }

  double distance = (end - pos).Mag();

  return distance;

}

void emshower::EMPi0Energy::reset() {
  trueEnergyPi0 = 0;
  trueEnergyHighEPhoton = 0;
  trueEnergyLowEPhoton = 0;
  trueTrackIDPi0 = 0;
  trueTrackIDHighEPhoton = 0;
  trueTrackIDLowEPhoton = 0;
  depositHighEPhoton = 0;
  depositLowEPhoton = 0;
  vertexDetectorDistHighEPhoton = 0;
  vertexDetectorDistLowEPhoton = 0;
  nclusters = 0;
  for (int cluster = 0; cluster < kMaxClusters; ++cluster) {
    cluster_plane      [cluster] = 0;
    cluster_size       [cluster] = 0;
    cluster_truetrackid[cluster] = 0;
    cluster_charge     [cluster] = 0;
    cluster_energy     [cluster] = 0;
  }
}

DEFINE_ART_MODULE(emshower::EMPi0Energy)
