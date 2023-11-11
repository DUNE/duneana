//////////////////////////////////////////////////////////////////////////
// Class:       ShowerAnalysis
// Module type: analyser
// File:        ShowerAnalysis_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), April 2016
//
// Analyser module to evaluate the shower reconstruction performance.
// Produces histograms and information within a tree for further
// analysis.
//////////////////////////////////////////////////////////////////////////

// framework includes
#include "art/Framework/Core/EDAnalyzer.h"
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
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT
#include "TTree.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

// c++
#include <string>

namespace showerAna {
  class ShowerAnalysis;
  class ShowerParticle;
}

class showerAna::ShowerParticle {
public:

  ShowerParticle();                               // default constructor
  ShowerParticle(int id);                         // standard constructor
  ShowerParticle(const ShowerParticle& particle); // copy constructor
  ~ShowerParticle();

  // Setters
  void SetEnergy(double energy);
  void SetDepositedEnergy(std::map<int,double> depositedEnergy);
  void SetDirection(TVector3 direction);
  void SetStart(TVector3 start);
  void SetEnd(TVector3 end);
  void SetPDG(int pdg);

  void AddAssociatedHit(const art::Ptr<recob::Hit>& hit);
  void AddAssociatedCluster(const art::Ptr<recob::Cluster>& cluster, const std::vector<art::Ptr<recob::Hit> >& hits, const std::vector<art::Ptr<recob::Hit> >& trueHits);
  void AddAssociatedShower(const art::Ptr<recob::Shower>& shower, const std::vector<art::Ptr<recob::Hit> >& hits, const std::vector<art::Ptr<recob::Hit> >& trueHits);

  // Getters
  int ID() const { return fID; }
  int PDG() const { return fPDG; }

  TVector3 Start() const { return fStart; }
  TVector3 End() const { return fEnd; }
  TVector3 Direction() const { return fDirection; }
  double Energy() const { return fEnergy; }
  double DepositedEnergy() const { art::Ptr<recob::Shower> s = fShowers.at(fLargestShower); return DepositedEnergy(s->best_plane()); }
  double DepositedEnergy(int plane) const { return fDepositedEnergy.at(plane); }

  TVector3 ShowerStart() const { return ShowerStart(fLargestShower); }
  TVector3 ShowerStart(int shower) const { return fShowers.at(shower)->ShowerStart(); }
  TVector3 ShowerDirection() const { return ShowerDirection(fLargestShower); }
  TVector3 ShowerDirection(int shower) const { return fShowers.at(shower)->Direction(); }
  double ShowerEnergy() const { return ShowerEnergy(fLargestShower); }
  double ShowerEnergy(int shower) const { art::Ptr<recob::Shower> s = fShowers.at(shower); return s->Energy().at(s->best_plane()); }
  double ShowerdEdx() const { return ShowerdEdx(fLargestShower); }
  double ShowerdEdx(int shower) const { art::Ptr<recob::Shower> s = fShowers.at(shower); return s->dEdx().at(s->best_plane()); }

  int NumHits() const { return fHits.size(); }
  int NumClusters() const { return fClusters.size(); }
  int NumShowers() const { return fShowers.size(); }

  int LargestCluster() const { return fLargestCluster; }
  bool LargestCluster(int cluster) const { return fLargestCluster == cluster; }
  int LargestShower() const { return fLargestShower; }
  bool LargestShower(int shower) const { return fLargestShower == shower; }

  double ClusterCompleteness(int cluster) const { return fClusterTrueHits.at(cluster).size() / (double)fHits.size(); }
  double ClusterPurity(int cluster) const { return fClusterTrueHits.at(cluster).size() / (double)fClusterHits.at(cluster).size(); }
  double ShowerCompleteness(int shower) const { return fShowerTrueHits.at(shower).size() / (double)fHits.size(); }
  double ShowerPurity(int shower) const { return fShowerTrueHits.at(shower).size() / (double)fShowerHits.at(shower).size(); }

private:

  int fID;
  int fPDG;

  TVector3 fStart, fEnd, fDirection;
  double fEnergy;
  std::map<int,double> fDepositedEnergy;

  std::vector<art::Ptr<recob::Hit> > fHits;
  std::vector<art::Ptr<recob::Cluster> > fClusters;
  std::vector<std::vector<art::Ptr<recob::Hit> > > fClusterHits, fClusterTrueHits;
  std::vector<art::Ptr<recob::Shower> > fShowers;
  std::vector<std::vector<art::Ptr<recob::Hit> > > fShowerHits, fShowerTrueHits;

  int fLargestCluster, fLargestShower;

};

showerAna::ShowerParticle::ShowerParticle(int id) {
  fID = id;
  fLargestCluster = 0;
  fLargestShower = 0;
}

showerAna::ShowerParticle::ShowerParticle(const ShowerParticle& particle) {
  fID = particle.ID();
  fLargestCluster = 0;
  fLargestShower = 0;
}

showerAna::ShowerParticle::~ShowerParticle() {
}

void showerAna::ShowerParticle::SetEnergy(double energy) {
  fEnergy = energy;
}

void showerAna::ShowerParticle::SetDepositedEnergy(std::map<int,double> depositedEnergy) {
  fDepositedEnergy = depositedEnergy;
}

void showerAna::ShowerParticle::SetDirection(TVector3 direction) {
  fDirection = direction;
}

void showerAna::ShowerParticle::SetStart(TVector3 start) {
  fStart = start;
}

void showerAna::ShowerParticle::SetEnd(TVector3 end) {
  fEnd = end;
}

void showerAna::ShowerParticle::SetPDG(int pdg) {
  fPDG = pdg;
}

void showerAna::ShowerParticle::AddAssociatedHit(const art::Ptr<recob::Hit>& hit) {
  fHits.push_back(hit);
}

void showerAna::ShowerParticle::AddAssociatedCluster(const art::Ptr<recob::Cluster>& cluster,
						     const std::vector<art::Ptr<recob::Hit> >& hits,
						     const std::vector<art::Ptr<recob::Hit> >& trueHits) {
  fClusters.push_back(cluster);
  fClusterHits.push_back(hits);
  fClusterTrueHits.push_back(trueHits);
  if (hits.size() > fClusterHits[fLargestCluster].size())
    fLargestCluster = fClusters.size() - 1;
}

void showerAna::ShowerParticle::AddAssociatedShower(const art::Ptr<recob::Shower>& shower,
						    const std::vector<art::Ptr<recob::Hit> >& hits,
						    const std::vector<art::Ptr<recob::Hit> >& trueHits) {
  fShowers.push_back(shower);
  fShowerHits.push_back(hits);
  fShowerTrueHits.push_back(trueHits);
  if (hits.size() > fShowerHits[fLargestShower].size())
    fLargestShower = fShowers.size() - 1;
}

class showerAna::ShowerAnalysis : public art::EDAnalyzer {
 public:

  ShowerAnalysis(const fhicl::ParameterSet& pset);
  ~ShowerAnalysis();

  void analyze(const art::Event& evt);

  void MakeDataProducts();
  void FillData(const std::map<int,std::shared_ptr<ShowerParticle> >& particles);
  void FillPi0Data(const std::map<int,std::shared_ptr<ShowerParticle> >& particles, const std::vector<int>& pi0Decays);
  int FindTrueParticle(detinfo::DetectorClocksData const& clockData,
                       const std::vector<art::Ptr<recob::Hit> >& hits);
  int FindParticleID(detinfo::DetectorClocksData const& clockData,
                     const art::Ptr<recob::Hit>& hit);
  std::vector<art::Ptr<recob::Hit> > FindTrueHits(detinfo::DetectorClocksData const& clockData,
                                                  const std::vector<art::Ptr<recob::Hit> >& hits, int trueParticle);

 private:

  std::string fShowerModuleLabel, fClusterModuleLabel, fHitsModuleLabel;

  double fElectrondEdx, fPhotondEdx, fElectrondEdxWidth, fPhotondEdxWidth;

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<art::TFileService> tfs;
  geo::WireReadoutGeom const* wireReadout = &art::ServiceHandle<geo::WireReadout const>{}->Get();

  TTree* fTree;

  // Showers
  TH1D *hClusterCompleteness, *hLargestClusterCompleteness, *hClusterPurity, *hLargestClusterPurity;
  TH2D *hClusterCompletenessEnergy, *hLargestClusterCompletenessEnergy, *hClusterCompletenessDirection, *hLargestClusterCompletenessDirection;
  TH1D *hShowerCompleteness, *hLargestShowerCompleteness, *hShowerPurity, *hLargestShowerPurity;
  TH2D *hShowerCompletenessEnergy, *hLargestShowerCompletenessEnergy, *hShowerCompletenessDirection, *hLargestShowerCompletenessDirection;
  TH1D *hShowerEnergy, *hShowerDepositedEnergy, *hShowerDirection, *hShowerdEdx, *hShowerStart, *hShowerReconstructed, *hNumShowersReconstructed, *hNumShowersReconstructedNonZero;
  TH2D *hShowerdEdxEnergy, *hNumShowersReconstructedEnergy;
  TProfile *hShowerReconstructedEnergy, *hShowerdEdxEnergyProfile, *hNumShowersReconstructedEnergyProfile, *hShowerEnergyCompleteness, *hShowerStartEnergy;
  TH1D *hElectronPull, *hPhotonPull;

  // Pi0
  TH1D* hPi0MassPeakReconEnergyReconAngle;
  TH1D* hPi0MassPeakReconEnergyTrueAngle;
  TH1D* hPi0MassPeakTrueEnergyReconAngle;
  TH1D* hPi0MassPeakTrueEnergyTrueAngle;

};

showerAna::ShowerAnalysis::ShowerAnalysis(const fhicl::ParameterSet& pset) : EDAnalyzer(pset) {
  fShowerModuleLabel  = pset.get<std::string>("ShowerModuleLabel");
  fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
  fHitsModuleLabel    = pset.get<std::string>("HitsModuleLabel");
  fElectrondEdx       = pset.get<double>("ElectrondEdx",2.1);
  fPhotondEdx         = pset.get<double>("PhotondEdx",4.2);
  fElectrondEdxWidth  = pset.get<double>("ElectrondEdxWidth",0.410);
  fPhotondEdxWidth    = pset.get<double>("PhotondEdxWidth",1.217);
  this->MakeDataProducts();
}

showerAna::ShowerAnalysis::~ShowerAnalysis() {
}

void showerAna::ShowerAnalysis::analyze(const art::Event& evt) {

  // Get showers out of event
  std::vector<art::Ptr<recob::Shower> > showers;
  auto showerHandle = evt.getHandle<std::vector<recob::Shower> >(fShowerModuleLabel);
  if (showerHandle)
    art::fill_ptr_vector(showers, showerHandle);

  // Get clusters out of event
  std::vector<art::Ptr<recob::Cluster> > clusters;
  auto clusterHandle = evt.getHandle<std::vector<recob::Cluster> >(fClusterModuleLabel);
  if (clusterHandle)
    art::fill_ptr_vector(clusters, clusterHandle);

  // Get hits out of event
  std::vector<art::Ptr<recob::Hit> > hits;
  auto hitHandle = evt.getHandle<std::vector<recob::Hit> >(fHitsModuleLabel);
  if (hitHandle)
    art::fill_ptr_vector(hits, hitHandle);

  // Get space points out of event
  std::vector<art::Ptr<recob::SpacePoint> > spacePoints;
  auto spacePointHandle = evt.getHandle<std::vector<recob::SpacePoint> >(fShowerModuleLabel);
  if (spacePointHandle)
    art::fill_ptr_vector(spacePoints, spacePointHandle);

  // Get associations out of event
  art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::Hit> fmhs(showerHandle, evt, fShowerModuleLabel);
  art::FindManyP<recob::Hit> fmhsp(spacePointHandle, evt, fShowerModuleLabel);

  // Map all the true and reconstructed information for each particle
  std::map<int,std::shared_ptr<ShowerParticle> > particles;

  // Keep an eye out for pi0s!
  bool isPi0 = false;
  std::vector<int> pi0Decays;

  // Fill true properties
  const sim::ParticleList& trueParticles = pi_serv->ParticleList();
  for (sim::ParticleList::const_iterator particleIt = trueParticles.begin(); particleIt != trueParticles.end(); ++particleIt) {

    const simb::MCParticle* trueParticle = particleIt->second;
    int mother = trueParticle->Mother();
    if (mother != 0 and trueParticles.at(mother)->PdgCode() == 111) {
      isPi0 = true;
      pi0Decays.push_back(particleIt->first);
    }

    std::map<int,double> depositedEnergy;
    const std::vector<art::Ptr< sim::SimChannel >>& simChannels = bt_serv->SimChannels();
    for (std::vector<art::Ptr< sim::SimChannel >>::const_iterator channelIt = simChannels.begin(); channelIt != simChannels.end(); ++channelIt) {
      int plane = wireReadout->View((*channelIt)->Channel());
      auto const & tdcidemap = (*channelIt)->TDCIDEMap();
      for (auto const& tdcIt : tdcidemap) {
        auto const& idevec = tdcIt.second;
        for (auto const& ideIt : idevec) {
          if (TMath::Abs(ideIt.trackID) != trueParticle->TrackId())
            continue;
          depositedEnergy[plane] += ideIt.energy / 1000;
        }
      }
    }

    std::shared_ptr<ShowerParticle> particle = std::make_shared<ShowerParticle>(trueParticle->TrackId());
    particle->SetEnergy(trueParticle->E());
    particle->SetDepositedEnergy(depositedEnergy);
    particle->SetDirection(trueParticle->Momentum().Vect().Unit());
    particle->SetStart(trueParticle->Position().Vect());
    particle->SetEnd(trueParticle->EndPosition().Vect());
    particle->SetPDG(trueParticle->PdgCode());

    particles[particleIt->first] = std::move(particle);
  }

  // Fill recon properties
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    int trueParticle = FindParticleID(clockData, *hitIt);
    if (particles.count(trueParticle))
      particles[trueParticle]->AddAssociatedHit(*hitIt);
  }

  for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = clusters.begin(); clusterIt != clusters.end(); ++clusterIt) {
    std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(clusterIt->key());
    int trueParticle = FindTrueParticle(clockData, hits);
    std::vector<art::Ptr<recob::Hit> > trueHits = FindTrueHits(clockData, hits, trueParticle);
    if (particles.count(trueParticle))
      particles[trueParticle]->AddAssociatedCluster(*clusterIt, hits, trueHits);
  }

  for (std::vector<art::Ptr<recob::Shower> >::iterator showerIt = showers.begin(); showerIt != showers.end(); ++showerIt) {
    std::vector<art::Ptr<recob::Hit> > hits = fmhs.at(showerIt->key());
    int trueParticle = FindTrueParticle(clockData, hits);
    std::vector<art::Ptr<recob::Hit> > trueHits = FindTrueHits(clockData, hits, trueParticle);
    if (particles.count(trueParticle))
      particles[trueParticle]->AddAssociatedShower(*showerIt, hits, trueHits);
  }

  // Fill output data products
  this->FillData(particles);
  if (isPi0 and pi0Decays.size() == 2)
    this->FillPi0Data(particles, pi0Decays);

  return;

}

void showerAna::ShowerAnalysis::FillData(const std::map<int,std::shared_ptr<ShowerParticle> >& particles) {

  // Look at each particle
  for (std::map<int,std::shared_ptr<ShowerParticle> >::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt) {

    std::shared_ptr<ShowerParticle> particle = particleIt->second;

    // No point making these plots for particles which don't shower!
    if (TMath::Abs(particle->PDG()) != 11 and TMath::Abs(particle->PDG()) != 22)
      continue;

    // Cluster plots
    for (int cluster = 0; cluster < particle->NumClusters(); ++cluster) {
      hClusterCompleteness		->Fill(particle->ClusterCompleteness(cluster));
      hClusterPurity			->Fill(particle->ClusterPurity(cluster));
      hClusterCompletenessEnergy	->Fill(particle->Energy(), particle->ClusterCompleteness(cluster));
      hClusterCompletenessDirection	->Fill(particle->Direction().Angle(TVector3(0,1,0)), particle->ClusterCompleteness(cluster));
      if (particle->LargestCluster(cluster)) {
	hLargestClusterCompleteness		->Fill(particle->ClusterCompleteness(cluster));
	hLargestClusterPurity			->Fill(particle->ClusterPurity(cluster));
	hLargestClusterCompletenessEnergy	->Fill(particle->Energy(), particle->ClusterCompleteness(cluster));
	hLargestClusterCompletenessDirection	->Fill(particle->Direction().Angle(TVector3(0,1,0)), particle->ClusterCompleteness(cluster));
      }
    }

    // Shower plots
    for (int shower = 0; shower < particle->NumShowers(); ++shower) {
      hShowerCompleteness		->Fill(particle->ShowerCompleteness(shower));
      hShowerPurity			->Fill(particle->ShowerPurity(shower));
      hShowerCompletenessEnergy		->Fill(particle->Energy(), particle->ShowerCompleteness(shower));
      hShowerCompletenessDirection	->Fill(particle->Direction().Angle(TVector3(0,1,0)), particle->ShowerCompleteness(shower));
      if (particle->LargestShower(shower)) {
	hLargestShowerCompleteness		->Fill(particle->ShowerCompleteness(shower));
	hLargestShowerPurity			->Fill(particle->ShowerPurity(shower));
	hLargestShowerCompletenessEnergy	->Fill(particle->Energy(), particle->ShowerCompleteness(shower));
	hLargestShowerCompletenessDirection	->Fill(particle->Direction().Angle(TVector3(0,1,0)), particle->ShowerCompleteness(shower));
	hShowerEnergy				->Fill(particle->ShowerEnergy() / particle->Energy());
	hShowerDepositedEnergy                  ->Fill(particle->ShowerEnergy() / particle->DepositedEnergy());
	hShowerEnergyCompleteness               ->Fill(particle->Energy(), particle->ShowerEnergy()/particle->DepositedEnergy());
	hShowerDirection                        ->Fill(particle->Direction().Dot(particle->ShowerDirection()));
	hShowerdEdx                             ->Fill(particle->ShowerdEdx());
	if (particle->PDG() == 22) {
	  hShowerStart                          ->Fill((particle->End() - particle->ShowerStart()).Mag());
	  hShowerStartEnergy                    ->Fill(particle->Energy(), (particle->End() - particle->ShowerStart()).Mag());
	}
	else {
	  hShowerStart                          ->Fill((particle->Start() - particle->ShowerStart()).Mag());
	  hShowerStartEnergy                    ->Fill(particle->Energy(), (particle->Start() - particle->ShowerStart()).Mag());
	}
      }
    }

    hShowerReconstructed->Fill(particle->NumShowers() != 0);
    hNumShowersReconstructed->Fill(particle->NumShowers());
    hShowerReconstructedEnergy->Fill(particle->Energy(), particle->NumShowers() != 0);
    hNumShowersReconstructedEnergy->Fill(particle->Energy(), particle->NumShowers());
    hNumShowersReconstructedEnergyProfile->Fill(particle->Energy(), particle->NumShowers());
    if (particle->NumShowers()) {
      hNumShowersReconstructedNonZero->Fill(particle->NumShowers());
      hShowerdEdxEnergy->Fill(particle->Energy(), particle->ShowerdEdx());
      hShowerdEdxEnergyProfile->Fill(particle->Energy(), particle->ShowerdEdx());
      hElectronPull->Fill((particle->ShowerdEdx() - fElectrondEdx)/fElectrondEdxWidth);
      hPhotonPull->Fill((particle->ShowerdEdx() - fPhotondEdx)/fPhotondEdxWidth);
    }

  }

  return;

}

void showerAna::ShowerAnalysis::FillPi0Data(const std::map<int,std::shared_ptr<ShowerParticle> >& particles, const std::vector<int>& pi0Decays) {

  /// Fill data products related to pi0s
  /// We only call this function if there are exactly two decay products, so we can assume the decay was pi0 -> gamma gamma

  std::shared_ptr<ShowerParticle> photon1 = particles.at(pi0Decays.at(0));
  std::shared_ptr<ShowerParticle> photon2 = particles.at(pi0Decays.at(1));

  // Make sure there is at least one shower and it was fully reconstructed
  if (photon1->NumShowers() == 0 or photon2->NumShowers() == 0 or
      photon1->ShowerStart() == TVector3(0,0,0) or photon2->ShowerStart() == TVector3(0,0,0)) {
    std::cout << "Event doesn't have enough info for pi0 mass determination" << std::endl;
    return;
  }

  double reconOpeningAngle = TMath::ASin(TMath::Sin(photon1->ShowerDirection().Angle(photon2->ShowerDirection())));
  double trueOpeningAngle = photon1->Direction().Angle(photon2->Direction());

  // Reconstruct the pi0 mass
  double massReconEnergyReconAngle = TMath::Sqrt(4 * photon1->ShowerEnergy() * photon2->ShowerEnergy() * TMath::Power(TMath::Sin(reconOpeningAngle/2),2));
  double massReconEnergyTrueAngle = TMath::Sqrt(4 * photon1->ShowerEnergy() * photon2->ShowerEnergy() * TMath::Power(TMath::Sin(trueOpeningAngle/2),2));
  double massTrueEnergyReconAngle = TMath::Sqrt(4 * photon1->Energy() * photon2->Energy() * TMath::Power(TMath::Sin(reconOpeningAngle/2),2));
  double massTrueEnergyTrueAngle = TMath::Sqrt(4 * photon1->Energy() * photon2->Energy() * TMath::Power(TMath::Sin(trueOpeningAngle/2),2));

  hPi0MassPeakReconEnergyReconAngle->Fill(massReconEnergyReconAngle);
  hPi0MassPeakReconEnergyTrueAngle->Fill(massReconEnergyTrueAngle);
  hPi0MassPeakTrueEnergyReconAngle->Fill(massTrueEnergyReconAngle);
  hPi0MassPeakTrueEnergyTrueAngle->Fill(massTrueEnergyTrueAngle);

  return;
  
}

int showerAna::ShowerAnalysis::FindTrueParticle(detinfo::DetectorClocksData const& clockData,
                                                const std::vector<art::Ptr<recob::Hit> >& showerHits) {

  /// Returns the true particle most likely associated with this shower

  // Make a map of the tracks which are associated with this shower and the charge each contributes
  std::map<int,double> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitIt = showerHits.begin(); showerHitIt != showerHits.end(); ++showerHitIt) {
    art::Ptr<recob::Hit> hit = *showerHitIt;
    int trackID = FindParticleID(clockData, hit);
    trackMap[trackID] += hit->Integral();
  }

  // Pick the track with the highest charge as the 'true track'
  double highestCharge = 0;
  int showerTrack = 0;
  for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt) {
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      showerTrack  = trackIt->first;
    }
  }

  return showerTrack;

}

int showerAna::ShowerAnalysis::FindParticleID(detinfo::DetectorClocksData const& clockData,
                                              const art::Ptr<recob::Hit>& hit) {

  /// Returns the true track ID associated with this hit (if more than one, returns the one with highest energy)

  double particleEnergy = 0;
  int likelyTrackID = 0;
  std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = TMath::Abs(trackIDs.at(idIt).trackID);
    }
  }

  return likelyTrackID;

}

std::vector<art::Ptr<recob::Hit> > showerAna::ShowerAnalysis::FindTrueHits(detinfo::DetectorClocksData const& clockData,
                                                                           const std::vector<art::Ptr<recob::Hit> >& hits, int trueParticle) {

  std::vector<art::Ptr<recob::Hit> > trueHits;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
    if (FindParticleID(clockData, *hitIt) == trueParticle)
      trueHits.push_back(*hitIt);

  return trueHits;

}

void showerAna::ShowerAnalysis::MakeDataProducts() {

  fTree = tfs->make<TTree>("ShowerAnalysis","ShowerAnalysis");

  // Cluster
  hClusterCompleteness = tfs->make<TH1D>("ClusterCompleteness","Completeness of all clusters;Completeness;",101,0,1.01);
  hLargestClusterCompleteness = tfs->make<TH1D>("LargestClusterCompleteness","Completeness of largest cluster;Completeness;",101,0,1.01);
  hClusterPurity = tfs->make<TH1D>("ClusterPurity","Purity of all clusters;Purity;",101,0,1.01);
  hLargestClusterPurity = tfs->make<TH1D>("LargestClusterPurity","Purity of largest cluster;Purity;",101,0,1.01);
  hClusterCompletenessEnergy = tfs->make<TH2D>("ClusterCompletenessEnergy","Completeness of all clusters vs Energy;Energy (GeV);Completeness;",100,0,10,101,0,1.01);
  hLargestClusterCompletenessEnergy = tfs->make<TH2D>("LargestClusterCompletenessEnergy","Completeness of largest cluster vs Energy;Energy (GeV);Completeness;",100,0,10,101,0,1.01);
  hClusterCompletenessDirection = tfs->make<TH2D>("ClusterCompletenessDirection","Completeness of all clusters vs Direction;Direction;Completeness;",100,0,5,101,0,1.01);
  hLargestClusterCompletenessDirection = tfs->make<TH2D>("LargestClusterCompletenessDirection","Completeness of largest cluster vs Direction;Direction;Completeness;",100,0,5,101,0,1.01);

  // Shower
  hShowerCompleteness = tfs->make<TH1D>("ShowerCompleteness","Completeness of all showers;Completeness;",101,0,1.01);
  hLargestShowerCompleteness = tfs->make<TH1D>("LargestShowerCompleteness","Completeness of largest shower;Completeness;",101,0,1.01);
  hShowerPurity = tfs->make<TH1D>("ShowerPurity","Purity of all showers;Purity;",101,0,1.01);
  hLargestShowerPurity = tfs->make<TH1D>("LargestShowerPurity","Purity of largest shower;Purity;",101,0,1.01);
  hShowerCompletenessEnergy = tfs->make<TH2D>("ShowerCompletenessEnergy","Completeness of all showers vs Energy;Energy (GeV);Completeness;",100,0,10,101,0,1.01);
  hLargestShowerCompletenessEnergy = tfs->make<TH2D>("LargestShowerCompletenessEnergy","Completeness of largest shower vs Energy;Energy (GeV);Completeness;",100,0,10,101,0,1.01);
  hShowerCompletenessDirection = tfs->make<TH2D>("ShowerCompletenessDirection","Completeness of all showers vs Direction;Direction;Completeness;",100,0,5,101,0,1.01);
  hLargestShowerCompletenessDirection = tfs->make<TH2D>("LargestShowerCompletenessDirection","Completeness of largest shower vs Direction;Direction;Completeness;",100,0,5,101,0,1.01);
  hShowerEnergy = tfs->make<TH1D>("ShowerEnergy","Shower energy;Recon Energy/True Energy;",120,0,1.2);
  hShowerDepositedEnergy = tfs->make<TH1D>("ShowerDepositedEnergy","Shower deposited energy;Recon Energy/True (Deposited) Energy;",120,0,1.2);
  hShowerDirection = tfs->make<TH1D>("ShowerDirection","Shower direction;True Direction.(Recon Direction);",202,-1.01,1.01);
  hShowerdEdx = tfs->make<TH1D>("ShowerdEdx","dEdx of Shower;dE/dx (MeV/cm)",50,0,10);
  hShowerStart = tfs->make<TH1D>("ShowerStart","Distance from reconstructed shower start to true shower start;True Start - Recon Start (cm);",100,0,20);
  hShowerReconstructed = tfs->make<TH1D>("ShowerReconstructed","% of showering particles with reconstructed shower",101,0,1.01);
  hShowerEnergyCompleteness = tfs->make<TProfile>("ShowerEnergyCompleteness","Shower energy completeness vs true deposited energy;True (Deposited) Energy (GeV);Energy Completeness;",100,0,5);
  hShowerStartEnergy = tfs->make<TProfile>("ShowerStartEnergy","Shower start distance vs true deposited energy;True (Deposited) Energy (GeV);True Start - Recon Start (cm);",100,0,5);
  hNumShowersReconstructed = tfs->make<TH1D>("NumShowersReconstructed","Number of showers reconstructed for each showering particle;Number of Showers",10,0,10);
  hNumShowersReconstructedNonZero = tfs->make<TH1D>("NumShowersReconstructedNonZero","Number of showers reconstructed for each showering particle;Number of Showers",9,1,10);
  for (int n_showers = 0; n_showers < 9; ++n_showers) {
    hNumShowersReconstructed->GetXaxis()->SetBinLabel(n_showers+1,std::to_string(n_showers).c_str());
    hNumShowersReconstructedNonZero->GetXaxis()->SetBinLabel(n_showers+1,std::to_string(n_showers+1).c_str());
  }
  hNumShowersReconstructed->GetXaxis()->SetBinLabel(10,"9 or more");
  hNumShowersReconstructedNonZero->GetXaxis()->SetBinLabel(9,"9 or more");
  hShowerReconstructedEnergy = tfs->make<TProfile>("ShowerReconstructedEnergyProfile","% of showering particles with reconstructed shower vs true energy;True Energy (GeV);Fraction of particles with reconstructed shower;",100,0,5);
  hNumShowersReconstructedEnergy = tfs->make<TH2D>("NumShowersReconstructedEnergy","Number of showers reconstructed per showering particle vs true energy;True Energy (GeV);Average Number of Showers;",100,0,5,10,0,10);
  hNumShowersReconstructedEnergyProfile = tfs->make<TProfile>("NumShowersReconstructedEnergyProfile",";True Energy (GeV);Number of Showers;",100,0,5);
  hShowerdEdxEnergy = tfs->make<TH2D>("ShowerdEdxEnergy","Shower dE/dx vs true energy;True Energy (GeV);dE/dx (MeV/cm);",100,0,5,50,0,10);
  hShowerdEdxEnergyProfile = tfs->make<TProfile>("ShowerdEdxEnergyProfile",";True Energy (GeV);dE/dx (MeV/cm);",100,0,5);
  hElectronPull = tfs->make<TH1D>("ElectronPull","Electron pull;Electron pull;",100,-10,10);
  hPhotonPull = tfs->make<TH1D>("PhotonPull","Photon pull;Photon pull;",100,-10,10);

  // pi0
  hPi0MassPeakReconEnergyReconAngle = tfs->make<TH1D>("Pi0MassPeakReconEnergyReconAngle",";Invariant Mass (GeV);",40,0,0.5);
  hPi0MassPeakTrueEnergyReconAngle = tfs->make<TH1D>("Pi0MassPeakTrueEnergyReconAngle",";Invariant Mass (GeV);",40,0,0.5);
  hPi0MassPeakReconEnergyTrueAngle = tfs->make<TH1D>("Pi0MassPeakReconEnergyTrueAngle",";Invariant Mass (GeV);",40,0,0.5);
  hPi0MassPeakTrueEnergyTrueAngle = tfs->make<TH1D>("Pi0MassPeakTrueEnergyTrueAngle",";Invariant Mass (GeV);",40,0,0.5);

}

DEFINE_ART_MODULE(showerAna::ShowerAnalysis)
