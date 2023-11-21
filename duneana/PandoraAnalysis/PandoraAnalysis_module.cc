////////////////////////////////////////////////////////////////////////
// Class:       pandoraAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        pandoraAnalysis_module.cc
//
// Generated at Fri Aug  7 15:01:35 2020 by Maria Brigida Brunetti using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <TTree.h>
#include <vector>
#include <string>
#include <TH1D.h>
#include <TLorentzVector.h>


namespace test {
  class pandoraAnalysis;
}


class test::pandoraAnalysis : public art::EDAnalyzer {
public:
  explicit pandoraAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  pandoraAnalysis(pandoraAnalysis const&) = delete;
  pandoraAnalysis(pandoraAnalysis&&) = delete;
  pandoraAnalysis& operator=(pandoraAnalysis const&) = delete;
  pandoraAnalysis& operator=(pandoraAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void reset(bool deepClean=false);

private:

  // Declare member data here.

  TTree *fTree;

  unsigned int fEventID;
  unsigned int fRunID;
  unsigned int fSubRunID;

  unsigned int fNMCParticles;
  unsigned int fNPFParticles;

/*  static const int kNMaxMCParticles = 200;
  static const int kNMaxPFParticles = 200;
  static const int kNMaxPFPClusters = 3;*/

  static const int kNMaxMCParticles = 20000;
  static const int kNMaxPFParticles = 2000;
  static const int kNMaxPFPClusters = 100;
  static const int kNViews = 3;


  bool fMCIsPrimary[kNMaxMCParticles];
  int fMCParticlePdgCode[kNMaxMCParticles];
  double fMCParticleTrueEnergy [kNMaxMCParticles]; 
  int fMCParticleTrackID[kNMaxMCParticles]; 
  int fMCParticleParentTrackID[kNMaxMCParticles]; 
  std::string fMCParticleStartProcess[kNMaxMCParticles]; 
  std::string fMCParticleEndProcess[kNMaxMCParticles]; 
  int fMCParticleNTrajectoryPoint[kNMaxMCParticles]; 
  double fMCParticleStartPositionX[kNMaxMCParticles];
  double fMCParticleStartPositionY[kNMaxMCParticles];
  double fMCParticleStartPositionZ[kNMaxMCParticles];
  double fMCParticleStartPositionT[kNMaxMCParticles];
  double fMCParticleStartMomentumX[kNMaxMCParticles];
  double fMCParticleStartMomentumY[kNMaxMCParticles];
  double fMCParticleStartMomentumZ[kNMaxMCParticles];
  double fMCParticleStartMomentumE[kNMaxMCParticles];
  double fMCParticleEndPositionX[kNMaxMCParticles];
  double fMCParticleEndPositionY[kNMaxMCParticles];
  double fMCParticleEndPositionZ[kNMaxMCParticles];
  double fMCParticleEndPositionT[kNMaxMCParticles];
  double fMCParticleEndMomentumX[kNMaxMCParticles];
  double fMCParticleEndMomentumY[kNMaxMCParticles];
  double fMCParticleEndMomentumZ[kNMaxMCParticles];
  double fMCParticleEndMomentumE[kNMaxMCParticles];
  double fMCParticleVertexTime[kNMaxMCParticles];
  double fMCParticleEndTime[kNMaxMCParticles];
  int fMCParticleNHits[kNMaxMCParticles];
  int fMCParticleNHitsView[kNMaxMCParticles][kNViews];
  //int fMCPfoMatchedPosition[kNMaxMCParticles];

  int fPFPID[kNMaxPFParticles];
  bool fPFPIsPrimary[kNMaxPFParticles];
  int fPFPTrueParticleMatchedID[kNMaxPFParticles];
  int fPFPTrueParticleMatchedPosition[kNMaxPFParticles];
  int fPFPParentID[kNMaxPFParticles];
  int fPFPPdgCode[kNMaxPFParticles];
  int fPFPNChildren[kNMaxPFParticles];
  int fPFPNClusters[kNMaxPFParticles];
  int fPFPNHits[kNMaxPFParticles];
  int fPFPNHitsView[kNMaxPFParticles][kNViews];
  int fPFPNSharedTrueParticleHits[kNMaxPFParticles];
  int fPFPNSharedTrueParticleHitsView[kNMaxPFParticles][kNViews];
  int fPFPTrueParticleMatchedIDView[kNMaxPFParticles][kNViews];
  int fPFPTrueParticleMatchedPositionView[kNMaxPFParticles][kNViews];
  bool fPFPIsTrack[kNMaxPFParticles];
  bool fPFPIsShower[kNMaxPFParticles];
  int fPFPCluPlane[kNMaxPFParticles][kNMaxPFPClusters];
  int fPFPCluView[kNMaxPFParticles][kNMaxPFPClusters];
  int fPFPCluNHits[kNMaxPFParticles][kNMaxPFPClusters];
  double fPFPCluIntegral[kNMaxPFParticles][kNMaxPFPClusters];


  int fPFPTrackID[kNMaxPFParticles];
  double fPFPTrackLength[kNMaxPFParticles];
  double fPFPTrackStartX[kNMaxPFParticles];
  double fPFPTrackStartY[kNMaxPFParticles];
  double fPFPTrackStartZ[kNMaxPFParticles];
  double fPFPTrackVertexX[kNMaxPFParticles];
  double fPFPTrackVertexY[kNMaxPFParticles];
  double fPFPTrackVertexZ[kNMaxPFParticles];
  double fPFPTrackEndX[kNMaxPFParticles];
  double fPFPTrackEndY[kNMaxPFParticles];
  double fPFPTrackEndZ[kNMaxPFParticles];
  double fPFPTrackTheta[kNMaxPFParticles];
  double fPFPTrackPhi[kNMaxPFParticles];
  double fPFPTrackZenithAngle[kNMaxPFParticles];
  double fPFPTrackAzimuthAngle[kNMaxPFParticles];
  double fPFPTrackStartDirectionX[kNMaxPFParticles];
  double fPFPTrackStartDirectionY[kNMaxPFParticles];
  double fPFPTrackStartDirectionZ[kNMaxPFParticles];
  double fPFPTrackVertexDirectionX[kNMaxPFParticles];
  double fPFPTrackVertexDirectionY[kNMaxPFParticles];
  double fPFPTrackVertexDirectionZ[kNMaxPFParticles];
  double fPFPTrackEndDirectionX[kNMaxPFParticles];
  double fPFPTrackEndDirectionY[kNMaxPFParticles];
  double fPFPTrackEndDirectionZ[kNMaxPFParticles];
  float fPFPTrackChi2[kNMaxPFParticles];
  int fPFPTrackNdof[kNMaxPFParticles];
 
  int    fPFPShowerID[kNMaxPFParticles];
  int    fPFPShowerBestPlane[kNMaxPFParticles];
  double fPFPShowerDirectionX[kNMaxPFParticles];
  double fPFPShowerDirectionY[kNMaxPFParticles];
  double fPFPShowerDirectionZ[kNMaxPFParticles];
  double fPFPShowerDirectionErrX[kNMaxPFParticles];
  double fPFPShowerDirectionErrY[kNMaxPFParticles];
  double fPFPShowerDirectionErrZ[kNMaxPFParticles];
  double fPFPShowerStartX[kNMaxPFParticles];
  double fPFPShowerStartY[kNMaxPFParticles];
  double fPFPShowerStartZ[kNMaxPFParticles];
  double fPFPShowerStartErrX[kNMaxPFParticles];
  double fPFPShowerStartErrY[kNMaxPFParticles];
  double fPFPShowerStartErrZ[kNMaxPFParticles];
  double fPFPShowerLength[kNMaxPFParticles];
  double fPFPShowerOpenAngle[kNMaxPFParticles];

  double fPFPCompleteness[kNMaxMCParticles];
  double fPFPCompletenessView[kNMaxMCParticles][kNViews];
  double fPFPPurity[kNMaxMCParticles];
  double fPFPPurityView[kNMaxMCParticles][kNViews];

  std::string fTruthLabel;
  std::string fHitLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPFParticleLabel;
  bool fRollUpUnsavedIDs;

  const geo::Geometry* fGeom;
};


test::pandoraAnalysis::pandoraAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fTruthLabel = p.get<std::string>("TruthLabel");
  fHitLabel = p.get<std::string>("HitLabel");
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fShowerLabel = p.get<std::string>("ShowerLabel");
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
  fRollUpUnsavedIDs = p.get<bool>("RollUpUnsavedIDs"); 
} 

void test::pandoraAnalysis::analyze(art::Event const& e)
{
  reset(); //Don't deep clean
  const art::ServiceHandle<cheat::BackTrackerService> btServ;
  fEventID = e.id().event();
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  std::cout << "=============== EVENT ID " << fEventID << " == RUN ID " << fRunID << " == SUBRUN ID " << fSubRunID << " ================" << std::endl;


  //Get all hits
  std::vector<art::Ptr<recob::Hit> > allHits;
  auto hitHandle = e.getHandle<std::vector<recob::Hit>>(fHitLabel);
  if (hitHandle)
  {art::fill_ptr_vector(allHits, hitHandle);}

  //Fill MC particle to hits map
  std::map<int,int> trueParticleHits, trueParticleHitsView0, trueParticleHitsView1, trueParticleHitsView2;
  for (const auto& hit: allHits){
      TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(clockData, hit, fRollUpUnsavedIDs));
      if (TruthMatchUtils::Valid(g4ID)){
          trueParticleHits[g4ID]++;
          if(hit->View()==0)trueParticleHitsView0[g4ID]++;
          if(hit->View()==1)trueParticleHitsView1[g4ID]++;
          if(hit->View()==2)trueParticleHitsView2[g4ID]++;
      }
  }

  //Access the truth information
  if(!e.isRealData()){
    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    if(mcParticles.isValid()){
 
     fNMCParticles=mcParticles->size();
      bool isMCPrimary(false);
      for(unsigned int iMc=0; iMc< mcParticles->size(); iMc++){
        const simb::MCParticle trueParticle = mcParticles->at(iMc);
        fMCParticleTrueEnergy[iMc]=trueParticle.E();
        fMCParticlePdgCode[iMc]=trueParticle.PdgCode();
        fMCParticleTrackID[iMc]=trueParticle.TrackId();
        fMCParticleVertexTime[iMc]=trueParticle.T();
        fMCParticleEndTime[iMc]=trueParticle.EndT();
        fMCParticleParentTrackID[iMc]=trueParticle.Mother();
        fMCParticleStartProcess[iMc]=trueParticle.Process();
        fMCParticleEndProcess[iMc]=trueParticle.EndProcess();
        fMCParticleNTrajectoryPoint[iMc]=trueParticle.NumberTrajectoryPoints();
        trueParticle.Process()=="primary"? isMCPrimary=true:isMCPrimary=false;
        fMCIsPrimary[iMc] = isMCPrimary;
        fMCParticleNHits[iMc]=trueParticleHits[trueParticle.TrackId()];
        fMCParticleNHitsView[iMc][0]=trueParticleHitsView0[trueParticle.TrackId()];
        fMCParticleNHitsView[iMc][1]=trueParticleHitsView1[trueParticle.TrackId()];
        fMCParticleNHitsView[iMc][2]=trueParticleHitsView2[trueParticle.TrackId()];
        fMCParticleStartPositionX[iMc] = trueParticle.Position().X();
        fMCParticleStartPositionY[iMc] = trueParticle.Position().Y();
        fMCParticleStartPositionZ[iMc] = trueParticle.Position().Z();
        fMCParticleStartPositionT[iMc] = trueParticle.Position().T();
        fMCParticleEndPositionX[iMc] = trueParticle.EndPosition().X();
        fMCParticleEndPositionY[iMc] = trueParticle.EndPosition().Y();
        fMCParticleEndPositionZ[iMc] = trueParticle.EndPosition().Z();
        fMCParticleEndPositionT[iMc] = trueParticle.EndPosition().T();

        fMCParticleStartMomentumX[iMc] = trueParticle.Momentum().X();
        fMCParticleStartMomentumY[iMc] = trueParticle.Momentum().Y();
        fMCParticleStartMomentumZ[iMc] = trueParticle.Momentum().Z();
        fMCParticleStartMomentumE[iMc] = trueParticle.Momentum().E();
        fMCParticleEndMomentumX[iMc] = trueParticle.EndMomentum().X();
        fMCParticleEndMomentumY[iMc] = trueParticle.EndMomentum().Y();
        fMCParticleEndMomentumZ[iMc] = trueParticle.EndMomentum().Z();
        fMCParticleEndMomentumE[iMc] = trueParticle.EndMomentum().E();

      }
    }
  }

  //Access the PFParticles from Pandora
  const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(e,fPFParticleLabel);
  fNPFParticles = pfparticleVect.size();
  if(!fNPFParticles) {
    std::cout << "No PFParticles found!" << std::endl;
    return;
  }

  //Access the Clusters
  std::vector<art::Ptr<recob::Cluster>> clusterVect;
  auto clusterHandle = e.getHandle<std::vector<recob::Cluster> >(fPFParticleLabel);
  if (clusterHandle)
    art::fill_ptr_vector(clusterVect,clusterHandle);

  art::FindManyP<recob::Cluster> clusterParticleAssoc(pfparticleVect, e, fPFParticleLabel);

  auto trackHandle = e.getHandle<std::vector<recob::Track> >(fTrackLabel);
  if (!trackHandle){
    std::cout<<"Unable to find std::vector<recob::Track> with module label: " << fTrackLabel << std::endl;
    return;
  }
  std::vector<art::Ptr<recob::Track> > trackList;
  art::fill_ptr_vector(trackList, trackHandle);

  //std::map<int,int> pfpToMcMap;
  int iPfp(0);
  for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){

    fPFPID[iPfp]=iPfp;
    fPFPIsPrimary[iPfp]=pfp->IsPrimary();
    fPFPPdgCode[iPfp]=pfp->PdgCode();
    fPFPNChildren[iPfp]=pfp->NumDaughters();
    (pfp->IsPrimary())? fPFPParentID[iPfp]=-1 : fPFPParentID[iPfp]=pfp->Parent();

    std::vector<art::Ptr<recob::Cluster>> pfpClusters = clusterParticleAssoc.at(pfp.key());
    fPFPNClusters[iPfp]=pfpClusters.size();
    if(!pfpClusters.empty()){
      int iClu(0);
      for(const art::Ptr<recob::Cluster> &clu:pfpClusters){
	fPFPCluPlane[iPfp][iClu]=clu->Plane().asPlaneID().Plane;
	fPFPCluView[iPfp][iClu]=clu->View();
        fPFPCluNHits[iPfp][iClu]=clu->NHits();
	fPFPCluIntegral[iPfp][iClu]=clu->Integral();
        iClu++;
        if (iClu == kNMaxPFPClusters)
            break;
      }
    }

    std::vector<art::Ptr<recob::Hit>> pfpHits;

    if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp,e,fPFParticleLabel,fTrackLabel)){
      fPFPIsTrack[iPfp]=true;
      art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp,e,fPFParticleLabel, fTrackLabel); 

      fPFPTrackID[iPfp]=track->ID();
      fPFPTrackLength[iPfp]=track->Length();
      fPFPTrackStartX[iPfp]=track->Start().X();
      fPFPTrackStartY[iPfp]=track->Start().Y();
      fPFPTrackStartZ[iPfp]=track->Start().Z();
      fPFPTrackVertexX[iPfp]=track->Vertex().X();
      fPFPTrackVertexY[iPfp]=track->Vertex().Y();
      fPFPTrackVertexZ[iPfp]=track->Vertex().Z();
      fPFPTrackEndX[iPfp]=track->End().X();
      fPFPTrackEndY[iPfp]=track->End().Y();
      fPFPTrackEndZ[iPfp]=track->End().Z();
      fPFPTrackTheta[iPfp]=track->Theta();
      fPFPTrackPhi[iPfp]=track->Phi();
      fPFPTrackZenithAngle[iPfp]=track->ZenithAngle();
      fPFPTrackAzimuthAngle[iPfp]=track->AzimuthAngle();
      fPFPTrackStartDirectionX[iPfp]=track->StartDirection().X();
      fPFPTrackStartDirectionY[iPfp]=track->StartDirection().Y();
      fPFPTrackStartDirectionZ[iPfp]=track->StartDirection().Z();
      fPFPTrackVertexDirectionX[iPfp]=track->VertexDirection().X();
      fPFPTrackVertexDirectionY[iPfp]=track->VertexDirection().Y();
      fPFPTrackVertexDirectionZ[iPfp]=track->VertexDirection().Z();
      fPFPTrackEndDirectionX[iPfp]=track->EndDirection().X();
      fPFPTrackEndDirectionY[iPfp]=track->EndDirection().Y();
      fPFPTrackEndDirectionZ[iPfp]=track->EndDirection().Z();
      fPFPTrackChi2[iPfp]=track->Chi2();
      fPFPTrackNdof[iPfp]=track->Ndof();

      pfpHits = dune_ana::DUNEAnaTrackUtils::GetHits(track,e,fTrackLabel);
      
      fPFPNHits[iPfp]=pfpHits.size();
      std::vector<art::Ptr<recob::Hit> > pfpHitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp,e,fPFParticleLabel, 0);
      fPFPNHitsView[iPfp][0] = pfpHitsView0.size();
      std::vector<art::Ptr<recob::Hit> > pfpHitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp,e,fPFParticleLabel, 1);
      fPFPNHitsView[iPfp][1] = pfpHitsView1.size();
      std::vector<art::Ptr<recob::Hit> > pfpHitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp,e,fPFParticleLabel, 2);
      fPFPNHitsView[iPfp][2] = pfpHitsView2.size();

      if(!e.isRealData()) {
        TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHits,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID)){
          fPFPTrueParticleMatchedID[iPfp] = g4ID;

	  int pos(999999); for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) {
	    if(fMCParticleTrackID[ipos]==g4ID)pos=ipos;
 	  }
          fPFPTrueParticleMatchedPosition[iPfp] = pos;

        }
        TruthMatchUtils::G4ID g4IDView0(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView0,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID))
        {
            fPFPTrueParticleMatchedIDView[iPfp][0] = g4IDView0;
         	  for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
                if (fMCParticleTrackID[ipos] != g4IDView0)
                    continue;
                fPFPTrueParticleMatchedPositionView[iPfp][0] = ipos;
                break;
            }
        }
        TruthMatchUtils::G4ID g4IDView1(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView1,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID))
        {
            fPFPTrueParticleMatchedIDView[iPfp][1] = g4IDView1;
         	  for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
                if (fMCParticleTrackID[ipos] != g4IDView1)
                    continue;
                fPFPTrueParticleMatchedPositionView[iPfp][1] = ipos;
                break;
            }
        }
        TruthMatchUtils::G4ID g4IDView2(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView2,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID))
        {
            fPFPTrueParticleMatchedIDView[iPfp][2] = g4IDView2;
         	  for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
                if (fMCParticleTrackID[ipos] != g4IDView2)
                    continue;
                fPFPTrueParticleMatchedPositionView[iPfp][2] = ipos;
                break;
            }
        }

      }
    }

    if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp,e,fPFParticleLabel,fShowerLabel)){
      fPFPIsShower[iPfp]=true;
      art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp,e,fPFParticleLabel, fShowerLabel); 
      fPFPShowerID[iPfp]=shower->ID();
      fPFPShowerBestPlane[iPfp]=shower->best_plane();
      fPFPShowerDirectionX[iPfp]=shower->Direction().X();
      fPFPShowerDirectionY[iPfp]=shower->Direction().Y();
      fPFPShowerDirectionZ[iPfp]=shower->Direction().Z();
      fPFPShowerDirectionErrX[iPfp]=shower->DirectionErr().X();
      fPFPShowerDirectionErrY[iPfp]=shower->DirectionErr().Y();
      fPFPShowerDirectionErrZ[iPfp]=shower->DirectionErr().Z();
      fPFPShowerStartX[iPfp]=shower->ShowerStart().X();
      fPFPShowerStartY[iPfp]=shower->ShowerStart().Y();
      fPFPShowerStartZ[iPfp]=shower->ShowerStart().Z();
      fPFPShowerStartErrX[iPfp]=shower->ShowerStartErr().X();
      fPFPShowerStartErrY[iPfp]=shower->ShowerStartErr().Y();
      fPFPShowerStartErrZ[iPfp]=shower->ShowerStartErr().Z();
      fPFPShowerLength[iPfp]=shower->Length();
      fPFPShowerOpenAngle[iPfp]=shower->OpenAngle();

      pfpHits = dune_ana::DUNEAnaShowerUtils::GetHits(shower,e,fShowerLabel);
      fPFPNHits[iPfp]=pfpHits.size();
      std::vector<art::Ptr<recob::Hit> > pfpHitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp,e,fPFParticleLabel, 0);
      fPFPNHitsView[iPfp][0] = pfpHitsView0.size();
      std::vector<art::Ptr<recob::Hit> > pfpHitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp,e,fPFParticleLabel, 1);
      fPFPNHitsView[iPfp][1] = pfpHitsView1.size();
      std::vector<art::Ptr<recob::Hit> > pfpHitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp,e,fPFParticleLabel, 2);
      fPFPNHitsView[iPfp][2] = pfpHitsView2.size();

      if(!e.isRealData()) {
        TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHits,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID)){
          fPFPTrueParticleMatchedID[iPfp] = g4ID;
	  int pos(999999); for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) {
	    if(fMCParticleTrackID[ipos]==g4ID)pos=ipos;
  	  }
          fPFPTrueParticleMatchedPosition[iPfp] = pos;
        }
        TruthMatchUtils::G4ID g4IDView0(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView0,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID))
        {
            fPFPTrueParticleMatchedIDView[iPfp][0] = g4IDView0;
         	  for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
                if (fMCParticleTrackID[ipos] != g4IDView0)
                    continue;
                fPFPTrueParticleMatchedPositionView[iPfp][0] = ipos;
                break;
            }
        }
        TruthMatchUtils::G4ID g4IDView1(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView1,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID))
        {
            fPFPTrueParticleMatchedIDView[iPfp][1] = g4IDView1;
         	  for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
                if (fMCParticleTrackID[ipos] != g4IDView1)
                    continue;
                fPFPTrueParticleMatchedPositionView[iPfp][1] = ipos;
                break;
            }
        }
        TruthMatchUtils::G4ID g4IDView2(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView2,fRollUpUnsavedIDs));
        if (TruthMatchUtils::Valid(g4ID))
        {
            fPFPTrueParticleMatchedIDView[iPfp][2] = g4IDView2;
         	  for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
                if (fMCParticleTrackID[ipos] != g4IDView2)
                    continue;
                fPFPTrueParticleMatchedPositionView[iPfp][2] = ipos;
                break;
            }
        }

      }
    }

    if(!e.isRealData()) {
      //Fill shared MC particle to hits map for this specific Pfp
      for (const auto& hit: pfpHits){
          TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(clockData, hit, fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID)){
              if(g4ID==fPFPTrueParticleMatchedID[iPfp])++fPFPNSharedTrueParticleHits[iPfp];
              if(g4ID==fPFPTrueParticleMatchedID[iPfp] && hit->View()==0)++fPFPNSharedTrueParticleHitsView[iPfp][0];
              if(g4ID==fPFPTrueParticleMatchedID[iPfp] && hit->View()==1)++fPFPNSharedTrueParticleHitsView[iPfp][1];
              if(g4ID==fPFPTrueParticleMatchedID[iPfp] && hit->View()==2)++fPFPNSharedTrueParticleHitsView[iPfp][2];            
          }
      }

      if(fPFPNHits[iPfp] > 0 && fPFPNHits[iPfp] < 999999) fPFPPurity[iPfp]         = (float)fPFPNSharedTrueParticleHits[iPfp] / fPFPNHits[iPfp];
      if(fPFPNHitsView[iPfp][0] > 0 && fPFPNHitsView[iPfp][0] < 999999) fPFPPurityView[iPfp][0]        = (float)fPFPNSharedTrueParticleHitsView[iPfp][0] / fPFPNHitsView[iPfp][0];
      if(fPFPNHitsView[iPfp][1] > 0 && fPFPNHitsView[iPfp][1] < 999999) fPFPPurityView[iPfp][1]         = (float)fPFPNSharedTrueParticleHitsView[iPfp][1] / fPFPNHitsView[iPfp][1];
      if(fPFPNHitsView[iPfp][2] > 0 && fPFPNHitsView[iPfp][2] < 999999) fPFPPurityView[iPfp][2]         = (float)fPFPNSharedTrueParticleHitsView[iPfp][2] / fPFPNHitsView[iPfp][2];

      if(fPFPTrueParticleMatchedPosition[iPfp]<999999 && fMCParticleNHits[fPFPTrueParticleMatchedPosition[iPfp]] > 0 && fMCParticleNHits[fPFPTrueParticleMatchedPosition[iPfp]] < 999999) fPFPCompleteness[iPfp]  = (float)fPFPNSharedTrueParticleHits[iPfp] / fMCParticleNHits[fPFPTrueParticleMatchedPosition[iPfp]];
      if(fPFPTrueParticleMatchedPosition[iPfp]<999999 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][0] > 0 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][0] < 999999) fPFPCompletenessView[iPfp][0]  = (float)fPFPNSharedTrueParticleHitsView[iPfp][0] / fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][0];
      if(fPFPTrueParticleMatchedPosition[iPfp]<999999 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][1] > 0 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][1] < 999999) fPFPCompletenessView[iPfp][1]  = (float)fPFPNSharedTrueParticleHitsView[iPfp][1] / fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][1];
      if(fPFPTrueParticleMatchedPosition[iPfp]<999999 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][2] > 0 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][2] < 999999) fPFPCompletenessView[iPfp][2]  = (float)fPFPNSharedTrueParticleHitsView[iPfp][2] / fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][2];

    }
    iPfp++;
  }
  fTree->Fill();
}

void test::pandoraAnalysis::beginJob()
{

  //deep clean the variables
  reset(true);
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("pandoraOutput","Pandora Output Tree");

  //Event branches
  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("runID",&fRunID,"runID/i");
  fTree->Branch("subrunID",&fSubRunID,"subrunID/i");
  fTree->Branch("nMCParticles",&fNMCParticles,"nMCParticles/i");
  fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");

  //MC truth branches
  //fTree->Branch("mcIsMCPrimary",&fMCIsPrimary);
  fTree->Branch("mcIsMCPrimary",&fMCIsPrimary,"MCIsPrimary[nMCParticles]/O");
  fTree->Branch("mcParticlePdgCode",&fMCParticlePdgCode,"MCParticlePdgCode[nMCParticles]/I");
  fTree->Branch("mcParticleTrueEnergy",&fMCParticleTrueEnergy,"MCParticleTrueEnergy[nMCParticles]/D");
  fTree->Branch("mcParticleTrackID",&fMCParticleTrackID,"MCParticleTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleParentTrackID",&fMCParticleParentTrackID,"MCParticleParentTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleStartProcess",&fMCParticleStartProcess,"MCParticleStartProcess[nMCParticles]/C");
  fTree->Branch("mcParticleEndProcess",&fMCParticleEndProcess, "MCParticleEndProcess[nMCParticles]/C");
  fTree->Branch("mcParticleNTrajectoryPoints",&fMCParticleNTrajectoryPoint, "MCParticleNTrajectoryPoint[nMCParticles]/I");
  fTree->Branch("mcParticleStartPositionX",&fMCParticleStartPositionX,"MCParticleStartPositionX[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionY",&fMCParticleStartPositionY,"MCParticleStartPositionY[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionZ",&fMCParticleStartPositionZ,"MCParticleStartPositionZ[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionT",&fMCParticleStartPositionT,"MCParticleStartPositionT[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumX",&fMCParticleStartMomentumX, "MCParticleStartMomentumX[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumY",&fMCParticleStartMomentumY, "MCParticleStartMomentumY[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumZ",&fMCParticleStartMomentumZ, "MCParticleStartMomentumZ[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumE",&fMCParticleStartMomentumE, "MCParticleStartMomentumE[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionX",&fMCParticleEndPositionX, "MCParticleEndPositionX[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionY",&fMCParticleEndPositionY, "MCParticleEndPositionY[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionZ",&fMCParticleEndPositionZ, "MCParticleEndPositionZ[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionT",&fMCParticleEndPositionT, "MCParticleEndPositionT[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumX",&fMCParticleEndMomentumX, "MCParticleEndMomentumX[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumY",&fMCParticleEndMomentumY, "MCParticleEndMomentumY[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumZ",&fMCParticleEndMomentumZ, "MCParticleEndMomentumZ[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumE",&fMCParticleEndMomentumE, "MCParticleEndMomentumE[nMCParticles]/D");
  fTree->Branch("mcParticleVertexTime",&fMCParticleVertexTime, "MCParticleVertexTime[nMCParticles]/D");
  fTree->Branch("mcParticleEndTime",&fMCParticleEndTime, "MCParticleEndTime[nMCParticles]/D");
  fTree->Branch("mcParticleNHits", &fMCParticleNHits, "MCParticleNHits[nMCParticles]/I");
  fTree->Branch("mcParticleNHitsView", &fMCParticleNHitsView, "MCParticleNHitsView[nMCParticles][3]/I");
  //fTree->Branch("mcPfoMatchedPosition", &fMCPfoMatchedPosition, "fMCPfoMatchedPosition[nMCParticles]/I");

  //PFP branches
  fTree->Branch("pfpTrueParticleMatchedID",&fPFPTrueParticleMatchedID,"PFPTrueParticleMatchedID[nPFParticles]/I");
  fTree->Branch("pfpTrueParticleMatchedPosition",&fPFPTrueParticleMatchedPosition,"PFPTrueParticleMatchedPosition[nPFParticles]/I");
  fTree->Branch("pfpIsPrimary",&fPFPIsPrimary,"PFPIsPrimary[nPFParticles]/O");
  fTree->Branch("pfpID",&fPFPID, "PFPID[nPFParticles]/I");
  fTree->Branch("pfpParentID",&fPFPParentID, "PFPParentID[nPFParticles]/I");
  fTree->Branch("pfpPdgCode",&fPFPPdgCode, "PFPPdgCode[nPFParticles]/I");
  fTree->Branch("pfpNChildren",&fPFPNChildren,"PFPNChildren[nPFParticles]/I");
  fTree->Branch("pfpNClusters",&fPFPNClusters,"PFPNClusters[nPFParticles]/I");
  fTree->Branch("pfpNHits",&fPFPNHits,"PFPNHits[nPFParticles]/I");
  fTree->Branch("pfpNHitsView",&fPFPNHitsView,"PFPNHitsView[nPFParticles][3]/I");
  fTree->Branch("pfpNSharedTrueParticleHits",&fPFPNSharedTrueParticleHits,"PFPNSharedTrueParticleHits[nPFParticles]/I");
  fTree->Branch("pfpNSharedTrueParticleHitsView",&fPFPNSharedTrueParticleHitsView,"PFPNSharedTrueParticleHitsView[nPFParticles][3]/I");
  fTree->Branch("pfpTrueParticleMatchedIDView",&fPFPTrueParticleMatchedIDView,"PFPTrueParticleMatchedIDView[nPFParticles][3]/I");
  fTree->Branch("pfpTrueParticleMatchedPositionView",&fPFPTrueParticleMatchedPositionView,"PFPTrueParticleMatchedPositionView[nPFParticles][3]/I");
  fTree->Branch("pfpIsTrack",&fPFPIsTrack,"PFPIsTrack[nPFParticles]/O");
  fTree->Branch("pfpIsShower",&fPFPIsShower,"PFPIsShower[nPFParticles]/O");
  fTree->Branch("pfpTrackID", &fPFPTrackID,"PFPNClusters[nPFParticles]/I");
  fTree->Branch("pfpTrackLength",&fPFPTrackLength,"PFPTrackLength[nPFParticles]/D");
  fTree->Branch("pfpTrackStartX",&fPFPTrackStartX,"PFPTrackStartX[nPFParticles]/D");
  fTree->Branch("pfpTrackStartY",&fPFPTrackStartY,"PFPTrackStartY[nPFParticles]/D");
  fTree->Branch("pfpTrackStartZ",&fPFPTrackStartZ,"PFPTrackStartZ[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexX",&fPFPTrackVertexX,"PFPTrackVertexX[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexY",&fPFPTrackVertexY,"PFPTrackVertexY[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexZ",&fPFPTrackVertexZ,"PFPTrackVertexZ[nPFParticles]/D");
  fTree->Branch("pfpTrackEndX",&fPFPTrackEndX,"PFPTrackEndX[nPFParticles]/D");
  fTree->Branch("pfpTrackEndY",&fPFPTrackEndY,"PFPTrackEndY[nPFParticles]/D");
  fTree->Branch("pfpTrackEndZ",&fPFPTrackEndZ,"PFPTrackEndZ[nPFParticles]/D");
  fTree->Branch("pfpTrackTheta",&fPFPTrackTheta,"PFPTrackTheta[nPFParticles]/D");
  fTree->Branch("pfpTrackPhi",&fPFPTrackPhi,"PFPTrackPhi[nPFParticles]/D");
  fTree->Branch("pfpTrackZenithAngle",&fPFPTrackZenithAngle,"PFPTrackZenithAngle[nPFParticles]/D");
  fTree->Branch("pfpTrackAzimuthAngle",&fPFPTrackAzimuthAngle,"PFPTrackAzimuthAngle[nPFParticles]/D");
  fTree->Branch("pfpTrackStartDirectionX",&fPFPTrackStartDirectionX,"PFPTrackStartDirectionX[nPFParticles]/D");
  fTree->Branch("pfpTrackStartDirectionY",&fPFPTrackStartDirectionY,"PFPTrackStartDirectionY[nPFParticles]/D");
  fTree->Branch("pfpTrackStartDirectionZ",&fPFPTrackStartDirectionZ,"PFPTrackStartDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexDirectionX",&fPFPTrackVertexDirectionX,"PFPTrackVertexDirectionX[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexDirectionY",&fPFPTrackVertexDirectionY,"PFPTrackVertexDirectionY[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexDirectionZ",&fPFPTrackVertexDirectionZ,"PFPTrackVertexDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpTrackEndDirectionX",&fPFPTrackEndDirectionX,"PFPTrackEndDirectionX[nPFParticles]/D");
  fTree->Branch("pfpTrackEndDirectionY",&fPFPTrackEndDirectionY,"PFPTrackEndDirectionY[nPFParticles]/D");
  fTree->Branch("pfpTrackEndDirectionZ",&fPFPTrackEndDirectionZ,"PFPTrackEndDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpTrackChi2",&fPFPTrackChi2,"PFPTrackChi2[nPFParticles]/F");
  fTree->Branch("pfpTrackStartNdof",&fPFPTrackNdof,"PFPTrackNdof[nPFParticles]/I");

  fTree->Branch("pfpCluPlane",fPFPCluPlane,"PFPCluPlane[nPFParticles][100]/I");
  fTree->Branch("pfpCluView",fPFPCluView,"PFPCluView[nPFParticles][100]/I");
  fTree->Branch("pfpCluNHits",fPFPCluNHits,"PFPCluNHits[nPFParticles][100]/I");
  fTree->Branch("pfpCluIntegral",fPFPCluIntegral,"PFPCluIntegral[nPFParticles][100]/D");

  fTree->Branch("pfpShowerID",&fPFPShowerID,"PFPShowerID[nPFParticles]/I");
  fTree->Branch("pfpShowerBestPlane",&fPFPShowerBestPlane,"PFPShowerBestPlane[nPFParticles]/I");
  fTree->Branch("pfpShowerDirectionX",&fPFPShowerDirectionX,"PFPShowerDirectionX[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionY",&fPFPShowerDirectionY,"PFPShowerDirectionY[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionZ",&fPFPShowerDirectionZ,"PFPShowerDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionErrX",&fPFPShowerDirectionErrX,"PFPShowerDirectionErrX[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionErrY",&fPFPShowerDirectionErrY,"PFPShowerDirectionErrY[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionErrZ",&fPFPShowerDirectionErrZ,"PFPShowerDirectionErrZ[nPFParticles]/D");
  fTree->Branch("pfpShowerStartX",&fPFPShowerStartX,"PFPShowerStartX[nPFParticles]/D");
  fTree->Branch("pfpShowerStartY",&fPFPShowerStartY,"PFPShowerStartY[nPFParticles]/D");
  fTree->Branch("pfpShowerStartZ",&fPFPShowerStartZ,"PFPShowerStartZ[nPFParticles]/D");
  fTree->Branch("pfpShowerStartErrX",&fPFPShowerStartErrX,"PFPShowerStartErrX[nPFParticles]/D");
  fTree->Branch("pfpShowerStartErrY",&fPFPShowerStartErrY,"PFPShowerStartErrY[nPFParticles]/D");
  fTree->Branch("pfpShowerStartErrZ",&fPFPShowerStartErrZ,"PFPShowerStartErrZ[nPFParticles]/D");
  fTree->Branch("pfpShowerLength",&fPFPShowerLength,"PFPShowerLength[nPFParticles]/D");
  fTree->Branch("pfpShowerOpeningAngle",&fPFPShowerOpenAngle,"PFPShowerOpenAngle[nPFParticles]/D");

  fTree->Branch("pfpCompleteness", &fPFPCompleteness, "PFPCompleteness[nPFParticles]/D");
  fTree->Branch("pfpCompletenessView", &fPFPCompletenessView, "PFPCompletenessView[nPFParticles][3]/D");
  fTree->Branch("pfpPurity", &fPFPPurity, "PFPPurity[nPFParticles]/D");
  fTree->Branch("pfpPurityView", &fPFPPurityView, "PFPPurityView[nPFParticles][3]/D");
}
void test::pandoraAnalysis::endJob()
{
  // Implementation of optional member function here.
}

void test::pandoraAnalysis::reset(bool deepClean)
{
    for(unsigned int iMc=0; iMc<(deepClean ? kNMaxMCParticles : fNMCParticles); iMc++){
      fMCIsPrimary[iMc]=0;
      fMCParticlePdgCode[iMc]=0;
      fMCParticleTrueEnergy[iMc]=-999999;
      fMCParticleTrackID[iMc]=999999;
      fMCParticleParentTrackID[iMc]=999999;
      fMCParticleStartProcess[iMc]="";
      fMCParticleEndProcess[iMc]="";
      fMCParticleNTrajectoryPoint[iMc]=999999;
      fMCParticleStartPositionX[iMc]=999999;
      fMCParticleStartPositionY[iMc]=999999;
      fMCParticleStartPositionZ[iMc]=999999;
      fMCParticleStartPositionT[iMc]=999999;
      fMCParticleStartMomentumX[iMc]=999999;
      fMCParticleStartMomentumY[iMc]=999999;
      fMCParticleStartMomentumZ[iMc]=999999;
      fMCParticleStartMomentumE[iMc]=999999;
      fMCParticleEndPositionX[iMc]=999999;
      fMCParticleEndPositionY[iMc]=999999;
      fMCParticleEndPositionZ[iMc]=999999;
      fMCParticleEndPositionT[iMc]=999999;
      fMCParticleEndMomentumX[iMc]=999999;
      fMCParticleEndMomentumY[iMc]=999999;
      fMCParticleEndMomentumZ[iMc]=999999;
      fMCParticleEndMomentumE[iMc]=999999;
      fMCParticleVertexTime[iMc]=999999;
      fMCParticleEndTime[iMc]=999999;
      fMCParticleNHits[iMc]=999999;
      fMCParticleNHitsView[iMc][0]=999999;
      fMCParticleNHitsView[iMc][1]=999999;
      fMCParticleNHitsView[iMc][2]=999999;
      //fMCPfoMatchedPosition[iMc]=999999;
    }
    fNMCParticles = 0;

    for(unsigned int iPfp=0; iPfp<(deepClean ? kNMaxPFParticles : fNPFParticles); iPfp++){

     fPFPID[iPfp]=999999;
     fPFPTrueParticleMatchedID[iPfp]=999999;
     fPFPTrueParticleMatchedPosition[iPfp]=999999;

     fPFPNHits[iPfp]=999999;
     fPFPNSharedTrueParticleHits[iPfp]=0;

     fPFPNClusters[iPfp]=999999;
     fPFPIsTrack[iPfp]=0;
     fPFPIsShower[iPfp]=0;

     fPFPTrackID[iPfp]=999999;
     fPFPTrackLength[iPfp]=999999;
     fPFPTrackStartX[iPfp]=999999;
     fPFPTrackStartY[iPfp]=999999;
     fPFPTrackStartZ[iPfp]=999999;
     fPFPTrackVertexX[iPfp]=999999;
     fPFPTrackVertexY[iPfp]=999999;
     fPFPTrackVertexZ[iPfp]=999999;
     fPFPTrackEndX[iPfp]=999999;
     fPFPTrackEndY[iPfp]=999999;
     fPFPTrackEndZ[iPfp]=999999;
     fPFPTrackTheta[iPfp]=999999;
     fPFPTrackPhi[iPfp]=999999;
     fPFPTrackZenithAngle[iPfp]=999999;
     fPFPTrackAzimuthAngle[iPfp]=999999;
     fPFPTrackStartDirectionX[iPfp]=999999;
     fPFPTrackStartDirectionY[iPfp]=999999;
     fPFPTrackStartDirectionZ[iPfp]=999999;
     fPFPTrackVertexDirectionX[iPfp]=999999;
     fPFPTrackVertexDirectionY[iPfp]=999999;
     fPFPTrackVertexDirectionZ[iPfp]=999999;
     fPFPTrackEndDirectionX[iPfp]=999999;
     fPFPTrackEndDirectionY[iPfp]=999999;
     fPFPTrackEndDirectionZ[iPfp]=999999;
     fPFPTrackChi2[iPfp]=999999;
     fPFPTrackNdof[iPfp]=999999;
     
     fPFPShowerID[iPfp]=999999;
     fPFPShowerBestPlane[iPfp]=999999;
     fPFPShowerDirectionX[iPfp]=999999;
     fPFPShowerDirectionY[iPfp]=999999;
     fPFPShowerDirectionZ[iPfp]=999999;
     fPFPShowerDirectionErrX[iPfp]=999999;
     fPFPShowerDirectionErrY[iPfp]=999999;
     fPFPShowerDirectionErrZ[iPfp]=999999;
     fPFPShowerStartX[iPfp]=999999;
     fPFPShowerStartY[iPfp]=999999;
     fPFPShowerStartZ[iPfp]=999999;
     fPFPShowerStartErrX[iPfp]=999999;
     fPFPShowerStartErrY[iPfp]=999999;
     fPFPShowerStartErrZ[iPfp]=999999;
     fPFPShowerLength[iPfp]=999999;
     fPFPShowerOpenAngle[iPfp]=999999;

     for(int iClu=0; iClu<kNMaxPFPClusters; iClu++){
       fPFPCluPlane[iPfp][iClu]=999999; 
       fPFPCluView[iPfp][iClu]=999999; 
       fPFPCluNHits[iPfp][iClu]=999999; 
       fPFPCluIntegral[iPfp][iClu]=999999; 
     }

     fPFPCompleteness[iPfp]=999999;
     fPFPPurity[iPfp]=999999;

     for (unsigned int iView = 0; iView < kNViews; iView++)
     {
          fPFPNHitsView[iPfp][iView] = 999999;
          fPFPNSharedTrueParticleHitsView[iPfp][iView]=0;
          fPFPTrueParticleMatchedIDView[iPfp][iView] = 999999;
          fPFPTrueParticleMatchedPositionView[iPfp][iView] = 999999;
          fPFPPurityView[iPfp][iView]=999999;
          fPFPCompletenessView[iPfp][iView]=999999;
     }
    }
    fNPFParticles = 0;

    return;
}

DEFINE_ART_MODULE(test::pandoraAnalysis)
