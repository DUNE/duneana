////////////////////////////////////////////////////////////////////////
// Class:       CosmicEfficiency
// Plugin Type: analyzer (art v2_07_03)
// File:        CosmicEfficiency_module.cc
//
// Generated at Mon Sep  4 06:55:33 2017 by Leigh Whitehead using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TH1F.h"

namespace cosmic {
  class CosmicEfficiency;
}


class cosmic::CosmicEfficiency : public art::EDAnalyzer {
public:

  // Just make an enum for plotting in a histogram to avoid
  // histograms of the pdg code. These definitions are
  // valid for particle and antiparticle
  // Note the enum starts at one since this is the bin
  // number for the histograms
  enum particleType
  {
    kNotSet = 0,
    kElectron = 1,
    kMuon = 2,
    kPion = 3,
    kProton = 4,
    kKaon = 5,
    kOther = 6
  };

  explicit CosmicEfficiency(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicEfficiency(CosmicEfficiency const &) = delete;
  CosmicEfficiency(CosmicEfficiency &&) = delete;
  CosmicEfficiency & operator = (CosmicEfficiency const &) = delete;
  CosmicEfficiency & operator = (CosmicEfficiency &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  const simb::MCParticle* GetTruthParticle(detinfo::DetectorClocksData const& clockData,
                                           const std::vector< art::Ptr<recob::Hit> > & hits) const;

  const particleType ConvertPDGCode(int pdg) const;

  int ConvertCosmicType(anab::CosmicTagID_t tag) const;

  // fcl parameters
  art::InputTag fTrackerTag;
  bool fVerbose;

  // Useful counters
  unsigned int fTrueCosmicMuon;
  unsigned int fTaggedCosmicMuon;
  unsigned int fTrueOther;
  unsigned int fTaggedOther;
  // Also tag just anything from cosmic or beam origin
  unsigned int fTrueCosmicOrigin;
  unsigned int fTaggedCosmicOrigin;
  unsigned int fTrueOtherOrigin;
  unsigned int fTaggedOtherOrigin;

  // Histograms
  TH1F *hCosmicMuonEff;
  TH1F *hCosmicMuonPur;
  TH1F *hCosmicEff;
  TH1F *hCosmicPur;
  TH1F *hTaggedPDG;
  TH1F *hTrkLenAll;
  TH1F *hTrkLenTag;
  TH1F *hTrkLenEff;
  TH1F  *hTagType;
};


cosmic::CosmicEfficiency::CosmicEfficiency(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fTrackerTag(p.get<art::InputTag>("TrackerTag")),
  fVerbose(p.get<bool>("Verbose")),
  fTrueCosmicMuon(0), fTaggedCosmicMuon(0),
  fTrueOther(0), fTaggedOther(0),
  fTrueCosmicOrigin(0), fTaggedCosmicOrigin(0),
  fTrueOtherOrigin(0), fTaggedOtherOrigin(0)
{

}

void cosmic::CosmicEfficiency::beginJob()
{

  // Get a file service so we can make some histograms to write out
  art::ServiceHandle<art::TFileService> tfs;

  hCosmicMuonEff = tfs->make<TH1F>("hCosmicMuonEff",";Efficiency",50,0,1.0001);
  hCosmicMuonPur = tfs->make<TH1F>("hCosmicMuonPur",";Purity",50,0,1.0001);
  hCosmicEff = tfs->make<TH1F>("hCosmicEff",";Efficiency",50,0,1.0001);
  hCosmicPur = tfs->make<TH1F>("hCosmicPur",";Purity",50,0,1.0001);
  hTaggedPDG = tfs->make<TH1F>("hTaggedParticleType","",6,0.5,6.5);
  hTrkLenAll = tfs->make<TH1F>("hTrkLenAll","",10,0,600);
  hTrkLenTag = tfs->make<TH1F>("hTrkLenTag","",10,0,600);
  hTrkLenEff = tfs->make<TH1F>("hTrkLenEff","",10,0,600);
  hTagType   = tfs->make<TH1F>("hTagType","",12,0,12);

  // Let's give the particle type plot useful labels
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kElectron,"Electron");
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kMuon,"Muon");
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kPion,"Pion");
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kProton,"Proton");
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kKaon,"Kaon");
  hTaggedPDG->GetXaxis()->SetBinLabel(cosmic::CosmicEfficiency::kOther,"Other");

  // Tag type labels would be nice too.
  hTagType->GetXaxis()->SetBinLabel(1,"YY");
  hTagType->GetXaxis()->SetBinLabel(2,"YZ");
  hTagType->GetXaxis()->SetBinLabel(3,"ZZ");
  hTagType->GetXaxis()->SetBinLabel(4,"XX");
  hTagType->GetXaxis()->SetBinLabel(5,"XY");
  hTagType->GetXaxis()->SetBinLabel(6,"XZ");
  hTagType->GetXaxis()->SetBinLabel(7,"Y");
  hTagType->GetXaxis()->SetBinLabel(8,"Z");
  hTagType->GetXaxis()->SetBinLabel(9,"X");
  hTagType->GetXaxis()->SetBinLabel(10,"Outside Drift (Part)");
  hTagType->GetXaxis()->SetBinLabel(11,"Outside Drift (Full)");
  hTagType->GetXaxis()->SetBinLabel(12,"Beam Incompatible");

}

void cosmic::CosmicEfficiency::analyze(art::Event const & evt)
{

  // We must have MC for this module to make sense
  if(evt.isRealData()) return;

  // Get the reconstructed tracks
  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackerTag);

  // We need the association between the tracks and the hits
  const art::FindManyP<recob::Hit> findTrackHits(recoTracks, evt, fTrackerTag);

  // Also want the association between tracks and anab::CosmicTag objects
  const art::FindManyP<anab::CosmicTag> findCosmicTags(recoTracks,evt,fTrackerTag);

  // Finally, get the backtracker
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  // Since we will make plots of efficiency per event, we need local versions of the class variables
  unsigned int local_TrueCosmicMuon = 0; 
  unsigned int local_TaggedCosmicMuon = 0;
  unsigned int local_TrueOther = 0;
  unsigned int local_TaggedOther = 0;
  unsigned int local_TrueCosmicOrigin = 0;
  unsigned int local_TaggedCosmicOrigin = 0;
  unsigned int local_TrueOtherOrigin = 0;
  unsigned int local_TaggedOtherOrigin = 0;

  // Fill a vector with the true tracks matched
  std::vector<int> matchedTruths;

  // Loop over the tracks
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  for(unsigned int t = 0; t < recoTracks->size(); ++t){
  
    // Get the hits
    auto const& trackHits = findTrackHits.at(t);

    // Get the best match MC Particle for this track
    const simb::MCParticle* trueMatch = GetTruthParticle(clockData, trackHits);

    if(trueMatch == 0x0) continue;

    bool alreadyMatched = false;
    if(std::find(matchedTruths.begin(),matchedTruths.end(),trueMatch->TrackId())==matchedTruths.end()){
      matchedTruths.push_back(trueMatch->TrackId());
    }
    else{
      alreadyMatched = true;
    }

    // We need to use the backtracker to find if this true particle is cosmic or not
    bool isCosmic = false;
    bool isCosmicMuon = false;
    const art::Ptr<simb::MCTruth> mc = pi_serv->TrackIdToMCTruth_P(trueMatch->TrackId());

    // Check if we have a cosmic muon (or anti-muon) and set the true counting variables
    if(mc->Origin() == simb::kCosmicRay){
      if(fabs(trueMatch->PdgCode()) == 13){
        ++local_TrueCosmicMuon;
        isCosmicMuon = true;
        hTrkLenAll->Fill(((*recoTracks)[t].Vertex()-(*recoTracks)[t].End()).R());
      }
      else{
        ++local_TrueOther;
      }
      ++local_TrueCosmicOrigin;
      isCosmic = true;
    }
    else{
      // Definitely not cosmic!
      ++local_TrueOther;
      ++local_TrueOtherOrigin;
    }

    // Now look for a cosmic tag object
    auto const& trackCosmicTag = findCosmicTags.at(t);
   
    // I'm not sure GEANT keeps the names of the particle processes properly. Let's
    // get the MCParticle straight from the MCTruth and see what happens.
    std::string correctProcess = "";
    for(int tp = 0; tp < mc->NParticles(); ++tp){

      if(mc->GetParticle(tp).Position().X() == trueMatch->Position().X()){
        if(mc->GetParticle(tp).Position().Y() == trueMatch->Position().Y()){
          if(mc->GetParticle(tp).Position().Z() == trueMatch->Position().Z()){
            correctProcess = mc->GetParticle(tp).Process();
            break;
          }
        }
      }
    }
 
    // Did we find a cosmic tag for this track?
    if(trackCosmicTag.size() == 0){
      // If we don't have a tag, let's see if we can find out why.
      if(isCosmicMuon && fVerbose){
        std::cout << "We've found a cosmic muon but we don't have a tag for it" << std::endl;
        //std::cout << "Reco vertex = "; (*recoTracks)[t].Vertex().Print();
        //std::cout << "True vertex = "; trueMatch->Position().Vect().Print();
        //std::cout << "Reco end = "; (*recoTracks)[t].End().Print();
        //std::cout << "True end = "; trueMatch->EndPosition().Vect().Print();
        std::cout << " - Track length = " << ((*recoTracks)[t].Vertex()-(*recoTracks)[t].End()).R() << std::endl; 
        if((*recoTracks)[t].Vertex().Y() < (*recoTracks)[t].End().Y()){
          std::cout << " - Track reconstructed with the vertex lower than the end point" << std::endl;
        }
        if((*recoTracks)[t].Vertex().Y() < 550 && (*recoTracks)[t].End().Y() < 550){
          std::cout << " - This track was reconstructed away from the top of the detector... no way to tag it. " << std::endl;
        }
        if(alreadyMatched){
          std::cout << " - This true particle " << trueMatch->TrackId() << "was already matched" << std::endl;
        }
      }
      continue;
    }
    else{
      // Specific cosmic muon counters
      if(isCosmicMuon){
        ++local_TaggedCosmicMuon;
        hTrkLenTag->Fill(((*recoTracks)[t].Vertex()-(*recoTracks)[t].End()).R());
      }
      else{
        ++local_TaggedOther;
      }
      // Origin counters
      if(isCosmic){
        ++local_TaggedCosmicOrigin;
      }
      else{
        ++local_TaggedOtherOrigin;
        if(fVerbose){
          if(correctProcess == "primary"){
            std::cout << " - Oops! Tagged a good beam particle (" << trueMatch->PdgCode() << ") as a cosmic (" << ConvertCosmicType(trackCosmicTag[0]->CosmicType()) << ")" << std::endl;
          }
          else if(correctProcess == "primaryBackground"){
            std::cout << " - Oops! Tagged a background beam particle (" << trueMatch->PdgCode() << ") as a cosmic (" << ConvertCosmicType(trackCosmicTag[0]->CosmicType()) << ")" << std::endl;
          }
          else{
            std::cout << " - Oops! Tagged a non-primary beam particle (" << trueMatch->PdgCode() << ") as a cosmic (" << ConvertCosmicType(trackCosmicTag[0]->CosmicType()) << ")" << std::endl;
          }
          (*recoTracks)[t].Vertex<TVector3>().Print();
          (*recoTracks)[t].End<TVector3>().Print();  
          trueMatch->Position().Vect().Print();
          trueMatch->EndPosition().Vect().Print();
        } 
      }
      // Fill the tagged particle hisrogram
      hTaggedPDG->Fill(ConvertPDGCode(trueMatch->PdgCode()));
      // Type of tag histogram
      hTagType->Fill(ConvertCosmicType(trackCosmicTag[0]->CosmicType()));
    }

  } // End loop over reconstructed tracks

  // Fill the plots
  if(local_TrueCosmicMuon > 0){
    hCosmicMuonEff->Fill(local_TaggedCosmicMuon/static_cast<float>(local_TrueCosmicMuon));
  }
  if(local_TaggedCosmicMuon+local_TaggedOther > 0){
    hCosmicMuonPur->Fill(local_TaggedCosmicMuon/static_cast<float>(local_TaggedCosmicMuon + local_TaggedOther));
  }
  if(local_TrueCosmicOrigin > 0){
    hCosmicEff->Fill(local_TaggedCosmicOrigin/static_cast<float>(local_TrueCosmicOrigin));
  }
  if(local_TaggedCosmicOrigin+local_TaggedOtherOrigin > 0){
    hCosmicPur->Fill(local_TaggedCosmicOrigin/static_cast<float>(local_TaggedCosmicOrigin + local_TaggedOtherOrigin));
  } 
  
  // Make sure to increment the main member counters
  fTrueCosmicMuon += local_TrueCosmicMuon;
  fTaggedCosmicMuon += local_TaggedCosmicMuon;
  fTrueOther += local_TrueOther;
  fTaggedOther += local_TaggedOther;

  fTrueCosmicOrigin += local_TrueCosmicOrigin;
  fTaggedCosmicOrigin += local_TaggedCosmicOrigin;
  fTrueOtherOrigin += local_TrueOtherOrigin;
  fTaggedOtherOrigin += local_TaggedOtherOrigin;
}

void cosmic::CosmicEfficiency::endJob()
{
  if(fTrueCosmicMuon != 0 && (fTaggedCosmicMuon+fTaggedOther) !=0){
    std::cout << " == Cosmic Efficiency Summary == " << std::endl;
    std::cout << " - Looking for cosmic muons: " << std::endl;
    std::cout << " Tagged " << fTaggedCosmicMuon << " of " << fTrueCosmicMuon << " :: " << 100.*fTaggedCosmicMuon/static_cast<float>(fTrueCosmicMuon) << "% efficiency" << std::endl;
    std::cout << " Total tags = " << fTaggedCosmicMuon + fTaggedOther << " :: " 
              << 100.*fTaggedCosmicMuon/static_cast<float>(fTaggedCosmicMuon + fTaggedOther) << "% purity" << std::endl;
    std::cout << " - Looking for just cosmic origin: " << std::endl;
    std::cout << " Tagged " << fTaggedCosmicOrigin << " of " << fTrueCosmicOrigin << " :: " << 100.*fTaggedCosmicOrigin/static_cast<float>(fTrueCosmicOrigin) << "% efficiency" << std::endl;
    std::cout << " Total tags = " << fTaggedCosmicOrigin + fTaggedOtherOrigin << " :: " 
              << 100.*fTaggedCosmicOrigin/static_cast<float>(fTaggedCosmicOrigin + fTaggedOtherOrigin) << "% purity" << std::endl;
  }
  for(int b = 1; b <= hTrkLenAll->GetNbinsX(); ++b){
    if(hTrkLenAll->GetBinContent(b) != 0){
      float all = hTrkLenAll->GetBinContent(b);
      float tag = hTrkLenTag->GetBinContent(b);
      float eff = tag / all;
      // Use binomial error... not perfect, but hopefully ok.
      float err = sqrt((eff*(1-eff))/all);
      hTrkLenEff->SetBinContent(b,eff);
      hTrkLenEff->SetBinError(b,err);
    }
  }
}

// Function to find the best matched true particle to a reconstructed particle
const simb::MCParticle* cosmic::CosmicEfficiency::GetTruthParticle(detinfo::DetectorClocksData const& clockData,
                                                                   const std::vector< art::Ptr<recob::Hit> > & hits) const
{
  const simb::MCParticle* mcParticle = 0;

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::unordered_map<int, double> trkIDE;
  for (auto const & h : hits)
  {
    for (auto const & ide : bt_serv->HitToTrackIDEs(clockData, h)) // loop over std::vector<sim::TrackIDE>
    {
        trkIDE[ide.trackID] += ide.energy; // sum energy contribution by each track ID
    }
  }

  int best_id = 0;
  double tot_e = 0, max_e = 0;
  for (auto const & contrib : trkIDE)
  {
    tot_e += contrib.second;     // sum total energy in these hits
    if (contrib.second > max_e)  // find track ID corresponding to max energy
    {
        max_e = contrib.second;
        best_id = contrib.first;
    }
  }

  if ((max_e > 0) && (tot_e > 0)) // ok, found something reasonable
  {
    if (best_id < 0)            // NOTE: negative ID means this is EM activity
    {                           // caused by track with the same but positive ID
//        best_id = -best_id;     // --> we'll find mother MCParticle of these hits
      return mcParticle;
    }
    mcParticle = pi_serv->TrackIdToParticle_P(best_id); // MCParticle corresponding to track ID
  }

  return mcParticle;
}

// Function to convert the PDG code into the scheme we want to use for our histogram
const cosmic::CosmicEfficiency::particleType cosmic::CosmicEfficiency::ConvertPDGCode(int pdg) const
{
  cosmic::CosmicEfficiency::particleType particle = cosmic::CosmicEfficiency::kNotSet;

  switch(pdg)
  {
    case -11:
    case  11: 
      particle = cosmic::CosmicEfficiency::kElectron;
      break;
    case -13:
    case  13: 
      particle = cosmic::CosmicEfficiency::kMuon;
      break;
    case  211:
    case -211:
      particle = cosmic::CosmicEfficiency::kPion;
      break;
    case  2212:
    case -2212:
      particle = cosmic::CosmicEfficiency::kProton;
      break;
    case  321:
    case -321:
      particle = cosmic::CosmicEfficiency::kKaon;
      break;
    default:
      particle = cosmic::CosmicEfficiency::kOther;
  }

  return particle;
}

// Function to convert the cosmic tag enum into a value for the histogram 
int cosmic::CosmicEfficiency::ConvertCosmicType(anab::CosmicTagID_t tag) const
{
  int val = -1;

  switch(tag){
    case anab::CosmicTagID_t::kGeometry_YY : val = 0; break;
    case anab::CosmicTagID_t::kGeometry_YZ : val = 1; break;
    case anab::CosmicTagID_t::kGeometry_ZZ : val = 2; break;
    case anab::CosmicTagID_t::kGeometry_XX : val = 3; break;
    case anab::CosmicTagID_t::kGeometry_XY : val = 4; break;
    case anab::CosmicTagID_t::kGeometry_XZ : val = 5; break;
    case anab::CosmicTagID_t::kGeometry_Y  : val = 6; break;
    case anab::CosmicTagID_t::kGeometry_Z  : val = 7; break;
    case anab::CosmicTagID_t::kGeometry_X  : val = 8; break;
    case anab::CosmicTagID_t::kOutsideDrift_Partial   : val = 9; break;
    case anab::CosmicTagID_t::kOutsideDrift_Complete  : val = 10; break;
    case anab::CosmicTagID_t::kFlash_BeamIncompatible : val = 11; break;
    default: val = -1;
  }
  return val;
}

DEFINE_ART_MODULE(cosmic::CosmicEfficiency)
