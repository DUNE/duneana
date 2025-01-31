 ////////////////////////////////////////////////////////////////////////
// Class:       SupernovaAna
// Module Type: analyzer
// File:        SupernovaAna_module.cc
//
// Generated at Mon Jul 11 21:36:48 2016 by Michael Baird using the old
// copy and paste...
////////////////////////////////////////////////////////////////////////

// C++ includes

// ROOT includes
#include "TTree.h"
#include "TH1F.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes
#include "dunecore/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
#include "dunecore/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"

// should this be somewhere else?
const int kMaxNumHits=50000;
const int kMaxNumParticles=500;

class SupernovaAna;

class SupernovaAna : public art::EDAnalyzer {

public:

  explicit SupernovaAna(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  SupernovaAna(SupernovaAna const &) = delete;
  SupernovaAna(SupernovaAna &&) = delete;
  SupernovaAna & operator = (SupernovaAna const &) = delete;
  SupernovaAna & operator = (SupernovaAna &&) = delete;

  // The main guts...
  void analyze(art::Event const & e) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;



private:

  // label for the modules that made trigger objects
  std::string fTruthLabel;
  std::string fHitLabel;

  // Define an ntuple to be filled
  TTree* fNTuple;

  // Define variables to be used as branches of the nTuple
  int EventNum;
  int EventNumCut;
  int NHits;
  double avgRMS;
  double avgRMSloop[kMaxNumHits];
  double avgRMSloopsum;
  double SummedADC;
  int NMCTruths;
  int NParticles[kMaxNumParticles];
  double MomentumYTruths[kMaxNumParticles];
  double MomentumZTruths[kMaxNumParticles];
  double NuEnergyTruths[kMaxNumParticles];
  double LeptonEnergyTruths[kMaxNumParticles];
  double ThetaTruths[kMaxNumParticles];
  int OriginTruths[kMaxNumParticles];
  int CCNCTruths[kMaxNumParticles];
  int ModeTruths[kMaxNumParticles];
  int NuPDGCodeTruths[kMaxNumParticles];  
  int LeptonPDGCodeTruths[kMaxNumParticles];  
  double xTruths[kMaxNumParticles];
  double yTruths[kMaxNumParticles];
  double zTruths[kMaxNumParticles];
  double tTruths[kMaxNumParticles];
  double EndxTruths[kMaxNumParticles];
  double EndyTruths[kMaxNumParticles];
  double EndzTruths[kMaxNumParticles];
  double EndtTruths[kMaxNumParticles];
  double DeltaXTruths[kMaxNumParticles];
  double DeltaYTruths[kMaxNumParticles];
  
};



//......................................................
SupernovaAna::SupernovaAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}



//......................................................
void SupernovaAna::reconfigure(fhicl::ParameterSet const & p)
{
  fTruthLabel = p.get<std::string> ("TruthLabel");
  fHitLabel = p.get<std::string> ("HitLabel");
}



//......................................................
void SupernovaAna::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;

  // Make nTuple
  fNTuple = tfs->make<TTree>("SupernovaAna","Supernonva analysis tree");

  // Add branches to nTuple
  fNTuple->Branch("EventNum",&EventNum);
  fNTuple->Branch("EventNumCut",&EventNumCut);
  fNTuple->Branch("NHits",&NHits);
  fNTuple->Branch("avgRMS",&avgRMS);
  fNTuple->Branch("SummedADC",&SummedADC);
  fNTuple->Branch("NMCTruths",&NMCTruths);
  fNTuple->Branch("NParticles",NParticles,"NParticles[NMCTruths]/I");
  fNTuple->Branch("MomentumYTruths",MomentumYTruths,"MomentumYTruths[NMCTruths]/D");
  fNTuple->Branch("MomentumZTruths",MomentumZTruths,"MomentumZTruths[NMCTruths]/D");
  fNTuple->Branch("NuEnergyTruths",NuEnergyTruths,"NuEnergyTruths[NMCTruths]/D");
  fNTuple->Branch("LeptonEnergyTruths",LeptonEnergyTruths,"LeptonEnergyTruths[NMCTruths]/D");
  fNTuple->Branch("ThetaTruths",ThetaTruths,"ThetaTruths[NMCTruths]/D");
  fNTuple->Branch("OriginTruths",OriginTruths,"OriginTruths[NMCTruths]/I");
  fNTuple->Branch("CCNCTruths",CCNCTruths,"CCNCTruths[NMCTruths]/I");
  fNTuple->Branch("ModeTruths",ModeTruths,"ModeTruths[NMCTruths]/I");
  fNTuple->Branch("NuPDGCodeTruths",NuPDGCodeTruths,"NuPDGCodeTruths[NMCTruths]/I");
  fNTuple->Branch("LeptonPDGCodeTruths",LeptonPDGCodeTruths,"LeptonPDGCodeTruths[NMCTruths]/I");
  fNTuple->Branch("xTruths",xTruths,"xTruths[NMCTruths]/D");
  fNTuple->Branch("yTruths",yTruths,"yTruths[NMCTruths]/D");
  fNTuple->Branch("zTruths",zTruths,"zTruths[NMCTruths]/D");
  fNTuple->Branch("tTruths",tTruths,"tTruths[NMCTruths]/D");
  fNTuple->Branch("EndxTruths",EndxTruths,"EndxTruths[NMCTruths]/D");
  fNTuple->Branch("EndyTruths",EndyTruths,"EndyTruths[NMCTruths]/D");
  fNTuple->Branch("EndzTruths",EndzTruths,"EndzTruths[NMCTruths]/D");
  fNTuple->Branch("EndtTruths",EndtTruths,"EndtTruths[NMCTruths]/D");
  fNTuple->Branch("DeltaXTruths",DeltaXTruths,"DeltaXTruths[NMCTruths]/D");
  fNTuple->Branch("DeltaYTruths",DeltaYTruths,"DeltaYTruths[NMCTruths]/D");

}



//......................................................
void SupernovaAna::analyze(art::Event const & e)
{

  
  //==============================HITS================================

  // Get the Hits out of the event
  auto hits_list = e.getHandle< std::vector< recob::Hit > >(fHitLabel);

  // Gets the number of hits per event and fills into a histogram
  NHits = hits_list->size();
 
  // Get summed ADC and average RMS of hit slope
  SummedADC = 0.0;
  for(unsigned int i = 0; i < hits_list->size();i++){
    recob::Hit const& hit = hits_list->at(i);  
    avgRMSloop[i] = hit.RMS();
    avgRMSloopsum += avgRMSloop[i];
    SummedADC += hit.ROISummedADC();
   }
  avgRMS = avgRMSloopsum/(hits_list->size());

  //============================MCTRUTHS==============================

  // Get the MCTruths out of the event
  auto truths_list = e.getHandle< std::vector< simb::MCTruth > >(fTruthLabel);

  // Gets the number of MC truths per event and fills into a histogram
  NMCTruths = truths_list->size();
  
  // Loop through the truth arrays and sets the values of the arrays
  // to an arbitrary default (unfilled) value
  for(int i=0;i<kMaxNumParticles;i++){
    xTruths[i] = -1.0e9;
    yTruths[i] = -1.0e9;
    zTruths[i] = -1.0e9;
    tTruths[i] = -1.0e9;
    EndxTruths[i] = -1.0e9;
    EndyTruths[i] = -1.0e9;
    EndzTruths[i] = -1.0e9;
    EndtTruths[i] = -1.0e9;
    NParticles[i] = -999;
    MomentumYTruths[i] = -1.0e9;
    MomentumZTruths[i] = -1.0e9;
    NuEnergyTruths[i] = -999;
    LeptonEnergyTruths[i] = -999;
    ThetaTruths[i] = -999;
    OriginTruths[i] = -999;
    CCNCTruths[i] = -999;
    ModeTruths[i] = -999;
    NuPDGCodeTruths[i] = 0;
    LeptonPDGCodeTruths[i] = 0;
    DeltaXTruths[i] = -1.0e9;
    DeltaYTruths[i] = -1.0e9;

  }
  
  // Detect if the number of MC truths per event is greater than 
  // or equal to the maximum number of particles per event else 
  // we'll fall off the end of the array and cause a seg fault
  if(NMCTruths >= kMaxNumParticles) { 
    std::cerr << "ERROR: NMCTruths " << NMCTruths << 
    " >= kMaxNumParticles " << kMaxNumParticles << 
    " , will cause a segmentation fault" << std::endl;}
  // Detect if the number of hits per event is greater than 
  // or equal to the maximum number of particles per event else 
  // we'll fall off the end of the array and cause a seg fault
  if(NHits >= kMaxNumHits) { 
    std::cerr << "ERROR: NHits " << NHits << 
    " >= kMaxNumHits " << kMaxNumHits << 
    " , will cause a segmentation fault" << std::endl;}
  
  // Create nTuple data for any number of truths 
  for(unsigned int i = 0; i < truths_list->size();i++){
    simb::MCTruth const& truth = truths_list->at(i);

    NParticles[i] = truth.NParticles();

    // Get the ith particle's Origin truths
    OriginTruths[i] = truth.Origin();
           
    // Only get neutrino truths if they exist to avoid a segmentation fault
    if (truth.NeutrinoSet()){
      simb::MCNeutrino const& mc_neutrino = truth.GetNeutrino();
      // should change notation such that mc_neutrino = truth.GetNeutrino().Nu();     // and then remove .Nu() when called
      simb::MCParticle const& mc_lepton = mc_neutrino.Lepton();

      // Get data to be filled into nTuple       
      xTruths[i] = mc_lepton.Vx();
      yTruths[i] = mc_lepton.Vy();
      zTruths[i] = mc_lepton.Vz();
      tTruths[i] = mc_lepton.T();
      EndxTruths[i] = mc_lepton.EndX();
      EndyTruths[i] = mc_lepton.EndY();
      EndzTruths[i] = mc_lepton.EndZ();
      EndtTruths[i] = mc_lepton.EndT();
      MomentumYTruths[i] = mc_lepton.Py();
      MomentumZTruths[i] = mc_lepton.Pz();
      NuEnergyTruths[i] = mc_neutrino.Nu().E();
      LeptonEnergyTruths[i] = mc_lepton.E();
      ThetaTruths[i] = mc_neutrino.Theta();
      CCNCTruths[i] = mc_neutrino.CCNC();
      ModeTruths[i] = mc_neutrino.Mode();
      NuPDGCodeTruths[i] = mc_neutrino.Nu().PdgCode();
      LeptonPDGCodeTruths[i] = mc_lepton.PdgCode();
      DeltaXTruths[i] = EndxTruths[i] - xTruths[i];
      DeltaYTruths[i] = EndyTruths[i] - yTruths[i];
    }
  }
  
    art::EventNumber_t event = e.id().event();
    EventNum = event;
    if (NHits<2) {
      EventNumCut = event;
    }


  // Now we should have NMCTruths = truths_list->size() and 
  // PDGCodeTruths [0 - (NMCTruths - 1)] should contain PdgCodes, 
  // [ NMCTruths - 500 ] should be equal to -999, signifying that 
  // we did not fill this part
  
  // Fill nTuple
  fNTuple->Fill();
}
DEFINE_ART_MODULE(SupernovaAna)
