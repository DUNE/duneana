////////////////////////////////////////////////////////////////////////
// Class:       MakeSNeTProfileHistos
// Module Type: analyzer
// File:        MakeSNeTProfileHistos_module.cc
//
// Generated at Mon Jan 09 2017 by Michael Baird using the old
// copy and paste method...
//
// Purpose:   The purpose of this module is to generate the input histos
//            for my time-profile, stand alone macro.
//
////////////////////////////////////////////////////////////////////////

// C++ includes

// ROOT includes
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

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

class MakeSNeTProfileHistos;

class MakeSNeTProfileHistos : public art::EDAnalyzer {

public:

  explicit MakeSNeTProfileHistos(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  MakeSNeTProfileHistos(MakeSNeTProfileHistos const &) = delete;
  MakeSNeTProfileHistos(MakeSNeTProfileHistos &&) = delete;
  MakeSNeTProfileHistos & operator = (MakeSNeTProfileHistos const &) = delete;
  MakeSNeTProfileHistos & operator = (MakeSNeTProfileHistos &&) = delete;

  // The main guts...
  void analyze(art::Event const & e) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;

  void endJob() override;



private:

  // labels for the processes that made the data products
  std::string fTruthLabel;   ///< label for process that made the MCTruth info
  std::string fHitLabel;     ///< label for process that made the reco hits



  // required histograms for the Time-Profile simulator:
  TH1F *fEnergySpec;          ///< 1D event energy spectrum (used only for backgrounds)
  TH1F *fTrigEff_numerator;   ///< how many times a trigger was issued per true event energy
  TH1F *fTrigEff;             ///< trigger efficiency per true event energy (fTrigEff_numerator/fEnergySpec)
  TH2F *fNHitsVsEnergy;       ///< number of reconstructed hits per true event energy
  TH2F *fSumADCVsEnergy;      ///< summed ADC for all reco hits per true event energy
  TH1F *fTPSplitProb;         ///< probability of splitting a trigger primitive into two (or more) TPs per true event energy



  // a few extra check histos for double checking general performance:
  TH1F *fNMCTruths;   ///< number of MCTruths per art::event
  TH1F *fCCNC;        ///< CCNC for each MCTruth
  TH1F *fOrigin;      ///< mctruth->Origin()
  TH1F *fNuPDG;       ///< mc_neutrino->Nu()->PdgCode()
  
};



//......................................................
MakeSNeTProfileHistos::MakeSNeTProfileHistos(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}



//......................................................
void MakeSNeTProfileHistos::reconfigure(fhicl::ParameterSet const & p)
{
  fTruthLabel = p.get<std::string>("TruthLabel");
  fHitLabel   = p.get<std::string>("HitLabel");
}



//......................................................
void MakeSNeTProfileHistos::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;

  // binning for the energy axes:
  int NbinsE = 200;
  float minE = 0.0;   // in MeV
  float maxE = 100.0; // in MeV

  // binning for the NHits axes:
  int NbinsNH = 5000;
  float minNH = 0.0;
  float maxNH = 5000.0; // a reasonable number determined from running over single nue events

  // binning for the ADC axes:
  int NbinsADC = 500;
  float minADC = 0.0;
  float maxADC = 500.0e3; // a reasonable number determined from running over single nue events



  //
  // Book required histograms for the Time-Profile simulator:
  //

  // 2D neutrino energy spectrum as a function of time:
  //   (Note: generated separately by Kate Schol.)

  // 1D event energy spectrum: (used only for backgrounds)
  fEnergySpec = tfs->make<TH1F>("EnergySpec",
				"True energy spectra;true energy [MeV];",
				NbinsE, minE, maxE);

  // Overall efficency for triggering as a function of true event energy:
  fTrigEff_numerator = tfs->make<TH1F>("TrigEff_numerator",
				       "Total number of triggers issued;true energy [MeV];",
				       NbinsE, minE, maxE);
  fTrigEff           = tfs->make<TH1F>("TrigEff",
				       "Trigger efficiency;true energy [MeV];",
				       NbinsE, minE, maxE);

  // Number of reco hits vs. true event energy:
  fNHitsVsEnergy = tfs->make<TH2F>("NHitsVsEnergy",
				   "Number of reco hits per true energy;true energy [MeV];NHits",
				   NbinsE, minE, maxE,
				   NbinsNH, minNH, maxNH);

  // Summed ADC vs. true event energy:
  fSumADCVsEnergy = tfs->make<TH2F>("SumADCVsEnergy",
				    "Summed ADC over all reco hits per true energy;true energy [MeV];Summed ADC",
				    NbinsE, minE, maxE,
				    NbinsADC, minADC, maxADC);

  // Probability of splitting a TP into 2 TPs as a function of true event energy:
  fTPSplitProb = tfs->make<TH1F>("TPSplitProb",
				 "Probability of spliting a TP into 2 or more per true energy;true energy [MeV];",
				 NbinsE, minE, maxE);



  //
  // Book a few extra check histos:
  //

  // Number of MCTruths per art::event:
  fNMCTruths = tfs->make<TH1F>("NMCTruths",
			       "Number of MCTruths per art::event;# of MCTruths;",
			       101,-0.5,100.5);

  // CCNC:
  fCCNC = tfs->make<TH1F>("CCNC",
			  "CCNC;CCNC;",
			  11,-5.5,5.5);

  // Origin:
  fOrigin = tfs->make<TH1F>("Origin",
			    "Neutrino origin;origin type;",
			    16,-5.5,10.5);

  // Neutrino pdg code:
  fNuPDG = tfs->make<TH1F>("NuPDG",
			   "Neutrino PDG code;pdg code;",
			   41,-20.5,20.5);

}



//......................................................
void MakeSNeTProfileHistos::endJob()
{

}



//......................................................
void MakeSNeTProfileHistos::analyze(art::Event const & e)
{

  // TODO:
  // 
  // 1. implement dividing the TrigEff histos in the endJob function
  // 2. implement anything that has anything to do with the trigger data products
  //       (including the TrigEff histos, TP splitting histos, etc...???)
  // 3. figure out a good way to generalize this for multiple MCTruths...




  // ===== Notes and Assumptions ===== //
  //
  // This module assumes the following:
  //
  // 1. There is only one neutrino interaction per art::event. This is
  //    important because it assumes that all reco hits go with that one
  //    neutrino.





  //
  // Lift out the desired data products:
  //

  // Get the hits out of the event
  auto hits_list = e.getHandle< std::vector< recob::Hit > >(fHitLabel);

  // Get the MCTruths out of the event
  auto truths_list = e.getHandle< std::vector< simb::MCTruth > >(fTruthLabel);
  
  fNMCTruths->Fill(truths_list->size());

  // Assert assumption that there is one and only one MCTruth in the event
  // (note:   this needs to be fixed later so that it can handle multiple MCTruths.
  //          This only works for single generated SNe neutrinos for now...)
  if(truths_list->size() == 1) {

    // calculate Summed ADC
    float SummedADC = 0.0; // ...I know... this number is a float?! Who would have guessed that...
    for(unsigned int i = 0; i < hits_list->size();i++){
      recob::Hit const& hit = hits_list->at(i);  
      SummedADC += hit.SummedADC();
    }

    
    simb::MCTruth const& truth = truths_list->at(0);

    fOrigin->Fill(truth.Origin());

    if (truth.NeutrinoSet()){
      simb::MCNeutrino const& mc_nu = truth.GetNeutrino();
      double NuE = 1000.0*mc_nu.Nu().E(); // convert to MeV

      fEnergySpec->Fill(NuE);
      fNHitsVsEnergy->Fill(NuE,hits_list->size());
      fSumADCVsEnergy->Fill(NuE,SummedADC);

      fCCNC->Fill(mc_nu.CCNC());
      fNuPDG->Fill(mc_nu.Nu().PdgCode());

    }    

  }
  else {
    std::cerr << "\n\n\nWARNING: Number of MCTruths != 1\n\n\n\n";
  }



}

DEFINE_ART_MODULE(MakeSNeTProfileHistos)
