////////////////////////////////////////////////////////////////////////
// Class:       DAQSimAna
// Module Type: analyzer
// File:        DAQSimAna_module.cc
//
// Generated by Michael Baird using the old copy and paste...
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

// C++ includes

// ROOT includes
#include "TH1F.h"

// Framework includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/RawDigit.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes



class DAQSimAna;

class DAQSimAna : public art::EDAnalyzer {

public:

  explicit DAQSimAna(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  DAQSimAna(DAQSimAna const &) = delete;
  DAQSimAna(DAQSimAna &&) = delete;
  DAQSimAna & operator = (DAQSimAna const &) = delete;
  DAQSimAna & operator = (DAQSimAna &&) = delete;

  // The main guts...
  void analyze(art::Event const & e) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob();



private:

  // label for module that made raw digits
  std::string fRawDigitLabel;

  // a simple histo to be filled
  
};



//......................................................
DAQSimAna::DAQSimAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}



//......................................................
void DAQSimAna::reconfigure(fhicl::ParameterSet const & p)
{
  fRawDigitLabel = p.get<std::string> ("RawDigitLabel");
}



//......................................................
void DAQSimAna::beginJob()
{
  /*
  art::ServiceHandle<art::TFileService> tfs;

  fTrigTypes = tfs->make<TH1F>("fTrigTypes","Trigger Types;trig type;count",
			       101,-0.5,100.5);
  */
}



//......................................................
void DAQSimAna::analyze(art::Event const & e)
{

  //
  // Lift out the TPC raw digits:
  //
  art::Handle<std::vector<raw::RawDigit>> digitsHandle;
  e.getByLabel(fRawDigitLabel, digitsHandle);  

  std::cout << "\n\n\nraw_digits.size() = " << digitsHandle->size() << "\n\n\n";

}

DEFINE_ART_MODULE(DAQSimAna)
