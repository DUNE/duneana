// dunezsanalysis_module.cc
// A basic "skeleton" to read in art::Event records from a file,
// access their information, and do something with them. 

// See
// <https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Using_the_Framework>
// for a description of the ART classes used here.

// Almost everything you see in the code below may have to be changed
// by you to suit your task. The example task is to make histograms
// and n-tuples related to dE/dx of particle tracks in the detector.

// As you try to understand why things are done a certain way in this
// example ("What's all this stuff about 'auto const&'?"), it will help
// to read ADDITIONAL_NOTES.txt in the same directory as this file.

#ifndef dunezsanalysis_Module
#define dunezsanalysis_Module

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/raw.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <http://root.cern.ch/root/html532/ClassIndex.html>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

namespace dunezsanalysis {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class dunezsanalysis : public art::EDAnalyzer 
  {
  public:
 
    // Standard constructor and destructor for an ART module.
    explicit dunezsanalysis(fhicl::ParameterSet const& pset);
    virtual ~dunezsanalysis();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event. 
    void analyze (const art::Event& evt); 

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fRawProducerLabel; // The name of the producer that tracked simulated particles through the detector

    // The n-tuples we'll create.
    //TTree* fSimulationNtuple; // unused
    TTree* fReconstructionNtuple;

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;

    // array of sums of charges for different choices of zero-suppression thresholds.

    double fChargeSum[200];

    // Other variables that will be shared between different methods.
    art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

  }; // class dunezsanalysis


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // Constructor
  dunezsanalysis::dunezsanalysis(fhicl::ParameterSet const& parameterSet) : EDAnalyzer(parameterSet)
  {
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  // Destructor
  dunezsanalysis::~dunezsanalysis() 
  {}
   
  //-----------------------------------------------------------------------
  void dunezsanalysis::beginJob()
  {

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;
  
    // The arguments to 'make<whatever>' are the same as those passed
    // to the 'whatever' constructor, provided 'whatever' is a ROOT
    // class that TFileService recognizes. 

    // Define our n-tuples, which are limited forms of ROOT
    // TTrees. Start with the TTree itself.

    fReconstructionNtuple = tfs->make<TTree>("dunezsanalysisReconstruction","dunezsanalysisReconstruction");

    // Define the branches (columns) of our simulation n-tuple. When
    // we write a variable, we give the address of the variable to
    // TTree::Branch.
    fReconstructionNtuple->Branch("Event",       &fEvent,          "Event/I");
    fReconstructionNtuple->Branch("Run",         &fRun,            "Run/I");
    // When we write arrays, we give the address of the array to
    // TTree::Branch; in C++ this is simply the array name.
    fReconstructionNtuple->Branch("ChargeSum",   fChargeSum,       "ChargeSum[200]/D");
  }
   
  //-----------------------------------------------------------------------
  void dunezsanalysis::beginRun(const art::Run& /*run*/)
  {
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${SRT_PUBLIC_CONTEXT}/SimpleTypesAndConstants/PhysicalConstants.h.
    art::ServiceHandle<sim::LArG4Parameters> larParameters;
    //fElectronsToGeV = 1./larParameters->GeVToElectrons();
  }

  //-----------------------------------------------------------------------
  void dunezsanalysis::reconfigure(fhicl::ParameterSet const& /*p*/)
  {
    // no parameters for now.  Just read in raw data
    return;
  }

  //-----------------------------------------------------------------------
  void dunezsanalysis::analyze(const art::Event& event) 
  {
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();

    auto rawDigitHandle = event.getHandle< std::vector<raw::RawDigit> >("daq");

    for (int i=0;i<200;i++)
      {
	fChargeSum[i] = 0;
      }

    unsigned int nNeighbors = 0;

    // add up all the digits on the collection wires in the entire event, for each value of
    // the zero-suppression threshold

    // put it in a more easily usable form
    std::vector< art::Ptr<raw::RawDigit> >  Digits;
    art::fill_ptr_vector(Digits, rawDigitHandle);
    //loop through all RawDigits (over entire channels)
    for(size_t d = 0; d < Digits.size(); d++)
      {
	art::Ptr<raw::RawDigit> digit;
	digit=Digits.at(d);
	uint32_t chan = digit->Channel();
	std::vector<short> uncompressed(digit->Samples());
	if (fGeometry->View(chan) == geo::kZ)  // for now only do charge sums for collection hits
	  {
	    raw::Uncompress(digit->ADCs(), uncompressed, digit->Compression());
	    for(unsigned int tick=0;tick<uncompressed.size();tick++) 
	      {
		//std::cout << "trjadc: " << fEvent << " " << chan << " " << uncompressed.at(tick) << std::endl;
		unsigned int tlow = (tick < nNeighbors)? 0: tick - nNeighbors;
		unsigned int thigh = tick + nNeighbors;
		if (thigh>=uncompressed.size()) thigh = uncompressed.size()-1;
		for (short zscut=0;zscut<200;zscut++)
		  {
		    if (uncompressed.at(tick) >= zscut)
		      {
			fChargeSum[zscut] += uncompressed.at(tick);
		      }
		    else
		      {
			for (unsigned int k=tlow;k<=thigh;k++)
			  {
			    if (uncompressed.at(k) >= zscut) 
			      {
				fChargeSum[zscut] += uncompressed.at(tick);
				break;
			      }
			  }
		      }
		  }
	      }
	  }
      }


    fReconstructionNtuple->Fill();
    return;
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see dunezsanalysis.fcl for more information.
  DEFINE_ART_MODULE(dunezsanalysis)

    } // namespace dunezsanalysis

#endif // dunezsanalysis_Module
