#ifndef dEdx_Module
#define dEdx_Module

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

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

// ROOT includes.
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "THStack.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

///////////////////////////////////////////////////////////////////////
//  This module is used to calculate dE/dx from simulated (cheated)  //
//  information.  It was intended for a e/gamma separation study,    //
//  but can be adapted for similar purposes if desired.              //
//                                                                   //
//  Kevin R. Wood -- krwood214@gmail.com (07/23/2014)                //
///////////////////////////////////////////////////////////////////////


namespace AnalysisExample
{

  class dEdx : public art::EDAnalyzer
  {
  public:
 
    explicit dEdx(fhicl::ParameterSet const& pset);
    virtual ~dEdx();
    void beginJob();
    void beginRun(const art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt); 
    void endJob();
    void findDepositionVertex(std::vector<sim::SimChannel> SCHandle);
    void checkDepositionVertex(std::vector<sim::SimChannel> SCHandle, double x, double y, double z);

  private:

    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<geo::Geometry> fGeom;
    std::string fSimulationProducerLabel;

    // dE/dx "cheated" information
    TH1D* fMCdEdxHist;
    std::vector<double> fMCdEdxVec;
    double fMaxMCdEdx = 0.;
    double dE;

    // distance of energy deposition from conversion point
    double trackLength;

    // energy deposition 3d positions for conversion points, ~2.5cm of "track"
    // upstream of conversion point (for dE/dx), and shower respectively
    TH3D* fStartEdepHist;
    TH3D* fTrackEdepHist;
    TH3D* fShowerEdepHist;

    // list of event numbers corresponding to unusually low dEdx measurement
    std::vector<unsigned int> fNullList;

    // position of conversion point
    double xi = 0.;
    double yi = 0.;
    double zi = 1000.;

    // for finding conversion point
    double ziFalse = -500.;
    unsigned int nearDeps = 0;
    bool goodDepVertex = false;

    // for conversion distance
    double fConvDist;
    TH1D* fConvDistHist;
    double Vxi;           //
    double Vyi;           // primary vertex coordinates
    double Vzi;           //

  };

  //-----------------------------------------------------------------------

  dEdx::dEdx(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    this->reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------

  dEdx::~dEdx()
  {
  }
   
  //-----------------------------------------------------------------------

  void dEdx::beginJob()
  {

    // hard-coded fd4apa cryostat dimensions
    fStartEdepHist = tfs->make<TH3D>("fStartEdepHist","fStartEdepHist",50,-240.,240.,50,-750.,800.,50,-50.,560.);
    fTrackEdepHist = tfs->make<TH3D>("fTrackEdepHist","fTrackEdepHist",50,-240.,240.,50,-750.,800.,50,-50.,560.);
    fShowerEdepHist = tfs->make<TH3D>("fShowerEdepHist","fShowerEdepHist",50,-240.,240.,50,-750.,800.,50,-50.,560.);

    fStartEdepHist->GetXaxis()->SetTitle("x");
    fStartEdepHist->GetYaxis()->SetTitle("y");
    fStartEdepHist->GetZaxis()->SetTitle("z");

    fShowerEdepHist->GetXaxis()->SetTitle("x");
    fShowerEdepHist->GetYaxis()->SetTitle("y");
    fShowerEdepHist->GetZaxis()->SetTitle("z");

    fTrackEdepHist->GetXaxis()->SetTitle("x");
    fTrackEdepHist->GetYaxis()->SetTitle("y");
    fTrackEdepHist->GetZaxis()->SetTitle("z");

    fConvDistHist = tfs->make<TH1D>("fConvDistHist","fConvDistHist",100,0.,50.);

  }

  //-----------------------------------------------------------------------

  void dEdx::beginRun(const art::Run& run)
  {
  }

  //-----------------------------------------------------------------------

  void dEdx::reconfigure(fhicl::ParameterSet const& p)
  {
    fSimulationProducerLabel = p.get<std::string>("SimulationLabel");
  }

  //-----------------------------------------------------------------------

  void dEdx::analyze( const art::Event& event )
  {

    auto simChanHandle = event.getHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);

    // initialize at the start of every event
    dE = 0.;
    trackLength = 0.;
    goodDepVertex = false;
    unsigned int nTries = 0;
    ziFalse = -500.;

    // find conversion point
    while(!goodDepVertex && nTries < 300)
    {
      this->findDepositionVertex( (*simChanHandle) );
      this->checkDepositionVertex( (*simChanHandle), xi, yi, zi );
      nTries++;
    }

    std::cout << "Number of attemps at finding deposition vertex = " << nTries << "." << std::endl;

    // find energy deposited into first 2.5 cm of "track"
    if(nTries < 300 && goodDepVertex)
    {
      fStartEdepHist->Fill(xi,yi,zi);
      for(auto const& channel : (*simChanHandle))
      {
	if(fGeom->SignalType(channel.Channel()) == geo::kCollection)
	{
	  auto const& timeSlices = channel.TDCIDEMap();
	  for(auto const& t : timeSlices)
	  {
	    auto const& eDeps = t.second;
	    for(auto const& eDep : eDeps)
	    {
	      fShowerEdepHist->Fill(eDep.x,eDep.y,eDep.z);
	      trackLength = std::sqrt( (eDep.x-xi)*(eDep.x-xi) + (eDep.y-yi)*(eDep.y-yi) + (eDep.z-zi)*(eDep.z-zi) );
	      if(trackLength < 2.5 /* && eDep.z >= zi*/)
	      {
		dE += eDep.energy;
		fTrackEdepHist->Fill(eDep.x,eDep.y,eDep.z);
     	      }
	    } // every energy deposition
	  } // in every time slice
	} // for the collection plane
      } // and for every wire
      dE /= 2.5;
      fMCdEdxVec.push_back(dE);
      if(dE > fMaxMCdEdx) fMaxMCdEdx = dE;
      if(dE < 0.5) fNullList.push_back(event.id().event());
    }

    auto particleHandle = event.getHandle< std::vector<simb::MCParticle> >(fSimulationProducerLabel);

    for(auto const& particle : (*particleHandle))
    {
      if(particle.Process() == "primary")
      {
	Vxi = particle.Vx();
	Vyi = particle.Vy();
	Vzi = particle.Vz();
      }
    }
      fConvDist =  std::sqrt( (Vxi-xi)*(Vxi-xi) + (Vyi-yi)*(Vyi-yi) + (Vzi-zi)*(Vzi-zi) );
      fConvDistHist->Fill(fConvDist);
  }

  //-----------------------------------------------------------------------

  void dEdx::endJob()
  {
    fMCdEdxHist = tfs->make<TH1D>("fMCdEdxHist",";dE/dx (MeV/cm)",200,0.0,10.);

    for(unsigned int e = 0; e < fMCdEdxVec.size(); e++)
    {
      fMCdEdxHist->Fill(fMCdEdxVec[e]);
    }

    std::cout << "Events with < 0.5 dE/dx entry: ";
    for(unsigned int a = 0; a < fNullList.size(); a++) std::cout << fNullList[a] << ", ";

  }

  void dEdx::findDepositionVertex(std::vector<sim::SimChannel> SCHandle)
  {
    zi = 1000.;
    std::cout << "Finding Deposition Vertex... " << std::endl;
    for(auto const& channel : SCHandle)
    {
      if(fGeom->SignalType(channel.Channel()) == geo::kCollection)
      {
	auto const& timeSlices = channel.TDCIDEMap();
	for(auto const& t : timeSlices)
        {
	  auto const& eDeps = t.second;
	  for(auto const& eDep : eDeps)
	  {
	    if( (eDep.z < zi) && (eDep.energy > 0.) && (eDep.z > ziFalse) )
            {
	      xi = eDep.x;
	      yi = eDep.y;
	      zi = eDep.z;
	    } // if z position is a "valid" minimum
	  } // energy deposition loop
	} // time slice loop 
      } // if collection wire
    } // sim channel loop
  } // findDepositionVertex

  void dEdx::checkDepositionVertex(std::vector<sim::SimChannel> SCHandle, double x,double y, double z)
  {
    std::cout << "Checking for Deposition Vertex" << std::endl;
    nearDeps = 0;
    double track = 0;
    for(auto const& channel : SCHandle)
    {
      if(fGeom->SignalType(channel.Channel()) == geo::kCollection)
      {
	auto const& timeSlices = channel.TDCIDEMap();
	for(auto const& t : timeSlices)
	{
	  auto const& eDeps = t.second;
	  for(auto const& eDep : eDeps)
	  {
	    track = std::sqrt( (eDep.x-x)*(eDep.x-x) + (eDep.y-y)*(eDep.y-y) + (eDep.z-z)*(eDep.z-z) );
	    if(track < 2.5 && eDep.energy > 0.)
	    {
	      nearDeps++;
	    } // if z position is minimum
	  } // energy deposition loop
	} // time slice loop 
      } // if collection wire
    } // sim channel loop
    if( nearDeps > 30) goodDepVertex = true;
    else ziFalse = z + 0.05;
  } // checkDepositionVertex

  DEFINE_ART_MODULE(dEdx)

} // namespace AnalysisExample

#endif // dEdx_Module
