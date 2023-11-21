////////////////////////////////////////////////////////////////////////
// Class:       TimeDist
// Module Type: analyzer
// File:        TimeDist_module.cc
//
// Calculate differences between flash times and hit times
//
// Celio Moura camj@fnal.gov celio.moura@ufabc.edu.br
//
////////////////////////////////////////////////////////////////////////

#ifndef TimeDist_Module
#define TimeDist_Module

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TH1.h"
//#include "TH2.h"
#include "TTree.h"
#include "TVector3.h"

// C++ Includes
#include <vector>
#include <string>

namespace TimeDist {

  class TimeDist : public art::EDAnalyzer
  {
  public:
 
    explicit TimeDist(fhicl::ParameterSet const& parameterSet);

    virtual void beginJob() override;
    virtual void reconfigure(fhicl::ParameterSet const& parameterSet) ;
    virtual void analyze (const art::Event& event) override;

  private:

    std::string fHitProducerLabel;        ///< The name of the producer that created hits
    std::string fFlashProducerLabel;      ///< The name of the producer that created flashes
    TH1D* fTimeHist;     ///< Hit time of all particles
    TH1D* fFlashHist;    ///< Flash time of all particles
    TH1D* fTmFshHist;   ///< Hit times minus Flash times 
    TH1D* fTmFshHistU;   ///< Hit times minus Flash times 
    TH1D* fTmFshHistV;   ///< Hit times minus Flash times 
    TH1D* fTmFshHistW;   ///< Hit times minus Flash times 

    double frequency;
    double hittime;

  }; // class TimeDist

  TimeDist::TimeDist(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  void TimeDist::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;

    fTimeHist     = tfs->make<TH1D>("timehist",";Histogram of Hit Times;",4500, -1000, 17000);
    fFlashHist    = tfs->make<TH1D>("flashhist",";Histogram of Flash Times;",180, -1000, 17000);
    fTmFshHist    = tfs->make<TH1D>("hit-flash_Times"  ,";Histogram of Hit-Flash Times;",7000, -3000, 4000);
    fTmFshHistU   = tfs->make<TH1D>("hit-flash_Times_U",";Histogram of Hit-Flash Times;",7000, -3000, 4000);
    fTmFshHistV   = tfs->make<TH1D>("hit-flash_Times_V",";Histogram of Hit-Flash Times;",7000, -3000, 4000);
    fTmFshHistW   = tfs->make<TH1D>("hit-flash_Times_W",";Histogram of Hit-Flash Times;",7000, -3000, 4000);
  }

  //-----------------------------------------------------------------------
  void TimeDist::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    fHitProducerLabel        = parameterSet.get< std::string >("HitLabel");
    fFlashProducerLabel      = parameterSet.get< std::string >("FlashLabel");
  }

  //-----------------------------------------------------------------------
  void TimeDist::analyze(const art::Event& event) 
  {
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);

    auto hitHandle = event.getHandle< std::vector<recob::Hit> >(fHitProducerLabel);

    // For every Hit:
    for ( auto const& hit : (*hitHandle) )
      {
        frequency = clockData.TPCClock().Frequency();
	hittime = hit.PeakTime()/frequency;

	fTimeHist->Fill(hittime);  // filling the historgram with hit times

      } // for each Hit

    auto flashHandle = event.getHandle< std::vector<recob::OpFlash> >(fFlashProducerLabel);

    // For every Flash:
    for ( auto const& opflash : (*flashHandle) )
      {
	// The channel associated with this flash.
	fFlashHist->Fill(opflash.Time());
      } // for each Flash

    for ( auto const& hit : (*hitHandle) )
      {
        double frequency = clockData.TPCClock().Frequency();
	double hittime = hit.PeakTime()/frequency;
	for ( auto const& opflash : (*flashHandle) )
	  {
	    fTmFshHist->Fill(hittime-opflash.Time());	    
	  }
	if (hit.View()==geo::kU)
	  {
	    for ( auto const& opflash : (*flashHandle) )
	      {
		fTmFshHistU->Fill(hittime-opflash.Time());	    
	      }
	  }
	else if (hit.View()==geo::kV)
	  {
	    for ( auto const& opflash : (*flashHandle) )
	      {
		fTmFshHistV->Fill(hittime-opflash.Time());	    
	      }
	  }
	else
	  {
	    for ( auto const& opflash : (*flashHandle) )
	      {
		fTmFshHistW->Fill(hittime-opflash.Time());	    
	      }
	  }
      } // hit for
  } // TimeDist::analyze()

  DEFINE_ART_MODULE(TimeDist)

} // namespace TimeDist

#endif // TimeDist_Module
