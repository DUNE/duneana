////////////////////////////////////////////////////////////////////////
// Class:       TimeDist
// Module Type: analyzer
// File:        RecoTrack_module.cc
//
// Under development: aims to calculate electron lifetime
//
// Celio Moura camj@fnal.gov celio.moura@ufabc.edu.br
//
////////////////////////////////////////////////////////////////////////
// Run like this:
// lar -c RecoTrack.fcl -s /dune/data/users/camj/input/prodcosmics_dune35t_milliblock_0_20150827T232050_merged.root -T recotrack.root -n 4

#ifndef RecoTrack_Module
#define RecoTrack_Module

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TVector3.h"

#include "TF1.h"
//#include "TCanvas.h"

// C++ Includes
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#define PI 3.14159265

namespace RecoTrack {

  class RecoTrack : public art::EDAnalyzer
  {
  public:
 
    explicit RecoTrack(fhicl::ParameterSet const& parameterSet);

    virtual void beginJob() override;
    virtual void reconfigure(fhicl::ParameterSet const& parameterSet) ;
    virtual void analyze (const art::Event& event) override;

  private:

    std::string fTrackProducerLabel;
    std::string fHitProducerLabel;        ///< The name of the producer that created hits
    TH1D* fTimeHist;     ///< Hit time of all particles
    TH1D* fTimeHist1;     ///< Hit time of all particles
    TH1D* fTimeHist2;     ///< Hit time of all particles
    TH1D* fTimeHist3;     ///< Hit time of all particles
    TH1D* fTimeHist4;     ///< Hit time of all particles
    TH1D* fChargeADCHist1;
    TH1D* fChargeADCHist2;
    TH1D* fChargeADCHist3;
    TH1D* fChargeADCHist31;
    TH1D* fChargeADCHist32;
    TH1D* fChargeADCHist33;
    TH1D* fChargeADCHist34;
    TH1D* fChargeADCHist35;
    TH1D* fChargeADCHist36;
    TH1D* fChargeADCHist37;
    TH1D* fChargeADCHist38;
    TH1D* fChargeADCHist4;
    TH1D* fTrackLengthHist;
    TH1D* fTrackThetaHist;
    TH1D* fTrackPhiHist;

    TH2D* fChADCvsTM;
    TH2D* fLogChADCvsTM;
    TH2D* fLogChADCvsTMkZ;
    TH2D* fLogChADCvsTM1;
    TH2D* fLogChADCvsTM2;
    TH2D* fLogChADCvsTM3;
    TH2D* fLogChADCvsTM31;
    TH2D* fLogChADCvsTM32;
    TH2D* fLogChADCvsTM33;
    TH2D* fLogChADCvsTM34;
    TH2D* fLogChADCvsTM35;
    TH2D* fLogChADCvsTM36;
    TH2D* fLogChADCvsTM37;
    TH2D* fLogChADCvsTM38;
    TH2D* fLogChADCvsTM4;

    //////    TProfile* prof;
    TF1 *f1;
    TF1 *f2;
    TF1 *f3;
    TF1 *f31;
    TF1 *f32;
    TF1 *f33;
    TF1 *f34;
    TF1 *f35;
    TF1 *f36;
    TF1 *f37;
    TF1 *f38;
    TF1 *f4;

    // for c2: remove unused variables
    //TF1 *g31,*g32,*g33,*g34,*g35,*g36,*g37;
    TF1 *g31,*g32;

    TH1D* fXHist;
    TH1D* fYHist;
    TH1D* fZHist;

    TH1D* fNofTrackHits;
    TH1D* fNofTrackWires;

    TTree* fSimulationNtuple;

    double frequency;
    double hittime;
    double hittime2;
    double chargeADC2;
    double chargeADC2Log;

    double chargeADCmax1,chargeADCmin1;
    double chargeADCmax2,chargeADCmin2;
    double chargeADCmax3,chargeADCmin3;
    double chargeADCmax4,chargeADCmin4;
    double chargeADCmax5,chargeADCmin5;
    double chargeADCmax6,chargeADCmin6;
    double chargeADCmax7,chargeADCmin7;

    double trackTheta;
    double trackPhi;

    double x;
    double y;
    double z;

    std::vector<double> hittest,timetest;
    std::vector<double> hittest2,timetest2;
    std::vector<double> hittest3,timetest3;
    std::vector<double> hittest4,timetest4;
    std::vector<double> hittest5,timetest5;
    std::vector<double> hittest6,timetest6;
    std::vector<double> hittest7,timetest7;
    std::vector<double> hittest8,timetest8;
    std::vector<int> wiretest,timevec;

    // for c2: these variables are not currently used
    //double hittimemin;
    //double hittimemax;

    // for c2: remove unused variables
    //double min1,min2,min3,min4,min5,min6,min7;
    //double max1,max2,max3,max4,max5,max6,max7;
    double min3;
    double max3;

    std::ofstream myfile;

    int wiretmp,wirehit,wirecount;

  }; // class RecoTrack

  RecoTrack::RecoTrack(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  void RecoTrack::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;

    fTrackLengthHist = tfs->make<TH1D>("trackVecLength",  ";particle track length (cm);", 500, -20, 480);
    fNofTrackWires = tfs->make<TH1D>("wires",";number of wires hit in a track;",500,-20,480);
    fNofTrackHits = tfs->make<TH1D>("hits",";number of hits in a track;",1000,-10,990);
    fTimeHist     = tfs->make<TH1D>("timehist",";Histogram of Hit Times;",1000, -10, 990);
    fTimeHist1     = tfs->make<TH1D>("timehist1",";Histogram2 of Hit Times;",500, -1000, 17000);
    fTimeHist2     = tfs->make<TH1D>("timehist2",";Histogram2 of Hit Times;",500, -1000, 17000);
    fTimeHist3     = tfs->make<TH1D>("timehist3",";Histogram of Hit Times;",500, -1000, 17000);
    fTimeHist4     = tfs->make<TH1D>("timehist4",";Histogram of Hit Times;",500, -1000, 17000);
    fChargeADCHist1     = tfs->make<TH1D>("ch-adc-hist1",";Histogram of Charge;",100, 100, 1100);
    fChargeADCHist2     = tfs->make<TH1D>("ch-adc-hist2",";Histogram of Charge;",30, 0, 1200);
    fChargeADCHist3     = tfs->make<TH1D>("ch-adc-hist3",";Histogram of Charge;",40, 0, 1400);
    fChargeADCHist31     = tfs->make<TH1D>("ch-adc-hist31",";Histogram of Charge;",30, 200, 800);
    fChargeADCHist32     = tfs->make<TH1D>("ch-adc-hist32",";Histogram of Charge;",30, 200, 800);
    fChargeADCHist33     = tfs->make<TH1D>("ch-adc-hist33",";Histogram of Charge;",30, 200, 800);
    fChargeADCHist34     = tfs->make<TH1D>("ch-adc-hist34",";Histogram of Charge;",30, 200, 800);
    fChargeADCHist35     = tfs->make<TH1D>("ch-adc-hist35",";Histogram of Charge;",30, 200, 800);
    fChargeADCHist36     = tfs->make<TH1D>("ch-adc-hist36",";Histogram of Charge;",30, 200, 800);
    fChargeADCHist37     = tfs->make<TH1D>("ch-adc-hist37",";Histogram of Charge;",30, 200, 800);
    fChargeADCHist38     = tfs->make<TH1D>("ch-adc-hist38",";Histogram of Charge;",30, 200, 800);
    fChargeADCHist4     = tfs->make<TH1D>("ch-adc-hist4",";Histogram of Charge;",20, 100, 700);
    fChADCvsTM    = tfs->make<TH2D>("chargeADC.vs.time",";ChargeADC vs Time 2D Plot;",7000,-500,6500,2000,0,2000);
    fLogChADCvsTM = tfs->make<TH2D>("chargeADCLog.vs.time",";Log(ChargeADC) vs Time 2D Plot;",7000,-500,6500,120,3,9);
    fLogChADCvsTMkZ = tfs->make<TH2D>("chargeADCLogkZ.vs.time",";Log(ChargeADC) vs Time 2D Plot;",7000,-500,6500,120,3,9);
    fLogChADCvsTM1 = tfs->make<TH2D>("chargeADCLog_1.vs.time",";Log(ChargeADC) vs Time 2D Plot;",300,1500,1800,120,3,9);
    fLogChADCvsTM2 = tfs->make<TH2D>("chargeADCLog_2.vs.time",";Log(ChargeADC) vs Time 2D Plot;",500,6500,7000,120,3,9);
    fLogChADCvsTM3 = tfs->make<TH2D>("chargeADCLog_3.vs.time",";Log(ChargeADC) vs Time 2D Plot;",1000,11600,12600,120,3,9);
    fLogChADCvsTM31 = tfs->make<TH2D>("chargeADCLog_31.vs.time",";Log(ChargeADC) vs Time 2D Plot;",200,11650,11850,60,3,9);
    fLogChADCvsTM32 = tfs->make<TH2D>("chargeADCLog_32.vs.time",";Log(ChargeADC) vs Time 2D Plot;",200,11750,11950,60,3,9);
    fLogChADCvsTM33 = tfs->make<TH2D>("chargeADCLog_33.vs.time",";Log(ChargeADC) vs Time 2D Plot;",7000,-500,1500,120,3,9);
    fLogChADCvsTM34 = tfs->make<TH2D>("chargeADCLog_34.vs.time",";Log(ChargeADC) vs Time 2D Plot;",200,11950,12150,60,3,9);
    fLogChADCvsTM35 = tfs->make<TH2D>("chargeADCLog_35.vs.time",";Log(ChargeADC) vs Time 2D Plot;",200,12050,12250,60,3,9);
    fLogChADCvsTM36 = tfs->make<TH2D>("chargeADCLog_36.vs.time",";Log(ChargeADC) vs Time 2D Plot;",200,12150,12350,60,3,9);
    fLogChADCvsTM37 = tfs->make<TH2D>("chargeADCLog_37.vs.time",";Log(ChargeADC) vs Time 2D Plot;",200,12250,12450,60,3,9);
    fLogChADCvsTM38 = tfs->make<TH2D>("chargeADCLog_38.vs.time",";Log(ChargeADC) vs Time 2D Plot;",1000,11600,12600,60,3,9);
    fLogChADCvsTM4 = tfs->make<TH2D>("chargeADCLog_4.vs.time",";Log(ChargeADC) vs Time 2D Plot;",200,12750,12950,120,3,9);
    f1 = tfs->make<TF1>("f1","[0]+[1]*x",1500,1800);
    f2 = tfs->make<TF1>("f2","[0]+[1]*x",6500,7000);
    f3 = tfs->make<TF1>("f3","[0]+[1]*x",11600,12600);
    f31 = tfs->make<TF1>("f31","[0]+[1]*x",11650,11850);
    f32 = tfs->make<TF1>("f32","[0]+[1]*x",11750,11950);
    f33 = tfs->make<TF1>("f33","[0]+[1]*x",11850,12050);
    f34 = tfs->make<TF1>("f34","[0]+[1]*x",11950,12150);
    f35 = tfs->make<TF1>("f35","[0]+[1]*x",12050,12250);
    f36 = tfs->make<TF1>("f36","[0]+[1]*x",12150,12350);
    f37 = tfs->make<TF1>("f37","[0]+[1]*x",12250,12450);
    f38 = tfs->make<TF1>("f38","[0]+[1]*x",11600,12600);
    f4 = tfs->make<TF1>("f4","[0]+[1]*x",12750,12950);
    g31 = tfs->make<TF1>("g31","[0]/([1]*sqrt(2.0*3.14159265))*exp(-pow(x-[2],2.0)/(2.0*[1]*[1]))",200,800);
    g32 = tfs->make<TF1>("g32","[0]/([1]*sqrt(2.0*3.14159265))*exp(-pow(x-[2],2.0)/(2.0*[1]*[1]))",200,800);

    fTrackThetaHist = tfs->make<TH1D>("trackTheta",  ";particle track theta angle (deg);", 440, -20, 200);
    fTrackPhiHist = tfs->make<TH1D>("trackPhi",  ";particle track phi angle (deg);", 800, -200, 200);

    fXHist = tfs->make<TH1D>("x",";projection on x;",20,-1,1);
    fYHist = tfs->make<TH1D>("y",";projection on y;",20,-1,1);
    fZHist = tfs->make<TH1D>("z",";projection on z;",20,-1,1);

    fSimulationNtuple     = tfs->make<TTree>("ElectronLifetime",    "ElectronLifetime");

    fSimulationNtuple->Branch("hittest", "vector<double>", &hittest);
    //    fSimulationNtuple->Branch("chargeSummedADC", "vector<double>", &chargeADCtest);
  }

  //-----------------------------------------------------------------------
  void RecoTrack::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    fHitProducerLabel        = parameterSet.get< std::string >("HitLabel");
    fTrackProducerLabel      = parameterSet.get< std::string >("TrackLabel");
  }

  //-----------------------------------------------------------------------
  void RecoTrack::analyze(const art::Event& event) 
  {
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event); // to get TPC clock and frequency
    frequency = clockData.TPCClock().Frequency();

    // get information about the hits
    auto hitHandle = event.getHandle< std::vector<recob::Hit> >(fHitProducerLabel);

    // get track information
    auto trackHandle = event.getHandle< std::vector<recob::Track> >(fTrackProducerLabel);
    std::vector<art::Ptr<recob::Track> > trackVec;
    art::fill_ptr_vector(trackVec,trackHandle);

    art::FindManyP<recob::Hit> fmtht(trackHandle, event, fTrackProducerLabel); // to associate tracks and hits

    chargeADCmin1 = 10000.0;
    chargeADCmax1 = 0.0;
    chargeADCmin2 = 10000.0;
    chargeADCmax2 = 0.0;
    chargeADCmin3 = 10000.0;
    chargeADCmax3 = 0.0;
    chargeADCmin4 = 10000.0;
    chargeADCmax4 = 0.0;
    chargeADCmin5 = 10000.0;
    chargeADCmax5 = 0.0;
    chargeADCmin6 = 10000.0;
    chargeADCmax6 = 0.0;
    chargeADCmin7 = 10000.0;
    chargeADCmax7 = 0.0;

    for ( size_t itrack=0 ; itrack < trackVec.size() ; ++itrack ) // Starts looping over tracks
      {
	/******************/
	/* Tracks' Angles */
	/******************/
	art::Ptr<recob::Track> ptrack(trackHandle, (int)itrack);
	const recob::Track& track = *ptrack;

	int ntraj = track.NumberTrajectoryPoints();
	if(ntraj > 0)
	  {
	    auto dir = track.VertexDirection();
	    trackTheta=dir.Theta()*180.0/PI;
	    fTrackThetaHist->Fill(trackTheta);
	    trackPhi=dir.Phi()*180.0/PI;
	    fTrackPhiHist->Fill(trackPhi);
	  }

	x = sin(trackTheta)*cos(trackPhi); // horizontal drift axis
	z = sin(trackTheta)*sin(trackPhi); // horizontal perpendicular to drift axis
	y = cos(trackTheta); // vertical

	/******************/
	/* End tracks' angles */
	/******************/

	/**********************/
	/* Filling Histograms */
	/**********************/
	fXHist->Fill(x);
	fYHist->Fill(y);
	fZHist->Fill(z);
	/**************************/
	/* End filling Histograms */
	/**************************/


	std::vector< art::Ptr<recob::Hit> > allHits = fmtht.at(itrack);

	/**************************************************************************/
	/* Getting hit times ******************************************************/
	/* Getting wire numbers, ordering, and counting number of different wires */
	/**************************************************************************/
	timevec.clear();       // clean the vector before every time the hits loop starts
	wiretest.clear();      // clean the vector before every time the hits loop starts

	// hittimemin = 20000.0;
	// hittimemax = 0.0;

	for (size_t hits=0; hits < allHits.size(); ++hits) // Loop over hits in a track
	  {
	    hittime = allHits[hits]->PeakTime()/frequency;
	    timevec.push_back(hittime); 

	    // if (hittime > hittimemax) hittimemax = hittime;
	    // if (hititme < hittimemin) hittimemin = hittime;


	    wirehit = allHits[hits]->WireID().Wire;
	    wiretest.push_back(int(wirehit));
	  } // end loop over track hits

	std::sort (timevec.begin(), timevec.end());
	double timeLength = timevec[allHits.size()-1]-timevec[0];

	std::sort (wiretest.begin(), wiretest.end());

	wirecount = 0;
	wiretmp = 60000;

	for (size_t i = 0; i < allHits.size(); ++i) // counting number of different wires
	  {
	    if (wiretest.at(i) != wiretmp)
	      ++wirecount;
	    wiretmp = wiretest.at(i);
	  } // end counting number of different wires

	/*********************************************************************************/
	/* End getting hit times *********************************************************/
	/* Finish getting wire numbers, ordering, and counting number of different wires */
	/*********************************************************************************/

	/**********************/
	/* Filling Histograms */
	/**********************/
	fTimeHist->Fill(timeLength); // Track time length
	fNofTrackWires->Fill(wirecount); // How many different wires in a track
	fNofTrackHits->Fill( allHits.size() ); // How many hits in a track
	fTrackLengthHist->Fill( trackVec[itrack]->Length() ); // Track length
	/**************************/
	/* End filling Histograms */
	/**************************/


	//	if ((wirecount < 300) && (timeLength < 600)) // cuts used for prodcosmics_dune35t_milliblock_0_20150827T232050_merged.root
	// if ((wirecount < 200) && (timeLength < 400))
	//   continue; // If condition is true, skip the rest of the loop in the current iteration

	hittest.clear();
	timetest.clear();
	hittest2.clear();
	timetest2.clear();
	hittest3.clear();
	timetest3.clear();
	hittest4.clear();
	timetest4.clear();
	hittest5.clear();
	timetest5.clear();
	hittest6.clear();
	timetest6.clear();
	hittest7.clear();
	timetest7.clear();
	hittest8.clear();
	timetest8.clear();

	hittime2 = 0;

	for(size_t h=0; h < allHits.size(); ++h)
	  {
	    art::Ptr<recob::Hit> itrack_hit = allHits[h];

	    hittime2 = itrack_hit->PeakTime()/frequency;
	    hittime2 = hittime2 + itrack*2000.0;

	    chargeADC2 = itrack_hit->ROISummedADC();
	    chargeADC2Log = log(chargeADC2);

	    fChADCvsTM->Fill(hittime2,chargeADC2);
	    fLogChADCvsTM->Fill(hittime2,chargeADC2Log);


	    if (itrack_hit->View() == geo::kZ)
	      {

		fLogChADCvsTMkZ->Fill(hittime2,chargeADC2Log);


	    // /* SELECTING TRACKS BY TIME 1 */
	    // 	//if ((hittime2 < 2000) && (hittime2 > 1000))
	    // 	if (itrack == 0)
	    // 	  {
	    // 	//chargeADC2 = itrack_hit->SummedADC();
	    // 	//		if (itrack_hit->View() == geo::kZ)
	    // 	//{
	    // 	    std::cout << "testing0" <<std::endl;
	    // 	    std::cout << "Track " << itrack << "; Hit charge " << h << " = " << chargeADC2 << std::endl;
	    // 	    fTimeHist1->Fill(hittime2);
	    // 	    fChargeADCHist1->Fill(chargeADC2);
	    // 	    if (chargeADC2 < 500) {
	    // 	    fLogChADCvsTM1->Fill(hittime2,chargeADC2Log);
	    // 	    }
	    // 	    //} // end if
	    // 	  } // end itrack == 0 // end if hittime2
	    // /* END SELECTING TRACKS BY TIME 1 */

	    // /* SELECTING TRACKS BY TIME 2 */
	    // if ((hittime2 < 8000) && (hittime2 > 6000))
	    //   {
	    // 	    fTimeHist2->Fill(hittime2);
	    // 	    fChargeADCHist2->Fill(chargeADC2);

	    // 	    //		    std::cout << "Track " << itrack << "; Hit charge " << h << " = " << chargeADC2 << std::endl;
	    // 	    //		    hittest2.push_back(chargeADC2);

	    // 	    //if ((chargeADC2Log < lcharge_arb*1.5) && (chargeADC2Log > lcharge_arb*0.5)) {
	    // 	    if ((chargeADC2 < 1500) && (chargeADC2 > 500)) {
	    // 	    fLogChADCvsTM2->Fill(hittime2,chargeADC2Log);
	    // 	    }

	    //   } // end if hittime2
	    // /* END SELECTING TRACKS BY TIME 2 */

	    // /* SELECTING TRACKS BY TIME 3 */
	    // /*>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<*/
	    // /*>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<*/
	    // /*>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<*/



	    // if ((hittime2 < 11800) && (hittime2 > 11700))
	    //   {
	    // 	// fChargeADCHist31->Fill(chargeADC2);
	    // 	// fLogChADCvsTM31->Fill(hittime2,chargeADC2Log);
	    // 	hittest.push_back(chargeADC2);
	    // 	timetest.push_back(hittime2);
	    //   }

	    // if ((hittime2 < 11900) && (hittime2 > 11800))
	    //   {
	    // 	// fChargeADCHist32->Fill(chargeADC2);
	    // 	// fLogChADCvsTM32->Fill(hittime2,chargeADC2Log);
	    // 	hittest2.push_back(chargeADC2);
	    // 	timetest2.push_back(hittime2);
	    //   }

	    if (hittime2 < 2000)
	      {
		fChargeADCHist33->Fill(chargeADC2);
		fLogChADCvsTM33->Fill(hittime2,chargeADC2Log);

		if (chargeADC2 > chargeADCmax3) chargeADCmax3 = chargeADC2;
		if (chargeADC2 < chargeADCmin3) chargeADCmin3 = chargeADC2;

		hittest3.push_back(chargeADC2);
		timetest3.push_back(hittime2);
	      }

	    // if ((hittime2 < 12100) && (hittime2 > 12000))
	    //   {
	    // 	// fChargeADCHist34->Fill(chargeADC2);
	    // 	// fLogChADCvsTM34->Fill(hittime2,chargeADC2Log);
	    // 	if (chargeADC2 > chargeADCmax4) chargeADCmax4 = chargeADC2;
	    // 	if (chargeADC2 < chargeADCmin4) chargeADCmin4 = chargeADC2;

	    // 	hittest4.push_back(chargeADC2);
	    // 	timetest4.push_back(hittime2);
	    //   }

	    // if ((hittime2 < 12200) && (hittime2 > 12100))
	    //   {
	    // 	// fChargeADCHist35->Fill(chargeADC2);
	    // 	// fLogChADCvsTM35->Fill(hittime2,chargeADC2Log);
	    // 	if (chargeADC2 > chargeADCmax5) chargeADCmax5 = chargeADC2;
	    // 	if (chargeADC2 < chargeADCmin5) chargeADCmin5 = chargeADC2;

	    // 	hittest5.push_back(chargeADC2);
	    // 	timetest5.push_back(hittime2);
	    //   }

	    // if ((hittime2 < 12300) && (hittime2 > 12200))
	    //   {
	    // 	// fChargeADCHist36->Fill(chargeADC2);
	    // 	// fLogChADCvsTM36->Fill(hittime2,chargeADC2Log);
	    // 	if (chargeADC2 > chargeADCmax6) chargeADCmax6 = chargeADC2;
	    // 	if (chargeADC2 < chargeADCmin6) chargeADCmin6 = chargeADC2;

	    // 	hittest6.push_back(chargeADC2);
	    // 	timetest6.push_back(hittime2);
	    //   }

	    // if ((hittime2 < 12400) && (hittime2 > 12300))
	    //   {
	    // 	// fChargeADCHist37->Fill(chargeADC2);
	    // 	// fLogChADCvsTM37->Fill(hittime2,chargeADC2Log);
	    // 	if (chargeADC2 > chargeADCmax7) chargeADCmax7 = chargeADC2;
	    // 	if (chargeADC2 < chargeADCmin7) chargeADCmin7 = chargeADC2;

	    // 	hittest7.push_back(chargeADC2);
	    // 	timetest7.push_back(hittime2);
	    //   }

	    // if ((hittime2 < 12600) && (hittime2 > 11600))
	    //   {
	    // 	std::cout << "testing1" <<std::endl;
	    // 	std::cout << "Track " << itrack << "; Hit charge " << h << " = " << chargeADC2 << std::endl;
	    // 	//chargeADC2 = itrack_hit->SummedADC();
	    // 	fTimeHist3->Fill(hittime2);
	    // 	fChargeADCHist3->Fill(chargeADC2);
	    // 	//if (chargeADC2 < 700) {
	    // 	  fLogChADCvsTM3->Fill(hittime2,chargeADC2Log);
	    // 	  //}

	    //   } // end if hittime2


	    // std::cout << "testing11" <<std::endl;




	    // /*>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<*/
	    // /*>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<*/
	    // /*>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<*/
	    // /* END SELECTING TRACKS BY TIME 3 */
	    
	    // /* SELECTING TRACKS BY TIME 4 */
	    // if ((hittime2 < 13000) && (hittime2 > 12700))
	    //   {
	    // 	std::cout << "testing2" <<std::endl;
	    // 	std::cout << "Track " << itrack << "; Hit charge " << h << " = " << chargeADC2 << std::endl;
	    // 	//chargeADC2 = itrack_hit->SummedADC();
	    // 	//if (itrack_hit->View() == geo::kZ)
	    // 	//{
	    // 	    fTimeHist4->Fill(hittime2);
	    // 	    fChargeADCHist4->Fill(chargeADC2);
	    // 	    if (chargeADC2 < 500) {
	    // 	    fLogChADCvsTM4->Fill(hittime2,chargeADC2Log);
	    // 	    }
	    // 	    //} // end if
	    //   } // end if hittime2
	    // /* END SELECTING TRACKS BY TIME 4 */
	    } // end if geo::kZ

	  } // end h for (all the hits for a track (itrack))


	max3 = chargeADCmax3;
	min3 = chargeADCmin3;

	std::cout<<" max3 = "<< max3 << "; min3 = " << min3 << std::endl;


	for(size_t j=0; j < hittest3.size(); ++j)
	  {
	    if ((hittest3[j] > ((max3-min3)*0.15+min3)) && (hittest3[j] < (max3-(max3-min3)*0.15)))
	      {

		fChargeADCHist33->Fill(hittest3[j]);
		fLogChADCvsTM33->Fill(timetest3[j],log(hittest3[j]));

		hittest8.push_back(hittest3[j]);
		timetest8.push_back(timetest3[j]);

	      }
	  }



      } // end itrack for

    myfile.open ("/dune/app/users/camj/dunelatest/work/lifetime.txt");

    // f1->SetParameters(6.,-0.000333);
    // f1->SetLineColor(kRed);
    // fLogChADCvsTM1->Fit(f1);
    // fLogChADCvsTM1->Draw();
    // f1->Draw("same");
    // //double p0 = f1->GetParameter(0);
    // double p01 = f1->GetParameter(0);
    // double Err01 = f1->GetParError(0);
    // double p1 = f1->GetParameter(1);
    // double Err1 = f1->GetParError(1);

    // //    std::cout << "p1+/-Err1 = (" << -0.001/p1 << " +/- " << (1./p1)*(1./p1)*Err1 <<") ms"<< std::endl;
    // myfile << "p1 = (" << -0.001/p1 << " +/- " << 1.e-6*(1./p1)*(1./p1)*Err1 <<") ms"<< std::endl;
    // myfile << "Fit = " << p1 << " +/- " << Err1 << "; " << p01 << " +/- " << Err01 << std::endl;

    // double maxbin = fChargeADCHist31->GetMaximumBin();
    // double binamp = fChargeADCHist31->GetBinContent(maxbin);
    // double bincent = fChargeADCHist31->GetBinCenter(maxbin);
    // g31->SetParameters(binamp,50.0,bincent);
    // g31->SetLineColor(kRed);
    // fChargeADCHist31->Fit(g31);
    // fChargeADCHist31->Draw();
    // g31->Draw("same");
    // double pamp1 = g31->GetParameter(0);
    // double psig1 = g31->GetParameter(1);
    // double pavr1 = g31->GetParameter(2);
    // myfile << "Amplitude = " << pamp1 << "; Sigma = " << psig1 << "; Average = " << pavr1 << std::endl;
    // myfile << "Bin of the max: " << maxbin << "; Amplitude of the max: " << binamp << "; Bin center: " << bincent << std::endl;

    // g32->SetParameters(4.0,100.0,420.0);
    // g32->SetLineColor(kRed);
    // fChargeADCHist32->Fit(g32);
    // fChargeADCHist32->Draw();
    // g32->Draw("same");
    // double pamp2 = g32->GetParameter(0);
    // double psig2 = g32->GetParameter(1);
    // double pavr2 = g32->GetParameter(2);
    // myfile << "Amplitude = " << pamp2 << "; Sigma = " << psig2 << "; Average = " << pavr2 << std::endl;

    // f2->SetParameters(6.,-0.000333);
    // f2->SetLineColor(kRed);
    // fLogChADCvsTM2->Fit(f2);
    // fLogChADCvsTM2->Draw();
    // f2->Draw("same");
    // double p2 = f2->GetParameter(1);
    // double Err2 = f2->GetParError(1);
    // double p02 = f2->GetParameter(0);
    // double Err02 = f2->GetParError(0);

    // myfile << "p2 = (" << -0.001/p2 << " +/- " << 1.e-6*(1./p2)*(1./p2)*Err2 <<") ms"<< std::endl;
    // myfile << "Fit = " << p2 << " +/- " << Err2 << "; " << p02 << " +/- " << Err02 << std::endl;

    // fLogChADCvsTM3->SetMarkerStyle(kCircle);
    // fLogChADCvsTM3->SetMarkerSize(0.5);
    // f3->SetParameters(6.,-0.000333);
    // fLogChADCvsTM3->Fit(f3);
    // fLogChADCvsTM3->Draw();
    // f3->SetLineColor(kRed);
    // f3->SetTitle("Electron Lifetime");
    // f3->Draw("same");
    // double p3 = f3->GetParameter(1);
    // double Err3 = f3->GetParError(1);
    // double p03 = f3->GetParameter(0);
    // double Err03 = f3->GetParError(0);
    
    // myfile << "p3 = (" << -0.001/p3 << " +/- " << 1.e-6*(1./p3)*(1./p3)*Err3 <<") ms"<< std::endl;
    // myfile << "Fit = " << p3 << " +/- " << Err3 << "; " << p03 << " +/- " << Err03 << std::endl;

    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    // fLogChADCvsTM31->SetMarkerStyle(kCircle);
    // fLogChADCvsTM31->SetMarkerSize(0.5);
    // f31->SetParameters(6.,-0.000333);
    // fLogChADCvsTM31->Fit(f31);
    // fLogChADCvsTM31->Draw();
    // f31->SetLineColor(kRed);
    // f31->SetTitle("Electron Lifetime");
    // f31->Draw("same");
    // double p31 = f31->GetParameter(1);
    // double Err31 = f31->GetParError(1);
    // double p031 = f31->GetParameter(0);
    // double Err031 = f31->GetParError(0);
    
    // myfile << "p31 = (" << -0.001/p31 << " +/- " << 1.e-6*(1./p31)*(1./p31)*Err31 <<") ms"<< std::endl;
    // myfile << "Fit = " << p31 << " +/- " << Err31 << "; " << p031 << " +/- " << Err031 << std::endl;
    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    // fLogChADCvsTM32->SetMarkerStyle(kCircle);
    // fLogChADCvsTM32->SetMarkerSize(0.5);
    // f32->SetParameters(6.,-0.000333);
    // fLogChADCvsTM32->Fit(f32);
    // fLogChADCvsTM32->Draw();
    // f32->SetLineColor(kRed);
    // f32->SetTitle("Electron Lifetime");
    // f32->Draw("same");
    // double p32 = f32->GetParameter(1);
    // double Err32 = f32->GetParError(1);
    // double p032 = f32->GetParameter(0);
    // double Err032 = f32->GetParError(0);
    
    // myfile << "p32 = (" << -0.001/p32 << " +/- " << 1.e-6*(1./p32)*(1./p32)*Err32 <<") ms"<< std::endl;
    // myfile << "Fit = " << p32 << " +/- " << Err32 << "; " << p032 << " +/- " << Err032 << std::endl;
    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    fLogChADCvsTM33->SetMarkerStyle(kCircle);
    fLogChADCvsTM33->SetMarkerSize(0.5);
    f33->SetParameters(6.,-0.000333);
    fLogChADCvsTM33->Fit(f33);
    fLogChADCvsTM33->Draw();
    f33->SetLineColor(kRed);
    f33->SetTitle("Electron Lifetime");
    f33->Draw("same");
    double p33 = f33->GetParameter(1);
    double Err33 = f33->GetParError(1);
    double p033 = f33->GetParameter(0);
    double Err033 = f33->GetParError(0);
    
    myfile << "p33 = (" << -0.001/p33 << " +/- " << 1.e-6*(1./p33)*(1./p33)*Err33 <<") ms"<< std::endl;
    myfile << "Fit = " << p33 << " +/- " << Err33 << "; " << p033 << " +/- " << Err033 << std::endl;
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    // fLogChADCvsTM34->SetMarkerStyle(kCircle);
    // fLogChADCvsTM34->SetMarkerSize(0.5);
    // f34->SetParameters(6.,-0.000333);
    // fLogChADCvsTM34->Fit(f34);
    // fLogChADCvsTM34->Draw();
    // f34->SetLineColor(kRed);
    // f34->SetTitle("Electron Lifetime");
    // f34->Draw("same");
    // double p34 = f34->GetParameter(1);
    // double Err34 = f34->GetParError(1);
    // double p034 = f34->GetParameter(0);
    // double Err034 = f34->GetParError(0);
    
    // myfile << "p34 = (" << -0.001/p34 << " +/- " << 1.e-6*(1./p34)*(1./p34)*Err34 <<") ms"<< std::endl;
    // myfile << "Fit = " << p34 << " +/- " << Err34 << "; " << p034 << " +/- " << Err034 << std::endl;
    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    // fLogChADCvsTM35->SetMarkerStyle(kCircle);
    // fLogChADCvsTM35->SetMarkerSize(0.5);
    // f35->SetParameters(6.,-0.000333);
    // fLogChADCvsTM35->Fit(f35);
    // fLogChADCvsTM35->Draw();
    // f35->SetLineColor(kRed);
    // f35->SetTitle("Electron Lifetime");
    // f35->Draw("same");
    // double p35 = f35->GetParameter(1);
    // double Err35 = f35->GetParError(1);
    // double p035 = f35->GetParameter(0);
    // double Err035 = f35->GetParError(0);
    
    // myfile << "p35 = (" << -0.001/p35 << " +/- " << 1.e-6*(1./p35)*(1./p35)*Err35 <<") ms"<< std::endl;
    // myfile << "Fit = " << p35 << " +/- " << Err35 << "; " << p035 << " +/- " << Err035 << std::endl;
    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    // fLogChADCvsTM36->SetMarkerStyle(kCircle);
    // fLogChADCvsTM36->SetMarkerSize(0.5);
    // f36->SetParameters(6.,-0.000333);
    // fLogChADCvsTM36->Fit(f36);
    // fLogChADCvsTM36->Draw();
    // f36->SetLineColor(kRed);
    // f36->SetTitle("Electron Lifetime");
    // f36->Draw("same");
    // double p36 = f36->GetParameter(1);
    // double Err36 = f36->GetParError(1);
    // double p036 = f36->GetParameter(0);
    // double Err036 = f36->GetParError(0);
    
    // myfile << "p36 = (" << -0.001/p36 << " +/- " << 1.e-6*(1./p36)*(1./p36)*Err36 <<") ms"<< std::endl;
    // myfile << "Fit = " << p36 << " +/- " << Err36 << "; " << p036 << " +/- " << Err036 << std::endl;
    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    // fLogChADCvsTM37->SetMarkerStyle(kCircle);
    // fLogChADCvsTM37->SetMarkerSize(0.5);
    // f37->SetParameters(6.,-0.000333);
    // fLogChADCvsTM37->Fit(f37);
    // fLogChADCvsTM37->Draw();
    // f37->SetLineColor(kRed);
    // f37->SetTitle("Electron Lifetime");
    // f37->Draw("same");
    // double p37 = f37->GetParameter(1);
    // double Err37 = f37->GetParError(1);
    // double p037 = f37->GetParameter(0);
    // double Err037 = f37->GetParError(0);
    
    // myfile << "p37 = (" << -0.001/p37 << " +/- " << 1.e-6*(1./p37)*(1./p37)*Err37 <<") ms"<< std::endl;
    // myfile << "Fit = " << p37 << " +/- " << Err37 << "; " << p037 << " +/- " << Err037 << std::endl;
    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    // fLogChADCvsTM38->SetMarkerStyle(kCircle);
    // fLogChADCvsTM38->SetMarkerSize(0.5);
    // f38->SetParameters(6.,-0.000333);
    // fLogChADCvsTM38->Fit(f38);
    // fLogChADCvsTM38->Draw();
    // f38->SetLineColor(kRed);
    // f38->SetTitle("Electron Lifetime");
    // f38->Draw("same");
    // double p38 = f38->GetParameter(1);
    // double Err38 = f38->GetParError(1);
    // double p038 = f38->GetParameter(0);
    // double Err038 = f38->GetParError(0);
    
    // myfile << "p38 = (" << -0.001/p38 << " +/- " << 1.e-6*(1./p38)*(1./p38)*Err38 <<") ms"<< std::endl;
    // myfile << "Fit = " << p38 << " +/- " << Err38 << "; " << p038 << " +/- " << Err038 << std::endl;
    // /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    // f4->SetParameters(6.,-0.000333);
    // f4->SetLineColor(kRed);
    // fLogChADCvsTM4->Fit(f4);
    // fLogChADCvsTM4->Draw();
    // f4->Draw("same");
    // double p4 = f4->GetParameter(1);
    // double Err4 = f4->GetParError(1);
    // double p04 = f4->GetParameter(0);
    // double Err04 = f4->GetParError(0);

    // myfile << "p4 = (" << -0.001/p4 << " +/- " << 1.e-6*(1./p4)*(1./p4)*Err4 <<") ms"<< std::endl;
    // myfile << "Fit = " << p4 << " +/- " << Err4 << "; " << p04 << " +/- " << Err04 << std::endl;

    myfile.close();

    //    TCanvas *c1 = new TCanvas();
    //    f22->Draw("Surf1");
    //    fLogChADCvsTM3->Draw("P0 Same");


    //    std::cout << lengthb << std::endl;

    // For every Hit:
    // for ( auto const& hit : (*hitHandle) )
    //   {
    //  frequency = clockData.TPCClock().Frequency();
    //  	hittime = hit.PeakTime()/frequency;

    //  	hittest.push_back(hittime); // vector with hit times

    //  	fTimeHist->Fill(hittime);  // filling the historgram with hit times


    //  	fChargeADCHist->Fill(chargeADC); // filling the histogram of charge frequency

    //  	if ((chargeADC > 10) && (chargeADC < 500))
    //  	  {
    //  	    fChADCvsTM->Fill(hittime,chargeADC);
    //  	    fLogChADCvsTM->Fill(hittime,chargeADCLog);
    //  	  } // end if for charge integral
    //   } // for each Hit

  //    fSimulationNtuple->Fill();

  } // RecoTrack::analyze()
    
  DEFINE_ART_MODULE(RecoTrack)

} // namespace RecoTrack

#endif // RecoTrack_Module

// 	chargeInteg = hit.Integral(); // charge per hit
// 	chargeIntegLog = log(chargeInteg); // calculating log of charge
// 	chargeInttest.push_back(chargeInteg); // vector with integrated charge for each hit

// 	fChargeIntegHist->Fill(chargeInteg); // filling the histogram of charge frequency

// 	    fCHvsTM->Fill(hittime,chargeInteg);
// 	    fLogCHvsTM->Fill(hittime,chargeIntegLog);

// trackTheta = trackVec[itrack]->Theta()*180.0/PI;
// fTrackThetaHist->Fill(trackTheta);
// trackPhi = trackVec[itrack]->Phi()*180.0/PI;
// fTrackPhiHist->Fill(trackPhi);

	    // chargeADC = itrack_hit->Integral(); // Charge per hit. Compare with ->Integral()
	    // chargeADCLog = log(chargeADC); // calculating log of charge
	    // chargeADCtest.push_back(chargeADC); // vector with summed ADC charge for each hit
	    // 	fChargeADCHist->Fill(chargeADCtest[h]); // filling the histogram of charge frequency
