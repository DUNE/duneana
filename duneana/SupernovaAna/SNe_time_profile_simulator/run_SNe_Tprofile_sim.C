//==================================================//
//
// Supernova Time-Profile Simulator
// 
// code:     run_SNe_Tprofile_sim.C
// author:   Michael Baird
//           (m.baird@sussex.ac.uk, mbaird42@fnal.gov)
//
// Description: This macro runs the supernova time 
//   profile simulator.
//
// This code must be run in compiled mode!
// Example:
//          root -l rune_SNe_Tprofile_sim.C+"(10.0,100,0.01)"
//
//==================================================//

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom1.h"

#include <iostream>



// TODO:
//
// * make sure that all of the code is generalized to not make any assumptions about the sizes or binnings of histograms...
// * make the input file name one of the input parameters
// *** code in some protection against the case that the time profile histo bin size is something other than the dt used in the loop over time
// * put some double checks (for things like out of range values) into the findValueFrom2DHisto function
// * add some assert statements (about the size of the nu vectors etc...)
// * code in max number of tries for the findValue functions



// ========== //
//
// ASSUMPTIONS:
//
// The following are a list of assumptions built into this code. The input histograms
// MUST adhere to these assumptions.
//
// 1. The SNe time-profile histogram assumes the number of expected events for a
//    distance of 10 kpc. This histo is renormalized by a factor of 1/r^2 where r is
//    the distance to the supernova (one of the macro inputs.)
//
// 2. The time profile histogram needs to have a linear binning for the X-axis.
//
// 3. The energy spectra (for both SNe and backgrounds) are interpreted as PDFs.
//    Therefore they must be normalized to have an area of 1.0.
//
// 4. The efficiency histograms MUST have the same X-axis binning as the energy
//    spectra, and they (obviously) must not exceed 1.0 for all bins.
//
// ========== //





// ================================================== //
//
// function to keep/reject an x-value from a 1D histo...
// (intended to be used for efficiencies)

bool keepValueFromEff(TH1F *h, float xValue, TRandom1 *rand) {

  // NOTE: It is assumed here that the efficiency is always a number between 0 and 1.0.

  bool keep = false;

  double xHeight = rand->Uniform(1.0);
  double content = h->GetBinContent(h->FindBin(xValue));
  if(xHeight <= content) {
    keep = true;
  }

  return keep;

}



// ================================================== //
//
// function to pick out an x-value from a 1D histo...

float findValueFrom1DHisto(TH1F *h, TRandom1 *rand) {

  float value = 0.0;

  // extract info about the 1D histo:
  double xMin = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst());
  double xMax = h->GetXaxis()->GetBinUpEdge (h->GetXaxis()->GetLast());
  double dX   = xMax-xMin;

  // find the max bin height:
  double maxHeightX = -1.0;
  for(int i = 1; i < h->GetNbinsX(); ++i) {
    double content = h->GetBinContent(i);
    if(content > maxHeightX) maxHeightX = content;
  }

  // pick a value:
  bool done = false;
  while(!done) {
    double xValue  = xMin + rand->Uniform(dX);
    double xHeight = rand->Uniform(maxHeightX);
    double content = h->GetBinContent(h->FindBin(xValue));
    if(xHeight <= content) {
      value = xValue;
      done  = true;
    }
  }

  return value;

}



// ================================================== //
//
// function to pick out a y-value from a 2D histo...

float findValueFrom2DHisto(TH2F *h, float xValue, TRandom1 *rand) {

  float value = 0.0;

  // extract info about the 2D histo:
  int    xBin = h->GetXaxis()->FindBin(xValue);
  double yMin = h->GetYaxis()->GetBinLowEdge(h->GetYaxis()->GetFirst());
  double yMax = h->GetYaxis()->GetBinUpEdge (h->GetYaxis()->GetLast());
  double dY   = yMax-yMin;

  // find the max bin height for the xValue column:
  double maxHeightY = -1.0;
  for(int j = 1; j < h->GetNbinsY(); ++j) {
    double content = h->GetBinContent(h->GetBin(xBin,j));
    if(content > maxHeightY) maxHeightY = content;
  }

  // pick a value:
  bool done = false;
  while(!done) {
    double yValue  = yMin + rand->Uniform(dY);
    double yHeight = rand->Uniform(maxHeightY);
    double content = h->GetBinContent(h->FindBin(xValue,yValue));
    if(yHeight <= content) {
      value = yValue;
      done  = true;
    }
  }

  return value;

}



// ================================================== //
//
// function to clear out old entries from the TP lists

void updateTPLists(std::vector<double> &times, std::vector<int> &nhits,
		   std::vector<int> &adc, std::vector<int> &truth,
		   double time, double winsize) {

  // NOTE: This function assumes that the elements in the list are
  //       in time order (which they should be unles the main code
  //       is drastically changed...)

  bool         useErase = false; // use the erase feature?
  unsigned int NErase   = 0;     // number of elements to erase
  for(unsigned int i = 0; i < times.size(); ++i) {
    if(times[i] < (time - winsize)) {
      useErase = true;
      NErase++;
    }
  }

  if(useErase) {
    times.erase(times.begin(),times.begin()+NErase);
    nhits.erase(nhits.begin(),nhits.begin()+NErase);
    adc  .erase(adc  .begin(),adc  .begin()+NErase);
    truth.erase(truth.begin(),truth.begin()+NErase);
  }

} 



// ================================================== //
//
// TRIGGER FUNCTION
//
// Simple trigger: just count # of TPs

bool TRIGGERcountTPs(std::vector<double> times, unsigned int Ntps) {

  bool trigger = false;

  if(times.size() >= Ntps) trigger = true;

  return trigger;

}





// ================================================== //
//
// main guts of the code...

void run_SNe_Tprofile_sim(double SNdist, double Nexpt, double WinSize) {

  std::cout << "\n\n\n"
	    << "Running supernova time-profile simulator with input parameters:\n"
	    << "\tSNe distance:         " << SNdist  << " kpc\n"
	    << "\tN experiments:        " << Nexpt   << "\n"
	    << "\tsliding window size:  " << WinSize << " sec\n"
	    << "\n\n";



  // open input file:
  TFile *inFile = new TFile("fake_SNe_input.root");



  // get histograms:
  TH1F *hTprofSN  = (TH1F*)inFile->FindObjectAny("hTprofSN");
  TH1F *hTprofBKG = (TH1F*)inFile->FindObjectAny("hTprofBKG");

  TH2F *hEspecSN  = (TH2F*)inFile->FindObjectAny("hEspecSN");
  TH1F *hEspecBKG = (TH1F*)inFile->FindObjectAny("hEspecBKG");

  TH2F *hNHitsSN  = (TH2F*)inFile->FindObjectAny("hNHitsSN");
  TH2F *hNHitsBKG = (TH2F*)inFile->FindObjectAny("hNHitsBKG");

  TH2F *hSumADCSN  = (TH2F*)inFile->FindObjectAny("hSumADCSN");
  TH2F *hSumADCBKG = (TH2F*)inFile->FindObjectAny("hSumADCBKG");

  TH1F *hEffSN  = (TH1F*)inFile->FindObjectAny("hEffSN");
  TH1F *hEffBKG = (TH1F*)inFile->FindObjectAny("hEffBKG");
  
  TH1F *hSplitTPprobSN  = (TH1F*)inFile->FindObjectAny("hSplitTPprobSN");
  TH1F *hSplitTPprobBKG = (TH1F*)inFile->FindObjectAny("hSplitTPprobBKG");
  
  
  
  // renormalize the SN Tprof histo to the input distance:
  // (assumes the input histogram was scalled to 10.0 kpc!)
  hTprofSN->Scale(100.0/SNdist/SNdist);
  
  
  
  // define a random number generator:
  TRandom1 *rand = new TRandom1();
  rand->SetSeed(0);



  // open an output file:
  TFile *outfile = TFile::Open("SNe_TProfile_output.root","RECREATE");

  // book some histos to be filled:
  TH1F *hDeltaT = new TH1F("hDeltaT","#Delta T btwn. Supernova start and trigger issued;#Delta T [sec];count",400,-0.2,0.2);
  
  
  // Generate some fake SNe...
  for(unsigned int n = 1; n <= Nexpt; ++n) {
    
    // only let a trigger be issued once per experiment
    bool trigIssued = false;
    
    // Generate a random start time for the supernova between 5-10 seconds into the expt.
    double SNstartTime = 5.0 + rand->Uniform(5.0);
    
    // background rate is assumed to be flat in time (thus the 1 bin histo)
    double NexpectedBKG = hTprofBKG->GetBinContent(1);
    
    
    
    // Define some vectors to keep running tallies of the TPs during the duration of the experiment
    std::vector<double> tpTimeList;
    std::vector<int>    tpNHitsList;
    std::vector<int>    tpSumADCList;
    std::vector<int>    tpTruthList;
    
    
    
    // Loop in time:  go for 20 seconds in 0.001 second ticks
    for(double t = 0.0; t <= 20.0; t += 0.001) {
      
      // determine number of SNe nu events in this time bin:
      unsigned int NevtsSN = 0;
      if(t > SNstartTime) {
	NevtsSN = rand->Poisson(hTprofSN->GetBinContent(hTprofSN->FindBin(t - SNstartTime)));
      }
      
      
      
      //
      // generate the info about the neutrinos (true E, NHits, Summed ADC)
      //
      
      // define vectors to hold the info about the individual neutrinos
      std::vector<float> nuE;
      std::vector<int>   nuNHits;
      std::vector<int>   nuSumADC;
      
      // loop over the number of SN events
      for(unsigned int i = 0; i < NevtsSN; ++i) {
	
	// pick the energy, number of hits, and the summed ADC for each true neutrino
	float energy = findValueFrom2DHisto(hEspecSN,t-SNstartTime,rand);
	int   nhits  = (int)findValueFrom2DHisto(hNHitsSN,energy,rand);
	int   sumadc = (int)findValueFrom2DHisto(hSumADCSN,energy,rand);
	
	bool keep = keepValueFromEff(hEffSN,energy,rand);

	// push the info into vectors if it passed the overall efficiency selector
	if(keep) {
	  nuE     .push_back(energy);
	  nuNHits .push_back(nhits);
	  nuSumADC.push_back(sumadc);
	}

      }
      
      
      
      // estimate pile up:
      // (loop over the list of true neutrinos, decide if some should be merged together, make new lists...)
      //
      // TODO:  This is not currently implemented!!!
      
      
      
      //
      // generate trigger primitives (TPs):
      //
      
      // define some vectors to store info about the TPs:
      std::vector<double> tpTime;
      std::vector<int>    tpNHits;
      std::vector<int>    tpSumADC;
      std::vector<int>    tpTruth; // use 1 for SN and 0 for background
      
      // loop over the new list of true nu info, decide to make one or two TPs, make final vectors of TPs
      for(unsigned int i = 0; i < nuE.size(); ++i) {
	double prob = rand->Uniform(1.0);
	if(prob <= hSplitTPprobSN->GetBinContent(hSplitTPprobSN->FindBin(nuE[i]))) {
	  // split this nu into two TPs:
	  // split the NHits and ADC the same between the two TPs
	  double splitFrac = 0.1 + rand->Uniform(0.8); // don't let the split fraction for one be less than 10%
	  int nhit1 = (int)(nuNHits[i]*splitFrac);
	  int nhit2 = nuNHits[i] - nhit1;
	  int adc1  = (int)(nuSumADC[i]*splitFrac);
	  int adc2  = nuSumADC[i] - adc1;
	  
	  tpTime  .push_back(t);
	  tpNHits .push_back(nhit1);
	  tpSumADC.push_back(adc1);
	  tpTruth .push_back(1);
	  
	  tpTime  .push_back(t);
	  tpNHits .push_back(nhit2);
	  tpSumADC.push_back(adc2);
	  tpTruth .push_back(1);
	}
	else {
	  tpTime  .push_back(t);
	  tpNHits .push_back(nuNHits[i]);
	  tpSumADC.push_back(nuSumADC[i]);
	  tpTruth .push_back(1);
	}
      } // end loop over i
      
      
      
      
      
      // ========== //
      //
      // Repeat all of the above for the background events...
      //
      // ========== //
      
      // determine number of noise events in this time bin
      unsigned int NevtsBKG = 0;
      NevtsBKG = rand->Poisson(NexpectedBKG);
      
      
      //
      // generate the info about the background events (true E, NHits, Summed ADC)
      //
      
      // define vectors to hold the info about the individual neutrinos
      std::vector<float> bkgE;
      std::vector<int>   bkgNHits;
      std::vector<int>   bkgSumADC;
      
      // loop over the number of BKG events
      for(unsigned int i = 0; i < NevtsBKG; ++i) {
	
	// pick the energy, number of hits, and the summed ADC for each true neutrino
	float energy = findValueFrom1DHisto(hEspecBKG,rand);
	int   nhits  = (int)findValueFrom2DHisto(hNHitsBKG,energy,rand);
	int   sumadc = (int)findValueFrom2DHisto(hSumADCBKG,energy,rand);

	bool keep = keepValueFromEff(hEffBKG,energy,rand);
	
	// push the info into vectors if it passes the overall efficiency selector
	if(keep) {
	  bkgE     .push_back(energy);
	  bkgNHits .push_back(nhits);
	  bkgSumADC.push_back(sumadc);
	}
      }
      
      
      
      // estimate pile up:
      // (loop over the list of true neutrinos, decide if some should be merged together, make new lists...)
      //
      // TODO:  This is not currently implemented!!!
      
      
      
      //
      // generate trigger primitives (TPs):
      //
      
      // loop over the new list of true nu info, decide to make one or two TPs, make final vectors of TPs
      for(unsigned int i = 0; i < bkgE.size(); ++i) {
	double prob = rand->Uniform(1.0);
	if(prob <= hSplitTPprobBKG->GetBinContent(hSplitTPprobBKG->FindBin(bkgE[i]))) {
	  // split this nu into two TPs:
	  // split the NHits and ADC the same between the two TPs
	  double splitFrac = 0.1 + rand->Uniform(0.8); // don't let the split fraction for one be less than 10%
	  int nhit1 = (int)(bkgNHits[i]*splitFrac);
	  int nhit2 = bkgNHits[i] - nhit1;
	  int adc1  = (int)(bkgSumADC[i]*splitFrac);
	  int adc2  = bkgSumADC[i] - adc1;
	  
	  tpTime  .push_back(t);
	  tpNHits .push_back(nhit1);
	  tpSumADC.push_back(adc1);
	  tpTruth .push_back(0);
	  
	  tpTime  .push_back(t);
	  tpNHits .push_back(nhit2);
	  tpSumADC.push_back(adc2);
	  tpTruth .push_back(0);
	}
	else {
	  tpTime  .push_back(t);
	  tpNHits .push_back(bkgNHits[i]);
	  tpSumADC.push_back(bkgSumADC[i]);
	  tpTruth .push_back(0);
	}
      } // end loop over i
      
      
      
      
      
      // Update the master list of TPs:
      //   * add new TPs to the list
      //   * wipe elements older than the look back time
      for(unsigned int i = 0; i < tpTime.size(); ++i) {
	tpTimeList  .push_back(tpTime[i]);
	tpNHitsList .push_back(tpNHits[i]);
	tpSumADCList.push_back(tpSumADC[i]);
	tpTruthList .push_back(tpTruth[i]);
      }
      updateTPLists(tpTimeList, tpNHitsList, tpSumADCList, tpTruthList, t, WinSize);
      
      
      
      // Call the trigger function
      bool trigger = false;
      if(!trigIssued) {
	trigger = TRIGGERcountTPs(tpTimeList,10);
	if(trigger) {
	  trigIssued = true;
	  /*
	  std::cout << "\n\n\nTrigger issued!\n"
		    << "at t = " << t << "\n"
		    << "SNe start time = " << SNstartTime
		    << "\n\n\n";
	  */
	  hDeltaT->Fill(t-SNstartTime);
	}
      }
      
      
      
    } // end loop over time t
    
    // Fill some plots with the experimental results
    
    
  } // end loop over n (number of experiments
  
  outfile->Write();

}
