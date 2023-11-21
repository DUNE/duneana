#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"

#include <iostream>



//==================================================//
//
// Author:   Michael Baird
// Description: This macro generates a set of fake
//   input hisotgrams for the macro run_SNe_Tprofile_sim.C.
//
// NOTE: These histograms are NOT representiative of
//   reality! They are just some place-holder best
//   guesses on my part...
//
//==================================================//
void make_fake_input_histos() {

  // create output file...
  TFile *outfile = TFile::Open("fake_SNe_input.root","RECREATE");

  // book histos (signal and background)
  TH1F *hTprofSN  = new TH1F("hTprofSN", "Supernova Time Profile;time [sec];N evts",10000,0.0,10.0);
  TH1F *hTprofBKG = new TH1F("hTprofBKG","Background Time Profile;time [sec];N evts",1,0.0,1.0);

  int NbinsE  = 100;
  double minE = 0.0;
  double maxE = 50.0;

  TH2F *hEspecSN  = new TH2F("hEspecSN", "Supernova Energy Spectrum;time [sec];energy [MeV]",10000,0.0,10.0,NbinsE,minE,maxE);
  TH1F *hEspecBKG = new TH1F("hEspecBKG","Background Energy Spectrum;energy [MeV];fraction of events",NbinsE,minE,maxE);

  TH1F *hEffSN  = new TH1F("hEffSN", "Supernova Efficiency vs. Energy;energy [MeV];eff",100,0.0,50.0);
  TH1F *hEffBKG = new TH1F("hEffBKG","Background Efficiency vs. Energy;energy [MeV];eff",100,0.0,50.0);

  int NbinsNH  = 500;
  double minNH = 0.0;
  double maxNH = 500.0;

  TH2F *hNHitsSN  = new TH2F("hNHitsSN", "Supernova - number of reco hits vs. true energy;energy [MeV];NHits",NbinsE,minE,maxE,NbinsNH,minNH,maxNH);
  TH2F *hNHitsBKG = new TH2F("hNHitsBKG","Background - number of reco hits vs. true energy;energy [MeV];NHits",NbinsE,minE,maxE,NbinsNH,minNH,maxNH);

  int NbinsADC  = 4096;
  double minADC = 0.0;
  double maxADC = 4095.0;

  TH2F *hSumADCSN  = new TH2F("hSumADCSN", "Supernova - summed ADC vs. true energy;energy [MeV];Summed ADC",NbinsE,minE,maxE,NbinsADC,minADC,maxADC);
  TH2F *hSumADCBKG = new TH2F("hSumADCBKG","Background - summed ADC vs. true energy;energy [MeV];Summed ADC",NbinsE,minE,maxE,NbinsADC,minADC,maxADC);

  TH1F *hSplitTPprobSN  = new TH1F("hSplitTPprobSN", "Supernova - Probability of splitting a TP into two;energy [MeV];split probability",NbinsE,minE,maxE);
  TH1F *hSplitTPprobBKG = new TH1F("hSplitTPprobBKG","Background - Probability of splitting a TP into two;energy [MeV];split probability",NbinsE,minE,maxE);

  //
  // Fill histograms:
  //

  // fill SN time profile histo:
  //
  // NOTE: this histo was eye-balled from a plot seen in ch.5 of the DUNE CDR
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.000),0);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.001),2);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.002),3);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.003),5);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.004),10);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.005),15);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.006),25);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.007),30);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.008),25);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.009),15);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.010),5);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.011),4);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.012),3);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.013),2);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.014),1);
  hTprofSN->SetBinContent(hTprofSN->FindBin(0.015),0.5);
  for(int i = 17; i < hTprofSN->GetNbinsX(); ++i) {
    hTprofSN->SetBinContent(i,0.05);
  }

  // fill SN energy spectra:
  // NOTE: for now we just use a constant energy spectra over time
  for(int i = 1; i <= hEspecSN->GetNbinsX(); ++i) {
    for(int j = 1; j <= hEspecSN->GetNbinsY(); ++j) {
      hEspecSN->SetBinContent(i,j,TMath::Poisson((maxE-minE)/(double)NbinsE/2.0+(j-1)*(maxE-minE)/(double)NbinsE,15.0));
    }
  }
  // This histo is a PDF, so normalize each vertical strip to 1.0
  for(int i = 1; i <= hEspecSN->GetNbinsX(); ++i) {
    double colSum = 0.0;
    for(int j = 1; j <= hEspecSN->GetNbinsY(); ++j) {
      colSum += hEspecSN->GetBinContent(i,j);
    }
    for(int j = 1; j <= hEspecSN->GetNbinsY(); ++j) {
      hEspecSN->SetBinContent(i,j,hEspecSN->GetBinContent(i,j)/colSum);
    }
  }

  // fill SN efficiency histo:
  for(int i = 1; i <= hEffSN->GetNbinsX(); ++i) {
    double value = 1.0 - TMath::Exp(-1.0*hEffSN->GetBinCenter(i));
    hEffSN->SetBinContent(i,value);
  }

  // fill SN NHits histo:
  for(int i = 1; i <= hNHitsSN->GetNbinsX(); ++i) {
    for(int j = 1; j <= hNHitsSN->GetNbinsY(); ++j) {
      hNHitsSN->SetBinContent(i,j,TMath::Poisson((maxNH-minNH)/(double)NbinsNH/2.0+(j-1)*(maxNH-minNH)/(double)NbinsNH,(maxE-minE)/(double)NbinsE/2.0+(i-1)*(maxE-minE)/(double)NbinsE));
    }
  }
  // This histo is a PDF, so normalize each vertical strip to 1.0
  for(int i = 1; i <= hNHitsSN->GetNbinsX(); ++i) {
    double colSum = 0.0;
    for(int j = 1; j <= hNHitsSN->GetNbinsY(); ++j) {
      colSum += hNHitsSN->GetBinContent(i,j);
    }
    for(int j = 1; j <= hNHitsSN->GetNbinsY(); ++j) {
      hNHitsSN->SetBinContent(i,j,hNHitsSN->GetBinContent(i,j)/colSum);
    }
  }

  // fill SN SumADC histo:
  for(int i = 1; i <= hSumADCSN->GetNbinsX(); ++i) {
    for(int j = 1; j <= hSumADCSN->GetNbinsY(); ++j) {
      hSumADCSN->SetBinContent(i,j,TMath::Poisson((maxADC-minADC)/(double)NbinsADC/2.0+(j-1)*(maxADC-minADC)/(double)NbinsADC,10.0*((maxE-minE)/(double)NbinsE/2.0+(i-1)*(maxE-minE)/(double)NbinsE)));
    }
  }
  // This histo is a PDF, so normalize each vertical strip to 1.0
  for(int i = 1; i <= hSumADCSN->GetNbinsX(); ++i) {
    double colSum = 0.0;
    for(int j = 1; j <= hSumADCSN->GetNbinsY(); ++j) {
      colSum += hSumADCSN->GetBinContent(i,j);
    }
    for(int j = 1; j <= hSumADCSN->GetNbinsY(); ++j) {
      hSumADCSN->SetBinContent(i,j,hSumADCSN->GetBinContent(i,j)/colSum);
    }
  }

  // fill SN split TP histo:
  for(int i = 1; i <= hSplitTPprobSN->GetNbinsX(); ++i) {
    double value = 0.5*(double)i/(double)NbinsE;
    hSplitTPprobSN->SetBinContent(i,value);
  }






  


  // fill BKG time profile histo:
  // NOTE: This is just one bin since the background rate is assumed here to be flat in time. This rate is not realistic...
  hTprofBKG->SetBinContent(1,0.005);

  // fill BKG energy spectra: 
  for(int i = 1; i <= hEspecBKG->GetNbinsX(); ++i) {
    double value = TMath::Exp(-1.0*hEffSN->GetBinCenter(i));
    hEspecBKG->SetBinContent(i,value);
  }
  // This histo is a PDF, so area normalize to 1.0.
  hEspecBKG->Scale(1.0/hEspecBKG->Integral());


  // fill BKG efficiency histo:
  for(int i = 1; i <= hEffBKG->GetNbinsX(); ++i) {
    double value = 1.0 - TMath::Exp(-1.0*hEffBKG->GetBinCenter(i));
    hEffBKG->SetBinContent(i,value);
  }

  // fill BKG NHits histo:
  for(int i = 1; i <= hNHitsBKG->GetNbinsX(); ++i) {
    for(int j = 1; j <= hNHitsBKG->GetNbinsY(); ++j) {
      hNHitsBKG->SetBinContent(i,j,TMath::Poisson((maxNH-minNH)/(double)NbinsNH/2.0+(j-1)*(maxNH-minNH)/(double)NbinsNH,(maxE-minE)/(double)NbinsE/2.0+(i-1)*(maxE-minE)/(double)NbinsE));
    }
  }
  // This histo is a PDF, so normalize each vertical strip to 1.0
  for(int i = 1; i <= hNHitsBKG->GetNbinsX(); ++i) {
    double colSum = 0.0;
    for(int j = 1; j <= hNHitsBKG->GetNbinsY(); ++j) {
      colSum += hNHitsBKG->GetBinContent(i,j);
    }
    for(int j = 1; j <= hNHitsBKG->GetNbinsY(); ++j) {
      hNHitsBKG->SetBinContent(i,j,hNHitsBKG->GetBinContent(i,j)/colSum);
    }
  }

  // fill BKG SumADC histo:
  for(int i = 1; i <= hSumADCBKG->GetNbinsX(); ++i) {
    for(int j = 1; j <= hSumADCBKG->GetNbinsY(); ++j) {
      hSumADCBKG->SetBinContent(i,j,TMath::Poisson((maxADC-minADC)/(double)NbinsADC/2.0+(j-1)*(maxADC-minADC)/(double)NbinsADC,10.0*((maxE-minE)/(double)NbinsE/2.0+(i-1)*(maxE-minE)/(double)NbinsE)));
    }
  }
  // This histo is a PDF, so normalize each vertical strip to 1.0
  for(int i = 1; i <= hSumADCBKG->GetNbinsX(); ++i) {
    double colSum = 0.0;
    for(int j = 1; j <= hSumADCBKG->GetNbinsY(); ++j) {
      colSum += hSumADCBKG->GetBinContent(i,j);
    }
    for(int j = 1; j <= hSumADCBKG->GetNbinsY(); ++j) {
      hSumADCBKG->SetBinContent(i,j,hSumADCBKG->GetBinContent(i,j)/colSum);
    }
  }

  // fill BKG split TP histo:
  for(int i = 1; i <= hSplitTPprobBKG->GetNbinsX(); ++i) {
    double value = 0.5*(double)i/(double)NbinsE;
    hSplitTPprobBKG->SetBinContent(i,value);
  }






  // write all histos to the output file:
  outfile->Write();

}
