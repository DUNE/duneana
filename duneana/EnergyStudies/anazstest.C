#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TPad.h"
#include "TAxis.h"
#include "TStyle.h"

//  This is a root script for analyzing the output of the
// dunezstest module -- it shows the energy resolution as a function
// of the zero suppression threshold after fitting Gaussians to
// the distribution of ADC sums for each value of the threshold
 
// 5E4 is the max for Babu's sample

//void anazstest(TString filebase="dunezsanalysis_kevin30mev", double hmax=5E5)
//void anazstest(TString filebase="dunezsanalysis_kevin30mev_2", double hmax=25E3)
void anazstest(TString filebase="dunezsanalysis", double hmax=5E4)
{
  TString filename;
  filename = filebase + ".root";
  TFile *f = new TFile(filename);
  TTree *tree1 = (TTree*) f->Get("dunezsanalysis/dunezsanalysisReconstruction");
  //f->Print();
  //tree1->Print();

  TH1F *hpt[200];
  for (int i=0;i<200;i++)
    {
      TString hname;
      hname = "chargesum";
      hname += i;
      hpt[i] = new TH1F(hname,hname,100,0,hmax);
    }

  double chargesum[200];
  tree1->SetBranchAddress("ChargeSum",chargesum);

  int nent = tree1->GetEntries();
  for (int i=0;i<nent;i++)
    {
      tree1->GetEntry(i);
      for (int j=0;j<200;j++)
	{
	  if (chargesum[j]<hmax)
	    {
	      hpt[j]->Fill(chargesum[j]);
	    }
	}
    }

  double hmean[200],hrms[200];
  double hcut[200];
  double hfr[200];

  for (int i=0;i<200;i++)
    {
      hmean[i] = hpt[i]->GetMean(); 
      hpt[i]->Fit("gaus");
      TF1 *myfunc = hpt[i]->GetFunction("gaus");
      double gwid = myfunc->GetParameter(2);
      //hrms[i] = hpt[i]->GetRMS();
      hrms[i] = gwid;
      hcut[i] = i;
      if (hmean[i]>0) { hfr[i] = hrms[i]/hmean[i]; }
      else {hfr[i] = 0;}
    }

  TGraph *gmean = new TGraph(200,hcut,hmean);
  gmean->GetXaxis()->SetTitle("Zero-Suppression Cut (ADC)");
  gmean->GetYaxis()->SetTitle("Mean ZPlane ADC Sum");
  gmean->SetTitle("");
  TGraph *grms = new TGraph(200,hcut,hrms);
  grms->GetXaxis()->SetTitle("Zero-Suppression Cut (ADC)");
  grms->GetYaxis()->SetTitle("RMS ZPlane ADC Sum");
  grms->SetTitle("");
  TGraph *gfrac = new TGraph(200,hcut,hfr);
  gfrac->GetXaxis()->SetTitle("Zero-Suppression Cut (ADC)");
  gfrac->GetYaxis()->SetTitle("RMS ZPlane ADC Sum/Mean");
  gfrac->SetTitle("");

  TCanvas *c1 = new TCanvas("c1");
  gStyle->SetOptStat(0);
  c1->Divide(2,2);
  c1->cd(1);
  gmean->Draw("AC");
  c1->cd(2);
  grms->Draw("AC");
  c1->cd(3);
  gfrac->SetMinimum(0);
  gfrac->Draw("AC");
  c1->cd(4);
  gfrac->SetMinimum(0.03);
  gPad->SetLogy(1);
  gfrac->Draw("AC");
  TString outpic1;
  outpic1 = "zstest1_" + filebase + ".eps";
  c1->Print(outpic1);

  TCanvas *c2 = new TCanvas("c2");
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);
  c2->Divide(2,2);
  gPad->SetLogy(0);
  c2->cd(1);
  hpt[0]->Draw();
  c2->cd(2);
  hpt[50]->Draw();
  c2->cd(3);
  hpt[100]->Draw();
  c2->cd(4);
  hpt[150]->Draw();
  TString outpic1;
  outpic2 = "chargedistribs_" + filebase + ".eps";
  c2->Print(outpic2);
}
