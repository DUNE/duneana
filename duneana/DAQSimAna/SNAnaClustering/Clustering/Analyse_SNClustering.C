#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TString.h>
#include <TROOT.h>
#include "Module_SNClustering_Config.h"


std::map<int,double> map_TypeToWeight;

double getWeight(int const &ClusterType, TString const &s_FileName)
{
  double weight = map_TypeToWeight[ClusterType];

  return weight;
}


int getClusterType(std::vector<int>* vec_GenType)
{
  int type = -1;

  std::map<int,int> map_GenTypeToCount;
  for(unsigned int i = 0; i < vec_GenType->size(); i++)
  {
    map_GenTypeToCount[vec_GenType->at(i)]++;  
  }

  std::map<int,int>::iterator it_GenTypeToCount;
  int largestCount = 0;
  for(it_GenTypeToCount=map_GenTypeToCount.begin(); it_GenTypeToCount!=map_GenTypeToCount.end(); it_GenTypeToCount++)
  {
    if(it_GenTypeToCount->second>largestCount && it_GenTypeToCount->first!=1)
    {
      largestCount = it_GenTypeToCount->second;
      type         = it_GenTypeToCount->first; 
    }
  }

  return type;
}


int main()
{
  gROOT->ProcessLine("#include <vector>");
  TString s_FileName = "GH_SNMC";

  //TFile *f_Input = new TFile("/dune/data/users/abooth/Module_"+s_FileName+".root", "READ");
  TFile *f_Input = new TFile("Module_"+s_FileName+".root", "READ");
  TTree *t_Input = (TTree*)f_Input->Get("t_Output");
  int nEventsOriginally = (int)(((TH1I*)f_Input->Get("hNEvents"))->GetMean());
  int nConfigs          = (int)(((TH1I*)f_Input->Get("hNConfigs"))->GetMean());
  std::cout << "THERE WERE " << nEventsOriginally << " EVENTS IN THE ORIGINAL SAMPLE AND " << nConfigs << " CONFIGURATIONS" << std::endl;
 
  TFile *f_Output = new TFile("Analyse_"+s_FileName+".root", "RECREATE");

  map_TypeToWeight[0] = 1.; //NOISE 
  map_TypeToWeight[1] = 1.; //MARLEY
  map_TypeToWeight[2] = 1.; //CO60
  map_TypeToWeight[3] = 1.; //K40
  map_TypeToWeight[4] = 1.; //AR39
  map_TypeToWeight[5] = 1.; //n
  map_TypeToWeight[6] = 1.; //KR
  map_TypeToWeight[7] = 1.; //PO
  map_TypeToWeight[8] = 1.; //RN
  map_TypeToWeight[9] = 1.; //AR42

  std::ofstream txt_WeightMap;
  txt_WeightMap.open("Analyse_WeightMap_"+s_FileName+".txt");
  std::map<int,double>::iterator it_TypeToWeight;
  for(it_TypeToWeight=map_TypeToWeight.begin();it_TypeToWeight!=map_TypeToWeight.end();it_TypeToWeight++)
  {
    txt_WeightMap << it_TypeToWeight->first << " " << it_TypeToWeight->second << std::endl;  
  }

  int Cluster;
  int Event;
  int Config;
  int StartChan;
  int EndChan;
  int ChanWidth;
  int NChan;
  int Type;
  int NHits;
  float SumADC;
  float FirstTimeHit;
  float LastTimeHit;
  float TimeWidth;
  double ENu;
  double ENu_Lep;
  std::vector<int> *GenType = 0;

  long int nClusters = t_Input->GetEntries();
  t_Input->SetBranchAddress("Cluster",      &Cluster     );
  t_Input->SetBranchAddress("Event",        &Event       );
  t_Input->SetBranchAddress("Config",       &Config      );
  t_Input->SetBranchAddress("StartChan",    &StartChan   );
  t_Input->SetBranchAddress("EndChan",      &EndChan     );
  t_Input->SetBranchAddress("ChanWidth",    &ChanWidth   );
  t_Input->SetBranchAddress("NChan",        &NChan       );
  t_Input->SetBranchAddress("Type",         &Type        );
  t_Input->SetBranchAddress("NHits",        &NHits       );
  t_Input->SetBranchAddress("SumADC",       &SumADC      );
  t_Input->SetBranchAddress("FirstTimeHit", &FirstTimeHit);
  t_Input->SetBranchAddress("LastTimeHit",  &LastTimeHit );
  t_Input->SetBranchAddress("TimeWidth",    &TimeWidth   );
  t_Input->SetBranchAddress("ENu",          &ENu         );
  t_Input->SetBranchAddress("ENu_Lep",      &ENu_Lep     );
  t_Input->SetBranchAddress("GenType",      &GenType     );

  //DEFINE ANALYSIS HISTOGRAMS.
  std::vector<TH1D*> vec_h_EfficiencyVEnergy;
  std::vector<TH1D*> vec_h_BkgdVEnergy;
  std::vector<TH1I*> vec_h_RemainingBkgdGenType;
  std::vector<TH1I*> vec_h_NMarlClustersPerEvent;
  std::vector<TH1I*> vec_h_NMarlHit_Marl;
  std::vector<TH1I*> vec_h_NBkgdHit_Marl;
  std::vector<TH1I*> vec_h_NMarlHit_Bkgd;
  std::vector<TH1D*> vec_h_FracMarlHit_Marl;
  std::vector<TH1D*> vec_h_FracBkgdHit_Marl;
  std::vector<TH1D*> vec_h_FracMarlHit_Bkgd;
  std::vector<TH1I*> vec_h_ClusterGenType;
  for(unsigned int i = 0; i < nConfigs; i++)
  {
    TString s_Config = Form("%i", i);
    TH1D *h_EfficiencyVEnergy = new TH1D("h_EfficiencyVEnergy"+s_Config, "h_EfficiencyVEnergy"+s_Config, 35, 0, 50);
    h_EfficiencyVEnergy->Sumw2();
    vec_h_EfficiencyVEnergy.push_back(h_EfficiencyVEnergy);
    TH1D *h_BkgdVEnergy = new TH1D("h_BkgdVEnergy"+s_Config, "h_BkgdVEnergy"+s_Config, 30, 0, 1500);
    vec_h_BkgdVEnergy.push_back(h_BkgdVEnergy);
    TH1I *h_RemainingBkgdGenType = new TH1I("h_RemainingBkgdGenType"+s_Config,"h_RemainingBkgdGenType"+s_Config, 11,-0.5,10.5);
    vec_h_RemainingBkgdGenType.push_back(h_RemainingBkgdGenType);
    TH1I *h_NMarlClustersPerEvent = new TH1I("h_NMarlClustersPerEvent"+s_Config,"h_NMarlClustersPerEvent"+s_Config,7,-0.5,6.5);
    vec_h_NMarlClustersPerEvent.push_back(h_NMarlClustersPerEvent);
    TH1I *h_NMarlHit_Marl = new TH1I("h_NMarlHit_Marl"+s_Config,"h_NMarlHit_Marl"+s_Config,31,-0.5,30.5);
    vec_h_NMarlHit_Marl.push_back(h_NMarlHit_Marl);
    TH1I *h_NMarlHit_Bkgd = new TH1I("h_NMarlHit_Bkgd"+s_Config,"h_NMarlHit_Bkgd"+s_Config,3,-0.5,2.5);
    vec_h_NMarlHit_Bkgd.push_back(h_NMarlHit_Bkgd);
    TH1I *h_NBkgdHit_Marl = new TH1I("h_NBkgdHit_Marl"+s_Config,"h_NBgkdHit_Marl"+s_Config,5,-0.5,4.5);
    vec_h_NBkgdHit_Marl.push_back(h_NBkgdHit_Marl);
    TH1D *h_FracMarlHit_Marl = new TH1D("h_FracMarlHit_Marl"+s_Config,"h_FracMarlHit_Marl"+s_Config,11,-0.05,1.05);
    vec_h_FracMarlHit_Marl.push_back(h_FracMarlHit_Marl);
    TH1D *h_FracMarlHit_Bkgd = new TH1D("h_FracMarlHit_Bkgd"+s_Config,"h_FracMarlHit_Bkgd"+s_Config,11,-0.05,1.05);
    vec_h_FracMarlHit_Bkgd.push_back(h_FracMarlHit_Bkgd);
    TH1D *h_FracBkgdHit_Marl = new TH1D("h_FracBkgdHit_Marl"+s_Config,"h_FracBkgdHit_Marl"+s_Config,11,-0.05,1.05);
    vec_h_FracBkgdHit_Marl.push_back(h_FracBkgdHit_Marl);
    TH1I *h_ClusterGenType = new TH1I("h_ClusterGenType"+s_Config,"h_ClusterGenType"+s_Config, 11,-0.5,10.5);
    vec_h_ClusterGenType.push_back(h_ClusterGenType);
  }

  //OVERALL EFFICIENCIES AND BACKGROUND RATES.
  std::map<int,int> map_ConfigToBkgdCount;
  std::map<std::pair<int,int>,int> map_ConfigAndEventToNClusters;
  for(long int i = 0; i < nClusters; i++)
  {
    t_Input->GetEntry(i);

    if(Type == 1)
    {
      map_ConfigAndEventToNClusters[{Config,Event}]++;
      if(map_ConfigAndEventToNClusters[{Config,Event}]==1)
      {
        vec_h_EfficiencyVEnergy.at(Config)->Fill(ENu*1000);
      }
      int marlCount = 0;
      int bkgdCount = 0;
      for(int j = 0; j < NHits; j++)
      {
        if(GenType->at(j)==1)
        {
          marlCount++;
        }
        else
        {
          bkgdCount++;
        }
      }
      vec_h_NMarlHit_Marl.at(Config)->Fill(marlCount);
      vec_h_FracMarlHit_Marl.at(Config)->Fill((double)marlCount/(double)NHits);
      vec_h_NBkgdHit_Marl.at(Config)->Fill(bkgdCount);
      vec_h_FracBkgdHit_Marl.at(Config)->Fill((double)bkgdCount/(double)NHits);
    }
    else
    {
      int    ClusterType = getClusterType(GenType);
      double weight = getWeight(ClusterType, s_FileName);
      map_ConfigToBkgdCount[Config] += weight; 
      vec_h_BkgdVEnergy.at(Config)->Fill(SumADC);
      vec_h_ClusterGenType.at(Config)->Fill(ClusterType);
      int marlCount = 0;
      for(int j = 0; j < NHits; j++)
      {
        vec_h_RemainingBkgdGenType.at(Config)->Fill(GenType->at(j));
        if(GenType->at(j)==1)
        {
          marlCount++;
        }
      }
      vec_h_NMarlHit_Bkgd.at(Config)->Fill(marlCount);
      vec_h_FracMarlHit_Bkgd.at(Config)->Fill((double)marlCount/(double)NHits);
    }
  }

  std::map<int,std::pair<double,double>> map_ConfigToEfficiencyAndBkgd;
  std::map<std::pair<int,int>,int>::iterator it_ConfigAndEventToNClusters;
  for(int i = 0; i < nConfigs; i++)
  {
    double eff = 0;
    for(it_ConfigAndEventToNClusters=map_ConfigAndEventToNClusters.begin();
        it_ConfigAndEventToNClusters!=map_ConfigAndEventToNClusters.end(); it_ConfigAndEventToNClusters++)
    {
      if(it_ConfigAndEventToNClusters->first.first==i)
      {
        vec_h_NMarlClustersPerEvent.at(i)->Fill(it_ConfigAndEventToNClusters->second);
        eff++;
      }
    }
    eff/=(double)nEventsOriginally;
    double bkgdRate = map_ConfigToBkgdCount[i]/((double)(nEventsOriginally)*2.246e-3);
    
    map_ConfigToEfficiencyAndBkgd[i] = {eff,bkgdRate/detectorScaling};
    vec_h_BkgdVEnergy.at(i)->Scale(1/((double)(nEventsOriginally)*2.246e-3*detectorScaling));

    vec_h_EfficiencyVEnergy.at(i)->Write();
    vec_h_BkgdVEnergy.at(i)->Write();
    vec_h_RemainingBkgdGenType.at(i)->Write();
    vec_h_NMarlHit_Marl.at(i)->Write();
    vec_h_NBkgdHit_Marl.at(i)->Write();
    vec_h_NMarlHit_Bkgd.at(i)->Write();
    vec_h_FracMarlHit_Marl.at(i)->Write();
    vec_h_FracBkgdHit_Marl.at(i)->Write();
    vec_h_FracMarlHit_Bkgd.at(i)->Write();
    vec_h_ClusterGenType.at(i)->Write();

    int zeroBin = vec_h_NMarlClustersPerEvent.at(i)->FindBin(0);
    vec_h_NMarlClustersPerEvent.at(i)->SetBinContent(zeroBin, nEventsOriginally-vec_h_NMarlClustersPerEvent.at(i)->GetEntries());
    vec_h_NMarlClustersPerEvent.at(i)->Write();
  }

  TGraph *g_ROC = new TGraph(nConfigs);
  g_ROC->SetNameTitle("g_BkgdVsEff", "g_BkgdVsEff");
  std::map<int,std::pair<double,double>>::iterator it_ConfigToEfficiencyAndBkgd;
  int c = 0;
  std::ofstream txt_Result;
  txt_Result.open("Analyse_"+s_FileName+".txt");
  for(it_ConfigToEfficiencyAndBkgd=map_ConfigToEfficiencyAndBkgd.begin();
      it_ConfigToEfficiencyAndBkgd!=map_ConfigToEfficiencyAndBkgd.end(); it_ConfigToEfficiencyAndBkgd++)
  {
    std::cout  << "CONFIG "           << it_ConfigToEfficiencyAndBkgd->first 
               << ", SN EFFICIENCY: " << it_ConfigToEfficiencyAndBkgd->second.first 
               << " BKGD RATE IN 10kt (Hz): " << it_ConfigToEfficiencyAndBkgd->second.second << std::endl;
    g_ROC->SetPoint(c, it_ConfigToEfficiencyAndBkgd->second.first, it_ConfigToEfficiencyAndBkgd->second.second);
    txt_Result << it_ConfigToEfficiencyAndBkgd->first << " " << it_ConfigToEfficiencyAndBkgd->second.first 
                                                      << " " << it_ConfigToEfficiencyAndBkgd->second.second << std::endl; 
    c++;
  }

  g_ROC->Write();

  return 0;
}
