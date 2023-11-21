#include <iostream>
#include <algorithm>
#include <map>
#include <TStopwatch.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TSystem.h>
#include "Module_SNClustering_Config.h"
#include "class_RecoHit.C"
#include "class_Cluster.C"
#include "class_ChannelCluster.C"


void clusterChannels(std::vector<recoHit> &vec_Hits, std::vector<channelCluster> &vec_ChannelCluster, unsigned int const &config)
{
  //HERE IT IS ASSUMED THAT THE HITS APPEAR SEQUENTIALLY BY CHANNEL.

  for(unsigned int i = 0; i < vec_Hits.size()-1; i++)
  {
    std::vector<recoHit> vec_TempHits;
    if(std::abs(vec_Hits.at(i).getHitChan()-vec_Hits.at(i+1).getHitChan())<=cut_AdjChanTolerance.at(config))
    {
      int  channelCount = 1;

      vec_TempHits.push_back(vec_Hits.at(i));
      vec_TempHits.push_back(vec_Hits.at(i+1));

      while((i+channelCount+1)<vec_Hits.size() 
          && std::abs(vec_Hits.at(i+channelCount).getHitChan()-vec_Hits.at(i+channelCount+1).getHitChan())
          <=cut_AdjChanTolerance.at(config))
      {
        vec_TempHits.push_back(vec_Hits.at(i + channelCount + 1));
        channelCount++;
      }

      i = i + channelCount;
      channelCluster temp(vec_TempHits.at(0).getEvent(), vec_TempHits, config);
      vec_ChannelCluster.push_back(temp);
    }
  }

  return;
}


void clusterCut(std::vector<cluster> &vec_Clusters, unsigned int const &config)
{
  //REMEMBER WE NEED BOTH A MAXIMUM AND MINIMUM CHANNEL WIDTH DUE TO COSMICS.

  for(std::vector<cluster>::iterator it_Clusters=vec_Clusters.begin(); it_Clusters!=vec_Clusters.end();)
  {
    bool fail_MinChannels(false), fail_ChanWidth(false), fail_TotalADC(false);
    if(it_Clusters->getNChan()<cut_MinChannels.at(config))
    {
      fail_MinChannels = true;
    }
    if(it_Clusters->getChanWidth()<cut_MinChanWidth.at(config))
    {
      fail_ChanWidth = true;
    }
    if(it_Clusters->getHitSADC()<cut_TotalADC.at(config))
    {
      fail_TotalADC = true;
    }

    if(fail_MinChannels == true || fail_ChanWidth == true || fail_TotalADC == true)
    {
      it_Clusters = vec_Clusters.erase(it_Clusters);
    }
    else
    {
      it_Clusters++;
    }
  }
  return;
}


void trigger(std::vector<cluster> &vec_Clusters, unsigned int const &config)
{
  for(unsigned int i = 0; i < vec_Clusters.size(); i++)
  {
    if(vec_Clusters.at(i).getNHits()>=cut_HitsInWindow.at(config))
    {
      vec_Clusters.at(i).setTriggerFlag(1);
    }
  }

  return;
}


std::vector<TGraph*> makeConfigGraph()
{
  std::vector<TGraph*> vec_gr_Config;
  for(unsigned int i = 0; i < NConfigs; i++)
  {
    TGraph *gr = new TGraph(NCuts);
    gr->SetPoint(0, 1, cut_AdjChanTolerance.at(i));
    gr->SetPoint(1, 2, cut_HitsInWindow    .at(i));
    gr->SetPoint(2, 3, cut_MinChannels     .at(i));
    gr->SetPoint(3, 4, cut_MinChanWidth    .at(i));
    gr->SetPoint(4, 5, cut_TimeWindowSize  .at(i));
    gr->SetPoint(5, 6, cut_TotalADC        .at(i));
    gr->SetName(Form("g_ConfigDefinitions%i",i));
    vec_gr_Config.push_back(gr);
  }
 
  return vec_gr_Config;
}

int main()
{
  std::vector<TGraph*> vec_gr_Config = makeConfigGraph();
  TString s_FileName = "GH_SNMC";
  //TFile *f_Input = new TFile("/pnfs/dune/persistent/users/abooth/SNTrigger/GH_SNMC/snb_timedep_radio_dune10kt_1x2x6/"+s_FileName+".root", "READ");
  //TFile *f_Input = new TFile(s_FileName+".root", "READ");
  TFile *f_Input = new TFile("SNAna_hist.root", "READ");
  TTree *t_Input = (TTree*)f_Input->Get("snanagaushit/SNSimTree");

  TFile *f_Output = new TFile("Module_"+s_FileName+".root", "RECREATE");
  TTree *t_Output = new TTree("t_Output", "Output Clusters");

  int nMaxHits = 1000;
  std::vector<int> vec_ClusterCount(NConfigs);
  struct MCContainer{
    double fENu;
    double fENu_Lep;
    std::vector<double> *fMarlTime;
  };
  std::map<int,MCContainer> map_EventToMC;

  //INPUT VARIABLES
  int   Event;
  int   NColHits;
  int   GenType[nMaxHits];
  int   HitView[nMaxHits];
  int   HitChan[nMaxHits];
  float HitTime[nMaxHits];
  float HitSADC[nMaxHits];
  float HitRMS[nMaxHits];
  double ENu;
  double ENu_Lep;
  std::vector<double> *MarlTime = 0;

  //int nEvents = 5000;
  int nEvents = t_Input->GetEntries();
  TH1I *hNEvents = new TH1I("hNEvents", "hNEvents", 100,0,10e6);
  hNEvents->Fill(nEvents);
  hNEvents->Write();
  TH1I *hNConfigs = new TH1I("hNConfigs", "hNConfigs", 50,0,50);
  hNConfigs->Fill(NConfigs);
  hNConfigs->Write();

  t_Input->SetBranchAddress("Event",    &Event);
  t_Input->SetBranchAddress("NColHits", &NColHits);
  t_Input->SetBranchAddress("GenType",  &GenType);
  t_Input->SetBranchAddress("HitView",  &HitView);
  t_Input->SetBranchAddress("HitChan",  &HitChan);
  t_Input->SetBranchAddress("HitTime",  &HitTime);
  t_Input->SetBranchAddress("HitSADC",  &HitSADC);
  t_Input->SetBranchAddress("HitRMS",   &HitRMS);
  t_Input->SetBranchAddress("ENu",      &ENu);
  t_Input->SetBranchAddress("ENu_Lep",  &ENu_Lep);
  t_Input->SetBranchAddress("MarlTime", &MarlTime);

  //OUTPUT VARIABLES
  int out_Cluster;
  int out_Event;
  int out_Config;
  int out_StartChan;
  int out_EndChan;
  int out_ChanWidth;
  int out_NChan;
  int out_Type;
  int out_NHits;
  float out_SumADC;
  float out_FirstTimeHit;
  float out_LastTimeHit;
  float out_TimeWidth;
  double out_ENu;
  double out_ENu_Lep;
  std::vector<double> out_MarlTime;
  std::vector<int> out_HitView;
  std::vector<int> out_GenType;
  std::vector<int> out_HitChan;
  std::vector<double> out_HitTime;
  std::vector<double> out_HitSADC;
  std::vector<double> out_HitRMS;

  t_Output->Branch("Cluster",      &out_Cluster,     "Cluster/I");
  t_Output->Branch("Event",        &out_Event,       "Event/I");
  t_Output->Branch("Config",       &out_Config,      "Config/I");
  t_Output->Branch("StartChan",    &out_StartChan,   "StartChan/I");
  t_Output->Branch("EndChan",      &out_EndChan,     "EndChan/I");
  t_Output->Branch("ChanWidth",    &out_ChanWidth,   "ChanWidth/I");
  t_Output->Branch("NChan",        &out_NChan,       "NChan/I");
  t_Output->Branch("Type",         &out_Type,        "Type/I");
  t_Output->Branch("NHits",        &out_NHits,       "NHits/I");
  t_Output->Branch("SumADC",       &out_SumADC,      "SumADC/F");
  t_Output->Branch("FirstTimeHit", &out_FirstTimeHit,"FirstTimeHit/F");
  t_Output->Branch("LastTimeHit",  &out_LastTimeHit, "LastTimeHit/F");
  t_Output->Branch("TimeWidth",    &out_TimeWidth,   "TimeWidth/F");
  t_Output->Branch("ENu",          &out_ENu,         "ENu/D");
  t_Output->Branch("ENu_Lep",      &out_ENu_Lep,     "ENu_Lep/D");
  t_Output->Branch("MarlTime", &out_MarlTime);
  t_Output->Branch("HitView",  &out_HitView);
  t_Output->Branch("GenType",  &out_GenType);
  t_Output->Branch("HitChan",  &out_HitChan);
  t_Output->Branch("HitTime",  &out_HitTime);
  t_Output->Branch("HitSADC",  &out_HitSADC);
  t_Output->Branch("HitRMS",   &out_HitRMS);

  TH1D *h_ENu_MC = new TH1D("h_ENu_MC","h_ENu_MC",35,0,50);
  h_ENu_MC->Sumw2();
  TH1D* h_MarlTime_MC = new TH1D("h_MarlTime_MC","h_MarlTime_MC",100,-0.1,10.5);
  //TIME IS IN MILLISECONDS.
  TH1D *h_TimeElapsed = new TH1D("h_TimeElapsed", "h_TimeElapsed", 50,0,0.5);
  for(unsigned int i = 0; i < nEvents; i++)
  {
    if(i % 500 == 0)
    {
      std::cout << "WORKING ON EVENT: " << i << std::endl;
    }
    t_Input->GetEntry(i);

    h_ENu_MC->Fill(1000*ENu);
    for(unsigned int j = 0; j < MarlTime->size(); j++)
    {
      h_MarlTime_MC->Fill(MarlTime->at(j));
    }

    map_EventToMC[Event] = {ENu, ENu_Lep, MarlTime};

    //MAKE RECOHIT OBJECTS EVENTWISE FROM THE TREE.
    std::vector<recoHit> vec_Hits;
    for(unsigned int j = 0; j < NColHits; j++)
    {
      recoHit hit(Event, HitView[j], GenType[j], HitChan[j], HitTime[j], HitSADC[j], HitRMS[j]);
      vec_Hits.push_back(hit);
    }

    for(unsigned int j = 0; j < NConfigs; j++)
    {
      std::vector<channelCluster> vec_ChannelCluster;
      std::vector<cluster>        vec_Clusters;
      TStopwatch *timeElapsed = new TStopwatch();
      clusterChannels(vec_Hits, vec_ChannelCluster, j);
      for(unsigned int k = 0; k < vec_ChannelCluster.size(); k++)
      {
        std::vector<cluster> vec_Temp = vec_ChannelCluster.at(k).getClusterVector();
        for(unsigned int l = 0; l < vec_Temp.size(); l++)
        {
          vec_Clusters.push_back(vec_Temp.at(l));
        }
      }
      clusterCut(vec_Clusters, j);
      trigger(vec_Clusters, j);
      h_TimeElapsed->Fill(timeElapsed->RealTime()*1000);

      //FILL THE OUTPUT TREE.
      for(unsigned int k  = 0; k  < vec_Clusters.size(); k++)
      {
        if(vec_Clusters.at(k).getTriggerFlag()==1)
        {
          out_Config       = j;
          out_Cluster      = vec_ClusterCount.at(j);
          out_Event        = vec_Clusters.at(k).getEvent();
          out_StartChan    = vec_Clusters.at(k).getStartChan();
          out_EndChan      = vec_Clusters.at(k).getEndChan();
          out_ChanWidth    = vec_Clusters.at(k).getChanWidth();
          out_NChan        = vec_Clusters.at(k).getNChan();
          out_Type         = vec_Clusters.at(k).getType();
          out_NHits        = vec_Clusters.at(k).getNHits();
          out_SumADC       = vec_Clusters.at(k).getHitSADC();
          out_FirstTimeHit = vec_Clusters.at(k).getFirstTimeHit();
          out_LastTimeHit  = vec_Clusters.at(k).getLastTimeHit();
          out_TimeWidth    = vec_Clusters.at(k).getTimeWidth();
          out_ENu          = map_EventToMC[vec_Clusters.at(k).getEvent()].fENu;
          out_ENu_Lep      = map_EventToMC[vec_Clusters.at(k).getEvent()].fENu_Lep;
          for(unsigned int l = 0; l < map_EventToMC[vec_Clusters.at(k).getEvent()].fMarlTime->size(); l++)
          {
            out_MarlTime.push_back(map_EventToMC[vec_Clusters.at(k).getEvent()].fMarlTime->at(l));
          }
          for(unsigned int l = 0; l < vec_Clusters.at(k).getHits().size(); l++)
          {
            out_HitView.push_back(vec_Clusters.at(k).getHits().at(l).getHitView());  
            out_GenType.push_back(vec_Clusters.at(k).getHits().at(l).getGenType());  
            out_HitChan.push_back(vec_Clusters.at(k).getHits().at(l).getHitChan());  
            out_HitTime.push_back(vec_Clusters.at(k).getHits().at(l).getHitTime());  
            out_HitSADC.push_back(vec_Clusters.at(k).getHits().at(l).getHitSADC());  
            out_HitRMS .push_back(vec_Clusters.at(k).getHits().at(l).getHitRMS ());  
          }
          t_Output->Fill();
          out_HitView.clear(); out_GenType.clear(); out_HitChan.clear(); 
          out_HitTime.clear(); out_HitSADC.clear(); out_HitRMS.clear();

          vec_ClusterCount.at(j)++;
        }
      }
    }
  }

  t_Output->Write();
  h_ENu_MC->Write();
  h_MarlTime_MC->Write();
  h_TimeElapsed->Write();
  for(unsigned int i = 0; i < NConfigs; i++)
  {
    vec_gr_Config.at(i)->Write();
  }

  f_Output->Close();

  return 0;
}
