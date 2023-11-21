#include "DAQQuickClustering_module.h"
#include "larcore/CoreUtils/ServiceUtil.h"

void ClusterHitsInTime::DoIt(std::vector<recoHit> cHitVector)
{
  fHitVector.clear();
  fHitVector = cHitVector;

  //ORDER IN THE HIT VECTOR BY TIME.
  std::sort(fHitVector.begin(), fHitVector.end(),
            [](const recoHit& lhs, const recoHit& rhs){return lhs<rhs;});

  for(unsigned int i = 0; i < fHitVector.size()-1; i++)
  {
    //std::cout << "fHitVector.at(" << i << ").getHitTime() " << fHitVector.at(i).getHitTime()
    //         << std::endl;
    std::vector<recoHit> vec_TempHits;
    if(std::abs(fHitVector.at(i).getHitTime()-fHitVector.at(i+1).getHitTime())<=fTimeWindow)
    {
      int  timeCount = 1;

      vec_TempHits.push_back(fHitVector.at(i));
      vec_TempHits.push_back(fHitVector.at(i+1));

      while((i+timeCount+1)<fHitVector.size() &&
            std::abs(fHitVector.at(i+timeCount).getHitTime()-fHitVector.at(i+timeCount+1).getHitTime()) <= fTimeWindow)
      {
        vec_TempHits.push_back(fHitVector.at(i + timeCount + 1));
        timeCount++;
      }

      i = i + timeCount;
      cluster temp(fHitVector.at(0).getEvent(), vec_TempHits);
      fVecClusters.push_back(temp);
    }
  }

  fNClusters = fVecClusters.size();
}

//......................................................
cluster::cluster(int cEvent, std::vector<recoHit> cHitVector)
{
  
  int type(0);
  std::set<int> channels;
  fHitVector.clear();
  fStartChan    = 50000;
  fEndChan      = 0;
  fFirstHitTime = 50000;
  fLastHitTime  = 0;
  fEvent        = cEvent;
  fHitVector    = cHitVector;
  fNHits        = fHitVector.size();
  for(int i = 0; i < fNHits; i++)
  {
    fHitSADC  += fHitVector.at(i).getHitSADC(); 
    int   chan = fHitVector.at(i).getHitChan();
    float time = fHitVector.at(i).getHitTime(); 

    channels.insert(chan);

    if(chan < fStartChan)
      fStartChan = chan;

    if(chan > fEndChan)
      fEndChan = chan;

    if(time < fFirstHitTime)
      fFirstHitTime = time;

    if(time > fLastHitTime)
      fLastHitTime = time;

    if(fHitVector.at(i).getGenType()==1)
      type++;
  }

  fNChan = channels.size();

  //CALL THE CLUSTER MARLEY IF THERE ARE MORE THAN TWO MARLEY HITS IN IT.
  if(type>=2)
    fType = 1;
  else
    fType = 0;

  fChanWidth = fEndChan - fStartChan;
  fTimeWidth = fLastHitTime - fFirstHitTime;

  //ORDER IN THE HIT VECTOR BY TIME.
  std::sort(fHitVector.begin(), fHitVector.end(),
            [](const recoHit& lhs, const recoHit& rhs){return lhs.getHitTime() < rhs.getHitTime();});
  
}


//......................................................
DAQQuickClustering::DAQQuickClustering(fhicl::ParameterSet const & p):EDAnalyzer(p){

   this->reconfigure(p);

}


//......................................................
void DAQQuickClustering::reconfigure(fhicl::ParameterSet const & p)
{

  fRawDigitLabel = p.get<std::string>("RawDigitLabel");
  fHitLabel      = p.get<std::string>("HitLabel"     );
  
  fGEANTLabel    = p.get<std::string>("GEANT4Label"  );
  fMARLLabel     = p.get<std::string>("MARLEYLabel"  );
  fAPALabel      = p.get<std::string>("APALabel"     );
  fCPALabel      = p.get<std::string>("CPALabel"     );
  fAr39Label     = p.get<std::string>("Argon39Label" );
  fNeutLabel     = p.get<std::string>("NeutronLabel" );
  fKrypLabel     = p.get<std::string>("KryptonLabel" );
  fPlonLabel     = p.get<std::string>("PoloniumLabel");
  fRdonLabel     = p.get<std::string>("RadonLabel"   );
  fAr42Label     = p.get<std::string>("Argon42Label" );
                                      

  cut_AdjChanTolerance = p.get<std::vector<int>>  ("AdjChanTolerance");
  cut_HitsInWindow     = p.get<std::vector<int>>  ("HitsInWindow");
  cut_MinChannels      = p.get<std::vector<int>>  ("MinChannels");
  cut_MinChanWidth     = p.get<std::vector<int>>  ("MinChanWidth");
  cut_TimeWindowSize   = p.get<std::vector<float>>("TimeWindowSize");
  cut_TotalADC         = p.get<std::vector<float>>("TotalADC");

  NConfigs = cut_AdjChanTolerance.size();
  NCuts    = 6;

  detectorScaling = p.get<double>("detectorScaling");

} // Reconfigure


//......................................................
void DAQQuickClustering::ResetVariables()
{
  MarlParts.clear(); APAParts .clear(); CPAParts .clear(); Ar39Parts.clear();
  NeutParts.clear(); KrypParts.clear(); PlonParts.clear(); RdonParts.clear();
  Ar42Parts.clear();
  trkIDToPType.clear();
  Run = SubRun = Event = -1;

  TotGen_Marl = TotGen_APA  = TotGen_CPA  = TotGen_Ar39 = 0;
  TotGen_Neut = TotGen_Kryp = TotGen_Plon = TotGen_Rdon = 0;
  TotGen_Ar42 = 0;

  out_HitView.clear();
  out_GenType.clear();
  out_HitChan.clear();
  out_HitTime.clear();
  out_HitSADC.clear();
  out_HitRMS .clear(); 
   
  NTotHits = NColHits = NIndHits = 0; 
  nHitsNotBackTracked = 0;
  for (int hh=0; hh<nMaxHits; ++hh)
  {
    HitView[hh] = HitSize[hh] = HitChan[hh] = GenType[hh] = 0;
    HitTime[hh] = HitRMS [hh] = HitSADC[hh] = 0;
    HitInt [hh] = HitPeak[hh] = HitTPC [hh] = 0;
    NCorrespondingIDEs[hh] = 0;
    Hit_X[hh] = 0;                  
    Hit_Y[hh] = 0;                  
    Hit_Z[hh] = 0;                  
    Hit_Energy[hh] = 0;             
    Hit_NumElectrons[hh] = 0;
  } 

  MarlSample.clear();
  MarlTime  .clear();
  MarlWeight.clear();
  ENu = 0;
  ENu_Lep = 0;
  VertX = 0;
  VertY = 0;
  VertZ = 0;
  VertexT = 0;
  Px    = 0;
  Py    = 0;
  Pz    = 0;
  VertexChan = 0;
  map_EventToMC.clear();
   

}


//......................................................
void DAQQuickClustering::clusterChannels(std::vector<recoHit> &vec_Hits,
                                         std::vector<ClusterHitsInTime> &vec_ChannelCluster,
                                         unsigned int const &config)
{
  //HERE IT IS ASSUMED THAT THE HITS APPEAR SEQUENTIALLY BY CHANNEL.

  for(unsigned int i = 0; i < vec_Hits.size()-1; i++)
  {
    std::vector<recoHit> vec_TempHits;
    double width = std::abs(vec_Hits.at(i  ).getHitChan()-
                            vec_Hits.at(i+1).getHitChan());
    //std::cout << cut_AdjChanTolerance.at(config) << std::endl;
    if(width <= cut_AdjChanTolerance.at(config))
    {
      int  channelCount = 1;
      
      vec_TempHits.push_back(vec_Hits.at(i  ));
      vec_TempHits.push_back(vec_Hits.at(i+1));

      while((i+channelCount+1)<vec_Hits.size() &&
            std::abs(vec_Hits.at(i+channelCount).getHitChan()-vec_Hits.at(i+channelCount+1).getHitChan())<=cut_AdjChanTolerance.at(config))
      {
        vec_TempHits.push_back(vec_Hits.at(i + channelCount + 1));
        channelCount++;
      }

      i = i + channelCount;
      ClusterHitsInTime temp(cut_TimeWindowSize[config]);
      temp.DoIt(vec_TempHits);
      
      vec_ChannelCluster.push_back(temp);
      //std::cout << " vec_ChannelCluster.size() " << vec_ChannelCluster.size() << std::endl;
            
    }
  }

  return;
}


//......................................................
void DAQQuickClustering::clusterCut(std::vector<cluster> &vec_Clusters, unsigned int const &config)
{
  //REMEMBER WE NEED BOTH A MAXIMUM AND MINIMUM CHANNEL WIDTH DUE TO COSMICS.

  for(std::vector<cluster>::iterator it_Clusters=vec_Clusters.begin();
      it_Clusters!=vec_Clusters.end();)
  {
    bool fail_MinChannels(false), fail_ChanWidth(false), fail_TotalADC(false);
    if(it_Clusters->getNChan()<cut_MinChannels.at(config))
      fail_MinChannels = true;

    if(it_Clusters->getChanWidth()<cut_MinChanWidth.at(config))
      fail_ChanWidth = true;

    if(it_Clusters->getHitSADC()<cut_TotalADC.at(config))
      fail_TotalADC = true;

    if(fail_MinChannels == true || fail_ChanWidth == true || fail_TotalADC == true)
      it_Clusters = vec_Clusters.erase(it_Clusters);
    else
      it_Clusters++;
  }
  return;
}


//......................................................
void DAQQuickClustering::trigger(std::vector<cluster> &vec_Clusters, unsigned int const &config)
{
  
  for(unsigned int i = 0; i < vec_Clusters.size(); i++)
    if(vec_Clusters.at(i).getNHits()>=cut_HitsInWindow.at(config))
      vec_Clusters.at(i).setTriggerFlag(1);

  return;
}


//......................................................
void DAQQuickClustering::makeConfigGraph()
{
  art::ServiceHandle<art::TFileService> tfs;

  for(unsigned int i = 0; i < NConfigs; i++)
  {
    TGraph* g_config = tfs->make<TGraph>(NCuts);
    g_config->SetName(Form("g_ConfigDefinitions%i",i));
    g_config->SetPoint(0, 1, cut_AdjChanTolerance[i]);
    g_config->SetPoint(1, 2, cut_HitsInWindow    [i]);
    g_config->SetPoint(2, 3, cut_MinChannels     [i]);
    g_config->SetPoint(3, 4, cut_MinChanWidth    [i]);
    g_config->SetPoint(4, 5, cut_TimeWindowSize  [i]);
    g_config->SetPoint(5, 6, cut_TotalADC        [i]);
  }
 
}


//......................................................
void DAQQuickClustering::beginJob()
{
     
  art::ServiceHandle<art::TFileService> tfs;

  makeConfigGraph();

  t_Output_unusedhits = tfs->make<TTree>("DAQQuickClustering_unusedhits",
                              "DAQQuickClustering_unusedhits");

  t_Output_unusedhits->Branch("Cluster",      &out_Cluster,     "Cluster/I");
  t_Output_unusedhits->Branch("Event",        &out_Event,       "Event/I");
  t_Output_unusedhits->Branch("Config",       &out_Config,      "Config/I");
  t_Output_unusedhits->Branch("StartChan",    &out_StartChan,   "StartChan/I");
  t_Output_unusedhits->Branch("EndChan",      &out_EndChan,     "EndChan/I");
  t_Output_unusedhits->Branch("ChanWidth",    &out_ChanWidth,   "ChanWidth/I");
  t_Output_unusedhits->Branch("NChan",        &out_NChan,       "NChan/I");
  t_Output_unusedhits->Branch("Type",         &out_Type,        "Type/I");
  t_Output_unusedhits->Branch("NHits",        &out_NHits,       "NHits/I");
  t_Output_unusedhits->Branch("SumADC",       &out_SumADC,      "SumADC/F");
  t_Output_unusedhits->Branch("FirstTimeHit", &out_FirstTimeHit,"FirstTimeHit/F");
  t_Output_unusedhits->Branch("LastTimeHit",  &out_LastTimeHit, "LastTimeHit/F");
  t_Output_unusedhits->Branch("TimeWidth",    &out_TimeWidth,   "TimeWidth/F");
  t_Output_unusedhits->Branch("ENu",          &out_ENu,         "ENu/D");
  t_Output_unusedhits->Branch("ENu_Lep",      &out_ENu_Lep,     "ENu_Lep/D");
  t_Output_unusedhits->Branch("MarlTime",     &out_MarlTime,    "MarlTime/D");
  t_Output_unusedhits->Branch("HitView",      &out_HitView);
  t_Output_unusedhits->Branch("GenType",      &out_GenType);
  t_Output_unusedhits->Branch("HitChan",      &out_HitChan);
  t_Output_unusedhits->Branch("HitTime",      &out_HitTime);
  t_Output_unusedhits->Branch("HitSADC",      &out_HitSADC);
  t_Output_unusedhits->Branch("HitRMS",       &out_HitRMS);

  t_Output_clusteredhits = tfs->make<TTree>("DAQQuickClustering_clusteredhits",
                                            "DAQQuickClustering_clusteredhits");

  t_Output_clusteredhits->Branch("Cluster",      &out_Cluster,     "Cluster/I");
  t_Output_clusteredhits->Branch("Event",        &out_Event,       "Event/I");
  t_Output_clusteredhits->Branch("Config",       &out_Config,      "Config/I");
  t_Output_clusteredhits->Branch("StartChan",    &out_StartChan,   "StartChan/I");
  t_Output_clusteredhits->Branch("EndChan",      &out_EndChan,     "EndChan/I");
  t_Output_clusteredhits->Branch("ChanWidth",    &out_ChanWidth,   "ChanWidth/I");
  t_Output_clusteredhits->Branch("NChan",        &out_NChan,       "NChan/I");
  t_Output_clusteredhits->Branch("Type",         &out_Type,        "Type/I");
  t_Output_clusteredhits->Branch("NHits",        &out_NHits,       "NHits/I");
  t_Output_clusteredhits->Branch("SumADC",       &out_SumADC,      "SumADC/F");
  t_Output_clusteredhits->Branch("FirstTimeHit", &out_FirstTimeHit,"FirstTimeHit/F");
  t_Output_clusteredhits->Branch("LastTimeHit",  &out_LastTimeHit, "LastTimeHit/F");
  t_Output_clusteredhits->Branch("TimeWidth",    &out_TimeWidth,   "TimeWidth/F");
  t_Output_clusteredhits->Branch("ENu",          &out_ENu,         "ENu/D");
  t_Output_clusteredhits->Branch("ENu_Lep",      &out_ENu_Lep,     "ENu_Lep/D");
  t_Output_clusteredhits->Branch("MarlTime",     &out_MarlTime,    "MarlTime/D");
  t_Output_clusteredhits->Branch("HitView",      &out_HitView);
  t_Output_clusteredhits->Branch("GenType",      &out_GenType);
  t_Output_clusteredhits->Branch("HitChan",      &out_HitChan);
  t_Output_clusteredhits->Branch("HitTime",      &out_HitTime);
  t_Output_clusteredhits->Branch("HitSADC",      &out_HitSADC);
  t_Output_clusteredhits->Branch("HitRMS",       &out_HitRMS);
  h_ENu_MC      = tfs->make<TH1D>("h_ENu_MC","h_ENu_MC",35,0,50);
  h_MarlTime_MC = tfs->make<TH1D>("h_MarlTime_MC","h_MarlTime_MC",100,-0.1,10.5);
  h_TimeElapsed = tfs->make<TH1D>("h_TimeElapsed", "h_TimeElapsed", 50,0,0.5);
  h_ENu_MC->Sumw2();
  std::cout << "Finished beginJob" << std::endl;
}


PType DAQQuickClustering::WhichParType(int TrID)
{
  PType ThisPType = kUnknown;
  auto const& it=trkIDToPType.find(TrID);
  if(it!=trkIDToPType.end()){
    ThisPType=it->second;
  }
  
  return ThisPType;
}


bool DAQQuickClustering::InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap )
{
     
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find(TrID);

  if(ParIt != ParMap.end()) return true;
  else                       return false;

}



void DAQQuickClustering::endJob()
{
  std::cout << "Job ended." << std::endl; 
  std::cerr << "firstCatch  " << firstCatch  << std::endl; 
  std::cerr << "secondCatch " << secondCatch << std::endl; 
  std::cerr << "thirdCatch  " << thirdCatch  << std::endl; 
}


void DAQQuickClustering::FillMyMaps(std::map<int, simb::MCParticle> &MyMap, 
                                    art::FindManyP<simb::MCParticle> Assn,
                                    art::ValidHandle< std::vector<simb::MCTruth> > Hand)
{
  for ( size_t L1=0; L1 < Hand->size(); ++L1 ) {
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
      MyMap[ThisPar.TrackId()] = ThisPar;
    }
  }
  return;
}

//......................................................
void DAQQuickClustering::analyze(art::Event const & evt)
{
  ResetVariables();
  
  Run    = evt.run();
  SubRun = evt.subRun();
  Event  = evt.event();

  //GET INFORMATION ABOUT THE DETECTOR'S GEOMETRY.
  auto const* geo = lar::providerFrom<geo::Geometry>();

// GET THE RECO HITS.
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(fHitLabel);
  
  //LIFT OUT THE MARLEY PARTICLES.
  auto MarlTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fMARLLabel);
  art::FindManyP<simb::MCParticle> MarlAssn(MarlTrue,evt,fGEANTLabel);
  FillMyMaps(MarlParts, MarlAssn, MarlTrue);
  TotGen_Marl = MarlParts.size();

  //SUPERNOVA TRUTH.
  art::FindManyP<sim::SupernovaTruth> SNTruth(MarlTrue, evt, fMARLLabel);

  double Px_(0), Py_(0), Pz_(0), Pnorm(1);
  for(unsigned int i = 0; i < MarlTrue->size(); i++)
  {
    Nu_Type     = MarlTrue->at(i).GetNeutrino().Nu().PdgCode();
    ENu         = MarlTrue->at(i).GetNeutrino().Nu().E();
    Mode        = MarlTrue->at(i).GetNeutrino().Mode();
    CCNC        = MarlTrue->at(i).GetNeutrino().CCNC();
    Target      = MarlTrue->at(i).GetNeutrino().Target();
    HitNucleon  = MarlTrue->at(i).GetNeutrino().HitNuc();
    Nu_Lep_Type = MarlTrue->at(i).GetNeutrino().Lepton().PdgCode(); 
    VertX       = MarlTrue->at(i).GetNeutrino().Lepton().Vx(); 
    VertY       = MarlTrue->at(i).GetNeutrino().Lepton().Vy(); 
    VertZ       = MarlTrue->at(i).GetNeutrino().Lepton().Vz(); 
    Px_         = MarlTrue->at(i).GetNeutrino().Lepton().Px();
    Py_         = MarlTrue->at(i).GetNeutrino().Lepton().Py();
    Pz_         = MarlTrue->at(i).GetNeutrino().Lepton().Pz();
    ENu_Lep     = MarlTrue->at(i).GetNeutrino().Lepton().E();
    for (unsigned int j = 0; j < SNTruth.at(i).size(); j++) 
    {
      const sim::SupernovaTruth ThisTr = (*SNTruth.at(i).at(j));
      MarlTime  .push_back(ThisTr.SupernovaTime);
      MarlWeight.push_back(ThisTr.Weight);
      MarlSample.push_back(ThisTr.SamplingMode);
    }
  }
  Pnorm = std::sqrt(Px_*Px_+Py_*Py_+Pz_*Pz_);
  Px = Px_/Pnorm;
  Py = Py_/Pnorm;
  Pz = Pz_/Pnorm;
  bool caught = false;
  geo::Point_t const Vertex{VertX, VertY, VertZ};
  geo::WireID WireID;
  geo::PlaneID Plane(geo->FindTPCAtPosition(Vertex),geo::kZ);
  try
  {
    WireID = geo->NearestWireID(Vertex, Plane);
  }
  catch(...)
  {
    caught = true;
  }

  if(caught)
    VertexChan = -1; 
  else
    VertexChan = geo->PlaneWireToChannel(WireID);
  
  //CM/MICROSECOND.
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
  double drift_velocity = detProp.DriftVelocity(detProp.Efield(),detProp.Temperature());
  //CM/TICK
  drift_velocity = drift_velocity*0.5;
  VertexT = VertX/drift_velocity;

  auto APATrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fAPALabel);
  art::FindManyP<simb::MCParticle> APAAssn(APATrue,evt,fGEANTLabel);
  FillMyMaps( APAParts, APAAssn, APATrue );
  TotGen_APA = APAParts.size();

  auto CPATrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fCPALabel);
  art::FindManyP<simb::MCParticle> CPAAssn(CPATrue,evt,fGEANTLabel);
  FillMyMaps( CPAParts, CPAAssn, CPATrue );
  TotGen_CPA = CPAParts.size();

  auto Ar39True = evt.getValidHandle<std::vector<simb::MCTruth> >(fAr39Label);
  art::FindManyP<simb::MCParticle> Ar39Assn(Ar39True,evt,fGEANTLabel);
  FillMyMaps( Ar39Parts, Ar39Assn, Ar39True );
  TotGen_Ar39 = Ar39Parts.size();

  auto NeutTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fNeutLabel);
  art::FindManyP<simb::MCParticle> NeutAssn(NeutTrue,evt,fGEANTLabel);
  FillMyMaps( NeutParts, NeutAssn, NeutTrue );
  TotGen_Neut = NeutParts.size();

  auto KrypTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fKrypLabel);
  art::FindManyP<simb::MCParticle> KrypAssn(KrypTrue,evt,fGEANTLabel);
  FillMyMaps( KrypParts, KrypAssn, KrypTrue );
  TotGen_Kryp = KrypParts.size();

  auto PlonTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fPlonLabel);
  art::FindManyP<simb::MCParticle> PlonAssn(PlonTrue,evt,fGEANTLabel);
  FillMyMaps( PlonParts, PlonAssn, PlonTrue );
  TotGen_Plon = PlonParts.size();

  auto RdonTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fRdonLabel);
  art::FindManyP<simb::MCParticle> RdonAssn(RdonTrue,evt,fGEANTLabel);
  FillMyMaps( RdonParts, RdonAssn, RdonTrue );
  TotGen_Rdon = RdonParts.size();

  auto Ar42True = evt.getValidHandle<std::vector<simb::MCTruth> >(fAr42Label);
  art::FindManyP<simb::MCParticle> Ar42Assn(Ar42True,evt,fGEANTLabel);
  FillMyMaps( Ar42Parts, Ar42Assn, Ar42True );
  TotGen_Ar42 = Ar42Parts.size();
  
  std::map<PType, std::map< int, simb::MCParticle >&> PTypeToMap{
      { kMarl, MarlParts},
      {kAPA, APAParts},
      {kCPA, CPAParts},
      {kAr39, Ar39Parts},
      {kNeut, NeutParts},
      {kKryp, KrypParts},
      {kPlon, PlonParts},
      {kRdon, RdonParts},
      {kAr42, Ar42Parts}
  };

  for(auto const& it : PTypeToMap){
      const PType p=it.first;
      auto const& m=it.second;
      for(auto const& it2 : m){
          trkIDToPType.insert(std::make_pair(it2.first, p));
      }
  }
  
  //std::cout << "THE EVENTS NUMBER IS: " << Event << std::endl;

  std::vector< recob::Hit > ColHits_Marl;
  std::vector< recob::Hit > ColHits_CPA;
  std::vector< recob::Hit > ColHits_APA;
  std::vector< recob::Hit > ColHits_Ar39;
  std::vector< recob::Hit > ColHits_Neut;
  std::vector< recob::Hit > ColHits_Kryp;
  std::vector< recob::Hit > ColHits_Plon;
  std::vector< recob::Hit > ColHits_Rdon;
  std::vector< recob::Hit > ColHits_Oth;
  std::vector< recob::Hit > ColHits_Ar42;

  NTotHits = reco_hits->size();
  int colHitCount(0);
  int LoopHits = std::min( NTotHits, nMaxHits );
  //if(LoopHits != NTotHits)
  // std::cerr << "---- There are " << NTotHits
  //           << " hits in the event, but array is of size " 
  //           << nMaxHits << ", so looping over first " << LoopHits
  //           << " hits." << std::endl;

  for(int hit = 0; hit < LoopHits; ++hit) 
  {
    recob::Hit const& ThisHit = reco_hits->at(hit);  
    if(ThisHit.View() == geo::kU || ThisHit.View() == geo::kV) ++NIndHits;
    else                                                       ++NColHits;
  }

  for(int hit = 0; hit < LoopHits; ++hit) 
  {
    recob::Hit const& ThisHit = reco_hits->at(hit);  

    if (ThisHit.View() == 2) {
    std::vector< sim::TrackIDE > ThisHitIDE; 
    //GETTING HOLD OF THE SIM::IDEs.
    std::vector<const sim::IDE*> ThisSimIDE;
    try
    {
      ThisHitIDE = bt_serv->HitToTrackIDEs(clockData,  ThisHit );
    }
    catch(...)
    {
      //std::cout << "FIRST CATCH" << std::endl;
      firstCatch++;
      try
      {
        ThisSimIDE = bt_serv->HitToSimIDEs_Ps(clockData, ThisHit);
      }
      catch(...)
      {
        //std::cout << "SECOND CATCH" << std::endl; 
        secondCatch++;
        continue;
      }
      continue;
    }
    try
    {
      ThisSimIDE = bt_serv->HitToSimIDEs_Ps(clockData, ThisHit);
    }
    catch(...)
    {
      //std::cout << "THIRD CATCH" << std::endl; 
      thirdCatch++;
      continue;
    }
    
    HitView[colHitCount] = ThisHit.View();
    HitSize[colHitCount] = ThisHit.EndTick() - ThisHit.StartTick();
    HitTPC [colHitCount] = ThisHit.WireID().TPC;
    HitChan[colHitCount] = ThisHit.Channel();
    HitTime[colHitCount] = ThisHit.PeakTime();
    HitRMS [colHitCount] = ThisHit.RMS();
    HitSADC[colHitCount] = ThisHit.SummedADC();
    HitInt [colHitCount] = ThisHit.Integral();
    HitPeak[colHitCount] = ThisHit.PeakAmplitude();
    NCorrespondingIDEs[colHitCount] = ThisHitIDE.size();
    
    if(ThisHitIDE.size()==0)
      nHitsNotBackTracked++;

    //WHICH PARTICLE CONTRIBUTED MOST TO THE HIT.
    int    MainTrID = -1;
    double TopEFrac = -DBL_MAX;
    for (size_t ideL=0; ideL < ThisHitIDE.size(); ++ideL) 
    {
      if ( ThisHitIDE[ideL].energyFrac > TopEFrac ) 
      {
        TopEFrac = ThisHitIDE[ideL].energyFrac;
        MainTrID = ThisHitIDE[ideL].trackID;
      }
    }

    PType ThisPType      = WhichParType( MainTrID );
    GenType[colHitCount] = ThisPType;

    if(MainTrID == -1)
    {
      Hit_X[colHitCount]            = -1;
      Hit_Y[colHitCount]            = -1;
      Hit_Z[colHitCount]            = -1;
      Hit_Energy[colHitCount]       = -1;
      Hit_NumElectrons[colHitCount] = -1;
    }
    else
    {
      for(unsigned int i = 0; i < ThisSimIDE.size(); i++)
      {
        if(ThisSimIDE.at(i)->trackID==MainTrID)
        {
          Hit_X[colHitCount]            = ThisSimIDE.at(i)->x;
          Hit_Y[colHitCount]            = ThisSimIDE.at(i)->y;
          Hit_Z[colHitCount]            = ThisSimIDE.at(i)->z;
          Hit_Energy[colHitCount]       = ThisSimIDE.at(i)->energy;
          Hit_NumElectrons[colHitCount] = ThisSimIDE.at(i)->numElectrons;
          break;
        }
      }
    }
    switch (ThisPType){
    case kUnknown: { ColHits_Oth .push_back( ThisHit ); break; }
    case kMarl:    { ColHits_Marl.push_back( ThisHit ); break; }
    case kAPA:     { ColHits_APA .push_back( ThisHit ); break; }
    case kCPA:     { ColHits_CPA .push_back( ThisHit ); break; }
    case kAr39:    { ColHits_Ar39.push_back( ThisHit ); break; }
    case kNeut:    { ColHits_Neut.push_back( ThisHit ); break; }
    case kKryp:    { ColHits_Kryp.push_back( ThisHit ); break; }
    case kPlon:    { ColHits_Plon.push_back( ThisHit ); break; }
    case kRdon:    { ColHits_Rdon.push_back( ThisHit ); break; }
    case kAr42:    { ColHits_Ar42.push_back( ThisHit ); break; }
    default: break;
    }

    colHitCount++;
    }
  } 

  std::vector<int> vec_ClusterCount(NConfigs);
  
  h_ENu_MC     ->Fill(1000*ENu);
  h_MarlTime_MC->Fill(MarlTime.back());
  
  map_EventToMC[Event] = {ENu, ENu_Lep, MarlTime.back()};
  
  //MAKE RECOHIT OBJECTS EVENTWISE FROM THE TREE.
  std::vector<recoHit> vec_Hits;
  for(int j = 0; j < NColHits; j++)
  {
    recoHit hit(Event,      HitView[j], GenType[j],
                HitChan[j], HitTime[j], HitSADC[j], HitRMS[j]);
    vec_Hits.push_back(hit);
  }

  for(unsigned int j = 0; j < 1; j++)
  {
    std::vector<ClusterHitsInTime> vec_ChannelCluster;
    std::vector<cluster>           vec_Clusters;
    TStopwatch *timeElapsed = new TStopwatch();
    clusterChannels(vec_Hits, vec_ChannelCluster, j);

    for(unsigned int k = 0; k < vec_ChannelCluster.size(); k++)
    {
      std::vector<cluster> vec_Temp = vec_ChannelCluster.at(k).getClusterVector();
      //std::cout << " vec_temp.size() " << vec_Temp.size() << std::endl;
      for(unsigned int l = 0; l < vec_Temp.size(); l++)
    	vec_Clusters.push_back(vec_Temp.at(l));
    }
    clusterCut(vec_Clusters, j);
    trigger   (vec_Clusters, j);
    h_TimeElapsed->Fill(timeElapsed->RealTime()*1000);
    delete timeElapsed;
    timeElapsed = NULL;
    //std::cout << " vec_Cluster.size() " << vec_Clusters.size() << std::endl;

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
	out_ENu          = map_EventToMC[vec_Clusters.at(k).getEvent()].at(0);
	out_ENu_Lep      = map_EventToMC[vec_Clusters.at(k).getEvent()].at(1);
	out_MarlTime     = map_EventToMC[vec_Clusters.at(k).getEvent()].at(2);
	for(unsigned int l = 0; l < vec_Clusters.at(k).getHits().size(); l++)
	{
	  out_HitView.push_back(vec_Clusters.at(k).getHits().at(l).getHitView());  
	  out_GenType.push_back(vec_Clusters.at(k).getHits().at(l).getGenType());  
	  out_HitChan.push_back(vec_Clusters.at(k).getHits().at(l).getHitChan());  
	  out_HitTime.push_back(vec_Clusters.at(k).getHits().at(l).getHitTime());  
	  out_HitSADC.push_back(vec_Clusters.at(k).getHits().at(l).getHitSADC());  
	  out_HitRMS .push_back(vec_Clusters.at(k).getHits().at(l).getHitRMS ());  
	}
	t_Output_clusteredhits->Fill();
	out_HitView.clear(); out_GenType.clear(); out_HitChan.clear(); 
	out_HitTime.clear(); out_HitSADC.clear(); out_HitRMS.clear();

	vec_ClusterCount.at(j)++;
      }
    }
    vec_ChannelCluster.clear();
    vec_Clusters.clear();
  }
  vec_Hits.clear();
}

DEFINE_ART_MODULE(DAQQuickClustering)
