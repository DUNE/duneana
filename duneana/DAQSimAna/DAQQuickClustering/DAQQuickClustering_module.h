#ifndef DAQQUICKCLUSTERING_H
#define DAQQUICKCLUSTERING_H
// C++ includes

// ROOT includes
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TGraph.h"
#include "TStopwatch.h"

// Framework includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SupernovaTruth.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

enum PType{ kUnknown, kMarl, kAPA, kCPA, kAr39, kNeut, kKryp, kPlon, kRdon , kAr42};

const int nMaxHits=100000;

class recoHit
{
 public:
  recoHit(int   cEvent,   int   cHitView, int   cGenType, int   cHitChan,
          float cHitTime, float cHitSADC, float cHitRMS):fEvent  (cEvent  ),
    fHitView(cHitView),
    fGenType(cGenType),
    fHitChan(cHitChan),
    fHitSADC(cHitSADC),
    fHitRMS (cHitRMS ),
    fHitTime(cHitTime)
    {
    };
  friend bool operator<(recoHit lhs, recoHit rhs);

  int   getEvent  () const { return fEvent  ; };
  int   getHitView() const { return fHitView; };
  int   getGenType() const { return fGenType; };
  int   getHitChan() const { return fHitChan; };
  float getHitTime() const { return fHitTime; };
  float getHitSADC() const { return fHitSADC; };
  float getHitRMS () const { return fHitRMS;  };
  void  Print  () const
  {
    std::cout << "fHitView " << fHitView
              << ",  fGenType " << fGenType
              << ",  fHitChan " << fHitChan
              << ",  fHitSADC " << fHitSADC
              << ",  fHitRMS  " << fHitRMS 
              << ",  fHitTime " << fHitTime
              << std::endl;
  };

 private:
  int   fEvent   = 0;
  int   fHitView = 0;
  int   fGenType = 0;
  int   fHitChan = 0;
  float fHitSADC = 0;
  float fHitRMS  = 0;
  float fHitTime = 0;

};

bool operator<(recoHit lhs, recoHit rhs){
  if(lhs.fHitTime<rhs.fHitTime)
    return true;
  return false;
}

class cluster
{
 public:
  cluster(int cEvent, std::vector<recoHit> cHitVector);
  cluster();

  int    getEvent       () const { return fEvent       ; };
  int    getStartChan   () const { return fStartChan   ; };
  int    getEndChan     () const { return fEndChan     ; };
  int    getNChan       () const { return fNChan       ; };
  int    getChanWidth   () const { return fChanWidth   ; };
  int    getNHits       () const { return fNHits       ; };
  int    getType        () const { return fType        ; };
  int    getTriggerFlag () const { return fTriggerFlag ; };
  float  getHitSADC     () const { return fHitSADC     ; };
  float  getFirstTimeHit() const { return fFirstHitTime; };
  float  getLastTimeHit () const { return fLastHitTime ; };
  float  getTimeWidth   () const { return fTimeWidth   ; };
  double getMC_EnergyNu () const { return fMC_EnergyNu ; };
  double getMC_EnergyLep() const { return fMC_EnergyLep; };
  double getMC_MarlTime () const { return fMC_MarlTime ; };

  std::vector<recoHit> getHits() const { return fHitVector; };

  void   setHitSADC     (float  cHitSADC     ) { fHitSADC      = cHitSADC     ; };
  void   setTriggerFlag (int    cTriggerFlag ) { fTriggerFlag  = cTriggerFlag ; };
  void   setMC_EnergyNu (double cMC_EnergyNu ) { fMC_EnergyNu  = cMC_EnergyNu ; };
  void   setMC_EnergyLep(double cMC_EnergyLep) { fMC_EnergyLep = cMC_EnergyLep; };
  void   setMC_MarlTime (double cMC_Marltime ) { fMC_MarlTime  = cMC_Marltime ; };

  void  printCluster() const
  {
    std::cout << "*********************"          << std::endl;
    std::cout << "Event        " << fEvent        << std::endl;
    std::cout << "StartChan    " << fStartChan    << std::endl;
    std::cout << "EndChan      " << fEndChan      << std::endl;
    std::cout << "NChan        " << fNChan        << std::endl;
    std::cout << "ChanWidth    " << fChanWidth    << std::endl;
    std::cout << "NHits        " << fNHits        << std::endl;
    std::cout << "Type         " << fType         << std::endl;
    std::cout << "TriggerFlag  " << fTriggerFlag  << std::endl;
    std::cout << "HitSADC      " << fHitSADC      << std::endl;
    std::cout << "FirstHitTime " << fFirstHitTime << std::endl;
    std::cout << "LastHitTime  " << fLastHitTime  << std::endl;
    std::cout << "TimeWidth    " << fTimeWidth    << std::endl;
    std::cout << "MC_EnergyNu  " << fMC_EnergyNu  << std::endl;
    std::cout << "MC_EnergyLep " << fMC_EnergyLep << std::endl;
    std::cout << "MC_MarlTime  " << fMC_MarlTime  << std::endl;
    for(int i = 0; i < fNHits; i++)
       fHitVector.at(i).Print();
    std::cout << "*********************"          << std::endl;

  };


 private:
  int    fEvent        = 0;
  int    fStartChan    = 0;
  int    fEndChan      = 0;
  int    fNChan        = 0;
  int    fChanWidth    = 0;
  int    fNHits        = 0;
  int    fType         = 0;
  int    fTriggerFlag  = 0;
  float  fHitSADC      = 0;
  float  fFirstHitTime = 0;
  float  fLastHitTime  = 0;
  float  fTimeWidth    = 0;
  double fMC_EnergyNu  = 0;
  double fMC_EnergyLep = 0;
  double fMC_MarlTime  = 0;
  std::vector<recoHit> fHitVector;
};

class ClusterHitsInTime
{
 public:
  ClusterHitsInTime(double TimeWindow): fNClusters(0),
                                           fTimeWindow(TimeWindow),
                                           fHitVector(),
                                           fVecClusters(){};
  
  void DoIt(std::vector<recoHit> cHitVector);
  double               GetTimeWindow()    const { return fTimeWindow;  };
  int                  getNClusters()     const { return fNClusters;   };
  std::vector<cluster> getClusterVector() const { return fVecClusters; };

  void   SetTimeWindow(double inTimeWindow=1) { fTimeWindow = inTimeWindow; };

private:
  int fNClusters       = 0;
  double fTimeWindow   = 0;
  std::vector<recoHit> fHitVector;
  std::vector<cluster> fVecClusters;
};
    

class DAQQuickClustering : public art::EDAnalyzer{
public:
  explicit DAQQuickClustering(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  DAQQuickClustering(DAQQuickClustering const &) = delete;
  DAQQuickClustering(DAQQuickClustering &&) = delete;
  DAQQuickClustering & operator = (DAQQuickClustering const &) = delete;
  DAQQuickClustering & operator = (DAQQuickClustering &&) = delete;

  // The main guts...
  void analyze(art::Event const & evt) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;

  void endJob() override;
  
private:

  void ResetVariables();

  void trigger(std::vector<cluster> &vec_Clusters, unsigned int const &config);
  void clusterCut(std::vector<cluster> &vec_Clusters, unsigned int const &config);
  void clusterChannels(std::vector<recoHit> &vec_Hits, std::vector<ClusterHitsInTime> &vec_ChannelCluster, unsigned int const &config);
  void makeConfigGraph();
  void FillMyMaps(std::map<int,simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn,
                  art::ValidHandle<std::vector<simb::MCTruth> > Hand);
  PType WhichParType( int TrID );
  bool InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap );
  
  std::vector<int>   cut_AdjChanTolerance;
  std::vector<int>   cut_HitsInWindow;
  std::vector<int>   cut_MinChannels;
  std::vector<int>   cut_MinChanWidth;
  std::vector<float> cut_TimeWindowSize;
  std::vector<float> cut_TotalADC;

  int firstCatch;
  int secondCatch;
  int thirdCatch;

  std::string fRawDigitLabel;
  std::string fHitLabel  ;
  std::string fGEANTLabel;
  std::string fMARLLabel ; std::map< int, simb::MCParticle > MarlParts;
  std::string fAPALabel  ; std::map< int, simb::MCParticle > APAParts ;
  std::string fCPALabel  ; std::map< int, simb::MCParticle > CPAParts ;
  std::string fAr39Label ; std::map< int, simb::MCParticle > Ar39Parts;
  std::string fNeutLabel ; std::map< int, simb::MCParticle > NeutParts;
  std::string fKrypLabel ; std::map< int, simb::MCParticle > KrypParts;
  std::string fPlonLabel ; std::map< int, simb::MCParticle > PlonParts;
  std::string fRdonLabel ; std::map< int, simb::MCParticle > RdonParts;
  std::string fAr42Label ; std::map< int, simb::MCParticle > Ar42Parts;

  // Mapping from track ID to particle type, for use in WhichParType()
  std::map<int, PType> trkIDToPType;

  unsigned int NConfigs;
  unsigned int NCuts;

  double detectorScaling;

  int Run   ;
  int SubRun;
  int Event ;

  std::map<int,std::vector<double>> map_EventToMC;


  
  TTree* t_Output_unusedhits;
  TTree* t_Output_clusteredhits;

  TH1D* h_ENu_MC;
  TH1D* h_MarlTime_MC;
  TH1D* h_TimeElapsed;
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
  double out_MarlTime;
  std::vector<int> out_HitView;
  std::vector<int> out_GenType;
  std::vector<int> out_HitChan;
  std::vector<double> out_HitTime;
  std::vector<double> out_HitSADC;
  std::vector<double> out_HitRMS;

  int   NTotHits;
  int   NColHits;
  int   NIndHits;
  int   nHitsNotBackTracked;
  int   HitView[nMaxHits];
  int   HitSize[nMaxHits];
  int   HitTPC [nMaxHits];
  int   HitChan[nMaxHits];
  float HitTime[nMaxHits];
  float HitRMS [nMaxHits];
  float HitSADC[nMaxHits];
  float HitInt [nMaxHits];
  float HitPeak[nMaxHits];
  int   GenType[nMaxHits];
  float Hit_X[nMaxHits];                  
  float Hit_Y[nMaxHits];                  
  float Hit_Z[nMaxHits];                  
  float Hit_Energy[nMaxHits];             
  float Hit_NumElectrons[nMaxHits];
  int   NCorrespondingIDEs[nMaxHits]; 
  int   TotGen_Marl;
  int   TotGen_APA;
  int   TotGen_CPA;
  int   TotGen_Ar39;
  int   TotGen_Neut;
  int   TotGen_Kryp;
  int   TotGen_Plon;
  int   TotGen_Rdon;
  int   TotGen_Ar42;

  int    VertexChan;
  int    Nu_Type;
  int    Nu_Lep_Type;
  int    Mode;
  int    CCNC;
  int    HitNucleon;
  int    Target;
  std::vector<int>    MarlSample;
  std::vector<double> MarlTime;
  std::vector<double> MarlWeight;
  double ENu;
  double ENu_Lep;
  double VertX;
  double VertY;
  double VertZ;
  double VertexT;
  double Px;
  double Py;
  double Pz;

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  

};
#endif
