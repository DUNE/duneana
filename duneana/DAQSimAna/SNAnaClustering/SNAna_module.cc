#include <functional>   // std::greater

//ROOT includes
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TMath.h"
#include "TGeoMatrix.h" // TGeoHMatrix

// LArSoft libraries

//LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SupernovaTruth.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"

//ART includes
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
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

enum PType{kUnknown=0, kMarl, kAPA, kCPA, kAr39, kNeut, kKryp, kPlon, kRdon , kAr42};
std::map<PType, std::string> PTypeString{
  {kUnknown,"Unknown"},
  {kMarl   ,"Marl"   },
  {kAPA    ,"APA"    },
  {kCPA    ,"CPA"    },
  {kAr39   ,"Ar39"   },
  {kNeut   ,"Neut"   },
  {kKryp   ,"Kryp"   },
  {kPlon   ,"Plon"   },
  {kRdon   ,"Rdon"   },
  {kAr42   ,"Ar42"   }};


class SNAna : public art::EDAnalyzer {

public:
  explicit SNAna(fhicl::ParameterSet const & p);

  SNAna(SNAna const &) = delete;
  SNAna(SNAna &&) = delete;
  SNAna & operator = (SNAna const &) = delete;
  SNAna & operator = (SNAna &&) = delete;

  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;
  void endJob() override;

private:

  void ResetVariables();
  void FillMyMaps  ( std::map< int, simb::MCParticle> &MyMap,
                     art::FindManyP<simb::MCParticle> Assn,
                     art::Handle< std::vector<simb::MCTruth> > Hand,
                     std::map<int, int>* indexMap=nullptr);
  PType WhichParType( int TrID );
  bool  InMyMap     ( int TrID, std::map< int, simb::MCParticle> ParMap );
  void FillTruth(const art::FindManyP<simb::MCParticle> Assn,
                 const art::Handle<std::vector<simb::MCTruth>>& Hand,
                 const PType type) ;
  void SaveNeighbourADC(int channel,
                        art::Handle< std::vector<raw::RawDigit> >rawDigitsVecHandle,
                        std::set<int> badChannels,
                        recob::Hit const& hits);

  void SaveIDEs(art::Event const & evt);
  std::string fname;
  
  int firstCatch;
  int secondCatch;
  int thirdCatch;

  std::string fRawDigitLabel;
  std::string fHitLabel;
  std::string fCalDataModuleLabel;
  std::string fOpHitModuleLabel;

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
  std::map<int, const simb::MCParticle*> truthmap;
  // Which MARLEY interaction (if any) caused this true track ID?
  std::map<int, int> trkIDToMarleyIndex;

  // Mapping from track ID to particle type, for use in WhichParType()
  std::map<int, PType> trkIDToPType;

  bool fSaveNeighbourADCs;
  bool fSaveIDEs;
  bool fSaveTruth;
  bool fSaveTPC;
  bool fSavePDS;

  TTree* fSNAnaTree;

  int Run;
  int SubRun;
  int Event;

  int NTotHit   ;
  int NColHit   ;
  int NIndHit   ;
  int NHitNoBT  ;

  std::vector<int>                  Hit_View                 ;
  std::vector<int>                  Hit_Size                 ;
  std::vector<int>                  Hit_TPC                  ;
  std::vector<int>                  Hit_Chan                 ;
  std::vector<double>               Hit_X_start              ;
  std::vector<double>               Hit_Y_start              ;
  std::vector<double>               Hit_Z_start              ;
  std::vector<double>               Hit_X_end                ;
  std::vector<double>               Hit_Y_end                ;
  std::vector<double>               Hit_Z_end                ;
  std::vector<float>                Hit_Time                 ;
  std::vector<float>                Hit_RMS                  ;
  std::vector<float>                Hit_SADC                 ;
  std::vector<float>                Hit_Int                  ;
  std::vector<float>                Hit_Peak                 ;
  std::vector<int>                  Hit_True_GenType         ;
  std::vector<int>                  Hit_True_MainTrID        ;
  std::vector<int>                  Hit_True_TrackID         ;
  std::vector<float>                Hit_True_EvEnergy        ;
  std::vector<int>                  Hit_True_MarleyIndex     ;
  std::vector<float>                Hit_True_X               ;
  std::vector<float>                Hit_True_Y               ;
  std::vector<float>                Hit_True_Z               ;
  std::vector<float>                Hit_True_Energy          ;
  std::vector<float>                Hit_True_nElec           ;
  std::vector<int>                  Hit_True_nIDEs           ;
  std::vector<int>                  Hit_AdjM5SADC            ;
  std::vector<int>                  Hit_AdjM2SADC            ;
  std::vector<int>                  Hit_AdjM1SADC            ;
  std::vector<int>                  Hit_AdjP1SADC            ;
  std::vector<int>                  Hit_AdjP2SADC            ;
  std::vector<int>                  Hit_AdjP5SADC            ;
  std::vector<int>                  Hit_AdjM5Chan            ;
  std::vector<int>                  Hit_AdjM2Chan            ;
  std::vector<int>                  Hit_AdjM1Chan            ;
  std::vector<int>                  Hit_AdjP1Chan            ;
  std::vector<int>                  Hit_AdjP2Chan            ;
  std::vector<int>                  Hit_AdjP5Chan            ;

  std::vector<int>                  PDS_OpHit_OpChannel      ;
  std::vector<double>               PDS_OpHit_X              ;
  std::vector<double>               PDS_OpHit_Y              ;
  std::vector<double>               PDS_OpHit_Z              ;
  std::vector<double>               PDS_OpHit_PeakTimeAbs    ;
  std::vector<double>               PDS_OpHit_PeakTime       ;
  std::vector<unsigned short>       PDS_OpHit_Frame          ;
  std::vector<double>               PDS_OpHit_Width          ;
  std::vector<double>               PDS_OpHit_Area           ;
  std::vector<double>               PDS_OpHit_Amplitude      ;
  std::vector<double>               PDS_OpHit_PE             ;
  std::vector<double>               PDS_OpHit_FastToTotal    ;
  std::vector<int>                  PDS_OpHit_True_GenType   ;
  std::vector<int>                  PDS_OpHit_True_Index     ;
  std::vector<double>               PDS_OpHit_True_Energy    ;
  std::vector<int>                  PDS_OpHit_True_TrackID   ;
  std::vector<int>                  PDS_OpHit_True_GenTypeAll;
  std::vector<double>               PDS_OpHit_True_EnergyAll ;
  std::vector<int>                  PDS_OpHit_True_TrackIDAll;
  std::vector<int>                  PDS_OpHit_True_IndexAll  ;

  std::vector<int>                  True_VertexChan          ;
  std::vector<int>                  True_Nu_Type             ;
  std::vector<int>                  True_Nu_Lep_Type         ;
  std::vector<int>                  True_Mode                ;
  std::vector<int>                  True_CCNC                ;
  std::vector<int>                  True_HitNucleon          ;
  std::vector<int>                  True_Target              ;
  std::vector<int>                  True_MarlSample          ;
  std::vector<float>                True_MarlTime            ;
  std::vector<float>                True_MarlWeight          ;
  std::vector<float>                True_ENu                 ;
  std::vector<float>                True_ENu_Lep             ;
  std::vector<float>                True_VertX               ;
  std::vector<float>                True_VertY               ;
  std::vector<float>                True_VertZ               ;
  std::vector<float>                True_VertexT             ;
  std::vector<float>                True_Px                  ;
  std::vector<float>                True_Py                  ;
  std::vector<float>                True_Pz                  ;
  std::vector<float>                True_Dirx                ;
  std::vector<float>                True_Diry                ;
  std::vector<float>                True_Dirz                ;
  std::vector<float>                True_Time                ;

  std::vector<int>                  True_Bck_Mode            ;   
  std::vector<int>                  True_Bck_PDG             ;
  std::vector<int>                  True_Bck_ID              ;
  std::vector<std::string>          True_Bck_Process         ; // str because why not
  std::vector<std::string>          True_Bck_EndProcess      ;
  std::vector<int>                  True_Bck_Mother          ;
  std::vector<double>               True_Bck_P               ;
  std::vector<double>               True_Bck_VertX           ;
  std::vector<double>               True_Bck_VertY           ;
  std::vector<double>               True_Bck_VertZ           ;
  std::vector<double>               True_Bck_Time            ;
  std::vector<double>               True_Bck_Energy          ;
  std::vector<double>               True_Bck_EndX            ;
  std::vector<double>               True_Bck_EndY            ;
  std::vector<double>               True_Bck_EndZ            ;
  std::vector<double>               True_Bck_EndT            ;
  std::vector<double>               True_Bck_EndE            ;   

  int   NTotIDEs;
  std::vector<int>                  IDEChannel               ;
  std::vector<int>                  IDEStartTime             ;
  std::vector<int>                  IDEEndTime               ;
  std::vector<float>                IDEEnergy                ;
  std::vector<float>                IDEElectrons             ;
  std::vector<int>                  IDEParticle              ;

  int   TotGen_Marl;
  int   TotGen_APA ;
  int   TotGen_CPA ;
  int   TotGen_Ar39;
  int   TotGen_Neut;
  int   TotGen_Kryp;
  int   TotGen_Plon;
  int   TotGen_Rdon;
  int   TotGen_Ar42;


  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;

};

SNAna::SNAna(fhicl::ParameterSet const & p):EDAnalyzer(p),
                                            fname("SNAna_module")
{
  this->reconfigure(p);
}


void SNAna::reconfigure(fhicl::ParameterSet const & p)
{
  
  fRawDigitLabel      = p.get<std::string>("RawDigitLabel"     );
  fHitLabel           = p.get<std::string>("HitLabel"          );
  fCalDataModuleLabel = p.get<std::string>("CalDataModuleLabel");
  fOpHitModuleLabel   = p.get<std::string>("OpHitModuleLabel"  );

  fGEANTLabel = p.get<std::string> ("GEANT4Label"  );
  fMARLLabel  = p.get<std::string> ("MARLEYLabel"  );
  fAPALabel   = p.get<std::string> ("APALabel"     );
  fCPALabel   = p.get<std::string> ("CPALabel"     );
  fAr39Label  = p.get<std::string> ("Argon39Label" );
  fNeutLabel  = p.get<std::string> ("NeutronLabel" );
  fKrypLabel  = p.get<std::string> ("KryptonLabel" );
  fPlonLabel  = p.get<std::string> ("PoloniumLabel");
  fRdonLabel  = p.get<std::string> ("RadonLabel"   );
  fAr42Label  = p.get<std::string> ("Argon42Label" );

  fSaveNeighbourADCs = p.get<bool> ("SaveNeighbourADCs",0);
  fSaveTruth = p.get<bool>("SaveTruth",0);
  fSaveIDEs  = p.get<bool>("SaveIDEs",0);
  fSaveTPC   = p.get<bool>("SaveTPC",1);
  fSavePDS   = p.get<bool>("SavePDS",1);

  mf::LogInfo(fname) << "Reconfigured " << this->processName() << " with "
                     << " SaveNeighbourADCs: " << std::boolalpha << fSaveNeighbourADCs
                     << " SaveTruth: " << std::boolalpha << fSaveTruth
                     << " SaveIDEs: " << std::boolalpha << fSaveIDEs << std::endl;
}


void SNAna::ResetVariables()
{
  trkIDToPType.clear();

  MarlParts.clear(); APAParts .clear(); CPAParts .clear(); Ar39Parts.clear();
  NeutParts.clear(); KrypParts.clear(); PlonParts.clear(); RdonParts.clear();
  Ar42Parts.clear();

  Run = SubRun = Event = -1;

  TotGen_Marl = TotGen_APA  = TotGen_CPA  = TotGen_Ar39 = 0;
  TotGen_Neut = TotGen_Kryp = TotGen_Plon = TotGen_Rdon = 0;
  TotGen_Ar42 = 0;

  NTotHit    = 0;
  NColHit    = 0;
  NIndHit    = 0;
  NHitNoBT   = 0;

  Hit_View                 .clear();
  Hit_Size                 .clear();
  Hit_TPC                  .clear();
  Hit_Chan                 .clear();
  Hit_X_start              .clear();
  Hit_Y_start              .clear();
  Hit_Z_start              .clear();
  Hit_X_end                .clear();
  Hit_Y_end                .clear();
  Hit_Z_end                .clear();
  Hit_Time                 .clear();
  Hit_RMS                  .clear();
  Hit_SADC                 .clear();
  Hit_Int                  .clear();
  Hit_Peak                 .clear();
  Hit_True_GenType         .clear();
  Hit_True_MainTrID        .clear();
  Hit_True_TrackID         .clear();
  Hit_True_EvEnergy        .clear();
  Hit_True_MarleyIndex     .clear();
  Hit_True_X               .clear();
  Hit_True_Y               .clear();
  Hit_True_Z               .clear();
  Hit_True_Energy          .clear();
  Hit_True_nElec           .clear();
  Hit_True_nIDEs           .clear();

  Hit_AdjM5SADC            .clear();
  Hit_AdjM2SADC            .clear();
  Hit_AdjM1SADC            .clear();
  Hit_AdjP1SADC            .clear();
  Hit_AdjP2SADC            .clear();
  Hit_AdjP5SADC            .clear();
  Hit_AdjM5Chan            .clear();
  Hit_AdjM2Chan            .clear();
  Hit_AdjM1Chan            .clear();
  Hit_AdjP1Chan            .clear();
  Hit_AdjP2Chan            .clear();
  Hit_AdjP5Chan            .clear();

  PDS_OpHit_OpChannel      .clear();
  PDS_OpHit_X              .clear();
  PDS_OpHit_Y              .clear();
  PDS_OpHit_Z              .clear();
  PDS_OpHit_PeakTimeAbs    .clear();
  PDS_OpHit_PeakTime       .clear();
  PDS_OpHit_Frame          .clear();
  PDS_OpHit_Width          .clear();
  PDS_OpHit_Area           .clear();
  PDS_OpHit_Amplitude      .clear();
  PDS_OpHit_PE             .clear();
  PDS_OpHit_FastToTotal    .clear();
  PDS_OpHit_True_GenType   .clear();
  PDS_OpHit_True_Index     .clear();
  PDS_OpHit_True_Energy    .clear();
  PDS_OpHit_True_TrackID   .clear();
  PDS_OpHit_True_GenTypeAll.clear();
  PDS_OpHit_True_EnergyAll .clear();
  PDS_OpHit_True_TrackIDAll.clear();
  PDS_OpHit_True_IndexAll  .clear();

  True_VertexChan          .clear();
  True_Nu_Type             .clear();
  True_Nu_Lep_Type         .clear();
  True_Mode                .clear();
  True_CCNC                .clear();
  True_HitNucleon          .clear();
  True_Target              .clear();
  True_MarlSample          .clear();
  True_MarlTime            .clear();
  True_MarlWeight          .clear();
  True_ENu                 .clear();
  True_ENu_Lep             .clear();
  True_VertX               .clear();
  True_VertY               .clear();
  True_VertZ               .clear();
  True_VertexT             .clear();
  True_Px                  .clear();
  True_Py                  .clear();
  True_Pz                  .clear();
  True_Dirx                .clear();
  True_Diry                .clear();
  True_Dirz                .clear();
  True_Time                .clear();

  True_Bck_Mode            .clear();
  True_Bck_PDG             .clear();
  True_Bck_ID              .clear();
  True_Bck_Process         .clear();
  True_Bck_EndProcess      .clear();
  True_Bck_Mother          .clear();
  True_Bck_P               .clear();
  True_Bck_VertX           .clear();
  True_Bck_VertY           .clear();
  True_Bck_VertZ           .clear();
  True_Bck_Time            .clear();
  True_Bck_Energy          .clear();
  True_Bck_EndX            .clear();
  True_Bck_EndY            .clear();
  True_Bck_EndZ            .clear();
  True_Bck_EndT            .clear();
  True_Bck_EndE            .clear();

  // IDEs
  NTotIDEs=0;
  IDEChannel               .clear();
  IDEStartTime             .clear();
  IDEEndTime               .clear();
  IDEEnergy                .clear();
  IDEElectrons             .clear();
  IDEParticle              .clear();

}


void SNAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;


  fSNAnaTree = tfs->make<TTree>("SNSimTree","SN simulation analysis tree");

  fSNAnaTree->Branch("Run"       , &Run       , "Run/I"       );
  fSNAnaTree->Branch("SubRun"    , &SubRun    , "SubRun/I"    );
  fSNAnaTree->Branch("Event"     , &Event     , "Event/I"     );
  fSNAnaTree->Branch("NTotHit"   , &NTotHit   , "NTotHits/I"  );
  fSNAnaTree->Branch("NColHit"   , &NColHit   , "NColHits/I"  );
  fSNAnaTree->Branch("NIndHit"   , &NIndHit   , "NIndHits/I"  );
  fSNAnaTree->Branch("NHitNoBT"  , &NHitNoBT  , "NHitNoBT/I"  );

  if (fSaveTPC) {
    fSNAnaTree->Branch("Hit_View"                 , &Hit_View                 );
    fSNAnaTree->Branch("Hit_Size"                 , &Hit_Size                 );
    fSNAnaTree->Branch("Hit_TPC"                  , &Hit_TPC                  );
    fSNAnaTree->Branch("Hit_Chan"                 , &Hit_Chan                 );
    fSNAnaTree->Branch("Hit_X_start"              , &Hit_X_start              );
    fSNAnaTree->Branch("Hit_Y_start"              , &Hit_Y_start              );
    fSNAnaTree->Branch("Hit_Z_start"              , &Hit_Z_start              );
    fSNAnaTree->Branch("Hit_X_end"                , &Hit_X_end                );
    fSNAnaTree->Branch("Hit_Y_end"                , &Hit_Y_end                );
    fSNAnaTree->Branch("Hit_Z_end"                , &Hit_Z_end                );
    fSNAnaTree->Branch("Hit_Time"                 , &Hit_Time                 );
    fSNAnaTree->Branch("Hit_RMS"                  , &Hit_RMS                  );
    fSNAnaTree->Branch("Hit_SADC"                 , &Hit_SADC                 );
    fSNAnaTree->Branch("Hit_Int"                  , &Hit_Int                  );
    fSNAnaTree->Branch("Hit_Peak"                 , &Hit_Peak                 );
    fSNAnaTree->Branch("Hit_True_GenType"         , &Hit_True_GenType         );
    fSNAnaTree->Branch("Hit_True_MainTrID"        , &Hit_True_MainTrID        );
    fSNAnaTree->Branch("Hit_True_TrackID"         , &Hit_True_TrackID         );
    fSNAnaTree->Branch("Hit_True_EvEnergy"        , &Hit_True_EvEnergy        );
    fSNAnaTree->Branch("Hit_True_MarleyIndex"     , &Hit_True_MarleyIndex     );
    fSNAnaTree->Branch("Hit_True_X"               , &Hit_True_X               );
    fSNAnaTree->Branch("Hit_True_Y"               , &Hit_True_Y               );
    fSNAnaTree->Branch("Hit_True_Z"               , &Hit_True_Z               );
    fSNAnaTree->Branch("Hit_True_Energy"          , &Hit_True_Energy          );
    fSNAnaTree->Branch("Hit_True_nElec"           , &Hit_True_nElec           );
    fSNAnaTree->Branch("Hit_True_nIDEs"           , &Hit_True_nIDEs           );
  }

  if (fSaveNeighbourADCs) {
    fSNAnaTree->Branch("Hit_AdjM5SADC"            , &Hit_AdjM5SADC            );
    fSNAnaTree->Branch("Hit_AdjM2SADC"            , &Hit_AdjM2SADC            );
    fSNAnaTree->Branch("Hit_AdjM1SADC"            , &Hit_AdjM1SADC            );
    fSNAnaTree->Branch("Hit_AdjP1SADC"            , &Hit_AdjP1SADC            );
    fSNAnaTree->Branch("Hit_AdjP2SADC"            , &Hit_AdjP2SADC            );
    fSNAnaTree->Branch("Hit_AdjP5SADC"            , &Hit_AdjP5SADC            );
    fSNAnaTree->Branch("Hit_AdjM5Chan"            , &Hit_AdjM5Chan            );
    fSNAnaTree->Branch("Hit_AdjM2Chan"            , &Hit_AdjM2Chan            );
    fSNAnaTree->Branch("Hit_AdjM1Chan"            , &Hit_AdjM1Chan            );
    fSNAnaTree->Branch("Hit_AdjP1Chan"            , &Hit_AdjP1Chan            );
    fSNAnaTree->Branch("Hit_AdjP2Chan"            , &Hit_AdjP2Chan            );
    fSNAnaTree->Branch("Hit_AdjP5Chan"            , &Hit_AdjP5Chan            );
  }

  if (fSavePDS) {
    fSNAnaTree->Branch("PDS_OpHit_OpChannel"      , &PDS_OpHit_OpChannel      );
    fSNAnaTree->Branch("PDS_OpHit_X"              , &PDS_OpHit_X              );
    fSNAnaTree->Branch("PDS_OpHit_Y"              , &PDS_OpHit_Y              );
    fSNAnaTree->Branch("PDS_OpHit_Z"              , &PDS_OpHit_Z              );
    fSNAnaTree->Branch("PDS_OpHit_PeakTimeAbs"    , &PDS_OpHit_PeakTimeAbs    );
    fSNAnaTree->Branch("PDS_OpHit_PeakTime"       , &PDS_OpHit_PeakTime       );
    fSNAnaTree->Branch("PDS_OpHit_Frame"          , &PDS_OpHit_Frame          );
    fSNAnaTree->Branch("PDS_OpHit_Width"          , &PDS_OpHit_Width          );
    fSNAnaTree->Branch("PDS_OpHit_Area"           , &PDS_OpHit_Area           );
    fSNAnaTree->Branch("PDS_OpHit_Amplitude"      , &PDS_OpHit_Amplitude      );
    fSNAnaTree->Branch("PDS_OpHit_PE"             , &PDS_OpHit_PE             );
    fSNAnaTree->Branch("PDS_OpHit_FastToTotal"    , &PDS_OpHit_FastToTotal    );
    fSNAnaTree->Branch("PDS_OpHit_True_GenType"   , &PDS_OpHit_True_GenType   );
    fSNAnaTree->Branch("PDS_OpHit_True_Energy"    , &PDS_OpHit_True_Energy    );
    fSNAnaTree->Branch("PDS_OpHit_True_TrackID"   , &PDS_OpHit_True_TrackID   );
    fSNAnaTree->Branch("PDS_OpHit_True_GenTypeAll", &PDS_OpHit_True_GenTypeAll);
    fSNAnaTree->Branch("PDS_OpHit_True_EnergyAll" , &PDS_OpHit_True_EnergyAll );
    fSNAnaTree->Branch("PDS_OpHit_True_TrackIDAll", &PDS_OpHit_True_TrackIDAll);
    fSNAnaTree->Branch("PDS_OpHit_True_IndexAll"  , &PDS_OpHit_True_IndexAll  );
  }

  fSNAnaTree->Branch("True_VertexChan"          , &True_VertexChan          );
  fSNAnaTree->Branch("True_Nu_Type"             , &True_Nu_Type             );
  fSNAnaTree->Branch("True_Nu_Lep_Type"         , &True_Nu_Lep_Type         );
  fSNAnaTree->Branch("True_Mode"                , &True_Mode                );
  fSNAnaTree->Branch("True_CCNC"                , &True_CCNC                );
  fSNAnaTree->Branch("True_HitNucleon"          , &True_HitNucleon          );
  fSNAnaTree->Branch("True_Target"              , &True_Target              );
  fSNAnaTree->Branch("True_MarlSample"          , &True_MarlSample          );
  fSNAnaTree->Branch("True_MarlTime"            , &True_MarlTime            );
  fSNAnaTree->Branch("True_MarlWeight"          , &True_MarlWeight          );
  fSNAnaTree->Branch("True_ENu"                 , &True_ENu                 );
  fSNAnaTree->Branch("True_ENu_Lep"             , &True_ENu_Lep             );
  fSNAnaTree->Branch("True_VertX"               , &True_VertX               );
  fSNAnaTree->Branch("True_VertY"               , &True_VertY               );
  fSNAnaTree->Branch("True_VertZ"               , &True_VertZ               );
  fSNAnaTree->Branch("True_VertexT"             , &True_VertexT             );
  fSNAnaTree->Branch("True_Px"                  , &True_Px                  );
  fSNAnaTree->Branch("True_Py"                  , &True_Py                  );
  fSNAnaTree->Branch("True_Pz"                  , &True_Pz                  );
  fSNAnaTree->Branch("True_Dirx"                , &True_Dirx                );
  fSNAnaTree->Branch("True_Diry"                , &True_Diry                );
  fSNAnaTree->Branch("True_Dirz"                , &True_Dirz                );
  fSNAnaTree->Branch("True_Time"                , &True_Time                );


  fSNAnaTree->Branch("True_Bck_Mode"            , &True_Bck_Mode            );
  fSNAnaTree->Branch("True_Bck_PDG"             , &True_Bck_PDG             );
  fSNAnaTree->Branch("True_Bck_ID"              , &True_Bck_ID              );
  fSNAnaTree->Branch("True_Bck_Process"         , &True_Bck_Process         );
  fSNAnaTree->Branch("True_Bck_EndProcess"      , &True_Bck_EndProcess      );
  fSNAnaTree->Branch("True_Bck_Mother"          , &True_Bck_Mother          );
  fSNAnaTree->Branch("True_Bck_P"               , &True_Bck_P               );
  fSNAnaTree->Branch("True_Bck_VertX"           , &True_Bck_VertX           );
  fSNAnaTree->Branch("True_Bck_VertY"           , &True_Bck_VertY           );
  fSNAnaTree->Branch("True_Bck_VertZ"           , &True_Bck_VertZ           );
  fSNAnaTree->Branch("True_Bck_Time"            , &True_Bck_Time            );
  fSNAnaTree->Branch("True_Bck_Energy"          , &True_Bck_Energy          );
  fSNAnaTree->Branch("True_Bck_EndX"            , &True_Bck_EndX            );
  fSNAnaTree->Branch("True_Bck_EndY"            , &True_Bck_EndY            );
  fSNAnaTree->Branch("True_Bck_EndZ"            , &True_Bck_EndZ            );
  fSNAnaTree->Branch("True_Bck_EndT"            , &True_Bck_EndT            );
  fSNAnaTree->Branch("True_Bck_EndE"            , &True_Bck_EndE            );  

  // IDEs
  if(fSaveIDEs) {
    fSNAnaTree->Branch("NTotIDEs"                 , &NTotIDEs  , "NTotIDEs/I" );
    fSNAnaTree->Branch("IDEChannel"               , &IDEChannel               );
    fSNAnaTree->Branch("IDEStartTime"             , &IDEStartTime             );
    fSNAnaTree->Branch("IDEEndTime"               , &IDEEndTime               );
    fSNAnaTree->Branch("IDEEnergy"                , &IDEEnergy                );
    fSNAnaTree->Branch("IDEElectrons"             , &IDEElectrons             );
    fSNAnaTree->Branch("IDEParticle"              , &IDEParticle              );
  }

  fSNAnaTree->Branch("TotGen_Marl", &TotGen_Marl, "TotGen_Marl/I");
  fSNAnaTree->Branch("TotGen_APA" , &TotGen_APA , "TotGen_APA/I" );
  fSNAnaTree->Branch("TotGen_CPA" , &TotGen_CPA , "TotGen_CPA/I" );
  fSNAnaTree->Branch("TotGen_Ar39", &TotGen_Ar39, "TotGen_Ar39/I");
  fSNAnaTree->Branch("TotGen_Neut", &TotGen_Neut, "TotGen_Neut/I");
  fSNAnaTree->Branch("TotGen_Kryp", &TotGen_Kryp, "TotGen_Kryp/I");
  fSNAnaTree->Branch("TotGen_Plon", &TotGen_Plon, "TotGen_Plon/I");
  fSNAnaTree->Branch("TotGen_Rdon", &TotGen_Rdon, "TotGen_Rdon/I");
  fSNAnaTree->Branch("TotGen_Ar42", &TotGen_Ar42, "TotGen_Ar42/I");


}


void SNAna::analyze(art::Event const & evt)
{
  ResetVariables();

  Run    = evt.run();
  SubRun = evt.subRun();
  Event  = evt.event();

  //GET INFORMATION ABOUT THE DETECTOR'S GEOMETRY.
  auto const* geo = lar::providerFrom<geo::Geometry>();

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

  //LIFT OUT THE MARLEY PARTICLES.
  //auto MarlTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fMARLLabel);
  auto MarlTrue = evt.getHandle< std::vector<simb::MCTruth> >(fMARLLabel);
  if (MarlTrue) {
    art::FindManyP<simb::MCParticle> MarlAssn(MarlTrue,evt,fGEANTLabel);
    FillMyMaps( MarlParts, MarlAssn, MarlTrue, &trkIDToMarleyIndex );
    TotGen_Marl = MarlParts.size();
    double Px_(0), Py_(0), Pz_(0), Pnorm(1);

    for(size_t i = 0; i < MarlTrue->size(); i++)
    {
      True_Nu_Type    .push_back(MarlTrue->at(i).GetNeutrino().Nu().PdgCode());
      True_Nu_Lep_Type.push_back(MarlTrue->at(i).GetNeutrino().Lepton().PdgCode());
      True_Mode       .push_back(MarlTrue->at(i).GetNeutrino().Mode());
      True_CCNC       .push_back(MarlTrue->at(i).GetNeutrino().CCNC());
      True_Target     .push_back(MarlTrue->at(i).GetNeutrino().Target());
      True_HitNucleon .push_back(MarlTrue->at(i).GetNeutrino().HitNuc());
      True_VertX      .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Vx());
      True_VertY      .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Vy());
      True_VertZ      .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Vz());
      True_ENu_Lep    .push_back(MarlTrue->at(i).GetNeutrino().Lepton().E());
      True_ENu        .push_back(MarlTrue->at(i).GetNeutrino().Nu().E());
      True_Px         .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Px());
      True_Py         .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Py());
      True_Pz         .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Pz());
      Pnorm = std::sqrt(Px_*Px_+Py_*Py_+Pz_*Pz_);
      double Px = Px_/Pnorm;
      double Py = Py_/Pnorm;
      double Pz = Pz_/Pnorm;
      True_Dirx.push_back(Px);
      True_Diry.push_back(Py);
      True_Dirz.push_back(Pz);
      True_Time.push_back(MarlTrue->at(i).GetNeutrino().Lepton().T());

      try {
        art::FindManyP<sim::SupernovaTruth> SNTruth(MarlTrue, evt, fMARLLabel);
        for (size_t j = 0; j < SNTruth.at(i).size(); j++) {
          const sim::SupernovaTruth ThisTr = (*SNTruth.at(i).at(j));
          True_MarlTime  .push_back(ThisTr.SupernovaTime);
          True_MarlWeight.push_back(ThisTr.Weight);
          True_MarlSample.push_back(ThisTr.SamplingMode);
        }
      } catch (...) {
        mf::LogDebug(fname) << "Didnt find SN truth (few things wont available):\n"
                            << " - MarlTime\n"
                            << " - MarlWeight\n"
                            << " - MarlSample\n";
      }
    }

    for(size_t i=0; i<True_VertX.size(); ++i) {
      bool caught = false;
      geo::Point_t const Vertex{True_VertX[i], True_VertY[i], True_VertZ[i]};
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
      if(caught==true)
      {
        True_VertexChan.push_back(-1);
      }
      else
      {
        True_VertexChan.push_back(geo->PlaneWireToChannel(WireID));
      }

      //CM/MICROSECOND.
      double drift_velocity = detProp.DriftVelocity(detProp.Efield(),detProp.Temperature());
      //CM/TICK
      drift_velocity = drift_velocity*0.5;
      True_VertexT.push_back(True_VertX.back()/drift_velocity);
    }
    if (fSaveTruth)
      FillTruth(MarlAssn, MarlTrue, kMarl);
  }

  auto APATrue = evt.getHandle< std::vector<simb::MCTruth> >(fAPALabel);
  if (APATrue) {
    art::FindManyP<simb::MCParticle> APAAssn(APATrue,evt,fGEANTLabel);
    FillMyMaps( APAParts, APAAssn, APATrue );
    TotGen_APA = APAParts.size();
    if(fSaveTruth) FillTruth(APAAssn , APATrue , kAPA );

  }

  auto CPATrue = evt.getHandle< std::vector<simb::MCTruth> >(fCPALabel);
  if (CPATrue) {
    art::FindManyP<simb::MCParticle> CPAAssn(CPATrue,evt,fGEANTLabel);
    FillMyMaps( CPAParts, CPAAssn, CPATrue );
    TotGen_CPA = CPAParts.size();
    if(fSaveTruth) FillTruth(CPAAssn , CPATrue , kCPA );
  }

  auto Ar39True = evt.getHandle< std::vector<simb::MCTruth> >(fAr39Label);
  if (Ar39True) {
    art::FindManyP<simb::MCParticle> Ar39Assn(Ar39True,evt,fGEANTLabel);
    FillMyMaps( Ar39Parts, Ar39Assn, Ar39True );
    TotGen_Ar39 = Ar39Parts.size();
    if(fSaveTruth) FillTruth(Ar39Assn, Ar39True, kAr39);
  }

  auto NeutTrue = evt.getHandle< std::vector<simb::MCTruth> >(fNeutLabel);
  if (NeutTrue) {
    art::FindManyP<simb::MCParticle> NeutAssn(NeutTrue,evt,fGEANTLabel);
    FillMyMaps( NeutParts, NeutAssn, NeutTrue );
    TotGen_Neut = NeutParts.size();
    if(fSaveTruth) FillTruth(NeutAssn, NeutTrue, kNeut);
  }

  auto KrypTrue = evt.getHandle< std::vector<simb::MCTruth> >(fKrypLabel);
  if (KrypTrue) {
    art::FindManyP<simb::MCParticle> KrypAssn(KrypTrue,evt,fGEANTLabel);
    FillMyMaps( KrypParts, KrypAssn, KrypTrue );
    TotGen_Kryp = KrypParts.size();
    if(fSaveTruth) FillTruth(KrypAssn, KrypTrue, kKryp);
  }

  auto PlonTrue = evt.getHandle< std::vector<simb::MCTruth> >(fPlonLabel);
  if (PlonTrue) {
    art::FindManyP<simb::MCParticle> PlonAssn(PlonTrue,evt,fGEANTLabel);
    FillMyMaps( PlonParts, PlonAssn, PlonTrue );
    TotGen_Plon = PlonParts.size();
    if(fSaveTruth) FillTruth(PlonAssn, PlonTrue, kPlon);
  }

  auto RdonTrue = evt.getHandle< std::vector<simb::MCTruth> >(fRdonLabel);
  if (RdonTrue) {
    art::FindManyP<simb::MCParticle> RdonAssn(RdonTrue,evt,fGEANTLabel);
    FillMyMaps( RdonParts, RdonAssn, RdonTrue );
    TotGen_Rdon = RdonParts.size();
    if(fSaveTruth) FillTruth(RdonAssn, RdonTrue, kRdon);
  }

  auto Ar42True = evt.getHandle< std::vector<simb::MCTruth> >(fAr42Label);
  if (Ar42True) {
    art::FindManyP<simb::MCParticle> Ar42Assn(Ar42True,evt,fGEANTLabel);
    FillMyMaps( Ar42Parts, Ar42Assn, Ar42True );
    TotGen_Ar42 = Ar42Parts.size();
    if(fSaveTruth) FillTruth(Ar42Assn, Ar42True, kAr42);
  }

  std::vector<simb::MCParticle> allTruthParts;
  for(auto& it: APAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: CPAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Ar39Parts)
    allTruthParts.push_back(it.second);
  for(auto& it: NeutParts)
    allTruthParts.push_back(it.second);
  for(auto& it: KrypParts)
    allTruthParts.push_back(it.second);
  for(auto& it: PlonParts)
    allTruthParts.push_back(it.second);
  for(auto& it: RdonParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Ar42Parts)
    allTruthParts.push_back(it.second);

  std::map<PType, std::map< int, simb::MCParticle >&> PTypeToMap{
    {kMarl, MarlParts},
    {kAPA,  APAParts },
    {kCPA,  CPAParts },
    {kAr39, Ar39Parts},
    {kAr42, Ar42Parts},
    {kNeut, NeutParts},
    {kKryp, KrypParts},
    {kPlon, PlonParts},
    {kRdon, RdonParts}
  };

  std::map<PType,std::set<int>> PTypeToTrackID;

  for(auto const& it : PTypeToMap){
    const PType p=it.first;
    auto const& m=it.second;
    for(auto const& it2 : m){
      trkIDToPType.insert(std::make_pair(it2.first, p));
      PTypeToTrackID[p].insert(it2.first);
    }
  }

  mf::LogInfo(fname) << "THE EVENTS NUMBER IS: " << Event << std::endl;

  if (fSaveTPC) {
    auto reco_hits = evt.getHandle< std::vector<recob::Hit> >(fHitLabel);
    auto rawDigitsVecHandle = evt.getHandle< std::vector<raw::RawDigit> >(fRawDigitLabel);

    if ( reco_hits && rawDigitsVecHandle ) {
      std::vector< recob::Hit > ColHits_Marl;
      std::vector< recob::Hit > ColHits_CPA ;
      std::vector< recob::Hit > ColHits_APA ;
      std::vector< recob::Hit > ColHits_Ar39;
      std::vector< recob::Hit > ColHits_Neut;
      std::vector< recob::Hit > ColHits_Kryp;
      std::vector< recob::Hit > ColHits_Plon;
      std::vector< recob::Hit > ColHits_Rdon;
      std::vector< recob::Hit > ColHits_Oth ;
      std::vector< recob::Hit > ColHits_Ar42;

      NTotHit = reco_hits->size();
      int colHitCount(0);
      int LoopHits = NTotHit;

      for(int hit = 0; hit < LoopHits; ++hit) {
        recob::Hit const& ThisHit = reco_hits->at(hit);
        if(ThisHit.View() == geo::kU || ThisHit.View() == geo::kV) {
          ++NIndHit;
        } else {
          ++NColHit;
        }
      }


      // Fill a set of the channels that we don't want to consider in the
      // adjacent-channel loop below. Do it here once so we don't have to
      // re-do it for every single hit
      std::set<int> badChannels;
      for(size_t i=0; i<rawDigitsVecHandle->size(); ++i) {
        int rawWireChannel=(*rawDigitsVecHandle)[i].Channel();
        std::vector<geo::WireID> adjacentwire = geo->ChannelToWire(rawWireChannel);

        if (adjacentwire.size() < 1 || adjacentwire[0].Plane == geo::kU ||
            adjacentwire[0].Plane == geo::kV){
          badChannels.insert(rawWireChannel);
        }
      }
      // std::cout << "Inserted " << badChannels.size() << " out of " << rawDigitsVecHandle->size() << " channels into set" << std::endl;

      auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
      for(int hit = 0; hit < LoopHits; ++hit) {
        recob::Hit const& ThisHit = reco_hits->at(hit);

        if (ThisHit.View() == 2) {
          std::vector<sim::TrackIDE> ThisHitIDE;
          //GETTING HOLD OF THE SIM::IDEs.

          std::vector<const sim::IDE*> ThisSimIDE;
          try {
            // HitToTrackIDEs opens a specific window around the hit. I want a
            // wider one, because the filtering can delay the hit. So this bit
            // is a copy of HitToTrackIDEs from the backtracker, with some
            // modification
            const double start = ThisHit.PeakTime()-20;
            const double end   = ThisHit.PeakTime()+ThisHit.RMS()+20;
            ThisHitIDE = bt_serv->ChannelToTrackIDEs(clockData, ThisHit.Channel(), start, end);

            // ThisHitIDE = bt_serv->HitToTrackIDEs(clockData,  ThisHit );
          } catch(...){
            // std::cout << "FIRST CATCH" << std::endl;
            firstCatch++;
            try {
              ThisSimIDE = bt_serv->HitToSimIDEs_Ps(clockData, ThisHit);
            } catch(...) {
               // std::cout << "SECOND CATCH" << std::endl;
              secondCatch++;
              // continue;
            }
            // continue;
          }

          // Get the simIDEs.
          try {
            ThisSimIDE = bt_serv->HitToSimIDEs_Ps(clockData, ThisHit);
          } catch(...) {
            // std::cout << "THIRD CATCH" << std::endl;
            thirdCatch++;
            // continue;
          }

          Hit_View.push_back(ThisHit.View());
          Hit_Size.push_back(ThisHit.EndTick() - ThisHit.StartTick());
          Hit_TPC .push_back(ThisHit.WireID().TPC);
          int channel = ThisHit.Channel();
          Hit_Chan.push_back(channel);

          if(fSaveNeighbourADCs)
            SaveNeighbourADC(channel,rawDigitsVecHandle, badChannels, ThisHit);

          auto& wgeo = geo->WireIDToWireGeo(ThisHit.WireID());
          auto const wire_start = wgeo.GetStart();
          auto const wire_end = wgeo.GetEnd();
          Hit_X_start.push_back(wire_start.X());
          Hit_Y_start.push_back(wire_start.Y());
          Hit_Z_start.push_back(wire_start.Z());
          Hit_X_end  .push_back(wire_end.X());
          Hit_Y_end  .push_back(wire_end.Y());
          Hit_Z_end  .push_back(wire_end.Z());
          Hit_Time   .push_back(ThisHit.PeakTime());
          Hit_RMS    .push_back(ThisHit.RMS());
          Hit_SADC   .push_back(ThisHit.SummedADC());
          Hit_Int    .push_back(ThisHit.Integral());
          Hit_Peak   .push_back(ThisHit.PeakAmplitude());
          Hit_True_nIDEs.push_back(ThisHitIDE.size());

          if(ThisHitIDE.size()==0)
            NHitNoBT++;

          //WHICH PARTICLE CONTRIBUTED MOST TO THE HIT.
          double TopEFrac = -DBL_MAX;

          Hit_True_EvEnergy.push_back(0);
          Hit_True_MainTrID.push_back(-1);
          for (size_t ideL=0; ideL < ThisHitIDE.size(); ++ideL) {
            Hit_True_TrackID.push_back(ThisHitIDE[ideL].trackID);
            for (size_t ipart=0; ipart<allTruthParts.size(); ++ipart) {

              if (allTruthParts[ipart].TrackId() == ThisHitIDE[ideL].trackID) {
                Hit_True_EvEnergy.at(colHitCount) += allTruthParts[ipart].E();
              }
            }
            if (ThisHitIDE[ideL].energyFrac > TopEFrac) {
              TopEFrac = ThisHitIDE[ideL].energyFrac;
              Hit_True_MainTrID.at(colHitCount) = ThisHitIDE[ideL].trackID;

            }
          }

          PType ThisPType = WhichParType(Hit_True_MainTrID.at(colHitCount));
          Hit_True_GenType.push_back(ThisPType);

          int thisMarleyIndex=-1;
          int MainTrID=Hit_True_MainTrID.at(colHitCount);
          if(ThisPType==kMarl && MainTrID!=0){
            auto const it=trkIDToMarleyIndex.find(MainTrID);
            if(it==trkIDToMarleyIndex.end()){
              mf::LogDebug(fname) << "Track ID " << MainTrID << " is not in Marley index map" << std::endl;
            }
            else{
              thisMarleyIndex=it->second;
            }
          }
          Hit_True_MarleyIndex.push_back(thisMarleyIndex);

          if(Hit_True_MainTrID[colHitCount] == -1)
          {
            Hit_True_X     .push_back(-1);
            Hit_True_Y     .push_back(-1);
            Hit_True_Z     .push_back(-1);
            Hit_True_Energy.push_back(-1);
            Hit_True_nElec .push_back(-1);
          }
          else
          {
            for(unsigned int i = 0; i < ThisSimIDE.size(); i++)
            {
              if(ThisSimIDE.at(i)->trackID==Hit_True_MainTrID[colHitCount])
              {
                Hit_True_X     .push_back(ThisSimIDE.at(i)->x           );
                Hit_True_Y     .push_back(ThisSimIDE.at(i)->y           );
                Hit_True_Z     .push_back(ThisSimIDE.at(i)->z           );
                Hit_True_Energy.push_back(ThisSimIDE.at(i)->energy      );
                Hit_True_nElec .push_back(ThisSimIDE.at(i)->numElectrons);
                break;
              }
            }
          }

          if      (ThisPType == 0) { ColHits_Oth .push_back( ThisHit ); }
          else if (ThisPType == 1) { ColHits_Marl.push_back( ThisHit ); }
          else if (ThisPType == 2) { ColHits_APA .push_back( ThisHit ); }
          else if (ThisPType == 3) { ColHits_CPA .push_back( ThisHit ); }
          else if (ThisPType == 4) { ColHits_Ar39.push_back( ThisHit ); }
          else if (ThisPType == 5) { ColHits_Neut.push_back( ThisHit ); }
          else if (ThisPType == 6) { ColHits_Kryp.push_back( ThisHit ); }
          else if (ThisPType == 7) { ColHits_Plon.push_back( ThisHit ); }
          else if (ThisPType == 8) { ColHits_Rdon.push_back( ThisHit ); }
          else if (ThisPType == 9) { ColHits_Ar42.push_back( ThisHit ); }

          colHitCount++;
        }
      }
      mf::LogInfo(fname) << "Total of:\n"
                         << " - Other: " << ColHits_Oth.size() << " col plane hits\n"
                         << " - Marl : " << TotGen_Marl << " true parts\t| " << ColHits_Marl.size() << " col plane hits\n"
                         << " - APA  : " << TotGen_APA  << " true parts\t| " << ColHits_APA .size() << " col plane hits\n"
                         << " - CPA  : " << TotGen_CPA  << " true parts\t| " << ColHits_CPA .size() << " col plane hits\n"
                         << " - Ar39 : " << TotGen_Ar39 << " true parts\t| " << ColHits_Ar39.size() << " col plane hits\n"
                         << " - Neut : " << TotGen_Neut << " true parts\t| " << ColHits_Neut.size() << " col plane hits\n"
                         << " - Kryp : " << TotGen_Kryp << " true parts\t| " << ColHits_Kryp.size() << " col plane hits\n"
                         << " - Plon : " << TotGen_Plon << " true parts\t| " << ColHits_Plon.size() << " col plane hits\n"
                         << " - Rdon : " << TotGen_Rdon << " true parts\t| " << ColHits_Rdon.size() << " col plane hits\n"
                         << " - Ar42 : " << TotGen_Ar42 << " true parts\t| " << ColHits_Ar42.size() << " col plane hits\n";
    } else {
      mf::LogError(fname) << "Requested to save wire hits, but cannot load any wire hits";
      throw art::Exception(art::errors::NotFound) << "Requested to save wire hits, but cannot load any wire hits\n";
    }
  }

  if (fSavePDS) {

    std::vector<art::Ptr<recob::OpHit> > ophitlist;
    std::map<PType, std::vector<art::Ptr<recob::OpHit> > > map_of_ophit;

    auto OpHitHandle = evt.getHandle< std::vector< recob::OpHit > >(fOpHitModuleLabel);
    if (OpHitHandle) {
      art::fill_ptr_vector(ophitlist, OpHitHandle);

      mf::LogDebug(fname) << "There are " << ophitlist.size() << " optical hits in the event." << std::endl;

      for(size_t i = 0; i < ophitlist.size(); ++i)
      {
        std::vector<sim::TrackSDP> vec_tracksdp = pbt_serv->OpHitToTrackSDPs(ophitlist.at(i));
        PType gen = kUnknown;

        std::sort(vec_tracksdp.begin(), vec_tracksdp.end(),
                  [](const sim::TrackSDP& a, const sim::TrackSDP& b) -> bool { return a.energyFrac > b.energyFrac; });

        for (size_t iSDP=0; iSDP<vec_tracksdp.size(); ++iSDP) {
          PDS_OpHit_True_TrackIDAll.push_back(vec_tracksdp[iSDP].trackID);
          PDS_OpHit_True_GenTypeAll.push_back(WhichParType(vec_tracksdp[iSDP].trackID));
          PDS_OpHit_True_EnergyAll .push_back(vec_tracksdp[iSDP].energy);
          PDS_OpHit_True_IndexAll  .push_back((int)i);
        }

        if (vec_tracksdp.size()>0){
          int MainTrID = vec_tracksdp[0].trackID;
          PDS_OpHit_True_TrackID.push_back(vec_tracksdp[0].trackID);
          gen = WhichParType(vec_tracksdp[0].trackID);
          PDS_OpHit_True_GenType.push_back(gen);
          int thisMarleyIndex;
          if(gen==kMarl && MainTrID!=0){
            auto const it=trkIDToMarleyIndex.find(MainTrID);
            if(it==trkIDToMarleyIndex.end()){
              mf::LogDebug(fname) << "Track ID " << MainTrID << " is not in Marley index map" << std::endl;
            }
            else{
              thisMarleyIndex=it->second;
            }
          }
          PDS_OpHit_True_Index  .push_back(thisMarleyIndex);
          PDS_OpHit_True_Energy .push_back(vec_tracksdp[0].energy);
        } else {
          PDS_OpHit_True_TrackID.push_back(-1);
          PDS_OpHit_True_GenType.push_back(kUnknown);
          PDS_OpHit_True_Energy .push_back(-1);
        }

        map_of_ophit[gen].push_back(ophitlist.at(i));

        auto const xyz_world = geo->OpDetGeoFromOpChannel(ophitlist[i]->OpChannel()).GetCenter();
        PDS_OpHit_OpChannel   .push_back(ophitlist[i]->OpChannel());
        PDS_OpHit_X           .push_back(xyz_world.X());
        PDS_OpHit_Y           .push_back(xyz_world.Y());
        PDS_OpHit_Z           .push_back(xyz_world.Z());
        PDS_OpHit_PeakTimeAbs .push_back(ophitlist[i]->PeakTimeAbs());
        PDS_OpHit_PeakTime    .push_back(ophitlist[i]->PeakTime());
        PDS_OpHit_Frame       .push_back(ophitlist[i]->Frame());
        PDS_OpHit_Width       .push_back(ophitlist[i]->Width());
        PDS_OpHit_Area        .push_back(ophitlist[i]->Area());
        PDS_OpHit_Amplitude   .push_back(ophitlist[i]->Amplitude());
        PDS_OpHit_PE          .push_back(ophitlist[i]->PE());
        PDS_OpHit_FastToTotal .push_back(ophitlist[i]->FastToTotal());
      }
      mf::LogInfo(fname) << "Total of :\n"
                         << " - Other: " << map_of_ophit[kUnknown].size() << " opt hits\n" 
                         << " - Marl : " << TotGen_Marl << " true parts\t| " << map_of_ophit[kMarl].size() << " opt hits\n"
                         << " - APA  : " << TotGen_APA  << " true parts\t| " << map_of_ophit[kAPA] .size() << " opt hits\n"
                         << " - CPA  : " << TotGen_CPA  << " true parts\t| " << map_of_ophit[kCPA] .size() << " opt hits\n"
                         << " - Ar39 : " << TotGen_Ar39 << " true parts\t| " << map_of_ophit[kAr39].size() << " opt hits\n"
                         << " - Neut : " << TotGen_Neut << " true parts\t| " << map_of_ophit[kNeut].size() << " opt hits\n"
                         << " - Kryp : " << TotGen_Kryp << " true parts\t| " << map_of_ophit[kKryp].size() << " opt hits\n"
                         << " - Plon : " << TotGen_Plon << " true parts\t| " << map_of_ophit[kPlon].size() << " opt hits\n"
                         << " - Rdon : " << TotGen_Rdon << " true parts\t| " << map_of_ophit[kRdon].size() << " opt hits\n"
                         << " - Ar42 : " << TotGen_Ar42 << " true parts\t| " << map_of_ophit[kAr42].size() << " opt hits\n";
    }
    else {
      mf::LogError(fname) << "Requested to save optical hits, but cannot load any ophits";
      throw art::Exception(art::errors::NotFound) << "Requested to save optical hits, but cannot load any optical hits\n";
    }
  }

  if(fSaveIDEs) SaveIDEs(evt);


  fSNAnaTree->Fill();
}

void SNAna::FillTruth(const art::FindManyP<simb::MCParticle> Assn,
                      const art::Handle<std::vector<simb::MCTruth>>& Hand,
                      const PType type) {

  for(size_t i1=0; i1<Hand->size(); ++i1) {
    for ( size_t i2=0; i2 < Assn.at(i1).size(); ++i2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(i1).at(i2));
      True_Bck_Mode      .push_back(type);
      True_Bck_PDG       .push_back(ThisPar.PdgCode   ());
      True_Bck_ID        .push_back(ThisPar.TrackId   ());
      True_Bck_Mother    .push_back(ThisPar.Mother    ());
      True_Bck_P         .push_back(ThisPar.P         ());
      True_Bck_VertX     .push_back(ThisPar.Vx        ());
      True_Bck_VertY     .push_back(ThisPar.Vy        ());
      True_Bck_VertZ     .push_back(ThisPar.Vz        ());
      True_Bck_Time      .push_back(ThisPar.T         ());
      True_Bck_Energy    .push_back(ThisPar.E         ());
      True_Bck_EndX      .push_back(ThisPar.EndX      ());
      True_Bck_EndY      .push_back(ThisPar.EndY      ());
      True_Bck_EndZ      .push_back(ThisPar.EndZ      ());
      True_Bck_EndT      .push_back(ThisPar.EndT      ());
      True_Bck_EndE      .push_back(ThisPar.EndE      ());
      True_Bck_EndProcess.push_back(ThisPar.EndProcess());
      True_Bck_Process   .push_back(ThisPar.Process   ());
    }
  }
}

void SNAna::FillMyMaps(std::map< int, simb::MCParticle> &MyMap,
                       art::FindManyP<simb::MCParticle> Assn,
                       art::Handle< std::vector<simb::MCTruth> > Hand,
                       std::map<int, int>* indexMap) {
  for ( size_t L1=0; L1 < Hand->size(); ++L1 ) {
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
      MyMap[ThisPar.TrackId()] = ThisPar;
      if(indexMap) indexMap->insert({ThisPar.TrackId(), L1});
    }
  }
  return;
}


PType SNAna::WhichParType( int TrID )
{
  PType ThisPType = kUnknown;
  auto const& it=trkIDToPType.find(TrID);
  if(it!=trkIDToPType.end()){
    ThisPType=it->second;
  }
  return ThisPType;

}


bool SNAna::InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap )
{
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find( TrID );
  if ( ParIt != ParMap.end() ) {
    return true;
  } else
    return false;
}

//......................................................
void SNAna::SaveIDEs(art::Event const & evt)
{
  auto allParticles = evt.getValidHandle<std::vector<simb::MCParticle> >(fGEANTLabel);
  art::FindMany<simb::MCTruth> assn(allParticles,evt,fGEANTLabel);
  std::map<int, const simb::MCTruth*> idToTruth;
  for(size_t i=0; i<allParticles->size(); ++i){
    const simb::MCParticle& particle=allParticles->at(i);
    const std::vector<const simb::MCTruth*> truths=assn.at(i);
    if(truths.size()==1){
      idToTruth[particle.TrackId()]=truths[0];
    }
    else{
      mf::LogDebug("DAQSimAna") << "Particle " << particle.TrackId() << " has " << truths.size() << " truths";
      idToTruth[particle.TrackId()]=nullptr;
    }
  }

  // Get the SimChannels so we can see where the actual energy depositions were
  auto& simchs=*evt.getValidHandle<std::vector<sim::SimChannel>>("largeant");

  for(auto&& simch: simchs){
    // We only care about collection channels
    if(geo->SignalType(simch.Channel())!=geo::kCollection) continue;

    // The IDEs record energy depositions at every tick, but
    // mostly we have several depositions at contiguous times. So
    // we're going to save just one output IDE for each contiguous
    // block of hits on a channel. Each item in vector is a list
    // of (TDC, IDE*) for contiguous-in-time IDEs
    std::vector<std::vector<std::pair<int, const sim::IDE*> > > contigIDEs;
    int prevTDC=0;
    for (const auto& TDCinfo: simch.TDCIDEMap()) {
      // Do we need to start a new group of IDEs? Yes if this is
      // the first IDE in this channel. Yes if this IDE is not
      // contiguous with the previous one
      if(contigIDEs.empty() || TDCinfo.first-prevTDC>5){
        contigIDEs.push_back(std::vector<std::pair<int, const sim::IDE*> >());
      }
      std::vector<std::pair<int, const sim::IDE*> >& currentIDEs=contigIDEs.back();

      // Add all the current tick's IDEs to the list
      for (const sim::IDE& ide: TDCinfo.second) {
        currentIDEs.push_back(std::make_pair(TDCinfo.first, &ide));
      }
      prevTDC=TDCinfo.first;
    }

    for(auto const& contigs : contigIDEs){
      float energy=0;
      float electrons=0;
      int startTime=99999;
      int endTime=0;
      std::map<PType, float> ptypeToEnergy;
      for(auto const& timeide : contigs){
        const int tdc=timeide.first;
        startTime=std::min(tdc, startTime);
        endTime=std::max(tdc, endTime);
        const sim::IDE& ide=*timeide.second;
        const float thisEnergy=ide.energy;
        const PType thisPType=WhichParType(std::abs(ide.trackID));
        energy+=thisEnergy;
        electrons+=ide.numElectrons;
        ptypeToEnergy[thisPType]+=thisEnergy;
      }
      float bestEnergy=0;
      PType bestPType=kUnknown;
      for(auto const& it : ptypeToEnergy){
        if(it.second>bestEnergy){
          bestEnergy=it.second;
          bestPType=it.first;
        }
      }
      // Ignore anything past the end of the readout window
      if(startTime<4492){
        IDEChannel.push_back(simch.Channel());
        IDEStartTime.push_back(startTime);
        IDEEndTime.push_back(endTime);
        IDEEnergy.push_back(energy);
        IDEElectrons.push_back(electrons);
        IDEParticle.push_back(bestPType);
      }
    } // loop over our compressed IDEs
  } // loop over SimChannels
}


void SNAna::SaveNeighbourADC(int channel,
                             art::Handle< std::vector<raw::RawDigit> >rawDigitsVecHandle,
                             std::set<int>badChannels,
                             recob::Hit const& ThisHit) {

  Hit_AdjM5SADC.push_back(0);
  Hit_AdjM2SADC.push_back(0);
  Hit_AdjM1SADC.push_back(0);
  Hit_AdjP1SADC.push_back(0);
  Hit_AdjP2SADC.push_back(0);
  Hit_AdjP5SADC.push_back(0);
  Hit_AdjM5Chan.push_back(0);
  Hit_AdjM2Chan.push_back(0);
  Hit_AdjM1Chan.push_back(0);
  Hit_AdjP1Chan.push_back(0);
  Hit_AdjP2Chan.push_back(0);
  Hit_AdjP5Chan.push_back(0);
  int colHitCount = Hit_AdjM1Chan.size()-1;

  // A vector for us to uncompress the rawdigits into. Create it
  // outside the loop here so that we only have to allocate it once
  raw::RawDigit::ADCvector_t ADCs((*rawDigitsVecHandle)[0].Samples());

  std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
  for(size_t i=0; i<rawDigitsVecHandle->size(); ++i)
  {
    int rawWireChannel=(*rawDigitsVecHandle)[i].Channel();
    const int chanDiff=rawWireChannel-channel;
    if(abs(chanDiff)!=1 && abs(chanDiff)!=2 && abs(chanDiff)!=5) continue;

    if(badChannels.find(rawWireChannel)!=badChannels.end()) continue;

    switch(chanDiff)
    {
    case -5:
      Hit_AdjM5SADC[colHitCount] = 0;
      Hit_AdjM5Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjM5SADC[colHitCount]+=ADCs[i];
      break;
    case -2:
      Hit_AdjM2SADC[colHitCount] = 0;
      Hit_AdjM2Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjM2SADC[colHitCount]+=ADCs[i];
      break;
    case -1:
      Hit_AdjM1SADC[colHitCount] = 0;
      Hit_AdjM1Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjM1SADC[colHitCount]+=ADCs[i];
      break;
    case  1:
      Hit_AdjP1SADC[colHitCount] = 0;
      Hit_AdjP1Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjP1SADC[colHitCount]+=ADCs[i];
      break;
    case  2:
      Hit_AdjP2SADC[colHitCount] = 0;
      Hit_AdjP2Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjP2SADC[colHitCount]+=ADCs[i];
      break;
    case  5:
      Hit_AdjP5SADC[colHitCount] = 0;
      Hit_AdjP5Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjP5SADC[colHitCount]+=ADCs[i];
      break;
    default:
      break;
    }
  }
}


void SNAna::endJob()
{
  mf::LogDebug(fname) << firstCatch << " " << secondCatch << " " << thirdCatch << std::endl;
}

DEFINE_ART_MODULE(SNAna)
