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

#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/sim.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

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
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


class DataHit : public art::EDAnalyzer {

public:
    explicit DataHit(fhicl::ParameterSet const & p);

    DataHit(DataHit const &) = delete;
    DataHit(DataHit &&) = delete;
    DataHit & operator = (DataHit const &) = delete;
    DataHit & operator = (DataHit &&) = delete;

    void analyze(art::Event const & evt) override;
    void reconfigure(fhicl::ParameterSet const & p);
    void beginJob() override;
    void endJob() override;

private:

    void ResetVariables();

    // void FillMyMaps  ( std::map< int, simb::MCParticle> &MyMap, 
    //                    art::FindManyP<simb::MCParticle> Assn,
    //                    art::ValidHandle< std::vector<simb::MCTruth> > Hand,
    //                    std::map<int, int>* indexMap=nullptr);
    // PType WhichParType( int TrID );
    // bool  InMyMap     ( int TrID, std::map< int, simb::MCParticle> ParMap );
    // void FillTruth(const art::FindManyP<simb::MCParticle> Assn,
    //                const art::ValidHandle<std::vector<simb::MCTruth>>& Hand,
    //                const PType type) ;

    // void SaveIDEs(art::Event const & evt);

    // int firstCatch;
    // int secondCatch;
    // int thirdCatch;

    // std::string fRawDigitLabel;
    std::string fHitLabel;
    // std::string fCalDataModuleLabel;
    // std::string fOpFlashModuleLabel;
    // std::string fOpHitModuleLabel;

    // std::string fGEANTLabel;
    // std::string fMARLLabel ; std::map< int, simb::MCParticle > MarlParts;
    // std::string fAPALabel  ; std::map< int, simb::MCParticle > APAParts ;
    // std::string fCPALabel  ; std::map< int, simb::MCParticle > CPAParts ;
    // std::string fAr39Label ; std::map< int, simb::MCParticle > Ar39Parts;
    // std::string fNeutLabel ; std::map< int, simb::MCParticle > NeutParts;
    // std::string fKrypLabel ; std::map< int, simb::MCParticle > KrypParts;
    // std::string fPlonLabel ; std::map< int, simb::MCParticle > PlonParts;
    // std::string fRdonLabel ; std::map< int, simb::MCParticle > RdonParts;
    // std::string fAr42Label ; std::map< int, simb::MCParticle > Ar42Parts;
    // std::map<int, const simb::MCParticle*> truthmap;
    // // Which MARLEY interaction (if any) caused this true track ID?
    // std::map<int, int> trkIDToMarleyIndex;

    // // Mapping from track ID to particle type, for use in WhichParType()
    // std::map<int, PType> trkIDToPType;

    // bool fSaveNeighbourADCs;
    // bool fSaveIDEs;
  
    TTree* fDataHitTree;

    int Run;
    int SubRun;
    int Event;
  
    int NTotHit   ;
    int NColHit   ;
    int NIndHit   ;

    std::vector<int>                  Hit_View                 ;
    std::vector<int>                  Hit_Size                 ;
    std::vector<int>                  Hit_TPC                  ;
    std::vector<int>                  Hit_Chan                 ;
    std::vector<float>                Hit_Time                 ;
    std::vector<float>                Hit_RMS                  ;
    std::vector<float>                Hit_SADC                 ;
    std::vector<float>                Hit_Int                  ;
    std::vector<float>                Hit_Peak                 ;

    art::ServiceHandle<geo::Geometry> geo;
};

DataHit::DataHit(fhicl::ParameterSet const & p):EDAnalyzer(p)
{
    this->reconfigure(p);
}


void DataHit::reconfigure(fhicl::ParameterSet const & p)
{
    fHitLabel           = p.get<std::string>("HitLabel"          );
} 


void DataHit::ResetVariables()
{
    Run = SubRun = Event = -1;

    NTotHit    = 0;
    NColHit    = 0;
    NIndHit    = 0;

    Hit_View                 .clear();
    Hit_Size                 .clear();
    Hit_TPC                  .clear();
    Hit_Chan                 .clear();
    Hit_Time                 .clear();
    Hit_RMS                  .clear();
    Hit_SADC                 .clear();
    Hit_Int                  .clear();
    Hit_Peak                 .clear();
}


void DataHit::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;


    fDataHitTree = tfs->make<TTree>("DataHitTree","A tree of data hits");
 
    fDataHitTree->Branch("Run"       , &Run       , "Run/I"       );
    fDataHitTree->Branch("SubRun"    , &SubRun    , "SubRun/I"    );
    fDataHitTree->Branch("Event"     , &Event     , "Event/I"     );
    fDataHitTree->Branch("NTotHit"   , &NTotHit   , "NTotHit/I"  );
    fDataHitTree->Branch("NColHit"   , &NColHit   , "NColHit/I"  );
    fDataHitTree->Branch("NIndHit"   , &NIndHit   , "NIndHit/I"  );

    fDataHitTree->Branch("Hit_View"                 , &Hit_View                 );
    fDataHitTree->Branch("Hit_Size"                 , &Hit_Size                 );
    fDataHitTree->Branch("Hit_TPC"                  , &Hit_TPC                  );
    fDataHitTree->Branch("Hit_Chan"                 , &Hit_Chan                 );
    fDataHitTree->Branch("Hit_Time"                 , &Hit_Time                 );
    fDataHitTree->Branch("Hit_RMS"                  , &Hit_RMS                  );
    fDataHitTree->Branch("Hit_SADC"                 , &Hit_SADC                 );
    fDataHitTree->Branch("Hit_Int"                  , &Hit_Int                  );
    fDataHitTree->Branch("Hit_Peak"                 , &Hit_Peak                 );
} 


void DataHit::analyze(art::Event const & evt)
{
    ResetVariables();
 
    Run    = evt.run();
    SubRun = evt.subRun();
    Event  = evt.event();

    //GET INFORMATION ABOUT THE DETECTOR'S GEOMETRY.
    // auto const* geo = lar::providerFrom<geo::Geometry>();

    //GET THE RECO HITS.
    auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(fHitLabel);

    NTotHit = reco_hits->size();
    
    for(int hit = 0; hit < NTotHit; ++hit) 
    {
        recob::Hit const& ThisHit = reco_hits->at(hit);  

        if(ThisHit.View() == geo::kU || ThisHit.View() == geo::kV) {
            ++NIndHit;
        } else { 
            ++NColHit;
        }

        Hit_View.push_back(ThisHit.View());
        Hit_Size.push_back(ThisHit.EndTick() - ThisHit.StartTick());
        Hit_TPC .push_back(ThisHit.WireID().TPC);
        int channel = ThisHit.Channel();
        Hit_Chan.push_back(channel);
        Hit_Time   .push_back(ThisHit.PeakTime());
        Hit_RMS    .push_back(ThisHit.RMS());
        Hit_SADC   .push_back(ThisHit.ROISummedADC());
        Hit_Int    .push_back(ThisHit.Integral());
        Hit_Peak   .push_back(ThisHit.PeakAmplitude());
    } 

    fDataHitTree->Fill();
} 

void DataHit::endJob()
{
}

DEFINE_ART_MODULE(DataHit)
