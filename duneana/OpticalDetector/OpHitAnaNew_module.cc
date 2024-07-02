// -*- mode: c++; c-basic-offset: 2; -*-
// This analyzer writes out a TTree containing the properties of
// each reconstructed hit
//

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

// C++ includes
#include "math.h"
#include <cstring>

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/OpHit.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

namespace opdet {

  class OpHitAnaNew : public art::EDAnalyzer {
  public:
    // Standard constructor and destructor for an ART module.
    OpHitAnaNew(const fhicl::ParameterSet&);

    // The analyzer routine, called once per event.
    void analyze(const art::Event&);

  private:
    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fOpHitModuleLabel;   // Input tag for OpHit collection
    float fSampleFreq;               // in MHz
    float fTimeBegin;                // in us
    float fTimeEnd;                  // in us

    float fYMin, fYMax, fZMin, fZMax;

    int PosHistYRes, PosHistZRes;

    TTree* fPerOpHitTree;

    Int_t fEventID;
    Int_t fHitID;
    Double_t fAbsTime;
    bool fInBeamFrame;
    int fOnBeamTime;
    Float_t fTotalPE;

    Float_t fNPe;
    Float_t fYCenter;
    Float_t fYWidth;
    Float_t fZCenter;
    Float_t fZWidth;

    Int_t fOpChannel;
    Double_t fPeakTimeAbs;
    Double_t fPeakTime;
    Int_t fFrame;
    Float_t fWidth;
    Float_t fArea;
    Float_t fAmplitude;
    Float_t fPE;
    Float_t fFastToTotal;
  };

}

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpHitAnaNew::OpHitAnaNew(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fOpHitModuleLabel = pset.get<std::string>("OpHitModuleLabel");

    //    art::ServiceHandle<OpDigiProperties const> odp;
    //    fTimeBegin  = odp->TimeBegin();
    //    fTimeEnd    = odp->TimeEnd();
    //    fSampleFreq = odp->SampleFreq();

    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fTimeBegin = clock_data.OpticalClock().Time();
    fTimeEnd = clock_data.OpticalClock().FramePeriod();
    fSampleFreq = clock_data.OpticalClock().Frequency();

    fYMin = pset.get<float>("YMin");
    fYMax = pset.get<float>("YMax");
    fZMin = pset.get<float>("ZMin");
    fZMax = pset.get<float>("ZMax");

    PosHistYRes = 100;
    PosHistZRes = 100;

    art::ServiceHandle<art::TFileService const> tfs;

    fPerOpHitTree = tfs->make<TTree>("PerOpHitTree", "PerOpHitTree");
    fPerOpHitTree->Branch("EventID", &fEventID, "EventID/I");
    fPerOpHitTree->Branch("HitID", &fHitID, "HitID/I");
    fPerOpHitTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
    fPerOpHitTree->Branch("PeakTimeAbs", &fPeakTimeAbs, "PeakTimeAbs/D");
    fPerOpHitTree->Branch("PeakTime", &fPeakTime, "PeakTime/D");
    fPerOpHitTree->Branch("Frame", &fFrame, "Frame/I");
    fPerOpHitTree->Branch("Width", &fWidth, "Width/F");
    fPerOpHitTree->Branch("Area", &fArea, "Area/F");
    fPerOpHitTree->Branch("Amplitude", &fAmplitude, "Amplitude/F");
    fPerOpHitTree->Branch("PE", &fPE, "PE/F");
    fPerOpHitTree->Branch("FastToTotal", &fFastToTotal, "FastToTotal/F");

  }

  //-----------------------------------------------------------------------
  void OpHitAnaNew::analyze(const art::Event& evt)
  {
    // Create string for histogram name
    // char HistName[50];

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us.
    art::ServiceHandle<art::TFileService const> tfs;

    fEventID = evt.id().event();

    art::ServiceHandle<geo::Geometry const> geom;
    // unsigned int NOpChannels = geom->NOpChannels();

    art::Handle<std::vector<recob::OpHit>> OpHitHandle;
    evt.getByLabel(fOpHitModuleLabel, OpHitHandle);

    for (size_t i = 0; i != OpHitHandle->size(); ++i) {
    fOpChannel = OpHitHandle->at(i).OpChannel();
    fPeakTimeAbs = OpHitHandle->at(i).PeakTimeAbs();
    fPeakTime = OpHitHandle->at(i).PeakTime();
    fFrame = OpHitHandle->at(i).Frame();
    fWidth = OpHitHandle->at(i).Width();
    fArea = OpHitHandle->at(i).Area();
    fAmplitude = OpHitHandle->at(i).Amplitude();
    fPE = OpHitHandle->at(i).PE();
    fFastToTotal = OpHitHandle->at(i).FastToTotal();
    fHitID = i;
    fPerOpHitTree->Fill();
    }
    
  }

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(OpHitAnaNew)
}
