// Christie Chiu and Ben Jones, MIT, 2012
//
// This is an analyzer module which writes the raw optical
// detector pulses for each PMT to an output file
//

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TTree.h"

// C++ Includes
#include "math.h"
#include <iostream>

// LArSoft includes
#include "lardataobj/RecoBase/OpFlash.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace opdet {

  class OpFlashMCTruthAna : public art::EDAnalyzer {
  public:
    // Standard constructor and destructor for an ART module.
    OpFlashMCTruthAna(const fhicl::ParameterSet&);

    // The analyzer routine, called once per event.
    void analyze(const art::Event&);

  private:
    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fFlashInputModule;
    std::string fTruthInputModule;

    TTree* fAnalysisTree;
    TTree* fPerEventTree;

    Int_t fEventID, fNFlashes, fNTruths;
    Float_t fFlashY, fFlashZ, fFlashU, fFlashV, fFlashT;
    Float_t fFlashPE, fFlashFastToTotal;
    Float_t fFlashWidthY, fFlashWidthZ, fFlashWidthU, fFlashWidthV;
    Float_t fVertexX, fVertexY, fVertexZ /* , fVertexU, fVertexV */, fVertexT;
    Float_t fTrueE;
    Int_t fTruePDG;
    Float_t fCenterX, fCenterY, fCenterZ;
    Float_t fDistFlashCenter, fDistFlashVertex;
    Float_t fDistFlashCenterNorm, fDistFlashVertexNorm;
  };

}

// OpFlashMCTruthAna_module.cc

// This is a short program required by the ART framework.  It enables
// a program (OpFlashMCTruthAna, in this case) to be called as a module
// from a .fcl file. It is unlikely that you'll have to make any
// changes to this file.

namespace opdet {
  DEFINE_ART_MODULE(OpFlashMCTruthAna)
}

// OpFlashMCTruthAna.cxx

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpFlashMCTruthAna::OpFlashMCTruthAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fFlashInputModule = pset.get<std::string>("FlashInputModule");
    fTruthInputModule = pset.get<std::string>("TruthInputModule");

    art::ServiceHandle<art::TFileService const> tfs;

    fPerEventTree = tfs->make<TTree>("PerEventTree", "PerEventTree");

    fPerEventTree->Branch("EventID", &fEventID, "EventID/I");
    fPerEventTree->Branch("NFlashes", &fNFlashes, "NFlashes/I");

    fPerEventTree->Branch("VertexX", &fVertexX, "VertexX/F");
    fPerEventTree->Branch("VertexY", &fVertexY, "VertexY/F");
    fPerEventTree->Branch("VertexZ", &fVertexZ, "VertexZ/F");

    fPerEventTree->Branch("TrueE", &fTrueE, "TrueE/F");
    fPerEventTree->Branch("TruePDG", &fTruePDG, "TruePDG/I");

    fPerEventTree->Branch("CenterX", &fCenterX, "CenterX/F");
    fPerEventTree->Branch("CenterY", &fCenterY, "CenterY/F");
    fPerEventTree->Branch("CenterZ", &fCenterZ, "CenterZ/F");

    fAnalysisTree = tfs->make<TTree>("AnalysisTree", "AnalysisTree");

    fAnalysisTree->Branch("EventID", &fEventID, "EventID/I");
    fAnalysisTree->Branch("NFlashes", &fNFlashes, "NFlashes/I");

    fAnalysisTree->Branch("FlashY", &fFlashY, "FlashY/F");
    fAnalysisTree->Branch("FlashZ", &fFlashZ, "FlashZ/F");
    fAnalysisTree->Branch("FlashU", &fFlashU, "FlashU/F");
    fAnalysisTree->Branch("FlashV", &fFlashV, "FlashV/F");

    fAnalysisTree->Branch("FlashWidthY", &fFlashWidthY, "FlashWidthY/F");
    fAnalysisTree->Branch("FlashWidthZ", &fFlashWidthZ, "FlashWidthZ/F");
    fAnalysisTree->Branch("FlashWidthU", &fFlashWidthU, "FlashWidthU/F");
    fAnalysisTree->Branch("FlashWidthV", &fFlashWidthV, "FlashWidthV/F");

    fAnalysisTree->Branch("FlashT", &fFlashT, "FlashT/F");
    fAnalysisTree->Branch("FlashPE", &fFlashPE, "FlashPE/F");
    fAnalysisTree->Branch("FlashFastToTotal", &fFlashFastToTotal, "FlashFastToTotal/F");

    fAnalysisTree->Branch("VertexX", &fVertexX, "VertexX/F");
    fAnalysisTree->Branch("VertexY", &fVertexY, "VertexY/F");
    fAnalysisTree->Branch("VertexZ", &fVertexZ, "VertexZ/F");

    fAnalysisTree->Branch("TrueE", &fTrueE, "TrueE/F");
    fAnalysisTree->Branch("TruePDG", &fTruePDG, "TruePDG/I");

    fAnalysisTree->Branch("CenterX", &fCenterX, "CenterX/F");
    fAnalysisTree->Branch("CenterY", &fCenterY, "CenterY/F");
    fAnalysisTree->Branch("CenterZ", &fCenterZ, "CenterZ/F");

    fAnalysisTree->Branch("DistFlashCenter", &fDistFlashCenter, "DistFlashCenter/F");
    fAnalysisTree->Branch("DistFlashVertex", &fDistFlashVertex, "DistFlashVertex/F");
  }

  //-----------------------------------------------------------------------
  void OpFlashMCTruthAna::analyze(const art::Event& evt)
  {

    // Create a handles
    art::Handle<std::vector<recob::OpFlash>> FlashHandle;
    art::Handle<std::vector<simb::MCTruth>> TruthHandle;

    // Read in data
    evt.getByLabel(fFlashInputModule, FlashHandle);
    evt.getByLabel(fTruthInputModule, TruthHandle);

    fEventID = evt.id().event();

    fNTruths = TruthHandle->at(0).NParticles();
    fNFlashes = FlashHandle->size();

    std::cout << "Size of truth collection : " << TruthHandle->size() << std::endl;
    std::cout << "We found " << fNTruths << " truth particles and " << fNFlashes << " flashes"
              << std::endl;

    for (int iPart = 0; iPart != fNTruths; ++iPart) {
      const simb::MCParticle ThisPart = TruthHandle->at(0).GetParticle(iPart);
      fTruePDG = ThisPart.PdgCode();
      fTrueE = ThisPart.E(0);

      fVertexX = ThisPart.Vx(0);
      fVertexY = ThisPart.Vy(0);
      fVertexZ = ThisPart.Vz(0);
      fVertexT = ThisPart.T(0);

      fCenterX = (ThisPart.Vx(0) + ThisPart.EndX()) / 2.;
      fCenterY = (ThisPart.Vy(0) + ThisPart.EndY()) / 2.;
      fCenterZ = (ThisPart.Vz(0) + ThisPart.EndZ()) / 2.;

      for (unsigned int i = 0; i < FlashHandle->size(); ++i) {

        // Get OpFlash
        art::Ptr<recob::OpFlash> TheFlashPtr(FlashHandle, i);

        fFlashT = TheFlashPtr->Time();
        fFlashPE = TheFlashPtr->TotalPE();
        fFlashFastToTotal = TheFlashPtr->FastToTotal();

        fFlashY = TheFlashPtr->YCenter();
        fFlashZ = TheFlashPtr->ZCenter();
        fFlashU = TheFlashPtr->WireCenters().at(0);
        fFlashV = TheFlashPtr->WireCenters().at(1);

        fFlashWidthY = TheFlashPtr->YWidth();
        fFlashWidthZ = TheFlashPtr->ZWidth();
        fFlashWidthU = TheFlashPtr->WireWidths().at(0);
        fFlashWidthV = TheFlashPtr->WireWidths().at(1);

        fDistFlashCenter = pow(pow(fCenterY - fFlashY, 2) + pow(fCenterZ - fFlashZ, 2), 0.5);

        fDistFlashVertex = pow(pow(fVertexY - fFlashY, 2) + pow(fVertexZ - fFlashZ, 2), 0.5);

        fDistFlashCenterNorm = pow(pow((fCenterY - fFlashY) / fFlashWidthY, 2) +
                                     pow((fCenterZ - fFlashZ) / fFlashWidthZ, 2),
                                   0.5);

        fDistFlashVertexNorm = pow(pow((fVertexY - fFlashY) / fFlashWidthY, 2) +
                                     pow((fVertexZ - fFlashZ) / fFlashWidthZ, 2),
                                   0.5);

        fAnalysisTree->Fill();
      }
    }

    fPerEventTree->Fill();
  }

} // namespace opdet
