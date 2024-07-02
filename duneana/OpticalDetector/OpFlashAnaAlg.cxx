#include "larana/OpticalDetector/OpFlashAnaAlg.h"

#include "TTree.h"

void opdet::OpFlashAnaAlg::SetOpFlashTree(TTree* tree, bool makeOpDetPEHist)
{
  fOpFlashTree = tree;
  fMakeOpDetPEHist = makeOpDetPEHist;

  fOpFlashTree->Branch(
    fOpFlashTree->GetName(), &fOpFlashAnaStruct, fOpFlashAnaStruct.LeafList.c_str());
  if (fMakeOpDetPEHist) {
    fOpFlashAnaStruct.FlashOpDetPEHist->SetName("hOpDetPEs");
    fOpFlashAnaStruct.FlashOpDetPEHist->SetTitle("PEs per Optical Detector; OpDet Number; PE");
    fOpFlashTree->Branch("flash_OpDetPEHist", "TH1D", &(fOpFlashAnaStruct.FlashOpDetPEHist));
  }
}

void opdet::OpFlashAnaAlg::SetOpHitTree(TTree* tree)
{
  fOpHitTree = tree;
  fOpHitTree->Branch(fOpHitTree->GetName(), &fOpHitAnaStruct, fOpHitAnaStruct.LeafList.c_str());
}

void opdet::OpFlashAnaAlg::FillOpFlashes(const std::vector<recob::OpFlash>& flashVector)
{
  if (fOpFlashTree) FillOpFlashTree(flashVector);
}

void opdet::OpFlashAnaAlg::FillOpHits(const std::vector<recob::OpHit>& hitVector)
{
  if (fOpHitTree) FillOpHitTree(hitVector);
}

void opdet::OpFlashAnaAlg::FillOpFlashTree(const std::vector<recob::OpFlash>& flashVector)
{
  for (auto const& flash : flashVector) {

    fOpFlashAnaStruct.FlashTime = flash.Time();
    fOpFlashAnaStruct.FlashTimeWidth = flash.TimeWidth();
    fOpFlashAnaStruct.FlashAbsTime = flash.AbsTime();
    fOpFlashAnaStruct.FlashFrame = flash.Frame();
    fOpFlashAnaStruct.FlashY = flash.YCenter();
    fOpFlashAnaStruct.FlashYWidth = flash.YWidth();
    fOpFlashAnaStruct.FlashZ = flash.ZCenter();
    fOpFlashAnaStruct.FlashZWidth = flash.ZWidth();
    fOpFlashAnaStruct.FlashInBeamFrame = flash.InBeamFrame();
    fOpFlashAnaStruct.FlashOnBeamTime = flash.OnBeamTime();
    fOpFlashAnaStruct.FlashFastToTotal = flash.FastToTotal();
    fOpFlashAnaStruct.FlashTotalPE = flash.TotalPE();

    if (fMakeOpDetPEHist) {
      int n_opdets = flash.WireCenters().size();
      fOpFlashAnaStruct.FlashOpDetPEHist->SetBins(n_opdets, -0.5, (float)n_opdets - 0.5);
      for (int i = 0; i < n_opdets; i++)
        fOpFlashAnaStruct.FlashOpDetPEHist->SetBinContent(i + 1, flash.PE(i));
    }

    fOpFlashTree->Fill();
  }
}

void opdet::OpFlashAnaAlg::FillOpHitTree(const std::vector<recob::OpHit>& hitVector)
{
  for (auto const& hit : hitVector) {

    fOpHitAnaStruct.HitPeakTime = hit.PeakTime();
    fOpHitAnaStruct.HitPeakTimeAbs = hit.PeakTimeAbs();
    fOpHitAnaStruct.HitWidth = hit.Width();
    fOpHitAnaStruct.HitArea = hit.Area();
    fOpHitAnaStruct.HitAmplitude = hit.Amplitude();
    fOpHitAnaStruct.HitFastToTotal = hit.FastToTotal();
    fOpHitAnaStruct.HitPE = hit.PE();
    fOpHitAnaStruct.HitFrame = hit.Frame();
    fOpHitAnaStruct.HitOpChannel = hit.OpChannel();

    fOpHitTree->Fill();
  }
}
