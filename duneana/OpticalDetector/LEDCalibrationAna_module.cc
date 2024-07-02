// -*- mode: c++; c-basic-offset: 2; -*-
// Ben Jones, MIT, 2013
//
//  This ana module extracts pedestals and gains from
//  LED calibration run data (incomplete)
//

// LArSoft includes
#include "larana/OpticalDetector/OpHitFinder/AlgoThreshold.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoEdges.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TF1.h"
#include "TH1.h"
#include <TTree.h>

// C++ Includes
#include "math.h"
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>

namespace opdet {

  class LEDCalibrationAna : public art::EDAnalyzer {
  public:
    // Standard constructor and destructor for an ART module.
    LEDCalibrationAna(const fhicl::ParameterSet&);

    void endJob();

    // The analyzer routine, called once per event.
    void analyze(const art::Event&);

  private:
    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpDet collection
    uint32_t fTriggerChannel;
    float fCoincThreshold;
    float fMaxTimeMean;
    float fMaxTimeThresh;
    float fAreaDivs;
    float fAreaMin;
    float fAreaMax;
    float fTriggerDelay;

    pmtana::PulseRecoManager fPulseRecoMgr;
    pmtana::AlgoThreshold fThreshAlg;
    pmtana::PedAlgoEdges fPedAlg;
    TTree* fPulseTree;
    TTree* fPulseTreeNonCoinc;
    Float_t fPeak;
    Float_t fPedRMS;
    Float_t fPedMean;
    Float_t fArea;
    Float_t fTBegin;
    Float_t fTEnd;
    Float_t fTMax;
    Float_t fOffset;
    Int_t fEventID;
    Int_t fRunID;
    Int_t fChannel;
    Int_t fShaper;

    uint32_t ShaperToChannel(uint32_t Shaper);

    bool fMakeNonCoincTree;

    std::map<uint32_t, std::vector<double>> fAreas;
  };

}

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  LEDCalibrationAna::LEDCalibrationAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset), fPulseRecoMgr(), fThreshAlg(), fPedAlg()
  {
    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");
    fTriggerChannel = pset.get<uint32_t>("TriggerChannel");
    fTriggerDelay = pset.get<uint32_t>("TriggerDelay");
    fCoincThreshold = pset.get<float>("CoincThreshold");
    fMaxTimeThresh = pset.get<float>("MaxTimeThresh");
    fMaxTimeMean = pset.get<float>("MaxTimeMean");
    fAreaMin = pset.get<float>("AreaMin");
    fAreaMax = pset.get<float>("AreaMax");
    fAreaDivs = pset.get<float>("AreaDivs");
    fMakeNonCoincTree = pset.get<bool>("MakeNonCoincTree");

    fPulseRecoMgr.AddRecoAlgo(&fThreshAlg);
    fPulseRecoMgr.SetDefaultPedAlgo(&fPedAlg);

    art::ServiceHandle<art::TFileService const> tfs;

    fPulseTree = tfs->make<TTree>("fPulseTree", "fPulseTree");
    fPulseTree->Branch("Area", &fArea, "Area/F");
    fPulseTree->Branch("Peak", &fPeak, "Peak/F");
    fPulseTree->Branch("TBegin", &fTBegin, "TBegin/F");
    fPulseTree->Branch("TEnd", &fTEnd, "TEnd/F");
    fPulseTree->Branch("TMax", &fTMax, "TMax/F");
    fPulseTree->Branch("Channel", &fChannel, "Channel/I");
    fPulseTree->Branch("Shaper", &fShaper, "Shaper/I");
    fPulseTree->Branch("PedMean", &fPedMean, "PedMean/F");
    fPulseTree->Branch("PedRMS", &fPedRMS, "PedRMS/F");
    fPulseTree->Branch("Offset", &fOffset, "Offset/F");
    fPulseTree->Branch("EventID", &fEventID, "EventID/I");
    fPulseTree->Branch("RunID", &fRunID, "RunID/I");

    fPulseTreeNonCoinc = tfs->make<TTree>("fPulseTreeNonCoinc", "fPulseTreeNonCoinc");
    fPulseTreeNonCoinc->Branch("Area", &fArea, "Area/F");
    fPulseTreeNonCoinc->Branch("Peak", &fPeak, "Peak/F");
    fPulseTreeNonCoinc->Branch("TBegin", &fTBegin, "TBegin/F");
    fPulseTreeNonCoinc->Branch("TEnd", &fTEnd, "TEnd/F");
    fPulseTreeNonCoinc->Branch("TMax", &fTMax, "TMax/F");
    fPulseTreeNonCoinc->Branch("Channel", &fChannel, "Channel/I");
    fPulseTreeNonCoinc->Branch("Shaper", &fShaper, "Shaper/I");
    fPulseTreeNonCoinc->Branch("PedMean", &fPedMean, "PedMean/F");
    fPulseTreeNonCoinc->Branch("PedRMS", &fPedRMS, "PedRMS/F");
    fPulseTreeNonCoinc->Branch("EventID", &fEventID, "EventID/I");
    fPulseTreeNonCoinc->Branch("RunID", &fRunID, "RunID/I");
  }

  //-----------------------------------------------------------------------
  void LEDCalibrationAna::endJob()
  {
    art::ServiceHandle<art::TFileService const> tfs;

    for (auto it = fAreas.begin(); it != fAreas.end(); ++it) {
      uint32_t Channel = it->first;

      std::stringstream histname;
      histname.flush();
      histname << "ch" << Channel << "area";

      TH1D* HistArea = tfs->make<TH1D>(
        histname.str().c_str(), histname.str().c_str(), fAreaDivs, fAreaMin, fAreaMax);

      for (size_t j = 0; j != it->second.size(); ++j) {
        HistArea->Fill(it->second.at(j));
      }

      std::stringstream fitname;
      fitname.flush();
      fitname << "ch" << Channel << "fit";

      double Max = HistArea->GetMaximum();
      double Mid = HistArea->GetBinContent(fAreaDivs / 2.);

      TF1* GausFit = new TF1(fitname.str().c_str(), "gaus(0)+gaus(3)+gaus(6)", fAreaMin, fAreaMax);

      GausFit->SetParameters(Mid,
                             (fAreaMin + fAreaMax) / 2.,
                             (fAreaMax - fAreaMin) / 2.,
                             Max,
                             0,
                             (fAreaMin + fAreaMax) / 8.,
                             Max / 5.,
                             0,
                             (fAreaMin + fAreaMax) / 4.);

      GausFit->SetParLimits(0, 0, 1.1 * Max);
      GausFit->SetParLimits(1, 0, fAreaMax);
      GausFit->SetParLimits(2, 0, fAreaMax);

      GausFit->SetParLimits(3, 0, 1.1 * Max);
      GausFit->FixParameter(4, 0);
      GausFit->SetParLimits(5, 0, (fAreaMin + fAreaMax) / 2.);

      GausFit->SetParLimits(6, 0, 1.1 * Max);
      GausFit->FixParameter(7, 0);
      GausFit->SetParLimits(8, 0, (fAreaMin + fAreaMax) / 2.);

      HistArea->Fit(GausFit);

      double Mean = GausFit->GetParameter(1);
      double Width = GausFit->GetParameter(2);

      double MeanErr = GausFit->GetParError(1);
      double WidthErr = GausFit->GetParError(2);

      double NPE = pow(Mean, 2) / pow(Width, 2);
      double SPEScale = Mean / NPE;

      double NPEError = NPE * pow(2. * (pow(MeanErr / Mean, 2) + pow(WidthErr / Width, 2)), 0.5);
      double SPEError = SPEScale * pow(2. * pow(WidthErr / Width, 2) + pow(MeanErr / Mean, 2), 0.5);

      std::cout << "Channel " << Channel << ":\tSPE Scale \t" << SPEScale << "\t +/- \t" << SPEError
                << ",\t NPE \t" << NPE << "\t +/- \t" << NPEError << std::endl;
    }
  }

  //-----------------------------------------------------------------------
  void LEDCalibrationAna::analyze(const art::Event& evt)
  {
    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    fRunID = evt.run();
    fEventID = evt.event();

    // Create a handle for our vector of pulses
    art::Handle<std::vector<raw::OpDetWaveform>> OpDetWaveformHandle;

    // Read in WaveformHandle
    evt.getByLabel(fInputModule, OpDetWaveformHandle);

    std::map<uint32_t, std::vector<int>> OrgOpDigitByChannel;

    for (size_t i = 0; i != OpDetWaveformHandle->size(); ++i) {
      OrgOpDigitByChannel[ShaperToChannel(OpDetWaveformHandle->at(i).ChannelNumber())].push_back(i);
    }

    std::vector<uint32_t> FrameNumbersForTrig;
    std::vector<uint32_t> TimeSlicesForTrig;

    for (size_t i = 0; i != OrgOpDigitByChannel[fTriggerChannel].size(); ++i) {
      double TimeStamp =
        OpDetWaveformHandle->at(OrgOpDigitByChannel[fTriggerChannel][i]).TimeStamp();
      uint32_t Frame = clock_data.OpticalClock().Frame(TimeStamp);
      uint32_t TimeSlice = clock_data.OpticalClock().Sample(TimeStamp);
      FrameNumbersForTrig.push_back(Frame);
      TimeSlicesForTrig.push_back(TimeSlice);
    }

    for (size_t i = 0; i != OpDetWaveformHandle->size(); ++i) {
      double TimeStamp = OpDetWaveformHandle->at(i).TimeStamp();
      uint32_t Frame = clock_data.OpticalClock().Frame(TimeStamp);
      uint32_t TimeSlice = clock_data.OpticalClock().Sample(TimeStamp);
      fShaper = OpDetWaveformHandle->at(i).ChannelNumber();
      fChannel = ShaperToChannel(fShaper);

      if (uint32_t(fChannel) != fTriggerChannel) {
        for (size_t j = 0; j != FrameNumbersForTrig.size(); ++j) {
          if ((Frame == FrameNumbersForTrig.at(j)) &&
              (fabs(TimeSlice - TimeSlicesForTrig.at(j) - fTriggerDelay) < fCoincThreshold)) {

            const raw::OpDetWaveform& wf = OpDetWaveformHandle->at(i);

            //fPulseRecoMgr.RecoPulse(wf);
            fPulseRecoMgr.Reconstruct(wf);

            size_t NPulses = fThreshAlg.GetNPulse();

            fOffset = TimeSlice - TimeSlicesForTrig.at(j);
            //fPedMean = fThreshAlg.PedMean();
            //fPedRMS  = fThreshAlg.PedRms();

            for (size_t k = 0; k != NPulses; ++k) {
              if (fabs(fMaxTimeMean - fThreshAlg.GetPulse(k).t_max) < fMaxTimeThresh) {

                fPeak = fThreshAlg.GetPulse(k).peak;
                fArea = fThreshAlg.GetPulse(k).area;
                fTBegin = fThreshAlg.GetPulse(k).t_start;
                fTEnd = fThreshAlg.GetPulse(k).t_end;
                fTMax = fThreshAlg.GetPulse(k).t_max;
                fPedMean = fThreshAlg.GetPulse(k).ped_mean;
                fPedRMS = fThreshAlg.GetPulse(k).ped_sigma;

                fPulseTree->Fill();

                fAreas[fChannel].push_back(fArea);
              }
              else if (fMakeNonCoincTree) {
                fPeak = fThreshAlg.GetPulse(k).peak;
                fArea = fThreshAlg.GetPulse(k).area;
                fTBegin = fThreshAlg.GetPulse(k).t_start;
                fTEnd = fThreshAlg.GetPulse(k).t_end;
                fTMax = fThreshAlg.GetPulse(k).t_max;

                fPulseTreeNonCoinc->Fill();
              }
            }
          }
        }
      }
    }
  }

  //---------------------------------
  uint32_t LEDCalibrationAna::ShaperToChannel(uint32_t Shaper)
  {
    static std::map<uint32_t, uint32_t> ShaperToChannelMap;
    if (ShaperToChannelMap.size() == 0) {

      // temporary
      for (size_t i = 0; i != 40; ++i) {
        ShaperToChannelMap[i] = i;
      }
    }

    return ShaperToChannelMap[Shaper];
  }

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(LEDCalibrationAna)
}
