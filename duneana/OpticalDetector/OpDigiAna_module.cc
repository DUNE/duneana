// Christie Chiu and Ben Jones, MIT, 2012
//
// This is an analyzer module which writes the raw optical
// detector pulses for each PMT to an output file
//

// LArSoft includes
#include "larana/OpticalDetector/OpDigiProperties.h"
#include "lardataobj/RawData/OpDetPulse.h"

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
#include "TH1.h"

// C++ Includes
#include <algorithm>

namespace opdet {

  class OpDigiAna : public art::EDAnalyzer {
  public:
    // Standard constructor and destructor for an ART module.
    OpDigiAna(const fhicl::ParameterSet&);

    // The analyzer routine, called once per event.
    void analyze(const art::Event&);

  private:
    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fInputModule;  // Input tag for OpDet collection
    std::string fInstanceName; // Input tag for OpDet collection
    float fSampleFreq;         // in MHz
    float fTimeBegin;          // in us
    float fTimeEnd;            // in us
                               //  short fPEheight;                       // in ADC counts
    float fZeroSupThresh;

    double fSPEAmp;

    bool fMakeBipolarHist;
    bool fMakeUnipolarHist;
  };

}

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpDigiAna::OpDigiAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {

    art::ServiceHandle<OpDigiProperties> odp;
    fSPEAmp = odp->GetSPECumulativeAmplitude();
    fTimeBegin = odp->TimeBegin();
    fTimeEnd = odp->TimeEnd();
    fSampleFreq = odp->SampleFreq();

    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");
    fInstanceName = pset.get<std::string>("InstanceName");
    fZeroSupThresh = pset.get<double>("ZeroSupThresh") * fSPEAmp;

    fMakeBipolarHist = pset.get<bool>("MakeBipolarHist");
    fMakeUnipolarHist = pset.get<bool>("MakeUnipolarHist");
  }

  //-----------------------------------------------------------------------
  void OpDigiAna::analyze(const art::Event& evt)
  {

    // Create a handle for our vector of pulses
    art::Handle<std::vector<raw::OpDetPulse>> WaveformHandle;

    // Create string for histogram name
    char HistName[50];

    // Read in WaveformHandle
    evt.getByLabel(fInputModule, fInstanceName, WaveformHandle);

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us.
    art::ServiceHandle<art::TFileService const> tfs;

    std::vector<std::string> histnames;

    // For every OpDet waveform in the vector given by WaveformHandle
    for (unsigned int i = 0; i < WaveformHandle->size(); ++i) {
      // Get OpDetPulse
      art::Ptr<raw::OpDetPulse> ThePulsePtr(WaveformHandle, i);
      raw::OpDetPulse ThePulse = *ThePulsePtr;

      // Make an instance of histogram and its pointer, changing all units to ns
      // Notice that histogram axis is in ns but binned by 1/64MHz;
      //   appropriate conversions are made from beginning and end time
      //   in us, and frequency in MHz.

      // Make sure the histogram name is unique since there can be multiple pulses
      // per event and photo detector
      unsigned int pulsenum = 0;
      while (true) {
        sprintf(
          HistName, "Event_%d_OpDet_%i_Pulse_%i", evt.id().event(), ThePulse.OpChannel(), pulsenum);
        auto p = std::find(histnames.begin(), histnames.end(), HistName);
        if (p != histnames.end()) {
          // Found a duplicate
          pulsenum++;
        }
        else {
          // Found a unique name
          histnames.push_back(HistName);
          break;
        }
      }

      TH1D* PulseHist = nullptr;
      if (fMakeBipolarHist) {
        PulseHist = tfs->make<TH1D>(HistName,
                                    ";t (ns);",
                                    int((fTimeEnd - fTimeBegin) * fSampleFreq),
                                    fTimeBegin * 1000.,
                                    fTimeEnd * 1000.);
      }

      sprintf(HistName,
              "Event_%d_uni_OpDet_%i_Pulse_%i",
              evt.id().event(),
              ThePulse.OpChannel(),
              pulsenum);
      TH1D* UnipolarHist = nullptr;
      if (fMakeUnipolarHist) {
        UnipolarHist = tfs->make<TH1D>(HistName,
                                       ";t (ns);",
                                       int((fTimeEnd - fTimeBegin) * fSampleFreq),
                                       fTimeBegin * 1000.,
                                       fTimeEnd * 1000.);
      }

      for (unsigned int binNum = 0; binNum < ThePulse.Waveform().size(); ++binNum) {
        // Fill histogram with waveform

        if (fMakeBipolarHist) PulseHist->SetBinContent(binNum, (double)ThePulse.Waveform()[binNum]);
        if ((binNum > 0) && (fMakeUnipolarHist)) {
          double BinContent =
            (double)ThePulse.Waveform()[binNum] + (double)UnipolarHist->GetBinContent(binNum - 1);
          if (BinContent > fZeroSupThresh)
            UnipolarHist->SetBinContent(binNum, BinContent);
          else
            UnipolarHist->SetBinContent(binNum, 0);
        }
      }
    }
  }

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(OpDigiAna)
}
