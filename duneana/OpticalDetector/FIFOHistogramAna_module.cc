// Ben Jones, MIT, 2013
//
//  This simple ana module writes a root file of pulse histograms,
//  organized by event and channel number
//

// LArSoft includes
#include "lardataobj/OpticalDetectorData/FIFOChannel.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TH1.h"

// C++ Includes
#include <cstring>
#include <map>
#include <sstream>

namespace opdet {

  class FIFOHistogramAna : public art::EDAnalyzer {
  public:
    // Standard constructor and destructor for an ART module.
    FIFOHistogramAna(const fhicl::ParameterSet&);

    // The analyzer routine, called once per event.
    void analyze(const art::Event&);

  private:
    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpDet collection
  };

}

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  FIFOHistogramAna::FIFOHistogramAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {
    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");
  }

  //-----------------------------------------------------------------------
  void FIFOHistogramAna::analyze(const art::Event& evt)
  {

    art::ServiceHandle<art::TFileService> tfs;

    // Create a handle for our vector of pulses
    art::Handle<std::vector<optdata::FIFOChannel>> FIFOChannelHandle;

    // Read in WaveformHandle
    evt.getByLabel(fInputModule, FIFOChannelHandle);

    int Run = evt.run();
    int EID = evt.event();

    std::stringstream FolderName;
    FolderName.flush();
    FolderName << "run" << Run << "_evt" << EID;

    art::TFileDirectory evtfolder = tfs->mkdir(FolderName.str().c_str());

    std::map<int, bool> ChanFolderMade;
    std::map<uint32_t, int> ChanFolderIndex;
    std::vector<art::TFileDirectory> ChanFolders;

    for (size_t i = 0; i != FIFOChannelHandle->size(); ++i) {
      uint32_t Frame = FIFOChannelHandle->at(i).Frame();
      uint32_t TimeSlice = FIFOChannelHandle->at(i).TimeSlice();
      uint32_t Channel = FIFOChannelHandle->at(i).ChannelNumber();

      if (!ChanFolderMade[Channel]) {
        std::stringstream ChannelLabel;
        ChannelLabel.flush();
        ChannelLabel << "chan" << Channel;
        ChanFolderIndex[Channel] = ChanFolders.size();
        ChanFolders.push_back(evtfolder.mkdir(ChannelLabel.str().c_str()));
        ChanFolderMade[Channel] = true;
      }

      std::stringstream HistName;
      HistName.flush();
      HistName << "frm" << Frame << "_"
               << "tsl" << TimeSlice;

      TH1D* ThisHist = ChanFolders[ChanFolderIndex[Channel]].make<TH1D>(
        HistName.str().c_str(),
        HistName.str().c_str(),
        FIFOChannelHandle->at(i).size(),
        float(TimeSlice) - 0.0001,
        float(FIFOChannelHandle->at(i).size()) - 0.0001 + TimeSlice);

      for (size_t j = 0; j != FIFOChannelHandle->at(i).size(); ++j) {
        ThisHist->Fill(TimeSlice + j, FIFOChannelHandle->at(i).at(j));
      }
    }
  }

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(FIFOHistogramAna)
}
