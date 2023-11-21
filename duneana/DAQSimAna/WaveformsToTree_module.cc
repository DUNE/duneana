////////////////////////////////////////////////////////////////////////
// Class:       WaveformsToTree
// Plugin Type: producer (art v2_10_03)
// File:        WaveformsToTree_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"

#include <memory>
#include <fstream>

#include "TTree.h"

class WaveformsToTree : public art::EDAnalyzer {
public:
    explicit WaveformsToTree(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    WaveformsToTree(WaveformsToTree const &) = delete;
    WaveformsToTree(WaveformsToTree &&) = delete;
    WaveformsToTree & operator = (WaveformsToTree const &) = delete;
    WaveformsToTree & operator = (WaveformsToTree &&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    void beginJob() override;

    void endJob() override { m_tree->Write(); }
private:
    // The module name of the raw digits we're reading in
    std::string m_inputTag;
    int m_maxChannels;
    TTree* m_tree;
    std::vector<std::vector<int> > m_waveforms;
    std::vector<int> m_chans;
};


WaveformsToTree::WaveformsToTree(fhicl::ParameterSet const & p)
    : EDAnalyzer(p),
      m_inputTag(p.get<std::string>("InputTag", "daq")),
      m_maxChannels(p.get<int>("MaxChannels"))
{
}

void WaveformsToTree::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  m_tree = tfs->make<TTree>("WaveformTree","title");
  m_tree->Branch<std::vector<std::vector<int> > >("waveforms", &m_waveforms);
  m_tree->Branch("chans", &m_chans);
}

void WaveformsToTree::analyze(art::Event const& e)
{
    m_waveforms.clear();
    m_chans.clear();

    auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(m_inputTag);
    auto& digits_in =*digits_handle;

    int nChan=0;
    art::ServiceHandle<geo::Geometry> geo;
    for(auto&& digit: digits_in){
        bool isCollection=geo->SignalType(digit.Channel())==geo::kCollection;
        if(!isCollection) continue;
        if(nChan++ > m_maxChannels) break;

        m_chans.push_back(digit.Channel());
        std::vector<int> waveform;
        for(auto const& adc: digit.ADCs()){
            waveform.push_back(adc);
        }
        m_waveforms.push_back(std::move(waveform));
    }
    m_tree->Fill();
}

DEFINE_ART_MODULE(WaveformsToTree)
