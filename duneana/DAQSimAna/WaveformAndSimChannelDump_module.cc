////////////////////////////////////////////////////////////////////////
// Class:       WaveformAndSimChannelDump
// Plugin Type: producer (art v2_10_03)
// File:        WaveformAndSimChannelDump_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include <memory>
#include <fstream>

class WaveformAndSimChannelDump;


class WaveformAndSimChannelDump : public art::EDAnalyzer {
public:
  explicit WaveformAndSimChannelDump(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WaveformAndSimChannelDump(WaveformAndSimChannelDump const &) = delete;
  WaveformAndSimChannelDump(WaveformAndSimChannelDump &&) = delete;
  WaveformAndSimChannelDump & operator = (WaveformAndSimChannelDump const &) = delete;
  WaveformAndSimChannelDump & operator = (WaveformAndSimChannelDump &&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  // The module name of the raw digits we're reading in
  std::string m_inputTagGEANT;
  std::string m_inputTagTPC;
  std::string m_outputFilename_tpc     ;
  std::string m_outputFilename_true_tpc;
  size_t m_max_channel;
  std::ofstream m_outputFile_tpc     ;
  std::ofstream m_outputFile_true_tpc;
};


WaveformAndSimChannelDump::WaveformAndSimChannelDump(fhicl::ParameterSet const & p)
  : EDAnalyzer(p),
    m_inputTagGEANT (p.get<std::string>("InputTagGEANT" , "largeant")),
    m_inputTagTPC   (p.get<std::string>("InputTagTPC"   , "daq"     )),
    m_outputFilename_tpc     (p.get<std::string>("OutputFileTPC"    ,"OutputFileTPC.txt"    )),
    m_outputFilename_true_tpc(p.get<std::string>("OutputFileTrueTPC","OutputFileTrueTPC.txt")),
    m_max_channel(p.get<size_t>("MaxChannels", 2560)),
    m_outputFile_tpc     (m_outputFilename_tpc     ),
    m_outputFile_true_tpc(m_outputFilename_true_tpc)
{
}

void WaveformAndSimChannelDump::analyze(art::Event const& e)
{
  art::ServiceHandle<geo::Geometry> geo;

  // TPC Waveforms
  size_t n_ticks_tpc = 0;
  auto const& digits_handle_tpc=e.getValidHandle<std::vector<raw::RawDigit>>(m_inputTagTPC);
  auto& digits_tpc_in =*digits_handle_tpc;
  for (auto&& digit: digits_tpc_in) {
    bool isCollection=geo->SignalType(digit.Channel())==geo::kCollection;
    if (digit.Channel() >= m_max_channel) continue;
    m_outputFile_tpc << e.event() << " "
                     << digit.Channel() << " "
                     << isCollection << " ";
    for(auto const& adc: digit.ADCs()){
      m_outputFile_tpc << adc << " ";
    }
    n_ticks_tpc=digits_tpc_in.at(0).ADCs().size();
    m_outputFile_tpc << "\n";
  }

  // TPC truth
  auto const& truth_handle_tpc=e.getValidHandle<std::vector<sim::SimChannel>>(m_inputTagGEANT);
  auto& truth_tpc_in =*truth_handle_tpc;
  for (auto&& truth: truth_tpc_in) {
    bool isCollection=geo->SignalType(truth.Channel())==geo::kCollection;
    if (truth.Channel() >= m_max_channel) continue;
    m_outputFile_true_tpc << e.event() << " "
                           << truth.Channel() <<" "
                          << isCollection << " ";
    for(size_t ichge=0; ichge<n_ticks_tpc; ++ichge){
      m_outputFile_true_tpc << truth.Charge(ichge) << " ";
    }
    m_outputFile_true_tpc << "\n";
  }
}

DEFINE_ART_MODULE(WaveformAndSimChannelDump)
