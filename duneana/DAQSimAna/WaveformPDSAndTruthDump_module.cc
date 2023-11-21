////////////////////////////////////////////////////////////////////////
// Class:       WaveformPDSAndTruthDump
// Plugin Type: producer (art v2_10_03)
// File:        WaveformPDSAndTruthDump_module.cc
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

class WaveformPDSAndTruthDump;


class WaveformPDSAndTruthDump : public art::EDAnalyzer {
public:
  explicit WaveformPDSAndTruthDump(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WaveformPDSAndTruthDump(WaveformPDSAndTruthDump const &) = delete;
  WaveformPDSAndTruthDump(WaveformPDSAndTruthDump &&) = delete;
  WaveformPDSAndTruthDump & operator = (WaveformPDSAndTruthDump const &) = delete;
  WaveformPDSAndTruthDump & operator = (WaveformPDSAndTruthDump &&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  // The module name of the raw digits we're reading in
  std::string m_inputTagGEANT;
  std::string m_inputTagPDS;
  std::string m_outputFilename_pds     ;
  std::string m_outputFilename_true_pds;
  std::ofstream m_outputFile_pds     ;
  std::ofstream m_outputFile_true_pds;
};


WaveformPDSAndTruthDump::WaveformPDSAndTruthDump(fhicl::ParameterSet const & p)
  : EDAnalyzer(p),
    m_inputTagGEANT (p.get<std::string>("InputTagGEANT" , "largeant")),
    m_inputTagPDS   (p.get<std::string>("InputTagPDS"   , "opdigi"  )), 
    m_outputFilename_pds     (p.get<std::string>("OutputFilePDS"    ,"OutputFilePDS.txt"    )),
    m_outputFilename_true_pds(p.get<std::string>("OutputFileTruePDS","OutputFileTruePDS.txt")),
    m_outputFile_pds     (m_outputFilename_pds     ),
    m_outputFile_true_pds(m_outputFilename_true_pds)
{
}

void WaveformPDSAndTruthDump::analyze(art::Event const& e)
{
  art::ServiceHandle<geo::Geometry> geo;
  auto const& digits_handle_pds=e.getValidHandle<std::vector<raw::OpDetWaveform>>(m_inputTagPDS);
  std::map<int, std::pair<int, int>> channel_timestamp;
  
  auto& digits_pds_in = *digits_handle_pds;
  for (auto&& digit: digits_pds_in) {
    m_outputFile_pds << e.event() << " "
                     << digit.ChannelNumber() << " ";
                     
    auto f = channel_timestamp.find(digit.ChannelNumber());
    if (f == channel_timestamp.end()) {
      channel_timestamp[digit.ChannelNumber()] = std::make_pair(digit.TimeStamp(), digit.size());
    } else {
      int begin_old = f->second.first;
      int begin_new = digit.TimeStamp();
      begin_new = (begin_old > begin_new) ? begin_new : begin_old;
      int end_old = f->second.first+f->second.second;
      int end_new = digit.TimeStamp()+digit.size();
      end_new = (end_old > end_new) ? end_old : end_new;
      f->second.first = begin_new;
      f->second.second = end_new - begin_new;
    }
    for(auto const& adc: digit) {
      m_outputFile_pds << adc << " ";
    }
    
    m_outputFile_pds << "\n";
  }

  // PDS truth
  auto const& truth_handle_pds=e.getValidHandle<std::vector<sim::OpDetBacktrackerRecord>>(m_inputTagGEANT);
  auto& truth_pds_in =*truth_handle_pds;

  for (auto const& interesting: channel_timestamp) {
    int channel = interesting.first;
    int timestamp = interesting.second.first;
    int nticks = interesting.second.second;
    bool found = false;
    for (auto&& truth: truth_pds_in) {
      if (truth.OpDetNum() != channel) continue;
      m_outputFile_true_pds << e.event() << " "
                            << truth.OpDetNum() << " ";
      for(int ichge=timestamp; ichge<nticks+timestamp; ++ichge){
        m_outputFile_true_pds << truth.Photons(ichge) << " ";
      }
      m_outputFile_true_pds << "\n";
      found=true;
      break;
    }
    if (!found) {
      std::cout << "didn't find this channel " << channel
                << " timestamp " << timestamp
                << " ticks " << nticks <<" in truth\n";
    }
  }
}


DEFINE_ART_MODULE(WaveformPDSAndTruthDump)
