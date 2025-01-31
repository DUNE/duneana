////////////////////////////////////////////////////////////////////////
// Class:       WaveformDump
// Plugin Type: producer (art v2_10_03)
// File:        WaveformDump_module.cc
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

#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/RawData/RawDigit.h"

#include <memory>
#include <fstream>

class WaveformDump;


class WaveformDump : public art::EDAnalyzer {
public:
  explicit WaveformDump(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WaveformDump(WaveformDump const &) = delete;
  WaveformDump(WaveformDump &&) = delete;
  WaveformDump & operator = (WaveformDump const &) = delete;
  WaveformDump & operator = (WaveformDump &&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
    // The module name of the raw digits we're reading in
    std::string m_inputTag;
    std::string m_outputFilename;
    std::ofstream m_outputFile;
};


WaveformDump::WaveformDump(fhicl::ParameterSet const & p)
    : EDAnalyzer(p),
      m_inputTag(p.get<std::string>("InputTag", "daq")), 
      m_outputFilename(p.get<std::string>("OutputFile")),
      m_outputFile(m_outputFilename)
{
}

void WaveformDump::analyze(art::Event const& e)
{
    auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(m_inputTag);
    auto& digits_in =*digits_handle;

    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    for(auto&& digit: digits_in){
        bool isCollection=wireReadout.SignalType(digit.Channel())==geo::kCollection;
        m_outputFile << e.event() << " "
                     << digit.Channel() << " "
                     << isCollection << " ";
        for(auto const& adc: digit.ADCs()){
            m_outputFile << adc << " ";
        }
        m_outputFile << std::endl;
    }
}

DEFINE_ART_MODULE(WaveformDump)
