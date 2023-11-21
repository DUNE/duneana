////////////////////////////////////////////////////////////////////////
// Class:       WaveformThreshDump
// Plugin Type: producer (art v2_10_03)
// File:        WaveformThreshDump_module.cc
// Author:      Klaudia Wawrowska
// Date:        October 2023
//
// Module to save waveforms which contained ADC signal
// above a user-specified threshold to a file.   
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

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"

#include <memory>
#include <fstream>
#include <numeric> 

class WaveformThreshDump;


class WaveformThreshDump : public art::EDAnalyzer {
public:
  explicit WaveformThreshDump(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  
  // Plugins should not be copied or assigned.
  WaveformThreshDump(WaveformThreshDump const &) = delete;
  WaveformThreshDump(WaveformThreshDump &&) = delete;
  WaveformThreshDump & operator = (WaveformThreshDump const &) = delete;
  WaveformThreshDump & operator = (WaveformThreshDump &&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  
private:
  //Raw digits label
  std::string m_inputTag;

  //File where the waveforms will be saved 
  std::string m_outputFilename;
  std::ofstream m_outputFile;

  int m_thresh; // static threshold [ADC] 

  std::vector<int> GetProcessedWaveform( std::vector<short> waveform_short);
  bool aboveThresh( const std::vector<short>& vec, int threshold);
};


WaveformThreshDump::WaveformThreshDump(fhicl::ParameterSet const & p)
    : EDAnalyzer(p),
      m_inputTag(p.get<std::string>("InputTag", "daq")), 
      m_outputFilename(p.get<std::string>("OutputFile")),
      m_outputFile(m_outputFilename),
      m_thresh(p.get<int>("Threshold", 30))
{
}

//Only output waveforms with potential signal to get rid of useless, empty data
bool WaveformThreshDump::aboveThresh( const std::vector<short>& vec, int threshold){
  double mean = std::accumulate(vec.begin(), vec.end(), 0) / vec.size(); 
  for (const short& value : vec) {
    if (value > threshold + mean) { return true; } // Potential signal found
  } 
  return false; //No potential signal found.
}
 

void WaveformThreshDump::analyze(art::Event const& e)
{     
    auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(m_inputTag);
    auto& digits_in =*digits_handle;

    art::ServiceHandle<geo::Geometry> geo;

    //plane-tagged ADC data 
    std::vector<std::vector<short>> adc_samples;
    std::vector<unsigned int> channel_numbers;
    std::vector<unsigned int> plane;
    
      
    //Select only channels from first TPC and sort into respective fields
    for ( auto&& digit: digits_in){
      
      unsigned int planeID = 999; 

      //determine channel view 
      if (geo->SignalType(digit.Channel())==geo::kCollection){ planeID = 2; }
      if (geo->View(digit.Channel())==geo::kU){ planeID =0; }
      if (geo->View(digit.Channel())==geo::kV){ planeID =1; }
      
      //save sorted ADC info
      channel_numbers.push_back(digit.Channel());
      adc_samples.push_back(digit.ADCs());
      plane.push_back(planeID);

    } // Run over ADC data 

    for (size_t ich = 0; ich < adc_samples.size(); ++ich){

      std::vector<short> waveform = adc_samples[ich]; 

      // Save waveform data if it potentially contains signal 
      if ( aboveThresh(waveform, m_thresh) ) {
	m_outputFile << e.event() << " "
		     << channel_numbers[ich] << " "
		     << plane[ich] << " ";
	for (size_t j = 0; j < waveform.size(); ++j){ 
	  m_outputFile << waveform[j] << " "; 
	}
	m_outputFile << std::endl;
      }
    }
}


DEFINE_ART_MODULE(WaveformThreshDump)
