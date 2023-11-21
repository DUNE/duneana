////////////////////////////////////////////////////////////////////////
// Class:       AbsRunningSumTPFinderPass1
// File:        AbsRunningSumTPFinderPass1_service.cc
////////////////////////////////////////////////////////////////////////

#include "duneana/DAQSimAna/AbsRunningSumHitFinder/AbsRunningSumTPFinderPass1.h"
#include "duneana/DAQSimAna/AlgParts.h"

#include <algorithm> // for std::transform
#include <numeric> // for std::accumulate


AbsRunningSumTPFinderPass1::AbsRunningSumTPFinderPass1(fhicl::ParameterSet const & p)
  : m_threshold          (p.get<unsigned int>      ("Threshold"            )),
    m_useSignalKill      (p.get<bool>              ("UseSignalKill"        ,                        true)),
    m_signalKillLookahead(p.get<short>             ("SignalKillLookahead"  ,                           5)),
    m_signalKillThreshold(p.get<short>             ("SignalKillThreshold"  ,                          15)),
    m_signalKillNContig  (p.get<short>             ("SignalKillNContig"    ,                           1)),
    m_frugalNContig      (p.get<short>             ("FrugalPedestalNContig",                          10)),
    m_R                  (p.get<float>             ("R",                                              0.9))
   
{

}


std::vector<short> AbsRunningSumTPFinderPass1::findPedestal(const std::vector<short>& waveform)
{
  //---------------------------------------------
  // Pedestal subtraction
  //---------------------------------------------
  const std::vector<short>& pedestal=m_useSignalKill ?
    frugal_pedestal_sigkill(waveform,
			    m_signalKillLookahead,
			    m_signalKillThreshold,
			    m_signalKillNContig) :
    frugal_pedestal(waveform, m_frugalNContig);
  return pedestal;
}

std::vector<short> AbsRunningSumTPFinderPass1::AbsRunningSum(const std::vector<short>& pedsub_waveform, float R) {

  //---------------------------------------------
  // Absolute Running Sum Algorithm
  //---------------------------------------------

  
  //scaling factor to keep average noise RMS O(10) ADC, just like on raw waveforms
  short s = 2;
  int adcMax = SHRT_MAX; 

  std::vector<short> absRS(pedsub_waveform.size(), 0); absRS[0] = pedsub_waveform[0]/s;
  for (size_t i=0; i<pedsub_waveform.size(); ++i) { 

    //guard the signal from overflowing
    absRS[i] = std::min( int(R*absRS[i-1] + std::abs(pedsub_waveform[i])/s), adcMax);

  }
  return absRS;
}


void
AbsRunningSumTPFinderPass1::hitFinding(const std::vector<short>& waveform,
				    std::vector<AbsRunningSumTPFinderTool::Hit>& hits,
				    int channel) {

  //---------------------------------------------
  // Hit finding
  //---------------------------------------------
  bool is_hit  = false;
  bool was_hit = false;
  std::vector<int> hit_charge; 
  //initialize the hit 
  AbsRunningSumTPFinderTool::Hit hit(channel, 0, 0, 0, 0, 0);
  for(size_t isample=0; isample<waveform.size()-1; ++isample){

    //Ignore the first ~100 ticks to let the pedestal stabilise. 
    if (isample > 100) { 
      short adc         = waveform[isample];
      
      is_hit = adc >  (short)m_threshold;
      if(is_hit && !was_hit) {
	hit_charge.push_back(adc); 
	hit.startTime         = isample;
	hit.SADC              = adc;
	hit.timeOverThreshold = 1;
      }
      if(is_hit && was_hit) {
	hit.SADC              += adc;
	hit.timeOverThreshold += 1;
	hit_charge.push_back(adc);
      }
      if(!is_hit && was_hit) {
	hit.peakCharge = *std::max_element(hit_charge.begin(), hit_charge.end()); 
	hit.peakTime = std::distance(hit_charge.begin(), std::max_element(hit_charge.begin(), hit_charge.end())) + hit.startTime; 
	hits.push_back(hit);
	hit_charge.clear();
      }    
      was_hit = is_hit; 
    }
  }
} 



std::vector<AbsRunningSumTPFinderTool::Hit>
AbsRunningSumTPFinderPass1::findHits(const std::vector<unsigned int>& channel_numbers, 
				  const std::vector<std::vector<short>>& adc_samples) {

  auto hits = std::vector<AbsRunningSumTPFinderTool::Hit>();

  for(size_t ich=0; ich<adc_samples.size(); ++ich){
    const std::vector<short>& waveform = adc_samples[ich];

    //First pedestal subtraction 
    std::vector<short> pedestal  = findPedestal(waveform);
    std::vector<short> pedsub(waveform.size(), 0);
    for(size_t i=0; i<pedsub.size(); ++i)
      pedsub[i]=waveform[i]-pedestal[i];

    //AbsRS filtering and second pedsub 
    std::vector<short> AbsRS = AbsRunningSum(pedsub, m_R);
    std::vector<short> AbsRS_ped = findPedestal(AbsRS);
    std::vector<short> AbsRS_pedsub(AbsRS.size(), 0); 
    for (size_t j = 0; j < AbsRS_pedsub.size(); ++j){
      AbsRS_pedsub[j] = AbsRS[j] - AbsRS_ped[j];
    }

    //hit finding 
    hitFinding(AbsRS_pedsub, hits, channel_numbers[ich]);
  }
  std::cout<< " ---- Using AbsRS hit finder --- " << std::endl;
  std::cout << "Returning " << hits.size() << " hits for a threshold of " << m_threshold << std::endl;
  return hits;
}

DEFINE_ART_CLASS_TOOL(AbsRunningSumTPFinderPass1)
