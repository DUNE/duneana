#ifndef ABSRUNNINGSUMTPFINDERPASS1_H
#define ABSRUNNINGSUMTPFINDERPASS1_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "duneana/DAQSimAna/AbsRunningSumHitFinder/AbsRunningSumTPFinderTool.h"

class AbsRunningSumTPFinderPass1 : public AbsRunningSumTPFinderTool {
 public:
  explicit AbsRunningSumTPFinderPass1(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  virtual std::vector<AbsRunningSumTPFinderTool::Hit> findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& adc_samples);
    
  protected:
  std::vector<short> findPedestal      (const std::vector<short>& orig);
  std::vector<short> AbsRunningSum     (const std::vector<short>& orig, float R);

  void hitFinding(const std::vector<short>& waveform, std::vector<AbsRunningSumTPFinderTool::Hit>& hits, int channel);

  unsigned int       m_threshold;
  bool               m_useSignalKill;
  short              m_signalKillLookahead;
  short              m_signalKillThreshold;
  short              m_signalKillNContig;
  short              m_frugalNContig;
  float              m_R;
};

#endif
