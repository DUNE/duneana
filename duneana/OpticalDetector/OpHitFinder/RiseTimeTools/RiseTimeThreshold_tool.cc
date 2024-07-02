/**
 * \file RiseTimeThreshold_tool.cc
 *
 * \brief Rise time is defined as the time slot in which the pulse
 * goes above a certain fraction of the maximum ADC peak value
 * given by the "PeakRatio" fhicl parameter
 *
 * @author Fran Nicolas, June 2022
 */

#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/types/Atom.h"

#include "RiseTimeCalculatorBase.h"

namespace pmtana {

  class RiseTimeThreshold : RiseTimeCalculatorBase {

  public:
    //Configuration parameters
    struct Config {

      fhicl::Atom<double> PeakRatio{fhicl::Name("PeakRatio")};
    };

    // Default constructor
    explicit RiseTimeThreshold(art::ToolConfigTable<Config> const& config);

    // Method to calculate the OpFlash t0
    double RiseTime(const pmtana::Waveform_t& wf_pulse,
                    const pmtana::PedestalMean_t& ped_pulse,
                    bool _positive) const override;

  private:
    double fPeakRatio;
  };

  RiseTimeThreshold::RiseTimeThreshold(art::ToolConfigTable<Config> const& config)
    : fPeakRatio{config().PeakRatio()}
  {}

  double RiseTimeThreshold::RiseTime(const pmtana::Waveform_t& wf_pulse,
                                     const pmtana::PedestalMean_t& ped_pulse,
                                     bool _positive) const
  {

    // Pedestal-subtracted pulse
    std::vector<double> wf_aux(ped_pulse);
    if (_positive) {
      for (size_t ix = 0; ix < wf_aux.size(); ix++) {
        wf_aux[ix] = ((double)wf_pulse[ix]) - wf_aux[ix];
      }
    }
    else {
      for (size_t ix = 0; ix < wf_aux.size(); ix++) {
        wf_aux[ix] = wf_aux[ix] - ((double)wf_pulse[ix]);
      }
    }

    auto it_max = max_element(wf_aux.begin(), wf_aux.end());
    size_t rise = std::lower_bound(wf_aux.begin(), it_max, fPeakRatio * (*it_max)) - wf_aux.begin();

    return rise;
  }

}

DEFINE_ART_CLASS_TOOL(pmtana::RiseTimeThreshold)
