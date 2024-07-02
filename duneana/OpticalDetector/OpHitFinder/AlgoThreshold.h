/**
 * \file AlgoThreshold.h
 *
 * \ingroup PulseReco
 *
 * \brief Class definition file of AlgoThreshold
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup PulseReco

@{*/

#ifndef ALGOTHRESHOLD_H
#define ALGOTHRESHOLD_H

#include "PMTPulseRecoBase.h"
namespace fhicl {
  class ParameterSet;
}
#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"

#include <string>

namespace pmtana {

  /**
   \class AlgoThreshold
   This class implements threshold algorithm to AlgoThreshold class.
   The algorithm defines a pulse in user-specified time window.
   A typical usage is to set the beginning of the window to be 0 (= start of the waveform)
   and integrate over the time of interest. By default, the ending is set to index=0, in
   which case it uses the ending index of the input waveform (i.e. full integration).
  */
  class AlgoThreshold : public PMTPulseRecoBase {

  public:
    /// Default constructor
    AlgoThreshold(const std::string name = "AlgoThreshold");

    /// Alternative constructor
    AlgoThreshold(const fhicl::ParameterSet& pset,
                  std::unique_ptr<pmtana::RiseTimeCalculatorBase> risetimecalculator = nullptr,
                  const std::string name = "AlgoThreshold");
    //AlgoThreshold(const ::fcllite::PSet &pset,const std::string name="AlgoThreshold");

    /// Implementation of AlgoThreshold::reset() method
    void Reset();

  protected:
    /// Implementation of AlgoThreshold::reco() method
    bool RecoPulse(const pmtana::Waveform_t& wf,
                   const pmtana::PedestalMean_t& mean_v,
                   const pmtana::PedestalSigma_t& sigma_v);

    /// A variable holder for a user-defined absolute ADC threshold value
    //double _adc_thres;
    double _start_adc_thres;
    double _end_adc_thres;

    /// A variable holder for a multiplicative factor for the pedestal standard deviation to define the threshold.
    //double _nsigma;
    double _nsigma_start;
    double _nsigma_end;
  };

}
#endif

/** @} */ // end of doxygen group
