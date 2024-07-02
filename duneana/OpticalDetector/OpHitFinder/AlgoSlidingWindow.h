/**
 * \file AlgoSlidingWindow.h
 *
 * \ingroup PulseReco
 *
 * \brief Class definition file of AlgoSlidingWindow
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup PulseReco

@{*/

#ifndef ALGOSLIDINGWINDOW_H
#define ALGOSLIDINGWINDOW_H

#include "PMTPulseRecoBase.h"
namespace fhicl {
  class ParameterSet;
}

#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"

#include <string>

namespace pmtana {

  /**
   \class AlgoSlidingWindow
   This class implements threshold algorithm to AlgoSlidingWindow class.
  */
  class AlgoSlidingWindow : public PMTPulseRecoBase {

  public:
    /// Default constructor
    AlgoSlidingWindow(const std::string name = "SlidingWindow");

    /// Alternative ctor
    AlgoSlidingWindow(const fhicl::ParameterSet& pset,
                      std::unique_ptr<pmtana::RiseTimeCalculatorBase> risetimecalculator = nullptr,
                      const std::string name = "SlidingWindow");
    //AlgoSlidingWindow(const ::fcllite::PSet &pset,const std::string name="SlidingWindow");

    /// Implementation of AlgoSlidingWindow::reset() method
    void Reset();

  protected:
    /// Implementation of AlgoSlidingWindow::reco() method
    bool RecoPulse(const pmtana::Waveform_t&,
                   const pmtana::PedestalMean_t&,
                   const pmtana::PedestalSigma_t&);

    /// A boolean to set waveform positive/negative polarity
    bool _positive;

    /// A variable holder for a user-defined absolute ADC threshold value
    float _adc_thres, _tail_adc_thres, _end_adc_thres;

    /// A variable holder to ensure the minimum pulse width
    size_t _min_width;

    /// A variable holder for a multiplicative factor for the pedestal standard deviation to define the threshold.
    float _nsigma, _tail_nsigma, _end_nsigma;
    bool _verbose;
    size_t _num_presample, _num_postsample;
  };

}
#endif

/** @} */ // end of doxygen group
