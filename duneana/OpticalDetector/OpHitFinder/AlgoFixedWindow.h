/**
 * \file AlgoFixedWindow.h
 *
 * \ingroup PulseReco
 *
 * \brief Class definition file of AlgoFixedWindow
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup PulseReco

@{*/

#ifndef ALGOFIXEDWINDOW_H
#define ALGOFIXEDWINDOW_H

#include "PMTPulseRecoBase.h"

namespace fhicl {
  class ParameterSet;
}

#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"

#include <string>

namespace pmtana {

  /**
   \class AlgoFixedWindow
   This class implements fixed_window algorithm to AlgoFixedWindow class.
   The algorithm defines a pulse in user-specified time window.
   A typical usage is to set the beginning of the window to be 0 (= start of the waveform)
   and integrate over the time of interest. By default, the ending is set to index=0, in
   which case it uses the ending index of the input waveform (i.e. full integration).
  */
  class AlgoFixedWindow : public PMTPulseRecoBase {

  public:
    /// Default ctor
    AlgoFixedWindow(const std::string name = "FixedWindow");

    /// Alternative ctor
    AlgoFixedWindow(const fhicl::ParameterSet& pset,
                    std::unique_ptr<pmtana::RiseTimeCalculatorBase> risetimecalculator = nullptr,
                    const std::string name = "FixedWindow");
    //AlgoFixedWindow(const ::fcllite::PSet &pset,const std::string name="FixedWindow");

    /// Implementation of AlgoFixedWindow::reset() method
    void Reset();

  protected:
    /// Implementation of AlgoFixedWindow::reco() method
    bool RecoPulse(const pmtana::Waveform_t&,
                   const pmtana::PedestalMean_t&,
                   const pmtana::PedestalSigma_t&);

    size_t _index_start; ///< index marker for the beginning of the pulse time window
    size_t _index_end;   ///< index marker for the end of pulse time window
  };

}
#endif

/** @} */ // end of doxygen group
