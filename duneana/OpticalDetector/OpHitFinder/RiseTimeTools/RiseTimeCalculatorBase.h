/**
 * \file RiseTimeCalculatorBase.h
 *
 *
 * \brief Interfacce class for a tool to calculate the pulse rise time
 *
 * @author Fran Nicolas, June 2022
 */

#ifndef RISETIMECALCULATORBASE_H
#define RISETIMECALCULATORBASE_H

#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"

namespace pmtana {
  class RiseTimeCalculatorBase {

  public:
    // Default destructor
    virtual ~RiseTimeCalculatorBase() noexcept = default;

    // Method to calculate the OpFlash t0
    virtual double RiseTime(const pmtana::Waveform_t& wf_pulse,
                            const pmtana::PedestalMean_t& ped_pulse,
                            bool _positive) const = 0;

  private:
  };
}

#endif
