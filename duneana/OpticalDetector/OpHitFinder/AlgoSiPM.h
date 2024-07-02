//=============================================================================
// AlgoSiPM.h
// This is a hit finding algorithm that makes an optical hit out of
// everything above threshold, but uses only the first peak to assign hit time
//
// Gleb Sinev, Duke, 2015
// Based on AlgoTreshold.h
//=============================================================================

#ifndef ALGOSIPM_H
#define ALGOSIPM_H

namespace fhicl {
  class ParameterSet;
}

#include "PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"

#include <string>

namespace pmtana {

  class AlgoSiPM : public PMTPulseRecoBase {

  public:
    AlgoSiPM(const fhicl::ParameterSet& pset,
             std::unique_ptr<pmtana::RiseTimeCalculatorBase> risetimecalculator = nullptr,
             const std::string name = "AlgoSiPM");

    // Implementation of PMTPulseRecoBase::Reset() method
    void Reset();

    // A method to set user-defined ADC threshold value
    //      void SetADCThreshold(double v) {_adc_thres = v;};

    // A method to set a multiplication factor to the pedestal standard deviation
    // which is used as one of two input values to define a threshold
    //      void SetNSigma(double v) {_nsigma = v;};

  protected:
    bool RecoPulse(const pmtana::Waveform_t&,
                   const pmtana::PedestalMean_t&,
                   const pmtana::PedestalSigma_t&);

    // A variable holder for a user-defined absolute ADC threshold value
    double _adc_thres;

    // Minimum width for a hit to be recorded
    int _min_width;

    // Start recording hit information after this threshold is reached
    double _2nd_thres;

    // Use this pedestal instead of the one given by the pedestal algorithm
    double _pedestal;

    // A variable holder for a multiplicative factor for the pedestal
    // standard deviation to define the threshold
    //      double _nsigma;
  };

}

#endif
