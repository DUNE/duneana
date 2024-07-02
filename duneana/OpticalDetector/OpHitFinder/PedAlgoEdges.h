/**
 * \file PedAlgoEdges.h
 *
 * \ingroup PulseReco
 *
 * \brief Class definition file of PedAlgoEdges
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup PulseReco

@{*/

#ifndef larana_OPTICALDETECTOR_PEDALGOEDGES_H
#define larana_OPTICALDETECTOR_PEDALGOEDGES_H

#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"

#include "PMTPedestalBase.h"
namespace fhicl {
  class ParameterSet;
}

#include <string>

namespace pmtana {

  /**
   \class PedAlgoEdges
   A class that calculates pedestal mean & standard deviation (here and elsewhere called as "RMS").
  */
  class PedAlgoEdges : public PMTPedestalBase {

  public:
    /// Default constructor
    PedAlgoEdges(const std::string name = "PedEdges");

    /// Alternative ctor
    PedAlgoEdges(const fhicl::ParameterSet& pset, const std::string name = "PedEdges");
    //PedAlgoEdges(const ::fcllite::PSet &pset,const std::string name="PedEdges");

    /// enum to define algorithm options
    enum PED_METHOD {
      kHEAD = 0, ///< Use first N samples
      kTAIL,     ///< Use last N samples
      kBOTH      ///< Calculate both and use the one with smaller RMS
    };

  protected:
    /// Method to compute a pedestal of the input waveform using "nsample" ADC samples from "start" index.
    bool ComputePedestal(const pmtana::Waveform_t& wf,
                         pmtana::PedestalMean_t& mean_v,
                         pmtana::PedestalSigma_t& sigma_v);

  private:
    size_t _nsample_front; ///< # ADC sample in front to be used
    size_t _nsample_tail;  ///< # ADC sample in tail to be used
    PED_METHOD _method;    ///< Methods
  };
}
#endif

/** @} */ // end of doxygen group
