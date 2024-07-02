/**
 * \file PedAlgoUB.h
 *
 * \ingroup PulseReco
 *
 * \brief Class definition file of PedAlgoUB
 *
 * @author Kazu, Vic - Nevis 2013
 */

/** \addtogroup PulseReco

@{*/

#ifndef larana_OPTICALDETECTOR_PEDALGOUB_H
#define larana_OPTICALDETECTOR_PEDALGOUB_H

#include "PMTPedestalBase.h"
namespace fhicl {
  class ParameterSet;
}

#include "PedAlgoRmsSlider.h"

#include <string>

namespace pmtana {

  /**
   \class PedAlgoUB
   A class that calculates pedestal mean & standard deviation (here and elsewhere called as "RMS").
  */
  class PedAlgoUB : public PMTPedestalBase {

  public:
    /// Default constructor
    PedAlgoUB(const std::string name = "PedCD");

    ///Alternative ctor
    PedAlgoUB(const fhicl::ParameterSet& pset,
              //PedAlgoUB(const ::fcllite::PSet &pset,
              const std::string name = "PedAlgoUB");

  protected:
    /// Method to compute a pedestal of the input waveform using "nsample" ADC samples from "start" index.
    bool ComputePedestal(const pmtana::Waveform_t& wf,
                         pmtana::PedestalMean_t& mean_v,
                         pmtana::PedestalSigma_t& sigma_v);

  private:
    //m    PedAlgoRollingMean _beamgatealgo;
    PedAlgoRmsSlider _beamgatealgo;
    unsigned int _beam_gate_samples;
  };
}
#endif

/** @} */ // end of doxygen group
