/**
 * \file PMTPedestalBase.h
 *
 * \ingroup PulseReco
 *
 * \brief Class definition file of PMTPedestalBase
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup PulseReco

@{*/
#ifndef larana_OPTICALDETECTOR_PMTPEDESTALBASE_H
#define larana_OPTICALDETECTOR_PMTPEDESTALBASE_H

// STL
#include "OpticalRecoTypes.h"
#include <string>

namespace pmtana {

  /**
   \class PMTPedestalBase
   A base class for pedestal calculation
  */
  class PMTPedestalBase {

  public:
    /// Default constructor
    PMTPedestalBase(std::string name = "noname");

    /// Default destructor
    virtual ~PMTPedestalBase();

    /// Name getter
    const std::string& Name() const;

    /// Method to compute a pedestal
    bool Evaluate(const pmtana::Waveform_t& wf);

    /// Getter of the pedestal mean value
    double Mean(size_t i) const;

    /// Getter of the pedestal standard deviation
    double Sigma(size_t i) const;

    /// Getter of the pedestal mean value
    const pmtana::PedestalMean_t& Mean() const;

    /// Getter of the pedestal standard deviation
    const pmtana::PedestalSigma_t& Sigma() const;

  protected:
    /**
       Method to compute pedestal: mean and sigma array should be filled per ADC.
       The length of each array is guaranteed to be same.
    */
    virtual bool ComputePedestal(const ::pmtana::Waveform_t& wf,
                                 pmtana::PedestalMean_t& mean_v,
                                 pmtana::PedestalSigma_t& sigma_v) = 0;

  private:
    /// Name
    std::string _name;

    /// A variable holder for pedestal mean value
    pmtana::PedestalMean_t _mean_v;

    /// A variable holder for pedestal standard deviation
    pmtana::PedestalSigma_t _sigma_v;
  };
}
#endif

/** @} */ // end of doxygen group
