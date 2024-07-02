/**
 * \file PulseRecoManager.h
 *
 * \ingroup PulseReco
 *
 * \brief Class definition file of PulseRecoManager
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup PulseReco

@{*/

#ifndef PULSERECOMANAGER_H
#define PULSERECOMANAGER_H

#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"

#include <vector>

namespace pmtana {

  class PMTPedestalBase;
  class PMTPulseRecoBase;

  /**
   \class PulseRecoManager
   A manager class of pulse reconstruction which acts as an analysis unit (inherits from ana_base).
   This class executes various pulse reconstruction algorithm which inherits from PMTPulseRecoBase
   Refer to analyze() function implementation to check how a pulse reconstruction algorithm is called.
  */
  class PulseRecoManager {

  public:
    /// Default constructor
    PulseRecoManager();

    /// Implementation of ana_base::analyze method
    bool Reconstruct(const pmtana::Waveform_t&) const;

    /// A method to set pulse reconstruction algorithm
    void AddRecoAlgo(pmtana::PMTPulseRecoBase* algo, PMTPedestalBase* ped_algo = nullptr);

    /// A method to set a choice of pedestal estimation method
    void SetDefaultPedAlgo(pmtana::PMTPedestalBase* algo);

  private:
    /// pulse reconstruction algorithm pointer
    std::vector<std::pair<pmtana::PMTPulseRecoBase*, pmtana::PMTPedestalBase*>> _reco_algo_v;

    /// ped_estimator object
    PMTPedestalBase* _ped_algo;
  };
}
#endif

/** @} */ // end of doxygen group
